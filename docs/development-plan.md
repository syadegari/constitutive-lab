# Development plan for autoregressive material models

The active checklist is in [the README](../readme.md#project-status-and-todo).
This document explains the research questions, implementation details, and
acceptance checks behind it.

## Objective and material update

Capture crystal-plasticity material response with an autoregressive model that
predicts stress and internal-state histories under prescribed loading. The
success criterion is reliable prediction on unseen loading histories, unseen
crystal orientations, and both together.

Write a material step schematically as

```text
state_(n+1) = M(state_n, F_n, F_(n+1), dt, orientation)
stress_(n+1) = stress(F_(n+1), state_(n+1), orientation)
```

The evolving physical state includes plastic deformation Fp, accumulated slip,
slip resistance, and material beta. The learned part can predict joint slip and
resistance updates, with mechanics updating the remaining state and stress.
During autoregressive evaluation, the predicted state is carried into the next
increment. Whether direct prediction or prediction followed by physical
correction works better remains an experimental question.

The local constitutive problem couples 24 directed slip increments and 24 slip
resistances. Expressing its equations in PyTorch supports physical training
losses, automatic-differentiation Jacobians, and possible Newton corrections.
Reference qualification and differentiation support the material-response goal.
A successful local solve or a low one-step error does not establish accurate
learned rollouts.

## Current foundation and immediate corrections

[umat/umat.py](../umat/umat.py) implements rotations, hardening, plastic
deformation, driving forces, Cauchy stress reconstruction, and autoregression
with an `alpha` switch from Fortran to learned updates.
[umat/model.py](../umat/model.py) predicts 24 accumulated slips and 24 slip
resistances from orientation, endpoint deformation gradients, and previous
slip/resistance. [umat/train.py](../umat/train.py) supplies data and physical
losses, penalties, logging, and a final model checkpoint. Formula and batching
tests and Fortran comparison tools exist under `tests/`.

Resolve the following before interpreting training comparisons:

- Use the effective non-Schmid activation threshold `g1 > s1 - non_schmid_stress`
  in `get_rII`, in place of the current `g1 > s1` condition.
- Constrain slip resistance against its previous-step value, rather than the
  initial material minimum, and keep a precision-aware margin below saturation.
  Apply compatible constraints in training and inference, checking their effect
  on gradients. Nonnegative slip increments include zero in inactive systems.
- Diagnose nonfinite physical losses instead of replacing NaN with zero. Save
  the offending inputs, outputs, and residual components.
- Standardize single-sample rotation and driving-force functions, then batch with
  `vmap`. Update callers and numerical comparisons together.

The inverse slip-rate residual is called during training; the supplied default
configuration sets the physics-loss coefficient to zero. Having an equation
implemented and demonstrating useful optimization with it are separate steps.
The existing tests cover specific formulas and cases; broader reference
agreement still needs explicit checks.

## Reference equations, parameters, and data

Use the Fortran model in [umat_ggcm/](../umat_ggcm/README.md) as the numerical
reference, with observable convergence and documented material parameters.

For identical orientations, old states, deformation histories, and increments,
compare `r_I`, `r_II`, and reconstructed Cauchy stress. Evaluate the trial residual
`(r_II)_p` with zero slip increment and unchanged resistance as well. Check both
sides of activation thresholds, constraint boundaries, and multiple loading
cases. Set explicit tolerances and make failures fail the test.

Record return-mapping iterations, residual norms, Newton-cap exits, active-set
checks, substeps, and timing. An accepted update is not proof of convergence.
Compare iteration caps at unchanged tolerance, and distinguish convergence of
the discrete update from temporal-refinement accuracy. Qualify the reference
before treating generated states as training labels.

Reconcile Python and Fortran initial states, phase constants, state layouts, and
parameter provenance. In particular, review the initial material beta against
the intended calibration, including the candidate value 0.0056. Resolve any
parameter-table/code differences before claiming reproduction of a calibration.
Keep datasets for different material variants distinct.

Create a reproducible trajectory dataset with a small runnable example and a
complete generation command. Store orientation, deformation gradients, times,
physical states, stresses, solver diagnostics, and material/source identifiers.
Check positive deformation-gradient determinant along the entire loading path.
HDF5 is a suitable container; the existing training entry point currently reads
a simulation pickle, so configuration and loading need a consistent interface.

Freeze independent loading-path and orientation partitions. Keep complete
trajectories and every window drawn from them within one partition. Check for
duplicate or physically equivalent orientations, and fit normalization on
training data only. Distinguish unseen endpoints within proportional loading
from genuinely new path shapes such as unloading and direction changes.

## Baseline training and physical losses

Start with a scaled, uniform-sampling offline baseline and fixed partitions.
First fit one trajectory, then a small set, while measuring full rollouts.
These are implementation diagnostics, not the generalization endpoint.

Retain joint slip and resistance unknowns. Compare absolute outputs with
increments, and supply a justified material state including Fp, material beta,
loading, orientation, and dt. Merely storing a field in a batch does not ensure
that the network uses it. Accumulated slip alone generally does not determine
Fp for noncommuting deformation histories.

A reduced resistance representation or signed slip-pair parameterization needs
physical justification. Opposite accumulated slips can both be positive after
reversal. Latent hardening can change resistance on a system with no active slip;
resistance is not an invertible substitute for the individual slip increments.

Form a combined physical residual `[r_I, r_II]` while logging its components
separately. Scale targets and residual components by justified physical scales
before tuning coefficients. Concatenating residuals changes the loss scale
relative to adding two separate component means.

Define weighted loss as a weighted mean over samples with a separate mean over
output components. Normalizing weights to sum one and then averaging over the
batch again introduces an unintended factor of `1/batch_size`. Verify uniform
weights, nonuniform weights, and gradient flow.

Compare barriers or alternative constraint parameterizations when projection
removes useful physical gradients. Study physics-loss warmup, per-loss gradient
clipping, and L2 regularization under controlled settings. Learning-rate
scheduling is separate from physics-loss warmup. Compare original and inverse
slip-rate equations through an explicit configuration choice, and examine
residuals along interpolated states between converged endpoints.

## Replay, scheduling, and diagnostics

Add prioritized replay as a controlled comparison with uniform sampling at a
matched data and optimizer-update budget. Store consecutive material states,
refill from reference trajectories, and update priorities from defined data,
physical, or combined errors. Sampling policy is an experimental factor.

Verify partial and full buffers, circular replacement, priority updates,
undersized-buffer behavior, and finite positive sampling probabilities. Choose
sampling with or without replacement explicitly and ensure the importance
weights match that policy.

Distinguish replay's importance-sampling beta from material beta. Bound replay
beta by one and specify a positive annealing horizon. Define how an open-ended
run advances the schedule; adaptive annealing is a separate candidate. Log
sampled weights as well as priority statistics and quantiles.

Add optional cosine learning-rate scheduling with documented initial rate,
restart boundaries, declining peaks, and final-rate behavior. Preserve continuous
logging steps across training and refill cycles, and close writers on exit.
Save model, optimizer, scheduler, replay contents/priorities/cursor, step counter,
and RNG states so checkpoint continuation can be verified.

Log total and component losses, physical residuals, gradient diagnostics,
learning rate, priorities, weights, replay beta, and validation errors. Plot
stress components and slip/resistance histories, including worst cases. Use
identity for an undeformed initial deformation gradient and define normalization
when reference norms vanish. Set validation frequency explicitly and implement
a separate final held-out evaluation.

## Rollouts, trial features, and multistep training

Measure reference-state one-step predictions and fully predicted-state rollouts
separately. Report trajectory errors for stress, slip, resistance, Fp, and beta,
constraint violations, residuals, nonfinite/failure rates, and computation time.
Include a zero-increment persistence baseline. Preserve the `alpha` switch for
controlled transitions from reference to learned inference, while evaluating
fully autoregressive trajectories for the main acceptance criterion.

Diagnose the first divergence, positive slip before activation, opposite-pair
activity, and amplification through state or trial-residual feedback. An elastic
trial branch can remove spurious elastic increments without fixing plastic
predictions. Treat it and trial-vector inputs as controlled comparisons rather
than assumed remedies. Preserve outliers and failed trajectories.

Compare teacher-forced one-step training with short contiguous autoregressive
windows. Initialize a window from one reference state and propagate every
predicted physical state, including Fp and beta. Replacing some fields with
reference values conceals part of the accumulated error. Independent replay
pairs need a trajectory-aware extension to supply contiguous windows.

Use all-step losses as an initial comparison, and distinguish differentiable
unrolling from detached feedback. Compare increment-only losses with scaled
state/stress objectives under matched initialization, horizon, and target
exposure. Increase horizons or capacity in response to diagnostic evidence.
Improved worst-case behavior need not imply improved median error or reliable
generalization; report both.

## Residual differentiation and correction

For fixed old state and loading, let `z = [delta_gamma, s_new]` and
`R(z) = [r_I, r_II]`. A Newton iteration uses

```text
J = dR/dz
solve J * delta_z = -R(z)
z_next = z + delta_z
```

PyTorch provides residual-Jacobian construction through
[`torch.func.jacrev`](https://docs.pytorch.org/docs/stable/generated/torch.func.jacrev.html)
and a linear-system solve through
[`torch.linalg.solve`](https://docs.pytorch.org/docs/stable/generated/torch.linalg.solve.html).
Solve the linear system directly. A scalar residual-loss gradient is not the
full residual Jacobian required by this update.

Assemble the corrected residual and compare its Jacobian against directional
finite differences away from activation boundaries. Test elastic and plastic
states, conditioning, residual decrease, constraints, and cost, starting from
reference and learned trial states. The Fortran active-set solver need not use
all 48 unknowns in every solve; active systems and substeps require attention.
Clipping and switching affect derivatives. Use damping or acceptance checks
where needed; automatic differentiation alone gives no convergence guarantee.

Compare a learned predictor followed by physical correction, or a learned
residual-guided correction, with direct prediction. Define acceptance checks,
iteration budgets, and a reference fallback; report fallback frequency and cost.
A learned fixed point must still satisfy the material equations. A verified
PyTorch Newton solver remains an experiment, rather than an existing capability
or a prerequisite for every learning comparison.

## Broader loading, precision, and tangents

Extend loading to unloading, direction changes, and several time increments.
The constitutive law depends explicitly on dt; a constant or unused dt feature
cannot establish rate generalization. Assess temporal refinement rather than
demand exact equality between one large step and several smaller steps.

Use float64 as the correctness baseline. Then assess float32 or mixed precision
against application-relevant errors and failure rates. Make device placement and
CPU Fortran/device training boundaries explicit. Increase model capacity only
after the baseline and rollout diagnostics justify it.

For consistent tangents, differentiate the complete material update with the
old state fixed; holding Fp fixed produces an elastic partial derivative.
Distinguish the residual Jacobian `dR/dz` from the derivative of stress with
respect to deformation. For a smooth converged solve `R(z,F)=0`,

```text
dz/dF = -solve(R_z, R_F)
dP/dF = partial(P)/partial(F) + partial(P)/partial(z) * dz/dF
```

Verify derivatives against independent stress perturbations and match stress
measures and objective-rate conversions to the intended interface. Review
active-set changes, substepping, unloading, direction changes, and geometric
conversion terms before FE use. Forward rollout accuracy and tangent accuracy
are separate acceptance questions.

## Experiment records and acceptance

Freeze protocols, input values or deterministic generation rules, partitions,
seeds, variants, and evaluation criteria before running comparisons. Record
source revisions and local patches, dependency/compiler versions, precision,
device/thread settings, normalization, checkpoints, failures, and a rerun command.
A seed alone does not establish reproducibility.

Establish fresh-process numerical repeatability and checkpoint continuation
before long runs. Compare canonical numerical arrays; timing and archive bytes
are separate records. Do not promise identical numerics across hardware or
library changes.

Select models with training and validation data. Reserve test trajectories for
predeclared evaluations and document any prior use. Report seen-path/seen-
orientation diagnostics, unseen-path results, unseen-orientation results, and
joint generalization separately. Application tolerances need justification;
training fit alone is not an acceptance criterion.

The proposed sequence is physical/reference checks, a reproducible offline
baseline, rollout diagnosis, then controlled comparisons of sampling, loss,
features, and correction. Tools can be developed earlier when they resolve a
specific question. Barrier formulations, model capacity, and optimizer choices
should be driven by observed failures and useful comparisons.
