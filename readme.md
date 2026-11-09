# Introduction and scope

Computational experiments in constitutive modelling and material response, using
PyTorch to capture the response of a crystal-plasticity material with
autoregressive models. At each loading increment, the model updates the material
state and predicts the response; its updated state becomes the input to the next
increment. The aim is reliable stress and internal-state histories across unseen
loading paths and crystal orientations.

The Fortran UMAT supplies the constitutive equations and reference trajectories.
Its local return-mapping problem couples 24 slip increments and 24 slip
resistances. The original work brought those equations into PyTorch to explore
automatic differentiation and Newton updates. That physics foundation supports
the autoregressive models, their training losses, and possible residual-based
corrections.

The current implementation provides:

- tensor operations for crystal rotations, hardening, plastic deformation, physical residuals, and Cauchy stress;
- a network predicting slip and slip resistance, with data loss, physics loss, and constraint penalties;
- autoregressive inference with a selectable switch from Fortran to learned updates;
- Fortran sources and Python bindings in `umat_ggcm/`, plus formula and batching tests.

The physics implementation and training workflow are a starting point. Reliable
rollout accuracy and generalization remain to be established.

## Project status and TODO

Checked items describe existing capabilities. Unchecked items describe the next
implementation and research tasks. See [the development plan](docs/development-plan.md)
for their rationale, numerical details, and acceptance checks.

### Material model and reference checks

- [x] Implement crystal rotations, hardening, plastic deformation, and Cauchy stress in PyTorch.
- [x] Add inverse slip-rate residuals, constraint penalties, and formula/batching tests.
- [x] Include the Fortran material model and Python bindings.
- [ ] Correct non-Schmid activation and enforce previous-step slip-resistance bounds below saturation.
- [ ] Standardize single-sample physics helpers and verify their `vmap` batching.
- [ ] Verify Fortran/Python residual and stress agreement across orientations and loading cases.
- [ ] Audit material parameters and initial state, including beta, and qualify the reference solver.
- [ ] Record return-mapping convergence, active-set checks, substeps, and timing.

### Data and training

- [x] Implement supervised/physics training, configurable loss logging, and model saving.
- [ ] Add reproducible trajectory datasets with independent path/orientation partitions.
- [ ] Establish a scaled, uniform-sampling baseline with joint slip/resistance updates and complete state inputs.
- [ ] Implement normalized weighted losses and a combined residual vector; diagnose nonfinite values.
- [ ] Add prioritized replay and trajectory refill; verify sampling, priorities, and beta annealing.
- [ ] Add cosine learning-rate scheduling and resumable training checkpoints.
- [ ] Extend diagnostics with loss components, priorities, sampled weights, replay beta, and validation plots.

### Rollout validation and generalization

- [x] Implement autoregression with an `alpha` switch between Fortran and learned updates.
- [ ] Evaluate full stress and internal-state histories, with explicit validation frequency and held-out tests.
- [ ] Pass one-trajectory and small-set rollout diagnostics before expanding model capacity.
- [ ] Compare absolute outputs, increments, trial-residual inputs, and an elastic trial branch.
- [ ] Compare one-step and multistep training, including scaled state/stress losses with full state feedback.
- [ ] Demonstrate generalization to unseen paths, unseen orientations, and both together.

### Further development

- [ ] Compare uniform and prioritized sampling under matched training budgets.
- [ ] Investigate residual-guided correction with acceptance checks and a reference fallback.
- [ ] Study constraint formulations and loss balancing, including barriers, warmup, gradient clipping, and L2 regularization.
- [ ] Compare original/inverse slip-rate residuals and residuals between converged states.
- [ ] Verify automatic-differentiation residual Jacobians and evaluate Newton corrections to predicted updates.
- [ ] Extend loading histories to unloading, direction changes, and variable time increments.
- [ ] Assess model capacity, float32 accuracy, and device portability after reproducing the baseline.
- [ ] Verify derivatives of the complete update and consistent tangents before FE integration.

## Run

Run commands from the repository root. The code uses PyTorch, NumPy, einops,
PyYAML, and TensorBoard; the current packaging does not declare a complete
dependency environment. Fortran-backed imports also require a compiled
`umat_ggcm.pyUMATLink` extension and a compatible Fortran/f2py toolchain.
The existing build recipe is in [umat_ggcm/Makefile](umat_ggcm/Makefile):

```bash
make -C umat_ggcm compile
```

Set `sims_path` in [config_file.yaml](config_file.yaml) to an existing simulation
pickle, choose the loss coefficients and training settings, then run:

```bash
python driver.py --config-file config_file.yaml
```

The sample configuration contains a machine-specific data path. A reproducible
dataset and environment are part of the planned work above.

## Testing

With the dependencies and Fortran extension available:

```bash
python -m unittest tests.test_umat -v
```

The [rotation](tests/rotation_matrix_test/README.md) and
[tensor-contraction](tests/einsum_test/README.md) comparisons have separate
Fortran build instructions. For coverage:

```bash
coverage run -m unittest discover
coverage report -m
```

## Documentation

The formulation notes are in [docs/main.tex](docs/main.tex) and require LaTeX to
render. The [development plan](docs/development-plan.md) explains the physical
checks, training experiments, and rollout criteria needed to develop reliable
autoregressive material models.
