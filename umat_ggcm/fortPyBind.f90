!
! initialize the internal variables
!
subroutine pySDVINI(statev, nstatv)
!
double precision statev(nstatv)
integer         nstatv
!
double precision coords(3)
integer         ncrds, noel, npt, layer, kspt
!
! assigning the variables for subroutine call
ncrds = 3
noel  = 1
npt   = 1
layer = 1
kspt  = 1

call SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,LAYER,KSPT)
end subroutine


!
! UMAT Subroutine 
! Updates the stress, tangent (dP/dF) and internal variables
!
subroutine fortUMAT(F1, F2, t1, t2, iniTemp, &
     sigma, C, XI, angles, nstatv)
! implicit none
! passed in variables
double precision F1(3,3), F2(3,3), XI(nstatv)
double precision sigma(6), C(6,6)
double precision t1, t2
double precision iniTemp
double precision angles(3)
integer nIter, nStep, nstatv, matID

! internal variables
integer ntens, ndi, nshr, noel, npt, layer
integer kspt, kstep, kinc
integer nprops
CHARACTER*80 CMNAME
double precision rpl, drpldt, ddsddt(6), drplde(6), stran(6), dstran(6)
double precision time(2), dtime
double precision sse, spd, scd
double precision temp, dtemp, pnewdt, celent
double precision predef(1), dpred(1), coords(3,3), drot(3,3), props(1)
!
!write(*,*), "fortUMAT begining"
!
ntens  = 6
ndi    = 1
nshr   = 1
noel   = 1
npt    = 1
layer  = 1
kspt   = 1
nprops = 1   ! this is hard-coded and should be changed if we want to
             ! do a different UMAT. props array holds the material-id 
             ! for the current implementation of the UAMT. So it can 
             ! be set to 1.
! internal heat generation and its variations, small strain and its 
! increment
rpl    = 0.d0
drpldt = 0.d0
ddsddt = 0.d0
drplde = 0.d0
stran  = 0.d0
dstran = 0.d0
! time and the begining and the end of the step and its increment
time(1) = t1
time(2) = t2
dtime   = t2 - t1
! temperature and its increment
! since there is no plan for doing termomechanical, increment of 
! temperature is hard-coded and set to zero
dtemp = 0.d0
temp  = iniTemp
! setting the material id through props array
! in this branch, we directly provide the angles and we don't need
! the matid in order to grab these angles (1, 2 and 3) from the
! include file. Because of this, this value is only a dummy value
props(1) = dble(1)
! specific elastic strain, plastic dissipation and creep dissipation
sse = 0.d0
spd = 0.d0
scd = 0.d0
predef(1) = 0.d0
dpred (1) = 0.d0
coords    = 0.d0
drot      = 0.d0
pnewdt    = 1.d0
celent    = 1.d0
! step and increment
kstep = 1
kinc  = 1
!
! sigma : Cauchy stress tensor
! C     : variation of first Piola Kirchhoff w.r.t. deformation gradient
! F1    : deformation gradient at the begining of step
! F2    : deformation gradient at the end of step 
! t1,t2 : time at the begining and the end of step

! debug begin
!write(*,*) 'linker information begin'
!write(*,*) 'nstatv = ', nstatv, nIter
!write(*,*) 'linker information end'
! debug end
! call to the umat
call   UMAT(sigma,XI,C,sse,spd,scd,rpl,ddsddt,                         &
       drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,dpred,  &
       cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,   &
       celent,F1,F2,noel,npt,layer,kspt,kstep,kinc,angles)


end subroutine      



