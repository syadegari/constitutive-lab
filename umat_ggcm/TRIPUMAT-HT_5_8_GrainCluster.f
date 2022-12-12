c .......1.........2.........3.........4.........5.........6.........7..
c
c **********************************************************************
c
c   subroutine to initialize the state-dependent variables (statev)
c
c **********************************************************************
c
      subroutine SDVINI(statev,coords,nstatv,ncrds,noel,npt,layer,kspt)
c
c  include external files
      include 'aba_param.inc'
      include 'constants.for'
c
c  declaration of variables and parameters
      double precision statev(nstatv),coords(ncrds)
      integer i1,j1,k1,i2,j2,k2,i,j,k,l,m,n
c
      statev = 0.d0
c      write (*,*) 'SDVINI is called'
c
c  define the initial values for FCC austenite
      if (nstatv .eq. 151) then
c
c    initial plastic deformation gradient
         do i1 = 1,3
            do j1 = 1,3
               k1 = j1 + 3*(i1 - 1)
               if (j1 .eq. i1) then
                  statev(k1) = 1.d0
               else
                  statev(k1) = 0.d0
               endif
            enddo
         enddo
c
c    initial transformation deformation gradient
         do i2 = 1,3
            do j2 = 1,3
               k2 = j2 + 3*(i2 - 1) + 9
               if (j2 .eq. i2) then
                  statev(k2) = 1.d0
               else
                  statev(k2) = 0.d0
               endif
            enddo
         enddo
c
c    accumulative plastic strain, slip resistance parameter, dislocation 
c    density parameter, martensite volume fraction, total volume fraction
c    initial slip resistance = 0.188 for 1.4 wt.%C (grain size 2 micron)
c                              0.144 for 1.0 wt.%C
c                              0.172 for 1.2 wt.%C
c                              0.222 for 1.8 wt.%C
c                              0.254 for 2.2 wt.%C
c                              0.140 for aust. grain size 8 micron
c                              0.164 for aust. grain size 4 micron
c                              0.218 for aust. grain size 1 micron
c                              0.273 for aust. grain size 0.5 micron
         do i = 1,24
            j = i + 18
	    k = i + 42
            statev(j) = 0.d0    !initial plastic strain
	    statev(k) = s0_A    ! 1.88d-1 !initial slip resistance parameter
         enddo 
         statev(67) = 0.d0      !initial dislocation density
         do i = 1,24
            l = i + 67
            statev(l) = 0.d0    !initial martensite volume fraction
         enddo
         statev(92) = 0.d0      !initial total martensite volume fraction
         do i = 1,24
	    m = i + 92
	    n = i + 116
	    statev(m) = 0.d0    !slip system, tag 0 = inactive, 1 = active
	    statev(n) = 0.d0    !transformation system tag, idem
         enddo
         statev(141) = 0.d0     !initial work
c
c  define the initial values for BCC ferrite
      else

c    initial plastic deformation gradient
         do i1 = 1,3
            do j1 = 1,3
               k1 = j1 + 3*(i1 - 1)
               if (j1 .eq. i1) then
                  statev(k1) = 1.d0
               else
                  statev(k1) = 0.d0
               endif
            enddo
         enddo
c
c    accumulative plastic strain, slip resistance parameter, dislocation 
c    density parameter
c    initial slip resistance = 0.158 for original matrix
c                              0.072 for soft matrix (low yield stress)
c                              0.216 for hard matrix (high yield stress)
         do i = 1,24
            j = i + 9
	    k = i + 33
            statev(j) = 0.d0    !initial plastic strain
	    statev(k) = s0_F    !1.58d-1 !initial slip resistance parameter
         enddo 
         statev(58) = 0.d0      !initial dislocation density
         do i = 1,24
	    m = i + 58
	    statev(m) = 0.d0    !slip system, tag 0 = inactive, 1 = active
         enddo
         statev(83) = 0.d0      !initial work
c     
c     return the initial values
      endif
      return
      end
c .......1.........2.........3.........4.........5.........6.........7..
c
c **********************************************************************
c
c  user-material subroutine for elasto-plasto-transformation model
c  in tranformation-induced plasticity steels
c
c  additional state-variable: total work                    (14.09.2006)
c
c                                                             version 04
c **********************************************************************
c
      subroutine UMAT(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     & drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,dpred,
     & cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     & celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc,angles)
c
c  header for standard user-subroutine for ABAQUS
      include 'aba_param.inc'
c
c  declation of variables
      character*80 cmname
      double precision stress(ntens),statev(nstatv),ddsdde(6,6),
     & ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),time(2),
     & predef(1),dpred(1),props(nprops),coords(3),drot(3,3),dfgrd0(3,3),
     & dfgrd1(3,3),sse,spd,scd,rpl,drpldt,temp,dtemp,pnewdt,celent,
     & statevA(nstatv),statevF(nstatv),dtime,angles(3)
      integer nstatvA,nstatvF,i
      logical debugoption, dbg
      Integer myTmpVar
	     
      debugoption=.FALSE.
      dbg        =.FALSE. 

c The following will create an infinite loop,
c giving sufficient time to attach the debugger to
c the running process once the program is running
      do While(dbg)
         myTmpVar = 1
      enddo
      if (debugoption) then
      write (*,*) '****************************'
      write (*,*) 'variables before going into the UMAT '
      write (*,*) 'nstatv'       
      write (*,*) nstatv
      write (*,*) 'STATEV'
      write (*,'(E23.15E3)') (statev(i),i=1,nstatv)
      write (*,*) 'STRESS'
      write (*,'(E23.15E3)') (stress(i),i=1,6)
      write (*,*) 'NPT, KSTEP, KINC, TIME and DTIME'            
      write (*,*) npt,kstep,kinc
      write (*,'(E23.15E3)') time
      write (*,'(E23.15E3)') dtime
      write (*,*) 'Temp and dTemp'
      write (*,'(E23.15E3)') temp
      write (*,'(E23.15E3)') dtemp
      write (*,*) 'F0 and F1' 
      write (*, *) 'end of input variables'
      endif
c
c  call to UMATFCC if austenite
c  
      nstatvF = nstatv
      do i = 1, nstatvF
         statevF(i) = statev(i)
      enddo
      call UMATBCC(stress,statevF,ddsdde,sse,spd,scd,rpl,ddsddt,
     &     drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,
     &     predef,dpred,cmname,ndi,nshr,ntens,nstatvF,props,
     &     nprops,coords,drot,pnewdt,celent,dfgrd0,dfgrd1,
     &     noel,npt,layer,kspt,kstep,kinc,angles)
      do i = 1, nstatvF
         statev(i) = statevF(i)
      enddo
c
c  return to main program ABAQUS
c
      return
      end


c
c .......1.........2.........3.........4.........5.........6.........7..
c
c **********************************************************************
c
c  subroutine to compute the state variable update (DeltaGamma, DeltaXi,
c  and SlipRes) using Newton-Raphson iteration (return-mapping algorithm
c  procedure)
c
c **********************************************************************
c
      subroutine RETURNMAPPINGBCC(DeltaGammaA,SlipResA,nGammaA,PDF,
     & PDFCrit,idGammaA,FP0,F1,BetaA,DeltaT,SubStepping,ConTangent,
     & theta,
     & ElasStif,SlipSys,ShearMod)
      include 'constants.for'
c
c  declaration of variables
      double precision DeltaGammaA(nGammaA),SlipResA(24),PDF(24),
     & PDFCrit(24),I2(3,3),I4(3,3,3,3),ElasStif(3,3,3,3),
     & SlipSys(24,3,3),ShearMod,FP0(3,3),FPA(3,3),F1(3,3),InvFPA(3,3),
     & EEA(3,3),FEA(3,3),Temp,TransFEA(3,3),CEA(3,3),PK2(3,3),
     & InTransFPA(3,3),MTemp(3,3),PK1(3,3),CPK2(3,3),CPK1(3,3),Beta,
     & BetaA,DeltaSlipRes(24),GMod(24),HMod(24,24),SchmidP(24),PDFCold,
     & Omega,PDFTherm,M2Temp(3,3,3,3),Res2(24),GammaDot0,mExp,DeltaT,
     & SlipRes(24),L2Res,locToler,L2ResBef,NSchmidP(24)
      integer nGammaA,idGammaA(nGammaA),i,j,k,l,nSize,errorstate,m,iP,
     & ierr
      logical PlotDebug,PlotResult,SubStepping,ConTangent
c
      double precision, allocatable::AVar(:),PDFRed(:),PDFCritR(:),
     &     Res1(:),MVar(:,:),ARes(:),AVarP(:),AResP(:),radum(:),
     &     dVardRes(:,:)
      integer, allocatable::rpdum(:),ipiv(:)
c        
c  variables and parameters for thermal deformation gradient
c
      double precision theta, det_th,Fth(3,3),Fth_inv(3,3),Fth_tr(3,3)
      double precision Fth_t(3,3),Fth_invt(3,3)
      double precision DFthDt(3,3)
      double precision Ftr_t(3,3), Fp_t(3,3)
      double precision eye(3,3)
      double precision Phi,Rho,ThetaT,LambdaT
c
c  values: sInf = 0.412 for original matrix
c                 0.178 for soft matrix1 (low max stress)
c                 0.442 for soft matrix2 (med max stress)
c                 0.404 for hard matrix1 (med max stress)
c                 0.758 for hard matrix2 (high max stress)
      double precision, parameter::Toler = 1.d-11
      double precision, parameter::Pertb = 1.d-8
      double precision::sInf  = sInf_F     ! 4.12d-1 ! (*)
      integer, parameter::iterMax = 8
c
c  parameters in kinetic model and thermal softening (calibrated for BCC 
c  ferrite)
c  values: Omega     = 7
c          PDFTherm  = 0.01
c          GammaDot0 = 0.001
c          mExp      = 0.02
      Omega     = W_F           !7.d0
c     PDFTherm  = 1.d-2
      GammaDot0 = GammaDot0_F   ! 1.d-3
      mExp      = mExp_F        ! 2.d-2
      Phi       = Phi_F         ! 4.27350427d-9
      Rho       = Rho_0         ! 78.d2 (kg/m^3)
c  construct the array of variables
      I2 = 0.d0
      do i = 1,3
         I2(i,i) = 1.d0
      enddo
      nSize = nGammaA + 24
      allocate(AVar(nSize),PDFRed(nGammaA),
     &     PDFCritR(nGammaA),Res1(nGammaA),MVar(nSize,nSize),
     &     ARes(nSize),AVarP(nSize),AResP(nSize),radum(nSize),
     &     rpdum(nSize),ipiv(nSize),dVardRes(nSize,nSize),
     &     stat = errorstate)
      if (errorstate .gt. 0) then
         write(*,*)'UMATBCC: FAILED TO ALLOCATE 3'
         stop
      endif
c   SlipRes is s_n  and is fixed
c   SlipResA is updated in each iteration
      SlipRes  = SlipResA 
      Beta     = BetaA
      do i = 1,nGammaA
         AVar(i) = DeltaGammaA(i)
      enddo
      do i = 1,24
         AVar(i + nGammaA) = SlipResA(i)
      enddo
c
c
c  starting the Newton-Raphson iteration                             (!)
      SFact = 1.0d0
      if (ConTangent) SFact = 2.0d0
      L2ResBef = 1.d2
      do iT = 1,iterMax
c
c  convert the array of variables
      do i = 1,nGammaA
         if (AVar(i) .lt. 0.d0) then
	    AVar(i) = 0.d0
         endif
         DeltaGammaA(i) = AVar(i)
      enddo
      do i = 1,24
         if (AVar(i + nGammaA) .gt. sInf) then
	    AVar(i + nGammaA) = sInf
         endif
         SlipResA(i) = AVar(i + nGammaA)
      enddo
c
c  calculate the plastic deformation gradients
      call PLASTICDEFGRADBCC(DeltaGammaA,idGammaA,SlipSys,nGammaA,FP0,
     &     FPA)
c
c  calculate the transformation, plastic and thermal deformation gradients
      call comp_Fth(theta,0.d0,'F',Fth,1.d0) ! Fth
      call INVERSEMATRIX3X3(Fth, Fth_inv) !Fth^-1
      call DETERMINANTMATRIX3X3(Fth,det_th) ! det(Fth)
      Fth_t=transpose(Fth)  ! Fth_t
      Fth_invt=transpose(Fth_inv) ! Fth_t^-1
      FPA=matmul(Fth, FPA) !  Fp <- Fth*Fp
c
c  calculate the second and first Piola-Kirchhoff stresses
      call INVERSEMATRIX3X3(FPA,InvFPA) ! (Fth*Fp)^-1
      call MULTIMATRIX3X3(F1,InvFPA,FEA) ! Fe
      call TRANSPOSEMATRIX3X3(FEA,TransFEA) ! Fe_t
      call MULTIMATRIX3X3(TransFEA,FEA,CEA) ! Ce = Fe_t*Fe
      call ADDMATRIX3X3(5.d-1,CEA,-5.d-1,I2,EEA) ! Ee
      call TENSOR4TENSOR2(ElasStif,EEA,PK2) ! S =C*Ee
c
c  evaluate the Schmid stress for plasticity
      call MULTIMATRIX3X3(CEA,PK2,CPK2) ! Ce*S
      CPK2=matmul(CPK2,Fth_invt) ! Ce*S*Fth^-t
      CPK2=matmul(Fth_t,CPK2) ! Fth_t*Ce*S*Fth^-t
      CPK2=det_th*CPK2  ! Jth*Fth_t*Ce*S*Fth^-t
      SchmidP = 0.d0
      do k = 1,24
         Temp = 0.d0
         do i = 1,3
            do j = 1,3
               Temp = Temp + SlipSys(k,i,j)*CPK2(i,j)
            enddo
         enddo
         SchmidP(k) = Temp
      enddo
c
c  compute the corresponding non-Schmid stress
      call NONSCHMIDSTRESSBCC(SchmidP,NSchmidP)
c
c  compute the hardening modulus and shear hardening modulus
      call HARDENINGMODULUSBCC(SlipResA,HMod,GMod,
     &     ShearMod)
c      
c  calculate the correction of slip resistance parameter and the
c  dislocation density parameter
      do i = 1,24
         Temp = 0.d0  
         do j = 1,nGammaA
	    Temp = Temp + HMod(i,idGammaA(j))*DeltaGammaA(j)
         enddo
         DeltaSlipRes(i) = Temp
      enddo
      BetaA = Beta
      do i = 1,nGammaA
         BetaA = BetaA + GMod(idGammaA(i))*DeltaGammaA(i)
      enddo
c
      PDFTherm = Phi*theta*Rho
      do i = 1,24
         PDFCold = Omega*ShearMod*BetaA*GMod(i)
         PDF(i)  = SchmidP(i) + PDFTherm - PDFCold
      enddo
c  calculate the effective slip resistance parameter by including non-
c  Schmid stress
      do i = 1,24
         PDFCrit(i) = SlipResA(i) - NSchmidP(i)
      enddo
c
c  define the resized array of driving forces and slip resistance
      PDFRed   = 0.d0
      PDFCritR = 0.d0
      do i = 1,nGammaA
         PDFRed(i)   = PDF(idGammaA(i))
         PDFCritR(i) = PDFCrit(idGammaA(i))
      enddo
c
c  compute the residual coming from the plastic kinetic relation
      Res1 = 0.d0
      do i = 1,nGammaA
         if (DeltaGammaA(i) .lt. 0.d0) DeltaGammaA(i) = 0.d0
         Temp    = (DeltaGammaA(i)/(GammaDot0*DeltaT) + 1.d0)**mExp
         Res1(i) = PDFRed(i) - PDFCritR(i)*Temp
      enddo
c
c  compute the residual related to the slip resistance parameter
      Res2 = 0.d0
      do i = 1,24
         Res2(i) = SlipResA(i) - SlipRes(i) - DeltaSlipRes(i)
      enddo
c     
c     construct the array of residuals
      ARes = 0.d0
      do i = 1,nGammaA
         ARes(i) = Res1(i)
      enddo
      do i = 1,24
         ARes(i + nGammaA) = Res2(i)
      enddo
c
c  set the convergent criteria using the L2-norm convergence
      L2Res    = 0.d0
      do i = 1,nSize
         L2Res = L2Res + ARes(i)*ARes(i)
      enddo      
      L2Res = sqrt(L2Res)
      locToler = Toler
      if (L2Res .lt. locToler) goto 67
c     if ((L2Res .gt. 1.d1) .or. (L2Res .gt. SFact*L2ResBef)) then
      if (L2Res .gt. 1.d1) then
         SubStepping = .true.
         goto 67
      endif
      if (iT .eq. iterMax) then
         goto 67
      endif
      L2ResBef = L2Res
c
*  remark:
*  if relative tolerance is used, copy the following:
*    if (iT .eq. 1) locToler = Toler*L2Res
*    if (locToler .eq. 0.d0) locToler = Toler
*  otherwise:
*    locToler = Toler
*
c  computing the Jacobian matrix for Newton-Raphson iteration using
c  numerical perturbation technique
c
c    define the perturbation on the array of variables
      do iP = 1,nSize
         AVarP = AVar
         AVarP(iP) = AVar(iP) + Pertb
c
c    convert the array of variables
         do i = 1,nGammaA
	    if (AVarP(i) .lt. 0.d0) then
               AVarP(i) = 0.d0
	    endif
            DeltaGammaA(i) = AVarP(i)
         enddo
         do i = 1,24
	    if (AVarP(i + nGammaA) .gt. sInf) then
               AVarP(i + nGammaA) = sInf
	    endif
	    SlipResA(i) = AVarP(i + nGammaA)
         enddo
c    calculate the plastic deformation gradients
         call PLASTICDEFGRADBCC(DeltaGammaA,idGammaA,
     &        SlipSys,nGammaA,FP0, FPA)
         call comp_Fth(theta,0.d0,'F',Fth,1.d0) ! Fth
         call INVERSEMATRIX3X3(Fth, Fth_inv) !Fth^-1
         call DETERMINANTMATRIX3X3(Fth,det_th) ! det(Fth)
         Fth_t=transpose(Fth)   ! Fth_t
         Fth_invt=transpose(Fth_inv) ! Fth_t^-1
         FPA=matmul(Fth, FPA)   !  Fp <- Fth*Fp
c
c    calculate the second and first Piola-Kirchhoff stresses
         call INVERSEMATRIX3X3(FPA,InvFPA)
         call MULTIMATRIX3X3(F1,InvFPA,FEA)
         call TRANSPOSEMATRIX3X3(FEA,TransFEA)
         call MULTIMATRIX3X3(TransFEA,FEA,CEA)
         call ADDMATRIX3X3(5.d-1,CEA,-5.d-1,I2,EEA)
         call TENSOR4TENSOR2(ElasStif,EEA,PK2)
c
c    evaluate the Schmid stress for plasticity
         call MULTIMATRIX3X3(CEA,PK2,CPK2)
         CPK2=matmul(CPK2,Fth_invt) ! Ce*S*Fth^-t
         CPK2=matmul(Fth_t,CPK2) ! Fth_t*Ce*S*Fth^-t
         CPK2=det_th*CPK2       ! Jth*Fth_t*Ce*S*Fth^-t
         SchmidP = 0.d0
         do k = 1,24
	    Temp = 0.d0
	    do i = 1,3
               do j = 1,3
                  Temp = Temp + SlipSys(k,i,j)*CPK2(i,j)
               enddo
            enddo
	    SchmidP(k) = Temp
         enddo
c    compute the corresponding non-Schmid stress
         call NONSCHMIDSTRESSBCC(SchmidP,NSchmidP)
c    compute the hardening modulus and shear hardening modulus
         call HARDENINGMODULUSBCC(SlipResA,HMod,GMod,
     &        ShearMod)
c    calculate the correction of slip resistance parameter and the
c    dislocation density parameter
         do i = 1,24
            Temp = 0.d0  
	    do j = 1,nGammaA
               Temp = Temp + HMod(i,idGammaA(j))*DeltaGammaA(j)
	    enddo
	    DeltaSlipRes(i) = Temp
         enddo
         BetaA = Beta
         do i = 1,nGammaA
	    BetaA = BetaA + GMod(idGammaA(i))*DeltaGammaA(i)
         enddo
c    evaluate the total driving force for plasticity and effective slip
c    resistance parameter
         PDFTherm = Phi*theta*Rho
         do i = 1,24
            PDFCold    = Omega*ShearMod*BetaA*GMod(i)
	    PDF(i)     = SchmidP(i) + PDFTherm - PDFCold
	    PDFCrit(i) = SlipResA(i) - NSchmidP(i)
         enddo
c    define the resized array of driving forces and slip resistance
         PDFRed   = 0.d0
         PDFCritR = 0.d0
         do i = 1,nGammaA
	    PDFRed(i)   = PDF(idGammaA(i))
	    PDFCritR(i) = PDFCrit(idGammaA(i))
         enddo
c     
c     compute the residual coming from the plastic kinetic relation
         Res1 = 0.d0
         do i = 1,nGammaA
	    if (DeltaGammaA(i) .lt. 0.d0) DeltaGammaA(i) = 0.d0
	    Temp    = (DeltaGammaA(i)/(GammaDot0*DeltaT) + 1.d0)**mExp
	    Res1(i) = PDFRed(i) - PDFCritR(i)*Temp
         enddo
c     
c     compute the residual related to the slip resistance parameter
         Res2 = 0.d0
         do i = 1,24
	    Res2(i) = SlipResA(i) - SlipRes(i) - DeltaSlipRes(i)
         enddo
c
c    construct the array of residuals
         AResP = 0.d0
         do i = 1,nGammaA
            AResP(i) = Res1(i)
         enddo
         do i = 1,24
	    AResP(i + nGammaA) = Res2(i)
         enddo
c     
c     contruct the matrix of Jacobian
         do i = 1,nSize
            MVar(i,iP) = (AResP(i) - ARes(i))/Pertb
         enddo
      enddo
c
*  remark:
*  the Jacobian matrix obtained from numerical perturbation is dRes/dVar
*  while the required one is dVar/dRes. inversion is done using standard
*  linpack subroutine
*
c  inverting the Jacobian matrix
      dVardRes = MVar
      call INVPGB(dVardRes,nSize,rpdum,radum,toler,ierr)
      if (ierr .gt. 0) then
         write(*,*)'UMATBCC: SINGULAR JACOBIAN MATRIX'
         stop
      endif
c     
c     calculate the update for the array of variables
      do i = 1,nSize
         do j = 1,nSize
	    AVar(i) = AVar(i) - dVardRes(i,j)*ARes(j)
         enddo
      enddo
c     
c     end of Newton-Raphson iteration                                   (!)
      enddo
c     
 67   continue
      deallocate(AVar,PDFRed,PDFCritR,Res1,MVar,ARes,AVarP,
     &     AResP,radum,rpdum,ipiv,dVardRes,stat = errorstate)
      if (errorstate .gt. 0) then
         write(*,*)'UMATBCC: FAILED TO DEALLOCATE'
         stop
      endif

c     
      return
      end      

      subroutine Comp_DFthDt(XiA, Flag, F, Jtr)
c     
c     .......1.........2.........3.........4.........5.........6.........7..
c     
c     **********************************************************************
c     
c     XiA:= 1-sum (Xi(i))
c     
c     **********************************************************************
c     
c     outputs
      double precision F(3,3)
c     inputs
      double precision XiA, Jtr
      double precision, parameter::DeltaVT = 3.2627d-2
      character Flag           
c     internals
      double precision a, eye(3,3),deltaT
      integer i
      include 'constants.for'        
c     
      eye=0.d0
      F = 0.d0
      do i=1,3
         eye(i,i)=1.d0
      end do
      select case(Flag)
      case('a','A')
         a=(alpha_A*XiA+(1.d0+DeltaVT)*alpha_M*(1.d0-XiA))/Jtr
         F=a*eye
      case('f','F')
         a=alpha_F
         F=a*eye
      end select
c
      end subroutine           

      subroutine Comp_Fth(theta, XiA, Flag, F, Jtr)
c     
c     .......1.........2.........3.........4.........5.........6.........7..
c     
c     **********************************************************************
c     
c     Xi:= 1-sum (Xi(i))
c     Flag : 'a' for Austenite and 'f' for Ferrite
c     
c     **********************************************************************
c     outputs
      double precision F(3,3)
c     inputs
      double precision theta, XiA, Jtr
      double precision, parameter::DeltaVT = 3.2627d-2       
      character Flag           
c     internals
      double precision a, eye(3,3),deltaT
      integer i
      include 'constants.for'        
      
      deltaT = theta-Theta0
      F = 0.d0
      eye=0.d0
      do i=1,3
         eye(i,i)=1.d0
      end do
      select case(Flag)
      case('a','A')
         a=(alpha_A*XiA+(1.d0+DeltaVT)*alpha_M*(1.d0-XiA))/Jtr
         a = a*deltaT
         F=(a+1.d0)*eye
      case('f','F')
         a=(alpha_F*deltaT+1.d0)
         F=a*eye
      end select
c
      end subroutine           
      
      subroutine Comp_Etam(ElasStiff,Ee,theta,Xi,Flag,Jtr,eta)
c     
c     .......1.........2.........3.........4.........5.........6.........7..
c     
c     **********************************************************************
c     
c     Xi:= 1-sum (Xi(i))
c     
c     **********************************************************************
c     
c     Outputs
      double precision eta        
c     Inputs
      double precision ElasStiff(3,3,3,3),Ee(3,3)
      double precision theta, Xi, Jtr
      character Flag
c     Internals
      double precision CEe(3,3),eye(3,3),Ce(3,3),CeCEe(3,3),CEe_dot_Ee
      double precision Fth(3,3),DFthDt(3,3),FthInv(3,3),Jth
      double precision DFDtFinv(3,3)
      double precision temp1, temp2
      double precision rho
      include 'constants.for'       
c
c     Body
      rho=Rho_0
      eye=0.d0
      do i=1,3
         eye(i,i)=1.d0
      end do 
      
      Ce=2.d0*Ee+eye
      call TENSOR4TENSOR2(ElasStiff,Ee,CEe)       
      CeCEe=matmul(Ce,CEe)
      call DOTMATRIX3X3(CEe,Ee,CEe_dot_Ee) ! C*Ee.Ee
c
      call Comp_DFthDt(Xi,Flag,DFthDt,Jtr)
      call Comp_Fth(theta, Xi, Flag, Fth,Jtr)
c
      call INVERSEMATRIX3X3(Fth,FthInv) ! Fth^-1
      call DETERMINANTMATRIX3X3(Fth,Jth) ! Jth
      DFDtFinv=matmul(DFthDt,FthInv)       
      call DOTMATRIX3X3(CeCEe,DFDtFinv,temp1)
      temp1=temp1*Jtr*Jth/rho
      call DOTMATRIX3X3(FthInv,DFthDt, temp2)
      temp2=temp2*CEe_dot_Ee*Jtr*Jth/(2.d0*rho)
      eta=temp1
      end subroutine
c
c .......1.........2.........3.........4.........5.........6.........7..
c
c **********************************************************************
c
c  subroutine to calculate the product of fourth-order tensor with 
c  second-order tensor (e.g. elastic stiffness with elastic strain):
c  M2B = M4 . M2A
c
c **********************************************************************
c
      subroutine TENSOR4TENSOR2(M4,M2A,M2B)
c     
c     declaration of variables
      double precision M4(3,3,3,3),M2A(3,3),M2B(3,3)
      integer i,j,k,l
c     
c     compute the product of fourth-order tensor with second-order tensor
      do i = 1,3
         do j = 1,3
            M2B(i,j) = 0.d0
            do k = 1,3
               do l = 1,3
                  M2B(i,j) = M2B(i,j) + M4(i,j,k,l)*M2A(k,l)
               enddo
            enddo
         enddo
      enddo
c     
      return
      end

c     
c     **********************************************************************
c     
c     subroutine to compute inverse of matrix: InvM = M^(-1)
c     
c     **********************************************************************
c     
      subroutine INVERSEMATRIX3X3(M,InvM)
c     
c     declaration of variables
      double precision DetM,M(3,3),InvM(3,3)
      integer i,j
c     
c     compute the determinant of the matrix
      call DETERMINANTMATRIX3X3(M,DetM)
c     
c     construct the inverse of the matrix
      InvM = 0.d0
      InvM(1,1) = M(2,3)*M(3,2) - M(2,2)*M(3,3)
      InvM(1,2) = M(1,2)*M(3,3) - M(1,3)*M(3,2)
      InvM(1,3) = M(1,3)*M(2,2) - M(1,2)*M(2,3)
      InvM(2,1) = M(2,1)*M(3,3) - M(2,3)*M(3,1)
      InvM(2,2) = M(1,3)*M(3,1) - M(1,1)*M(3,3)
      InvM(2,3) = M(1,1)*M(2,3) - M(1,3)*M(2,1)
      InvM(3,1) = M(2,2)*M(3,1) - M(2,1)*M(3,2)
      InvM(3,2) = M(1,1)*M(3,2) - M(1,2)*M(3,1)
      InvM(3,3) = M(1,2)*M(2,1) - M(1,1)*M(2,2)
      do i = 1,3
         do j = 1,3
            InvM(i,j) = (-1.d0/DetM)*InvM(i,j)
         enddo
      enddo
c     
      return
      end

c
c **********************************************************************
c
c  subroutine to compute the determinant of matrix: detM = det(M)
c
c **********************************************************************
c
      subroutine DETERMINANTMATRIX3X3(M,DetM)
c     
c     declaration of variables
      double precision DetM,M(3,3)
c     
c     compute the determinant of the matrix
      DetM = M(1,1)*M(2,2)*M(3,3) + M(1,2)*M(2,3)*M(3,1) +
     &       M(1,3)*M(2,1)*M(3,2) - M(1,3)*M(2,2)*M(3,1) -
     &       M(1,2)*M(2,1)*M(3,3) - M(1,1)*M(2,3)*M(3,2)
c     
      return
      end

c
c **********************************************************************
c
c  subroutine to compute dot product of matrices: c = Ma dot Mb
c
c **********************************************************************
c
      subroutine DOTMATRIX3X3(Ma,Mb,c)
c     
c     declaration of variables
      double precision Ma(3,3),Mb(3,3),c
      integer i,j
c     
c     compute the dot product of matrices
      c = 0.d0
      do i = 1,3
         do j = 1,3
            c = c + Ma(i,j)*Mb(i,j)
         enddo
      enddo
c     
      return
      end       

c
c
c **********************************************************************
c
c  subroutine to compute summation of matrices: Mc = a.Ma + b.Mb
c
c **********************************************************************
c
      subroutine ADDMATRIX3X3(a,Ma,b,Mb,Mc)
c     
c     declaration of variables
      double precision a,b,Ma(3,3),Mb(3,3),Mc(3,3)
      integer i,j
c     
c     compute the summation of matrices
      do i = 1,3
         do j = 1,3
            Mc(i,j) = a*Ma(i,j) + b*Mb(i,j)
         enddo
      enddo
c     
      return
      end

c
c **********************************************************************
c
c  subroutine to compute transpose of matrix: TransM = M^T
c
c **********************************************************************
c
      subroutine TRANSPOSEMATRIX3X3(M,TransM)
c     
c     declaration of variables
      double precision M(3,3),TransM(3,3)
      integer i,j
c     
c     compute the transpose of the matrix
      do i = 1,3
         do j = 1,3
            TransM(i,j) = M(j,i)
         enddo
      enddo
c     
      return
      end

c
c **********************************************************************
c
c  subroutine to compute multiplication of matrices: Mc = Ma.Mb
c
c **********************************************************************
c
      subroutine MULTIMATRIX3X3(Ma,Mb,Mc)
c     
c     declaration of variables
      double precision Ma(3,3),Mb(3,3),Mc(3,3)
      integer i,j,k
c     
c     compute the multiplication of matrices
      do i = 1,3
         do j = 1,3
            Mc(i,j) = 0.d0
            do k = 1,3
               Mc(i,j) = Mc(i,j) + Ma(i,k)*Mb(k,j)
            enddo
         enddo
      enddo
c     
      return
      end

c
c .......1.........2.........3.........4.........5.........6.........7..
c
c **********************************************************************
c
c  subroutine to compute the plastic deformation gradient
c
c **********************************************************************
c
      subroutine PLASTICDEFGRADBCC(DeltaGammaA,idGammaA,SlipSys,
     &     nGammaA, FP0,FP1)
c     
c     declaration of variables
      double precision DeltaGammaA(nGammaA),FP0(3,3),FP1(3,3),
     &     I2(3,3),LP(3,3),InvLP(3,3),DGammaFull(24),SlipSys(24,3,3),
     &     I4(3,3,3,3)
      integer i,j,k,nGammaA,idGammaA(nGammaA)
c     calculate the inverse plastic velocity gradient
      I2 = 0.d0
      do i = 1,3
         I2(i,i) = 1.d0
      enddo
      DGammaFull = 0.d0
      do k = 1,nGammaA
         DGammaFull(idGammaA(k)) = DeltaGammaA(k)
      enddo
      LP = 0.d0
      do k = 1,24
         do i = 1,3
            do j = 1,3
               LP(i,j) = LP(i,j) + DGammaFull(k)*SlipSys(k,i,j)
            enddo
         enddo
      enddo
      LP = I2 - LP
      call INVERSEMATRIX3X3(LP,InvLP)
c     
c     calculate the plastic deformation gradient
      call MULTIMATRIX3X3(InvLP,FP0,FP1)
c     
      return
      end

c
c .......1.........2.........3.........4.........5.........6.........7..
c
c **********************************************************************
c
c  subroutine to compute non-Schmid stress according to Bassani's
c  convention
c
c **********************************************************************
c
      subroutine NONSCHMIDSTRESSBCC(Schmid,NonSchmid)
c     
c     declaration of variables
      double precision Schmid(24),NonSchmid(24),NGlide
c     
c     define the effective non-glide parameter (calibrated for BCC ferrite,
c     using the data of Franciosi (1983) for alpha-iron)
c     values: NGlide = 0.12
      NGlide = 1.2d-1
c     
c     compute the non-glide stress (non-Schmid stress)
      NonSchmid(1)  = NGlide*Schmid(6) 
      NonSchmid(2)  = NGlide*Schmid(3) 
      NonSchmid(3)  = NGlide*Schmid(2) 
      NonSchmid(4)  = NGlide*Schmid(5) 
      NonSchmid(5)  = NGlide*Schmid(4) 
      NonSchmid(6)  = NGlide*Schmid(1)
c     
      NonSchmid(7)  = NGlide*Schmid(10) 
      NonSchmid(8)  = NGlide*Schmid(11) 
      NonSchmid(9)  = NGlide*Schmid(12) 
      NonSchmid(10) = NGlide*Schmid(7) 
      NonSchmid(11) = NGlide*Schmid(8) 
      NonSchmid(12) = NGlide*Schmid(9)
c     
      NonSchmid(13) = NGlide*Schmid(16) 
      NonSchmid(14) = NGlide*Schmid(17) 
      NonSchmid(15) = NGlide*Schmid(18) 
      NonSchmid(16) = NGlide*Schmid(13) 
      NonSchmid(17) = NGlide*Schmid(14) 
      NonSchmid(18) = NGlide*Schmid(15)
c     
      NonSchmid(19) = NGlide*Schmid(24) 
      NonSchmid(20) = NGlide*Schmid(21) 
      NonSchmid(21) = NGlide*Schmid(20) 
      NonSchmid(22) = NGlide*Schmid(23) 
      NonSchmid(23) = NGlide*Schmid(22) 
      NonSchmid(24) = NGlide*Schmid(19)
c     
      return
      end

c
c .......1.........2.........3.........4.........5.........6.........7..
c
c **********************************************************************
c
c  subroutine to compute hardening modulus and shear hardening modulus
c  using the exponential hardening model (Brown et al. 1992)
c
c **********************************************************************
c
      subroutine HARDENINGMODULUSBCC(SlipResA,HMod,GMod,
     &     ShearMod)
      include 'constants.for'
c     declaration of variables
      double precision SlipResA(24),ShearMod,HMod(24,24),GMod(24),
     &     Ks(24)
      integer i,j
      logical PlotDebug,PlotResult
c     
c     hardening model parameters (calibrated for BCC ferrite, taken into 
c     account the mixture with 25% bainite, data from Jacques et al (2001))
c     values: sInf = 0.412 for original matrix
c     0.178 for soft matrix1 (low max stress)
c     0.442 for soft matrix2 (med max stress)
c     0.404 for hard matrix1 (med max stress)
c     0.758 for hard matrix2 (high max stress)
c     k0   = 1.9
c     w0   = 2.8
c     q0   = 1
c     c0   = 0.5
      double precision::sInf = sInf_F ! 4.12d-1  ! consistent with (*)
      double precision::k0   = k0_F ! 1.9d0
      double precision::w0   = w0_F ! 2.8d0
      double precision::q0   = q0_F ! 1.d0
      double precision::c0   = c0_F ! 5.d-1
c     compute the single slip hardening model
      do i = 1,24
         Ks(i) = k0*(1.d0 - SlipResA(i)/sInf)**w0
      enddo
c     
c     calculate the hardening modulus
      HMod = 0.d0
      do i = 1,24
         do j = 1,24
            HMod(i,j) = q0*Ks(j)
         enddo
         HMod(i,i) = Ks(i)
      enddo
c     
c     calculate the shear hardening modulus
      do i = 1,24
         GMod(i) = 0.d0
         do j = 1,24
            GMod(i) = GMod(i) + HMod(j,i)/(24 * c0 * ShearMod)
         enddo
      enddo
c     
      return
      end

c
c .......1.........2.........3.........4.........5.........6.........7..
c
C=======================================================================
C
      SUBROUTINE INVPGB( A, NA, RP, RA, EPS, IERR )
C
C-----------------------------------------------------------------------
C...  PURPOSE
C...     MATRIX INVERSION WITH PARTIAL PIVOTING( GAUSS JORDAN).
C...     THE COMPLETE INVERSE OF A IS PROJECTED ON THE ORIGINAL MATRIX
C...        [ A ].
C...        < RP > AND < RA > ARE USED FOR ROW INTERCHANGING.
C...
C...  INPUT
C...     ARG.LIST  A      -R(NA*NA) MATRIX TO BE INVERTED.
C...               NA     -I(1)     ORDER OF [ A ]
C...               EPS    -R(1)     TOLERANCE FOR PIVOT
C...
C...  OUTPUT
C...     ARG.LIST  A      -R(NA*NA) INVERSE OF [ A ]
C...               IERR   -I(1)     ERROR CODE
C...                                = 0 INVERSION SUCCESFULL;
C...                                =-1 SINGULARITY ENCOUNTERED;
C...
C...  WORKSPACE
C...     ARG.LIST  RP     -I(NA)    ROW INDICES OF [ A ]
C...               RA     -R(NA)    WORKSPACE FOR ROW OF [ A ]
C...
C-----------------------------------------------------------------------
C
C...ARGUMENTS
      INTEGER          NA, RP(*), IERR
      DOUBLE PRECISION A(NA,*), RA(*), EPS
C     
C...  LOCAL
      INTEGER          I, IR, J, K, N
      DOUBLE PRECISION APIV, AR, EPIV, HPIV
C     
      IERR = 0
      N = NA
      EPIV = EPS
C     
C...  INITIALIZE ROW-INDICES
      DO 10, I = 1, N
         RP(I) = I
 10   CONTINUE
C     
C...  ROW LOOP
      DO 200, I = 1, N
C     
C...  FIND PIVOT
         HPIV = ABS( A(I,I) )
         IR = I
         IF( I .LT. N ) THEN
            DO 50, K = (I+1), N
               APIV = ABS( A(K,I) )
               IF( APIV .GT. HPIV ) THEN
                  HPIV = APIV
                  IR = K
               ENDIF
 50         CONTINUE
         ENDIF
         IF( HPIV .LT. EPIV ) THEN
            IERR = -RP(IR)
            RETURN
         ENDIF
C     
         IF( IR .GT. I ) THEN
C...  INTERCHANGE ROWS
            DO 60, J = 1, N
               AR      = A(I, J)
               A(I, J) = A(IR,J)
               A(IR,J) = AR
 60         CONTINUE
            K      = RP(I)
            RP(I)  = RP(IR)
            RP(IR) = K
         ENDIF
C     
         APIV = 1.D0/A(I,I)
         CALL UVS( A(1,I), N, APIV, A(1,I) )
         A(I,I) = APIV
C     
C...  COLUMN LOOP
         DO 100, J = 1, N
            IF( J .NE. I ) THEN
               APIV = -A(I,J)
               IF( APIV .NE. 0.D0 ) THEN
                  CALL UVPWS( A(1,J), A(1,I), N, APIV, A(1,J) )
                  A(I,J) = APIV * A(I,I)
               ELSE
                  A(I,J) = 0.D0
               ENDIF
            ENDIF
 100     CONTINUE
C     
 200  CONTINUE
C     
C...  INTERCHANGE COLUMNS
      DO 250, I = 1, N
C     
         DO 210, J = 1, N
            K  = RP(J)
            RA(K) = A(I,J)
 210     CONTINUE
C     
         DO 220, J = 1, N
            A(I,J) = RA(J)
 220     CONTINUE
C     
 250  CONTINUE
      RETURN
      END

C
C=======================================================================
C
      SUBROUTINE UVS( V, N, S, U )
C
C-----------------------------------------------------------------------
C...
C...  PURPOSE
C...     PRODUCT OF VECTOR AND SCALAR
C...
C...             N           N
C...        < U >   =  < V >   *  S
C...             1           1
C...
C...  INPUT
C...     ARG.LIST  V(1N), N, S
C...
C...  OUTPUT
C...     ARG.LIST  U(1N)
C...
C-----------------------------------------------------------------------
C
      INTEGER          N
      DOUBLE PRECISION S, U(N), V(N)
      INTEGER          I
C     
      DO 10 I = 1, N
         U(I) = V(I) * S
 10   CONTINUE
      END
C
C=======================================================================
C
      SUBROUTINE UVPWS( V, W, N, S, U )
C
C-----------------------------------------------------------------------
C...
C...  PURPOSE
C...          N           N          N
C...     < U >   =  < V >   +  < W >  * S
C...          1           1          1
C...
C...  INPUT
C...     ARG.LIST  V(1N), W(1N), N, S
C...
C...  OUTPUT
C...     ARG.LIST  U(1N)
C...
C-----------------------------------------------------------------------
C
      INTEGER          N
      DOUBLE PRECISION S, U(N), V(N), W(N)
      INTEGER          I
C
      IF ( N .LE. 0 ) RETURN
      DO 10 I = 1, N
         U(I)  =  V(I)  +  W(I) * S
   10 CONTINUE
      RETURN
      END

c
c .......1.........2.........3.........4.........5.........6.........7..
c
c **********************************************************************
c
c  subroutine to compute the consistent tangent by means of numerical
c  perturbation technique
c
c **********************************************************************
c
      subroutine CONSISTENTTANGENTBCC(statev,dfgrd0,dfgrd1,DeltaT0,
     &     DDPDDF,FirstPK,nstatv,ConTangent,theta,
     &     ElasStif,SlipSys,ShearMod)

	double precision statev(nstatv),dfgrd0(3,3),dfgrd1(3,3),DeltaT,
     & DDPDDF(3,3,3,3),FirstPK(3,3),F0(3,3),F1(3,3),FP0(3,3),DeltaT0,
     & FE0(3,3),FE1(3,3),Gamma0(24),Gamma1(24),SigmaXi0,SlipRes0(24),
     & SlipRes1(24),ElasStif(3,3,3,3),SlipSys(24,3,3),Beta0,ShearMod,
     & SlipResA(24),PDF(24),PDFCrit(24),SubStep,SubFraction,DeltaF(3,3),
     & I2(3,3),JE1,FP1(3,3),I4(3,3,3,3),EE1(3,3),CE1(3,3),TransFE1(3,3),
     & PK2(3,3),PK1(3,3),MTemp(3,3),TransF1(3,3),BetaA,DF(3,3),
     & CopySDV(83),InvFP1(3,3),InTransFP1(3,3)
      integer i,j,k,TagGamma(24),MaterialID,nGammaA,errorstate,iterP,
     & iPer,jPer,nstatv
      logical PlotDebug,PlotResult,UpgradeSlipSys,SubStepping,ConTangent
      double precision theta, det_th,Fth(3,3)
      double precision, allocatable::DeltaGammaA(:)
      integer, allocatable::idGammaA(:)
      double precision, parameter::FPert = 1.d-8
c
      I2 = 0.d0
      do i = 1,3
         I2(i,i) = 1.d0
      enddo
c     begin the numerical perturbation
      do iPer = 1,3
         do jPer = 1,3
c
c  initialization of stress-update algorithm
            CopySDV     = statev
            SubFraction = 0.d0
            SubStep     = 1.d0
            DeltaF      = dfgrd1 - dfgrd0
            DeltaF(iPer,jPer) = DeltaF(iPer,jPer) + FPert
c
c  substepping procedure begins
            SubStepping = .false.
            DeltaT  = SubStep*DeltaT0
            DF      = SubStep*DeltaF
            F0 = dfgrd0
            F1 = F0 + DF
            do i = 1,3
               do j = 1,3
                  FP0(i,j)  = CopySDV(j + (i - 1)*3)
               enddo
            enddo
            do i = 1,24
               Gamma0(i)   = CopySDV(i + 9)
               SlipRes0(i) = CopySDV(i + 33)
               TagGamma(i) = int(CopySDV(i + 58))
            enddo
            Beta0    = CopySDV(58)
c     
c     set the initially active slip systems
            iterP = 0
            nGammaA = 0
            do i = 1,24
               if (TagGamma(i) .eq. 1) then
                  nGammaA = nGammaA + 1
               endif
            enddo
            allocate(DeltaGammaA(nGammaA),idGammaA(nGammaA),
     &           stat = errorstate)
            if (errorstate .gt. 0) then
               write(*,*)'UMATBCC: FAILED TO ALLOCATE 42'
               stop
            endif
            UpgradeSlipSys  = .false.
c
c     set the initial values for iteration
            j = 0
            k = 0
            do i = 1,24
               if (TagGamma(i) .eq. 1) then
                  j = j + 1
                  DeltaGammaA(j) = 0.d0
                  idGammaA(j)    = i
               endif
            enddo
            BetaA    = Beta0
            SlipResA = SlipRes0
c
c     compute the stress-update using return-mapping algorithm
            call RETURNMAPPINGBCC(DeltaGammaA,SlipResA,nGammaA,
     &           PDF,PDFCrit,
     &           idGammaA,FP0,F1,BetaA,DeltaT,SubStepping,
     &           ConTangent,theta,
     &           ElasStif,SlipSys,ShearMod)
c     
c     update the parameters 
            SlipRes1 = SlipResA
            Gamma1   = Gamma0
            do i = 1,nGammaA
               Gamma1(idGammaA(i)) = Gamma0(idGammaA(i)) +
     &              DeltaGammaA(i)
            enddo
            Beta1    = BetaA
            SlipRes1 = SlipResA
c     
c     calculate the plastic deformation gradients
            call PLASTICDEFGRADBCC(DeltaGammaA,idGammaA,SlipSys,
     &           nGammaA, FP0, FP1)
            call comp_Fth(theta,0.d0,'F',Fth,1.d0) ! Fth
            call DETERMINANTMATRIX3X3(Fth,det_th) ! det(Fth)
            FP1=matmul(Fth, FP1) !  Fp <- Fth*Fp      
c     calculate the second Piola-Kirchhoff stress
            call INVERSEMATRIX3X3(FP1,InvFP1)
            call MULTIMATRIX3X3(F1,InvFP1,FE1)
            call TRANSPOSEMATRIX3X3(FE1,TransFE1)
            call MULTIMATRIX3X3(TransFE1,FE1,CE1)
            call ADDMATRIX3X3(5.d-1,CE1,-5.d-1,I2,EE1)
            call TENSOR4TENSOR2(ElasStif,EE1,PK2)
c  calculate the first Piola-Kirchhoff stress
            call TRANSPOSEMATRIX3X3(InvFP1,InTransFP1)
            call MULTIMATRIX3X3(PK2,InTransFP1,MTemp)
            call MULTIMATRIX3X3(FE1,MTemp,PK1)
            PK1=det_th*PK1
c     construct the consistent tangent
            do i = 1,3
               do j = 1,3
                  DDPDDF(i,j,iPer,jPer) = (PK1(i,j) -
     &             FirstPK(i,j)) / FPert
               enddo
            enddo
c     
            deallocate(DeltaGammaA,idGammaA,stat = errorstate)
            if (errorstate .gt. 0) then
               write(*,*)'UMAT: FAILED TO DEALLOCATE '
               stop
            endif
c     end of perturbation algorithm
         enddo
      enddo
c     
c     
      return
      end

c .......1.........2.........3.........4.........5.........6.........7..
c
c **********************************************************************
c
c  subroutine to compute internal energy as a function of elastic deformation, 
c  temperature and internal variables
c **********************************************************************
c    
      subroutine comp_InEnBCC(Fe,T,C,Beta,Flag,energy)
c     inputs
      double precision Fe(3,3)
      double precision T ! Temperature
      double precision Beta,C(3,3,3,3)
      character Flag
c     outputs
      double precision energy
c     internals
      double precision Fthptr(3,3),F_inv(3,3),Fth(3,3)
      double precision SigXi,SigXiA
      double precision Fe_t(3,3),Ce(3,3),Ee(3,3)
      double precision eye(3,3)
      double precision PK2(3,3)
      double precision Rho,Jtr,Jth,mu
      double precision En_m,En_mth,En_d,En_s,En_th,Eta_m,SHeat
      integer i
      include 'constants.for'      
c     Body
      Rho=Rho_0
      eye=0.d0
      do i=1,3
      eye(i,i)=1.d0
      end do
c     Fe computation     
      call Comp_Fth(T,1.d0,Flag,Fth,1.d0)
      call DETERMINANTMATRIX3X3(Fth,Jth)
      Fe_t=TRANSPOSE(Fe)
      Ce=matmul(Fe_t,Fe)
      Ee=0.5d0*(Ce-eye)
c mechanical energy      
      call TENSOR4TENSOR2(C,Ee,PK2)
      call DOTMATRIX3X3(PK2,Ee,En_m)
      En_m=(Jth)/(2.d0*Rho)*En_m
c coupled (mechanical/thermal) energy      
      call Comp_Etam(C,Ee,T,1.d0,'F',1.d0,Eta_m)
      En_mth=Eta_m*T
c defect energy
      mu=Mu_F 
      En_d=(1.d0/(2.d0*Rho))*W_F*mu*Beta**2.d0
c thermal energy
      SHeat=SHeatF
      En_th=SHeat*T
c total energy
      energy=En_m+En_mth+En_d+En_th      
      end subroutine

c
c .......1.........2.........3.........4.........5.........6.........7..
c
c ______________________________________________________________________
c
c
c
c                       *** BCC FERRITE BLOCK ***
c
c
c ______________________________________________________________________
c
c .......1.........2.........3.........4.........5.........6.........7..
c
c **********************************************************************
c
c  user-material subroutine for elasto-plastic model of BCC ferrite
c  in tranformation-induced plasticity steels
c
c **********************************************************************
c
      subroutine UMATBCC(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     & drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,dpred,
     & cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     & celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc,angles)
c
c  input from ABAQUS:
c    stress => old-Cauchy stress
c    statev => (user-defined) state variables
c    stran  => old-total strain (small-strain)
c    dstran => total strain increment
c    time   => old-value of step/total time
c    dtime  => time increment
c    predef => old-interpolated predefined field variables
c    dpred  => increment of predefined field varibles
c    cmname => given material name
c    ndi    => number of direct stress components
c    nshr   => number of shear stress components
c    ntens  => size of stress array (= ndi + nshr)
c    nstatv => size of statev array
c    props  => material constants
c    nprops => size of props array
c    coords => coordinates of the material point
c    drot   => rotation increment matrix
c    dfgrd0 => old-deformation gradient
c    dfgrd1 => new-deformation gradient
c
c  return to ABAQUS (updated values):
c    stress => updated-Cauchy stress
c    ddsdde => elastic/elasto-plastic stiffness tangent
c    statev => updated-state variables
c
c  header for standard user-subroutine for ABAQUS
c      include 'aba_param.inc'
c
      character*80 cmname
      double precision stress(ntens),statev(nstatv),ddsdde(6,6),
     & ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),time(2),
     & predef(1),dpred(1),props(nprops),coords(3),drot(3,3),dfgrd0(3,3),
     & dfgrd1(3,3),sse,spd,scd,rpl,drpldt,temp,dtemp,pnewdt,celent,
     & dtime, angles(3)
c
c  declaration of variables used in the subroutine
      double precision F0(3,3),DF(3,3),F1(3,3),DeltaT,FP0(3,3),FP1(3,3),
     & FE0(3,3),FE1(3,3),Gamma0(24),Gamma1(24),SlipRes0(24),
     & SlipRes1(24),ElasStif(3,3,3,3),SlipSys(24,3,3),Beta0,ShearMod,
     & SlipResA(24),PDF(24),PDFCrit(24),SubStep,SubFraction,DeltaF(3,3),
     & I2(3,3),JE1,I4(3,3,3,3),EE1(3,3),CE1(3,3),InvFP1(3,3),
     & TransFE1(3,3),PK2(3,3),PK1(3,3),InTransFP1(3,3),MTemp(3,3),
     & CauchyT(3,3),KirchS(3,3),TransF1(3,3),XiA(24),DDPDDF(3,3,3,3),
     & BTgt(3,3,3,3),HTgt(3,3,3,3),J1,StateR(83),DefGrad0(3,3),
     & DefGrad1(3,3),BetaA,WorkF0,WorkF1,dWorkF
      integer i,j,k,TagGamma(24),MaterialID,nGammaA,
     & errorstate,iterP,l,m,n
      logical PlotDebug,PlotResult,UpgradeSlipSys,SubStepping,ConTangent
      double precision, allocatable::DeltaGammaA(:)
      integer, allocatable::idGammaA(:)
      double precision InEn0,InEn1,InEnPrt,dudt
      double precision Fth(3,3),Jth,tmpprt,Fthp(3,3),Ftemp(3,3)
      double precision FinvTemp(3,3)      
      double precision temp11,PFDot,temptemp,Beta1
c
      include 'constants.for'
      I2 = 0.d0
      do i = 1,3
         I2(i,i) = 1.d0
      enddo      
c  get the material identification from ABAQUS
      MaterialID = int(props(1)) !int((dble(noel) - 1.d-2)/8.d0) + 1
c
c  get the input of material properties
      call MATERIALPROPERTIESBCC(MaterialID,
     &     ElasStif,SlipSys,ShearMod,angles)
c	call IDENTITYTENSOR(I2,I4)
c
c  initialization of stress-update algorithm
      SubFraction = 0.d0
      SubStep     = 1.d0
      DeltaF      = dfgrd1 - dfgrd0
      StateR      = statev
      DefGrad0    = dfgrd0
      DefGrad1    = dfgrd1
c cc internal energy at the beginning of the step
      call comp_Fth(temp,1.d0,'F',Fth,1.d0)
      do i = 1,3
         do j = 1,3
            FP0(i,j)  = statev(j + (i - 1)*3)
         enddo
      enddo      
      FTemp = matmul(FP0,Fth)
      call INVERSEMATRIX3X3(FTemp,FinvTemp)
      FE0 = matmul(dfgrd0,FinvTemp)
      call comp_InEnBCC(Fe0,Temp,ElasStif,statev(58),'F',InEn0)
c
c  substepping procedure begins
 60   continue
      SubStepping = .false.
      DeltaT  = SubStep*dtime
      DF      = SubStep*DeltaF
      F0 = DefGrad0
      F1 = F0 + DF
      do i = 1,3
         do j = 1,3
	    FP0(i,j)  = StateR(j + (i - 1)*3)
         enddo
      enddo
      do i = 1,24
         Gamma0(i)   = StateR(i + 9)
         SlipRes0(i) = StateR(i + 33)
         TagGamma(i) = int(StateR(i + 58))
      enddo
      Beta0    = StateR(58)
      WorkF0   = StateR(83)
c     
c
c  set the initially active slip systems
      iterP = 0
 61   continue
      nGammaA = 0
c counting the number of active systems
      do i = 1,24
         if (TagGamma(i) .eq. 1) then
	    nGammaA = nGammaA + 1
         endif
      enddo
c allocating according to nGamma
      allocate(DeltaGammaA(nGammaA),idGammaA(nGammaA),
     &     stat = errorstate)
      if (errorstate .gt. 0) then
         write(*,*)'UMATBCC: FAILED TO ALLOCATE 2'
         stop
      endif
      UpgradeSlipSys  = .false.
c
c  set the initial values for iteration
      j = 0
      k = 0
      do i = 1,24
         if (TagGamma(i) .eq. 1) then
	    j = j + 1
            DeltaGammaA(j) = 0.d0
	    idGammaA(j)    = i
         endif
      enddo
      BetaA    = Beta0
      SlipResA = SlipRes0
c
c  compute the stress-update using return-mapping algorithm
      ConTangent = .false.
      call RETURNMAPPINGBCC(DeltaGammaA,SlipResA,nGammaA,PDF,PDFCrit,
     &     idGammaA,FP0,F1,BetaA,DeltaT,SubStepping,
     &     ConTangent,temp+dtemp,
     &     ElasStif,SlipSys,ShearMod)
c     
      if (SubStepping) goto 68
      SubStepping = .false.
c
c  consistency check for slip systems
      call SLIPSYSTEMCHECKBCC(PDF,PDFCrit,DeltaGammaA,idGammaA,nGammaA,
     &     TagGamma,UpgradeSlipSys,iterP)
      if (UpgradeSlipSys) then
         iterP = iterP + 1
         deallocate(DeltaGammaA,idGammaA,stat = errorstate)
         if (errorstate .gt. 0) then
            write(*,*)'UMATBCC: FAILED TO DEALLOCATE 3'
            stop
         endif
         goto 61
      endif
      if (iterP .gt. 6) then
         SubStepping = .true.
         goto 68
      endif
c
c  update the parameters 
      SlipRes1 = SlipResA
      Gamma1   = Gamma0
      do i = 1,nGammaA
         Gamma1(idGammaA(i)) = Gamma0(idGammaA(i)) + DeltaGammaA(i)
      enddo
      Beta1    = BetaA
      SlipRes1 = SlipResA
c
c  calculate the plastic deformation gradients
      call PLASTICDEFGRADBCC(DeltaGammaA,idGammaA,SlipSys,nGammaA,FP0,
     &     FP1)
c
c  performing substepping algorithm
      if (SubStep .lt. 1.d-4) SubStepping = .false.
 68   continue
      if (SubStepping) then
         SubStep = SubStep/2.d0
         deallocate(DeltaGammaA,idGammaA,stat = errorstate)
         if (errorstate .gt. 0) then
            write(*,*)'UMATBCC: FAILED TO DEALLOCATE 2'
	    stop
         endif
         goto 60
      else
         SubFraction = SubFraction + DeltaT/dtime
         SubStep = 1.d0 - SubFraction
         if (SubStep .gt. 1.d-6) then
	    DefGrad0 = F1
	    do i = 1,3
               do j = 1,3
                  StateR(j + (i - 1)*3)     = FP1(i,j)
               enddo
	    enddo
	    do i = 1,24
               StateR(i + 9)   = Gamma1(i)
               StateR(i + 33)  = SlipRes1(i)
               StateR(i + 58)  = dble(TagGamma(i))
	    enddo
	    StateR(58) = Beta1
            deallocate(DeltaGammaA,idGammaA,stat = errorstate)
            if (errorstate .gt. 0) then
               write(*,*)'UMATBCC: FAILED TO DEALLOCATE 1'
               stop
            endif
	    goto 60
         endif
      endif
c
*  remark:
*  for simplicity reason, the main parameter after converging substep is
*  updated. this parameter will be used as the initial value for the
*  next substep
*
*  besides, the calculation of numerical consistent tangent will start
*  from the previously convergent state, and not from the intial value
*  in order to avoid the use of substepping procedure (again) during the
*  computation of numerical consistent tangent with perturbed value
*
c  calculate the second Piola-Kirchhoff stress
c
      call comp_Fth(temp+dtemp,1.d0,'F',Fth,1.d0)
      Fthp = matmul(Fth,FP1)
      call INVERSEMATRIX3X3(Fthp,InvFP1)            
      call MULTIMATRIX3X3(F1,InvFP1,FE1)
      call TRANSPOSEMATRIX3X3(FE1,TransFE1)
      call MULTIMATRIX3X3(TransFE1,FE1,CE1) 
      call DETERMINANTMATRIX3X3(Fth,Jth)
      tmpprt = 1.d-8
      temptemp = temp + dtemp
      call comp_InEnBCC(FE1,temptemp,ElasStif,Beta1,'F',InEn1)
      temptemp = temp + dtemp + tmpPrt
      call comp_InEnBCC(FE1,temptemp,
     &     ElasStif,Beta1,'F',InEnPrt)
      dudt=(InEnPrt - InEn1) / tmpprt
      call ADDMATRIX3X3(5.d-1,CE1,-5.d-1,I2,EE1)
      call TENSOR4TENSOR2(ElasStif,EE1,PK2)
c  calculate the first Piola-Kirchhoff stress
      call TRANSPOSEMATRIX3X3(InvFP1,InTransFP1)
      call MULTIMATRIX3X3(PK2,InTransFP1,MTemp)
      call MULTIMATRIX3X3(FE1,MTemp,PK1)
      PK1 = Jth * PK1
      call DOTMATRIX3X3(dfgrd1 - dfgrd0, PK1, PFDot)
      RPL = PFDot - Rho_0 * (InEn1 - InEn0 - dudt * dtemp)
      if (dtime .NE. 0.d0) then
         RPL = RPL / dtime
         statev(84) = PFDot / dtime
      endif
c cc
c
c  calculate the work
      call DOTMATRIX3X3(PK1,DeltaF,dWorkF)
      WorkF1 = WorkF0 + dWorkF*dtime
c
c  calculate Cauchy stress and Kirchhoff stress
      call DETERMINANTMATRIX3X3(FE1,JE1)
      call MULTIMATRIX3X3(PK2,TransFE1,MTemp)
      call MULTIMATRIX3X3(FE1,MTemp,CauchyT)
      CauchyT = (1.d0 / JE1) * CauchyT
      call TRANSPOSEMATRIX3X3(F1,TransF1)
      call MULTIMATRIX3X3(PK1,TransF1,KirchS)
c
*  remark:
*  upto this stage, the variable that will be return to ABAQUS is still
*  expressed in terms of local variable as follows:
*    CauchyT  => stress
*    FP1      => statev(1-9)
*    Gamma1   => statev(10-33)
*    SlipRes1 => statev(34-57)
*    Beta1    => statev(58)
*    TagGamma => statev(59-82)
*  this parameter is not coverted yet since the convergent value at the
*  one step before is needed for computation of numerical consistent
*  tangent
*
c  return the Cauchy stress to ABAQUS
      stress(1) = CauchyT(1,1)
      stress(2) = CauchyT(2,2)
      stress(3) = CauchyT(3,3)
      stress(4) = CauchyT(1,2)
      stress(5) = CauchyT(1,3)
      stress(6) = CauchyT(2,3)
c
c  update the tag for active systems before computing the consistent
c  tangent
      do i = 1,24
         StateR(i + 58)  = dble(TagGamma(i))
      enddo
c
c  compute the numerical consistent tangent in terms of first Piola-
c  Kirchhoff stress and total deformation gradient
      ConTangent = .true.
      call CONSISTENTTANGENTBCC(StateR,DefGrad0,DefGrad1,DeltaT,DDPDDF,
     &     PK1,nstatv,ConTangent,temp+dtemp,
     &     ElasStif,SlipSys,ShearMod)
c      ddsdde = DDPDDF
c
c  convert the consistent tangent to be in correspondence to Lie's 
c  derivative of Kirchhoff stress
      BTgt = 0.d0
      do i = 1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3
                  Temp11 = 0.d0
                  do m = 1,3
                     do n = 1,3
                        ! TODO: Check why there are two versions??
                        Temp11 = Temp11 + F1(j,m)*F1(l,n) *
     &                       DDPDDF(i,m,k,n)
                        Temp11 = Temp11 + F1(j,m)*F1(l,n) *
     &                       DDPDDF(i,m,k,n) - F1(i,m)*PK1(k,m)*I2(j,l)
                     enddo
                  enddo
                  BTgt(i,j,k,l) = Temp11-KirchS(i,k)*I2(j,l)
                  BTgt(i,j,k,l) = Temp11
               enddo
            enddo
         enddo
      enddo
c
c  consistent tangent that corresponds to Jaumann's derivative of 
c  Kirchhoff stress
      HTgt = 0.d0
      do i = 1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3
                  HTgt(i,j,k,l) = BTgt(i,j,k,l) +
     &                 5.d-1*(KirchS(j,l)*I2(i,k) +
     &                 KirchS(i,k)*I2(j,l) +
     &                 KirchS(j,k)*I2(i,l) +
     &                 KirchS(i,l)*I2(j,k))
               enddo
            enddo
         enddo
	enddo
c
c  return the state-dependent variables to ABAQUS
	do i = 1,3
           do j = 1,3
              statev(j + (i - 1)*3)     = FP1(i,j)
           enddo
	enddo
	do i = 1,24
           statev(i + 9)   = Gamma1(i)
           statev(i + 33)  = SlipRes1(i)
           statev(i + 58)  = dble(TagGamma(i))
	enddo
	statev(58) = Beta1
	statev(83) = WorkF1
	statev(93) = dudt
c
c  return the consistent tangent to ABAQUS
      call DETERMINANTMATRIX3X3(F1,J1)
      ddsdde(1,1) = (1.d0/J1)*HTgt(1,1,1,1)
      ddsdde(1,2) = (1.d0/J1)*HTgt(1,1,2,2)
      ddsdde(1,3) = (1.d0/J1)*HTgt(1,1,3,3)
      ddsdde(1,4) = (1.d0/J1)*HTgt(1,1,1,2)
      ddsdde(1,5) = (1.d0/J1)*HTgt(1,1,1,3)
      ddsdde(1,6) = (1.d0/J1)*HTgt(1,1,2,3)
      ddsdde(2,1) = (1.d0/J1)*HTgt(2,2,1,1)
      ddsdde(2,2) = (1.d0/J1)*HTgt(2,2,2,2)
      ddsdde(2,3) = (1.d0/J1)*HTgt(2,2,3,3)
      ddsdde(2,4) = (1.d0/J1)*HTgt(2,2,1,2)
      ddsdde(2,5) = (1.d0/J1)*HTgt(2,2,1,3)
      ddsdde(2,6) = (1.d0/J1)*HTgt(2,2,2,3)
      ddsdde(3,1) = (1.d0/J1)*HTgt(3,3,1,1)
      ddsdde(3,2) = (1.d0/J1)*HTgt(3,3,2,2)
      ddsdde(3,3) = (1.d0/J1)*HTgt(3,3,3,3)
      ddsdde(3,4) = (1.d0/J1)*HTgt(3,3,1,2)
      ddsdde(3,5) = (1.d0/J1)*HTgt(3,3,1,3)
      ddsdde(3,6) = (1.d0/J1)*HTgt(3,3,2,3)
      ddsdde(4,1) = (1.d0/J1)*HTgt(1,2,1,1)
      ddsdde(4,2) = (1.d0/J1)*HTgt(1,2,2,2)
      ddsdde(4,3) = (1.d0/J1)*HTgt(1,2,3,3)
      ddsdde(4,4) = (1.d0/J1)*HTgt(1,2,1,2)
      ddsdde(4,5) = (1.d0/J1)*HTgt(1,2,1,3)
      ddsdde(4,6) = (1.d0/J1)*HTgt(1,2,2,3)
      ddsdde(5,1) = (1.d0/J1)*HTgt(1,3,1,1)
      ddsdde(5,2) = (1.d0/J1)*HTgt(1,3,2,2)
      ddsdde(5,3) = (1.d0/J1)*HTgt(1,3,3,3)
      ddsdde(5,4) = (1.d0/J1)*HTgt(1,3,1,2)
      ddsdde(5,5) = (1.d0/J1)*HTgt(1,3,1,3)
      ddsdde(5,6) = (1.d0/J1)*HTgt(1,3,2,3)
      ddsdde(6,1) = (1.d0/J1)*HTgt(2,3,1,1)
      ddsdde(6,2) = (1.d0/J1)*HTgt(2,3,2,2)
      ddsdde(6,3) = (1.d0/J1)*HTgt(2,3,3,3)
      ddsdde(6,4) = (1.d0/J1)*HTgt(2,3,1,2)
      ddsdde(6,5) = (1.d0/J1)*HTgt(2,3,1,3)
      ddsdde(6,6) = (1.d0/J1)*HTgt(2,3,2,3)
c
      statev(90) = rpl / J1
c  end of subroutine UMAT
      deallocate(DeltaGammaA,idGammaA,stat = errorstate)
      if (errorstate .gt. 0) then
         write(*,*)'UMATBCC: FAILED TO DEALLOCATE '
         stop
      endif
c
      return
      end	
c

c
c .......1.........2.........3.........4.........5.........6.........7..
c
c **********************************************************************
c
c  subroutine to rotate the material data in accordance with austenite
c  grain orientation with respect to the loading direction
c
c **********************************************************************
c
      subroutine GRAINORIENTATIONBCC(MaterialID, ElasStif, SlipSys,
     &     angles)
c
c  declaration of variables
      double precision ElasStif(3,3,3,3)
      double precision SlipSys(24,3,3)
      double precision RM(3,3)
      double precision angles(3)
      double precision Angle1, Angle2, Angle3
      double precision X(3,3), Z(3,3,3,3), Temp
      integer i,j,k,l,m,a,b,c,d,MaterialID
      logical PlotDebug,PlotResult
c
c      
c  get the rotation matrix and calculate its transpose
      Angle1 = angles(1)
      Angle2 = angles(2)
      Angle3 = angles(3)
      call ROTATIONMATRIX(Angle1,Angle2,Angle3,RM)
c
c  calculate the rotated slip system tensors
      do k = 1,24
         do i = 1,3
            do j = 1,3
               X(i,j) = SlipSys(k,i,j)
            enddo
         enddo
         do i = 1,3
            do j = 1,3
               Temp = 0.d0
               do a = 1,3
                  do b = 1,3
                     Temp = Temp + X(a,b)*RM(i,a)*RM(j,b)
                  enddo
               enddo
               SlipSys(k,i,j)  = Temp
            enddo
         enddo
      enddo
c     
c  debugging the rotated slip systems
c
c  calculate the rotated elastic stiffnesses
      Z = ElasStif
      do i = 1,3
      do j = 1,3
      do k = 1,3
      do l = 1,3
         Temp = 0.d0
         do a = 1,3
         do b = 1,3
         do c = 1,3
         do d = 1,3
            Temp = Temp + Z(a,b,c,d)*RM(i,a)*RM(j,b)*RM(k,c)*RM(l,d)
         enddo
         enddo
         enddo
         enddo
         ElasStif(i,j,k,l) = Temp
      enddo
      enddo
      enddo
      enddo
c
c  debugging the elastic stiffnesses
c
c    ..return the rotated Schmid-tensors and elastic stiffness tensor
      return
      end

c
c .......1.........2.........3.........4.........5.........6.........7..
c
c **********************************************************************
c
c  subroutine to construct second-and fourth-order identity tensors
c
c **********************************************************************
c
      subroutine IDENTITYTENSOR(I2,I4)
c
c  declaration of variables
      double precision I2(3,3),I4(3,3,3,3)
      integer i
c
c  construct the second-order identity tensor
      I2 = 0.d0
      do i = 1,3
         I2(i,i) = 1.d0
      enddo
c
c  construct the fourth-order identity tensor
      do i = 1,3
      do j = 1,3
      do k = 1,3
      do l = 1,3
         I4(i,j,k,l) = 5.d-1*(I2(j,k)*I2(i,l) + I2(i,k)*I2(j,l))
      enddo
      enddo
      enddo
      enddo
c
      return
      end
c
c
c .......1.........2.........3.........4.........5.........6.........7..
c
c **********************************************************************
c
c  subroutine to read the material properties from files
c
c **********************************************************************
c
      subroutine MATERIALPROPERTIESBCC(MaterialID,
     &     ElasStif,SlipSys,ShearMod,angles)
c
c  declaration of variables
      double precision ElasStif(3,3,3,3),SlipSys(24,3,3),ShearMod
      double precision angles(3)
      integer i,j,k,MaterialID
      logical PlotDebug,PlotResult
c  get the elastic stiffness, slip systems
      include 'TRIP_FerriteData-v01.for'
c
c  calculate the rotated elastic stiffness, slip systems
      call GRAINORIENTATIONBCC(MaterialID,ElasStif,SlipSys,angles)
c     
      return
      end

c
c .......1.........2.........3.........4.........5.........6.........7..
c
c **********************************************************************
c
c  subroutine to define the rotation matrix using 323 euler rotations
c
c **********************************************************************
c
      subroutine ROTATIONMATRIX(Angle1,Angle2,Angle3,RM)
c
c  declaration of variables
      double precision RM(3,3),R1(3,3),R2(3,3),R3(3,3),Angle1,Angle2,
     & Angle3,Temp(3,3)
      integer i,j,k
c
c  rotation matrix of the first rotation (Angle1)
	R1 = 0.d0
	R1(1,1) =  cos(Angle1)
	R1(1,2) =  sin(Angle1)
	R1(2,1) = -sin(Angle1)
	R1(2,2) =  cos(Angle1)
	R1(3,3) =  1.d0
c
c  rotation matrix of the second rotation (Angle2)
	R2 = 0.d0
	R2(1,1) =  cos(Angle2)
	R2(1,3) = -sin(Angle2)
	R2(2,2) =  1.d0
	R2(3,1) =  sin(Angle2)
	R2(3,3) =  cos(Angle2)
c
c  rotation matrix of the third rotation (Angle3)
	R3 = 0.d0
	R3(1,1) =  cos(Angle3)
	R3(1,2) =  sin(Angle3)
	R3(2,1) = -sin(Angle3)
	R3(2,2) =  cos(Angle3)
	R3(3,3) =  1.d0
c
c  calculate the overall rotation matrix
      call MULTIMATRIX3X3(R3,R2,Temp)
	call MULTIMATRIX3X3(Temp,R1,RM)
c
      return
      end
c      
c
c .......1.........2.........3.........4.........5.........6.........7..
c
c **********************************************************************
c
c  subroutine to check the consistency condition of active slip systems
c
c **********************************************************************
c
      subroutine SLIPSYSTEMCHECKBCC(PDF,PDFCrit,DeltaGammaA,idGammaA,
     & nGammaA,TagGamma,UpgradeSlipSys,iterP)
c
c  declaration of variables
      double precision PDF(24),PDFCrit(24),DeltaGammaA(nGammaA)
      integer i,nGammaA,TagGamma(24),idGammaA(nGammaA),iterP
      logical UpgradeSlipSys,PlotDebug,PlotResult
c
c 
c
c  first check: if there are any active systems become inactive
      do i = 1,nGammaA
         if (DeltaGammaA(i) .le. 0.d0) then
	    TagGamma(idGammaA(i)) = 0
	    UpgradeSlipSys = .true.
         endif
      enddo
      if (UpgradeSlipSys) goto 64
c
c  second check: if there are any inactive systems become active
      do i = 1,24
         if (TagGamma(i) .eq. 0) then
	    if (iterP .le. 6) then
               if (PDF(i) .gt. PDFCrit(i)) then
                  TagGamma(i)    = 1
                  UpgradeSlipSys = .true.
               endif
	    endif
         endif
      enddo
c
64    continue
      return
      end

	SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1  STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2  CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3  NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
c
c **********************************************************************
c  STATEV(151 OR 93) = dU/dT_n+1
c
c
c **********************************************************************
c
c
      include 'aba_param.inc'
      include 'constants.for'
c
      CHARACTER*80 CMNAME
      DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     1  DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD),
     2  TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3)
c
c
*       if (NOEL==1) then
*       	if (NPT==1) then 
*       	  write (*,*) 'UMATHT'
*       	  write(*,*) 'Increment=',KINC
*!       	  write(*,*) 'Step Number=',KSTEP
*!         write(*,*) 'Material Props:',props(1),props(2),props(3)
*          write(*,*) 'Temp=',TEMP ,'DTemp=',DTEMP
*	  write(*,*) 'U0 =',U
*	  write(*,*) 'U1 =',STATEV(1)
*	  write(*,*) '    '
*        endif
*       endif

      COND=HCond
c      U = STATEV(1)
      if (nstatv.eq.151) DUDT = STATEV(151)
      if (nstatv.eq.93)  DUDT = STATEV(93)
c      if (debugoption) then
c        write (*,*) 'NSTATV = ', NSTATV
c        write (*,*) 'STATEV(93) = ', STATEV(93)
c     endif
      du = dudt*dtemp
      u = u+du
c      write(*,*) 'DUDT =',DUDT
c      write(*,*) 'NTGRD =',NTGRD
      DFDG=0.d0
      Do I=1, NTGRD
        FLUX(I)=-COND*DTEMDX(I)
        DFDG(I,I)=-COND
      End Do
      DUDG=0.d0
      DFDT=0.d0
c
c
      RETURN
      END      
