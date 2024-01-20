MODULE mynn_dmp_mf
use module_bl_mynn_common,only: cp, cpv, ep_2, ep_3, grav, gtr, onethird, r_v, tv0, xlvcp
use mynn_functions
IMPLICIT NONE

CONTAINS
! ==================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine is the Dynamic Multi-Plume (DMP) Mass-Flux Scheme.
!! 
!! dmp_mf() calculates the nonlocal turbulent transport from the dynamic
!! multiplume mass-flux scheme as well as the shallow-cumulus component of 
!! the subgrid clouds. Note that this mass-flux scheme is called when the
!! namelist paramter \p bl_mynn_edmf is set to 1 (recommended).
!!
!! Much thanks to Kay Suslj of NASA-JPL for contributing the original version
!! of this mass-flux scheme. Considerable changes have been made from it's
!! original form. Some additions include:
!!  -# scale-aware tapering as dx -> 0
!!  -# transport of TKE (extra namelist option)
!!  -# Chaboureau-Bechtold cloud fraction & coupling to radiation (when icloud_bl > 0)
!!  -# some extra limits for numerical stability
!!
!! This scheme remains under development, so consider it experimental code. 
!!
  SUBROUTINE DMP_mf(                            &
                 & kts,kte,dt,zw,dz,p,rho,      &
                 & momentum_opt,                &
                 & tke_opt,                     &
                 & scalar_opt,                  &
                 & u,v,w,th,thl,thv,tk,         &
                 & qt,qv,qc,qke,                &
                 & qnc,qni,qnwfa,qnifa,qnbca,   &
                 & exner,vt,vq,sgm,             &
                 & ust,flt,fltv,flq,flqv,       &
                 & pblh,kpbl,dx,landsea,ts,     &
            ! outputs - updraft properties   
                 & edmf_a,edmf_w,               &
                 & edmf_qt,edmf_thl,            &
                 & edmf_ent,edmf_qc,            &
            ! outputs - variables needed for solver 
                 & s_aw,s_awthl,s_awqt,         &
                 & s_awqv,s_awqc,               &
                 & s_awu,s_awv,s_awqke,         &
                 & s_awqnc,s_awqni,             &
                 & s_awqnwfa,s_awqnifa,         &
                 & s_awqnbca,                   &
                 & sub_thl,sub_sqv,             &
                 & sub_u,sub_v,                 &
                 & det_thl,det_sqv,det_sqc,     &
                 & det_u,det_v,                 &
            ! chem/smoke
                 & nchem,chem1,s_awchem,        &
                 & mix_chem,                    &
            ! in/outputs - subgrid scale clouds
                 & qc_bl1d,cldfra_bl1d,         &
                 & qc_bl1D_old,cldfra_bl1D_old, &
            ! inputs - flags for moist arrays
                 & F_QC,F_QI,                   &
                 & F_QNC,F_QNI,                 &
                 & F_QNWFA,F_QNIFA,F_QNBCA,     &
                 & Psig_shcu,                   &
            ! output info
                 &nup2,ktop,maxmf,ztop,         &
            ! unputs for stochastic perturbations
                 &spp_pbl,rstoch_col,           &
                 &env_subs                      ) 

  ! inputs:
     INTEGER, INTENT(IN) :: KTS,KTE,KPBL,momentum_opt,tke_opt,scalar_opt
     LOGICAL, INTENT(IN) :: env_subs
#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

! Stochastic 
     INTEGER,  INTENT(IN)          :: spp_pbl
     REAL, DIMENSION(KTS:KTE)      :: rstoch_col

     REAL,DIMENSION(KTS:KTE), INTENT(IN) ::                            &
                   u,v,w,th,thl,tk,qt,qv,qc,                           &
                   exner,dz,THV,P,rho,qke,qnc,qni,qnwfa,qnifa,qnbca
     REAL,DIMENSION(KTS:KTE+1), INTENT(IN) :: zw    !height at full-sigma
     REAL, INTENT(IN) :: dt,ust,flt,fltv,flq,flqv,pblh,                &
                         dx,psig_shcu,landsea,ts
     LOGICAL, OPTIONAL :: f_qc,f_qi,f_qnc,f_qni,                       &
                   f_qnwfa,f_qnifa,f_qnbca

  ! outputs - updraft properties
     REAL,DIMENSION(KTS:KTE), INTENT(OUT) :: edmf_a,edmf_w,            &
                      & edmf_qt,edmf_thl,edmf_ent,edmf_qc
     !add one local edmf variable:
     REAL,DIMENSION(KTS:KTE) :: edmf_th
  ! output
     INTEGER, INTENT(OUT) :: nup2,ktop
     REAL, INTENT(OUT) :: maxmf,ztop
  ! outputs - variables needed for solver - sum ai*rho*wis_awphi
     REAL,DIMENSION(KTS:KTE+1) :: s_aw,s_awthl,s_awqt,                 &
                         s_awqv,s_awqc,s_awqnc,s_awqni,                &
                         s_awqnwfa,s_awqnifa,s_awqnbca,                &
                         s_awu,s_awv,s_awqke,s_aw2

     REAL,DIMENSION(KTS:KTE), INTENT(INOUT) :: qc_bl1d,cldfra_bl1d,    &
                                       qc_bl1d_old,cldfra_bl1d_old

    INTEGER, PARAMETER :: nup=10, debug_mf=0

  !------------- local variables -------------------
  ! updraft properties defined on interfaces (k=1 is the top of the
  ! first model layer
     REAL,DIMENSION(KTS:KTE+1,1:NUP) :: UPW,UPTHL,UPQT,UPQC,UPQV,      &
                                        UPA,UPU,UPV,UPTHV,UPQKE,UPQNC, &
                                        UPQNI,UPQNWFA,UPQNIFA,UPQNBCA
  ! entrainment variables
     REAL,DIMENSION(KTS:KTE,1:NUP) :: ENT,ENTf
     INTEGER,DIMENSION(KTS:KTE,1:NUP) :: ENTi
  ! internal variables
     INTEGER :: K,I,k50
     REAL :: fltv2,wstar,qstar,thstar,sigmaW,sigmaQT,sigmaTH,z0,       &
             pwmin,pwmax,wmin,wmax,wlv,Psig_w,maxw,maxqc,wpbl
     REAL :: B,QTn,THLn,THVn,QCn,Un,Vn,QKEn,QNCn,QNIn,                 &
             QNWFAn,QNIFAn,QNBCAn,                                     &
             Wn2,Wn,EntEXP,EntEXM,EntW,BCOEFF,THVkm1,THVk,Pk,rho_int

  ! w parameters
     REAL,PARAMETER :: &
          &Wa=2./3.,   &
          &Wb=0.002,   &
          &Wc=1.5 
        
  ! Lateral entrainment parameters ( L0=100 and ENT0=0.1) were taken from
  ! Suselj et al (2013, jas). Note that Suselj et al (2014,waf) use L0=200 and ENT0=0.2.
     REAL,PARAMETER :: &
         & L0=100.,    &
         & ENT0=0.1

  ! Implement ideas from Neggers (2016, JAMES):
     REAL, PARAMETER :: Atot = 0.10 ! Maximum total fractional area of all updrafts
     REAL, PARAMETER :: lmax = 1000.! diameter of largest plume
     REAL, PARAMETER :: dl   = 100. ! diff size of each plume - the differential multiplied by the integrand
     REAL, PARAMETER :: dcut = 1.2  ! max diameter of plume to parameterize relative to dx (km)
     REAL ::  d            != -2.3 to -1.7  ;=-1.9 in Neggers paper; power law exponent for number density (N=Cl^d).
          ! Note that changing d to -2.0 makes each size plume equally contribute to the total coverage of all plumes.
          ! Note that changing d to -1.7 doubles the area coverage of the largest plumes relative to the smallest plumes.
     REAL :: cn,c,l,n,an2,hux,maxwidth,wspd_pbl,cloud_base,width_flx

  ! chem/smoke
     INTEGER, INTENT(IN) :: nchem
     REAL,DIMENSION(:, :) :: chem1
     REAL,DIMENSION(kts:kte+1, nchem) :: s_awchem
     REAL,DIMENSION(nchem) :: chemn
     REAL,DIMENSION(KTS:KTE+1,1:NUP, nchem) :: UPCHEM
     INTEGER :: ic
     REAL,DIMENSION(KTS:KTE+1, nchem) :: edmf_chem
     LOGICAL, INTENT(IN) :: mix_chem

  !JOE: add declaration of ERF
   REAL :: ERF

   LOGICAL :: superadiabatic

  ! VARIABLES FOR CHABOUREAU-BECHTOLD CLOUD FRACTION
   REAL,DIMENSION(KTS:KTE), INTENT(INOUT) :: vt, vq, sgm
   REAL :: sigq,xl,rsl,cpm,a,qmq,mf_cf,Aup,Q1,diffqt,qsat_tk,&
           Fng,qww,alpha,beta,bb,f,pt,t,q2p,b9,satvp,rhgrid, &
           Ac_mf,Ac_strat,qc_mf
   REAL, PARAMETER :: cf_thresh = 0.5 ! only overwrite stratus CF less than this value

  ! Variables for plume interpolation/saturation check
   REAL,DIMENSION(KTS:KTE) :: exneri,dzi
   REAL :: THp, QTp, QCp, QCs, esat, qsl
   REAL :: csigma,acfac,ac_wsp,ac_cld

   !plume overshoot
   INTEGER :: overshoot
   REAL :: bvf, Frz, dzp

   !Flux limiter: not let mass-flux of heat between k=1&2 exceed (fluxportion)*(surface heat flux).
   !This limiter makes adjustments to the entire column.
   REAL :: adjustment, flx1
   REAL, PARAMETER :: fluxportion=0.75 ! set liberally, so has minimal impact. 0.5 starts to have a noticeable impact
                                       ! over land (decrease maxMF by 10-20%), but no impact over water.

   !Subsidence
   REAL,DIMENSION(KTS:KTE) :: sub_thl,sub_sqv,sub_u,sub_v,    &  !tendencies due to subsidence
                      det_thl,det_sqv,det_sqc,det_u,det_v,    &  !tendencied due to detrainment
                 envm_a,envm_w,envm_thl,envm_sqv,envm_sqc,    &
                                       envm_u,envm_v  !environmental variables defined at middle of layer
   REAL,DIMENSION(KTS:KTE+1) ::  envi_a,envi_w        !environmental variables defined at model interface
   REAL :: temp,sublim,qc_ent,qv_ent,qt_ent,thl_ent,detrate,  &
           detrateUV,oow,exc_fac,aratio,detturb,qc_grid,qc_sgs,&
           qc_plume,exc_heat,exc_moist,tk_int
   REAL, PARAMETER :: Cdet   = 1./45.
   REAL, PARAMETER :: dzpmax = 300. !limit dz used in detrainment - can be excessing in thick layers
   !parameter "Csub" determines the propotion of upward vertical velocity that contributes to
   !environmenatal subsidence. Some portion is expected to be compensated by downdrafts instead of
   !gentle environmental subsidence. 1.0 assumes all upward vertical velocity in the mass-flux scheme
   !is compensated by "gentle" environmental subsidence. 
   REAL, PARAMETER :: Csub=0.25

   !Factor for the pressure gradient effects on momentum transport
   REAL, PARAMETER :: pgfac = 0.00  ! Zhang and Wu showed 0.4 is more appropriate for lower troposphere
   REAL :: Uk,Ukm1,Vk,Vkm1,dxsa

! check the inputs
!     print *,'dt',dt
!     print *,'dz',dz
!     print *,'u',u
!     print *,'v',v
!     print *,'thl',thl
!     print *,'qt',qt
!     print *,'ust',ust
!     print *,'flt',flt
!     print *,'flq',flq
!     print *,'pblh',pblh

! Initialize individual updraft properties
  UPW=0.
  UPTHL=0.
  UPTHV=0.
  UPQT=0.
  UPA=0.
  UPU=0.
  UPV=0.
  UPQC=0.
  UPQV=0.
  UPQKE=0.
  UPQNC=0.
  UPQNI=0.
  UPQNWFA=0.
  UPQNIFA=0.
  UPQNBCA=0.
  IF ( mix_chem ) THEN
     UPCHEM(KTS:KTE+1,1:NUP,1:nchem)=0.0
  ENDIF

  ENT=0.001
! Initialize mean updraft properties
  edmf_a  =0.
  edmf_w  =0.
  edmf_qt =0.
  edmf_thl=0.
  edmf_ent=0.
  edmf_qc =0.
  IF ( mix_chem ) THEN
     edmf_chem(kts:kte+1,1:nchem) = 0.0
  ENDIF

! Initialize the variables needed for implicit solver
  s_aw=0.
  s_awthl=0.
  s_awqt=0.
  s_awqv=0.
  s_awqc=0.
  s_awu=0.
  s_awv=0.
  s_awqke=0.
  s_awqnc=0.
  s_awqni=0.
  s_awqnwfa=0.
  s_awqnifa=0.
  s_awqnbca=0.
  IF ( mix_chem ) THEN
     s_awchem(kts:kte+1,1:nchem) = 0.0
  ENDIF

! Initialize explicit tendencies for subsidence & detrainment
  sub_thl = 0.
  sub_sqv = 0.
  sub_u = 0.
  sub_v = 0.
  det_thl = 0.
  det_sqv = 0.
  det_sqc = 0.
  det_u = 0.
  det_v = 0.

  ! Taper off MF scheme when significant resolved-scale motions
  ! are present This function needs to be asymetric...
  k      = 1
  maxw   = 0.0
  cloud_base  = 9000.0
!  DO WHILE (ZW(k) < pblh + 500.)
  DO k=1,kte-1
     IF(zw(k) > pblh + 500.) exit

     wpbl = w(k)
     IF(w(k) < 0.)wpbl = 2.*w(k)
     maxw = MAX(maxw,ABS(wpbl))

     !Find highest k-level below 50m AGL
     IF(ZW(k)<=50.)k50=k

     !Search for cloud base
     qc_sgs = MAX(qc(k), qc_bl1d(k)*cldfra_bl1d(k))
     IF(qc_sgs> 1E-5 .AND. cloud_base == 9000.0)THEN
       cloud_base = 0.5*(ZW(k)+ZW(k+1))
     ENDIF

     !k = k + 1
  ENDDO
  !print*," maxw before manipulation=", maxw
  maxw = MAX(0.,maxw - 1.0)     ! do nothing for small w (< 1 m/s), but
  Psig_w = MAX(0.0, 1.0 - maxw) ! linearly taper off for w > 1.0 m/s
  Psig_w = MIN(Psig_w, Psig_shcu)
  !print*," maxw=", maxw," Psig_w=",Psig_w," Psig_shcu=",Psig_shcu

  !Completely shut off MF scheme for strong resolved-scale vertical velocities.
  fltv2 = fltv
  IF(Psig_w == 0.0 .and. fltv > 0.0) fltv2 = -1.*fltv

  ! If surface buoyancy is positive we do integration, otherwise no.
  ! Also, ensure that it is at least slightly superadiabatic up through 50 m
  superadiabatic = .false.
  IF((landsea-1.5).GE.0)THEN
     hux = -0.001   ! WATER  ! dT/dz must be < - 0.1 K per 100 m.
  ELSE
     hux = -0.005  ! LAND    ! dT/dz must be < - 0.5 K per 100 m.
  ENDIF
  DO k=1,MAX(1,k50-1) !use "-1" because k50 used interface heights (zw). 
    IF (k == 1) then
      IF ((th(k)-ts)/(0.5*dz(k)) < hux) THEN
        superadiabatic = .true.
      ELSE
        superadiabatic = .false.
        exit
      ENDIF
    ELSE
      IF ((th(k)-th(k-1))/(0.5*(dz(k)+dz(k-1))) < hux) THEN
        superadiabatic = .true.
      ELSE
        superadiabatic = .false.
        exit
      ENDIF
    ENDIF
  ENDDO

  ! Determine the numer of updrafts/plumes in the grid column:
  ! Some of these criteria may be a little redundant but useful for bullet-proofing.
  !   (1) largest plume = 1.0 * dx.
  !   (2) Apply a scale-break, assuming no plumes with diameter larger than PBLH can exist.
  !   (3) max plume size beneath clouds deck approx = 0.5 * cloud_base.
  !   (4) add wspd-dependent limit, when plume model breaks down. (hurricanes)
  !   (5) limit to reduce max plume sizes in weakly forced conditions. This is only
  !       meant to "soften" the activation of the mass-flux scheme.
  ! Criteria (1)
    NUP2 = max(1,min(NUP,INT(dx*dcut/dl)))
  !Criteria (2)
    maxwidth = 1.1*PBLH 
  ! Criteria (3)
    maxwidth = MIN(maxwidth,0.5*cloud_base)
  ! Criteria (4)
    wspd_pbl=SQRT(MAX(u(kts)**2 + v(kts)**2, 0.01))
    !Note: area fraction (acfac) is modified below
  ! Criteria (5) - only a function of flt (not fltv)
    if ((landsea-1.5).LT.0) then  !land
      !width_flx = MAX(MIN(1000.*(0.6*tanh((flt - 0.050)/0.03) + .5),1000.), 0.)
      width_flx = MAX(MIN(1000.*(0.6*tanh((flt - 0.040)/0.03) + .5),1000.), 0.) 
    else                          !water
      width_flx = MAX(MIN(1000.*(0.6*tanh((flt - 0.003)/0.01) + .5),1000.), 0.)
    endif
    maxwidth = MIN(maxwidth,width_flx)
  ! Convert maxwidth to number of plumes
    NUP2 = MIN(MAX(INT((maxwidth - MOD(maxwidth,100.))/100), 0), NUP2)

  !Initialize values for 2d output fields:
  ktop = 0
  ztop = 0.0
  maxmf= 0.0

  IF ( fltv2 > 0.002 .AND. NUP2 .GE. 1 .AND. superadiabatic) then
    !PRINT*," Conditions met to run mass-flux scheme",fltv2,pblh

    ! Find coef C for number size density N
    cn = 0.
    d=-1.9  !set d to value suggested by Neggers 2015 (JAMES).
    !d=-1.9 + .2*tanh((fltv2 - 0.05)/0.15) 
    do I=1,NUP !NUP2
       IF(I > NUP2) exit
       l  = dl*I                            ! diameter of plume
       cn = cn + l**d * (l*l)/(dx*dx) * dl  ! sum fractional area of each plume
    enddo
    C = Atot/cn   !Normalize C according to the defined total fraction (Atot)

    ! Make updraft area (UPA) a function of the buoyancy flux
    if ((landsea-1.5).LT.0) then  !land
       !acfac = .5*tanh((fltv2 - 0.03)/0.09) + .5
       !acfac = .5*tanh((fltv2 - 0.02)/0.09) + .5
       acfac = .5*tanh((fltv2 - 0.02)/0.05) + .5
    else                          !water
       acfac = .5*tanh((fltv2 - 0.01)/0.03) + .5
    endif
    !add a windspeed-dependent adjustment to acfac that tapers off
    !the mass-flux scheme linearly above sfc wind speeds of 20 m/s:
    ac_wsp = 1.0 - min(max(wspd_pbl - 20.0, 0.0), 10.0)/10.0
    !reduce area fraction beneath cloud bases < 1200 m AGL
    ac_cld = min(cloud_base/1200., 1.0)
    acfac  = acfac * min(ac_wsp, ac_cld)

    ! Find the portion of the total fraction (Atot) of each plume size:
    An2 = 0.
    do I=1,NUP !NUP2
       IF(I > NUP2) exit
       l  = dl*I                            ! diameter of plume
       N  = C*l**d                          ! number density of plume n
       UPA(1,I) = N*l*l/(dx*dx) * dl        ! fractional area of plume n

       UPA(1,I) = UPA(1,I)*acfac
       An2 = An2 + UPA(1,I)                 ! total fractional area of all plumes
       !print*," plume size=",l,"; area=",UPA(1,I),"; total=",An2
    end do

    ! set initial conditions for updrafts
    z0=50.
    pwmin=0.1       ! was 0.5
    pwmax=0.4       ! was 3.0

    wstar=max(1.E-2,(gtr*fltv2*pblh)**(onethird))
    qstar=max(flq,1.0E-5)/wstar
    thstar=flt/wstar

    IF((landsea-1.5).GE.0)THEN
       csigma = 1.34   ! WATER
    ELSE
       csigma = 1.34   ! LAND
    ENDIF

    if (env_subs) then
       exc_fac = 0.0
    else
       if ((landsea-1.5).GE.0) then
         !water: increase factor to compensate for decreased pwmin/pwmax
         exc_fac = 0.58*4.0*min(cloud_base/1000., 1.0)
       else
         !land: no need to increase factor - already sufficiently large superadiabatic layers
         exc_fac = 0.58
       endif
    endif

    !Note: sigmaW is typically about 0.5*wstar
    sigmaW =csigma*wstar*(z0/pblh)**(onethird)*(1 - 0.8*z0/pblh)
    sigmaQT=csigma*qstar*(z0/pblh)**(onethird)
    sigmaTH=csigma*thstar*(z0/pblh)**(onethird)

    !Note: Given the pwmin & pwmax set above, these max/mins are
    !      rarely exceeded. 
    wmin=MIN(sigmaW*pwmin,0.1)
    wmax=MIN(sigmaW*pwmax,0.5)

    !SPECIFY SURFACE UPDRAFT PROPERTIES AT MODEL INTERFACE BETWEEN K = 1 & 2
    DO I=1,NUP !NUP2
       IF(I > NUP2) exit
       wlv=wmin+(wmax-wmin)/NUP2*(i-1)

       !SURFACE UPDRAFT VERTICAL VELOCITY
       UPW(1,I)=wmin + REAL(i)/REAL(NUP)*(wmax-wmin)
       !IF (UPW(1,I) > 0.5*ZW(2)/dt) UPW(1,I) = 0.5*ZW(2)/dt

       UPU(1,I)=(U(KTS)*DZ(KTS+1)+U(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPV(1,I)=(V(KTS)*DZ(KTS+1)+V(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPQC(1,I)=0.0
       !UPQC(1,I)=(QC(KTS)*DZ(KTS+1)+QC(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))

       exc_heat = exc_fac*UPW(1,I)*sigmaTH/sigmaW
       UPTHV(1,I)=(THV(KTS)*DZ(KTS+1)+THV(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1)) &
           &     + exc_heat
!was       UPTHL(1,I)= UPTHV(1,I)/(1.+svp1*UPQT(1,I))  !assume no saturated parcel at surface
       UPTHL(1,I)=(THL(KTS)*DZ(KTS+1)+THL(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1)) &
           &     + exc_heat

       !calculate exc_moist by use of surface fluxes
       exc_moist=exc_fac*UPW(1,I)*sigmaQT/sigmaW
       !calculate exc_moist by conserving rh:
!       tk_int  =(tk(kts)*dz(kts+1)+tk(kts+1)*dz(kts))/(dz(kts+1)+dz(kts))
!       pk      =(p(kts)*dz(kts+1)+p(kts+1)*dz(kts))/(dz(kts+1)+dz(kts))
!       qtk     =(qt(kts)*dz(kts+1)+qt(kts+1)*dz(kts))/(dz(kts)+dz(kts+1))
!       qsat_tk = qsat_blend(tk_int,  pk)    ! saturation water vapor mixing ratio at tk and p
!       rhgrid  =MAX(MIN(1.0,qtk/MAX(1.E-8,qsat_tk)),0.001)
!       tk_int  = tk_int + exc_heat
!       qsat_tk = qsat_blend(tk_int,  pk) 
!       exc_moist= max(rhgrid*qsat_tk - qtk, 0.0)
       UPQT(1,I)=(QT(KTS)*DZ(KTS+1)+QT(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))&
            &     +exc_moist

       UPQKE(1,I)=(QKE(KTS)*DZ(KTS+1)+QKE(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPQNC(1,I)=(QNC(KTS)*DZ(KTS+1)+QNC(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPQNI(1,I)=(QNI(KTS)*DZ(KTS+1)+QNI(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPQNWFA(1,I)=(QNWFA(KTS)*DZ(KTS+1)+QNWFA(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPQNIFA(1,I)=(QNIFA(KTS)*DZ(KTS+1)+QNIFA(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPQNBCA(1,I)=(QNBCA(KTS)*DZ(KTS+1)+QNBCA(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
    ENDDO

    IF ( mix_chem ) THEN
      DO I=1,NUP !NUP2
        IF(I > NUP2) exit
        do ic = 1,nchem
          UPCHEM(1,I,ic)=(chem1(KTS,ic)*DZ(KTS+1)+chem1(KTS+1,ic)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
        enddo
      ENDDO
    ENDIF

    !Initialize environmental variables which can be modified by detrainment
    DO k=kts,kte
       envm_thl(k)=THL(k)
       envm_sqv(k)=QV(k)
       envm_sqc(k)=QC(k)
       envm_u(k)=U(k)
       envm_v(k)=V(k)
    ENDDO

    !dxsa is scale-adaptive factor governing the pressure-gradient term of the momentum transport
    dxsa = 1. - MIN(MAX((12000.0-dx)/(12000.0-3000.0), 0.), 1.)

    ! do integration  updraft
    DO I=1,NUP !NUP2
       IF(I > NUP2) exit
       QCn = 0.
       overshoot = 0
       l  = dl*I                            ! diameter of plume
       DO k=KTS+1,KTE-1
          !Entrainment from Tian and Kuang (2016)
          !ENT(k,i) = 0.35/(MIN(MAX(UPW(K-1,I),0.75),1.9)*l)
          wmin = 0.3 + l*0.0005 !* MAX(pblh-ZW(k+1), 0.0)/pblh
          ENT(k,i) = 0.33/(MIN(MAX(UPW(K-1,I),wmin),0.9)*l)

          !Entrainment from Negggers (2015, JAMES)
          !ENT(k,i) = 0.02*l**-0.35 - 0.0009
          !ENT(k,i) = 0.04*l**-0.50 - 0.0009   !more plume diversity
          !ENT(k,i) = 0.04*l**-0.495 - 0.0009  !"neg1+"

          !Minimum background entrainment 
          ENT(k,i) = max(ENT(k,i),0.0003)
          !ENT(k,i) = max(ENT(k,i),0.05/ZW(k))  !not needed for Tian and Kuang

          !JOE - increase entrainment for plumes extending very high.
          IF(ZW(k) >= MIN(pblh+1500., 4000.))THEN
            ENT(k,i)=ENT(k,i) + (ZW(k)-MIN(pblh+1500.,4000.))*5.0E-6
          ENDIF

          !SPP
          ENT(k,i) = ENT(k,i) * (1.0 - rstoch_col(k))

          ENT(k,i) = min(ENT(k,i),0.9/(ZW(k+1)-ZW(k)))

          ! Define environment U & V at the model interface levels
          Uk  =(U(k)*DZ(k+1)+U(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
          Ukm1=(U(k-1)*DZ(k)+U(k)*DZ(k-1))/(DZ(k-1)+DZ(k))
          Vk  =(V(k)*DZ(k+1)+V(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
          Vkm1=(V(k-1)*DZ(k)+V(k)*DZ(k-1))/(DZ(k-1)+DZ(k))

          ! Linear entrainment:
          EntExp= ENT(K,I)*(ZW(k+1)-ZW(k))
          EntExm= EntExp*0.3333    !reduce entrainment for momentum
          QTn =UPQT(k-1,I) *(1.-EntExp) + QT(k)*EntExp
          THLn=UPTHL(k-1,I)*(1.-EntExp) + THL(k)*EntExp
          Un  =UPU(k-1,I)  *(1.-EntExm) + U(k)*EntExm + dxsa*pgfac*(Uk - Ukm1)
          Vn  =UPV(k-1,I)  *(1.-EntExm) + V(k)*EntExm + dxsa*pgfac*(Vk - Vkm1)
          QKEn=UPQKE(k-1,I)*(1.-EntExp) + QKE(k)*EntExp
          QNCn=UPQNC(k-1,I)*(1.-EntExp) + QNC(k)*EntExp
          QNIn=UPQNI(k-1,I)*(1.-EntExp) + QNI(k)*EntExp
          QNWFAn=UPQNWFA(k-1,I)*(1.-EntExp) + QNWFA(k)*EntExp
          QNIFAn=UPQNIFA(k-1,I)*(1.-EntExp) + QNIFA(k)*EntExp
          QNBCAn=UPQNBCA(k-1,I)*(1.-EntExp) + QNBCA(k)*EntExp

          !capture the updated qc, qt & thl modified by entranment alone,
          !since they will be modified later if condensation occurs.
          qc_ent  = QCn
          qt_ent  = QTn
          thl_ent = THLn

          ! Exponential Entrainment:
          !EntExp= exp(-ENT(K,I)*(ZW(k)-ZW(k-1)))
          !QTn =QT(K) *(1-EntExp)+UPQT(K-1,I)*EntExp
          !THLn=THL(K)*(1-EntExp)+UPTHL(K-1,I)*EntExp
          !Un  =U(K)  *(1-EntExp)+UPU(K-1,I)*EntExp
          !Vn  =V(K)  *(1-EntExp)+UPV(K-1,I)*EntExp
          !QKEn=QKE(k)*(1-EntExp)+UPQKE(K-1,I)*EntExp

          if ( mix_chem ) then
            do ic = 1,nchem
              ! Exponential Entrainment:
              !chemn(ic) = chem(k,ic)*(1-EntExp)+UPCHEM(K-1,I,ic)*EntExp
              ! Linear entrainment:
              chemn(ic)=UPCHEM(k-1,i,ic)*(1.-EntExp) + chem1(k,ic)*EntExp
            enddo
          endif

          ! Define pressure at model interface
          Pk    =(P(k)*DZ(k+1)+P(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
          ! Compute plume properties thvn and qcn
          call condensation_edmf(QTn,THLn,Pk,ZW(k+1),THVn,QCn)

          ! Define environment THV at the model interface levels
          THVk  =(THV(k)*DZ(k+1)+THV(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
          THVkm1=(THV(k-1)*DZ(k)+THV(k)*DZ(k-1))/(DZ(k-1)+DZ(k))

!          B=g*(0.5*(THVn+UPTHV(k-1,I))/THV(k-1) - 1.0)
          B=grav*(THVn/THVk - 1.0)
          IF(B>0.)THEN
            BCOEFF = 0.15        !w typically stays < 2.5, so doesnt hit the limits nearly as much
          ELSE
            BCOEFF = 0.2 !0.33
          ENDIF

          ! Original StEM with exponential entrainment
          !EntW=exp(-2.*(Wb+Wc*ENT(K,I))*(ZW(k)-ZW(k-1)))
          !Wn2=UPW(K-1,I)**2*EntW + (1.-EntW)*0.5*Wa*B/(Wb+Wc*ENT(K,I))
          ! Original StEM with linear entrainment
          !Wn2=UPW(K-1,I)**2*(1.-EntExp) + EntExp*0.5*Wa*B/(Wb+Wc*ENT(K,I))
          !Wn2=MAX(Wn2,0.0)
          !WA: TEMF form
!          IF (B>0.0 .AND. UPW(K-1,I) < 0.2 ) THEN
          IF (UPW(K-1,I) < 0.2 ) THEN
             Wn = UPW(K-1,I) + (-2. * ENT(K,I) * UPW(K-1,I) + BCOEFF*B / MAX(UPW(K-1,I),0.2)) * MIN(ZW(k)-ZW(k-1), 250.)
          ELSE
             Wn = UPW(K-1,I) + (-2. * ENT(K,I) * UPW(K-1,I) + BCOEFF*B / UPW(K-1,I)) * MIN(ZW(k)-ZW(k-1), 250.)
          ENDIF
          !Do not allow a parcel to accelerate more than 1.25 m/s over 200 m.
          !Add max increase of 2.0 m/s for coarse vertical resolution.
          IF(Wn > UPW(K-1,I) + MIN(1.25*(ZW(k)-ZW(k-1))/200., 2.0) ) THEN
             Wn = UPW(K-1,I) + MIN(1.25*(ZW(k)-ZW(k-1))/200., 2.0)
          ENDIF
          !Add symmetrical max decrease in w
          IF(Wn < UPW(K-1,I) - MIN(1.25*(ZW(k)-ZW(k-1))/200., 2.0) ) THEN
             Wn = UPW(K-1,I) - MIN(1.25*(ZW(k)-ZW(k-1))/200., 2.0)
          ENDIF
          Wn = MIN(MAX(Wn,0.0), 3.0)

          !Check to make sure that the plume made it up at least one level.
          !if it failed, then set nup2=0 and exit the mass-flux portion.
          IF (k==kts+1 .AND. Wn == 0.) THEN
             NUP2=0
             exit
          ENDIF

          IF (debug_mf == 1) THEN
            IF (Wn .GE. 3.0) THEN
              ! surface values
              print *," **** SUSPICIOUSLY LARGE W:"
              print *,' QCn:',QCn,' ENT=',ENT(k,i),' Nup2=',Nup2
              print *,'pblh:',pblh,' Wn:',Wn,' UPW(k-1)=',UPW(K-1,I)
              print *,'K=',k,' B=',B,' dz=',ZW(k)-ZW(k-1)
            ENDIF
          ENDIF

          !Allow strongly forced plumes to overshoot if KE is sufficient
          IF (Wn <= 0.0 .AND. overshoot == 0) THEN
             overshoot = 1
             IF ( THVk-THVkm1 .GT. 0.0 ) THEN
                bvf = SQRT( gtr*(THVk-THVkm1)/dz(k) )
                !vertical Froude number
                Frz = UPW(K-1,I)/(bvf*dz(k))
                !IF ( Frz >= 0.5 ) Wn =  MIN(Frz,1.0)*UPW(K-1,I)
                dzp = dz(k)*MAX(MIN(Frz,1.0),0.0) ! portion of highest layer the plume penetrates
             ENDIF
          ELSE
             dzp = dz(k)
          ENDIF

          !Limit very tall plumes
          Wn=Wn*EXP(-MAX(ZW(k+1)-MIN(pblh+2000.,3500.),0.0)/1000.)

          !JOE- minimize the plume penetratration in stratocu-topped PBL
   !       IF (fltv2 < 0.06) THEN
   !          IF(ZW(k+1) >= pblh-200. .AND. qc(k) > 1e-5 .AND. I > 4) Wn=0.
   !       ENDIF

          !Modify environment variables (representative of the model layer - envm*)
          !following the updraft dynamical detrainment of Asai and Kasahara (1967, JAS).
          !Reminder: w is limited to be non-negative (above)
          aratio   = MIN(UPA(K-1,I)/(1.-UPA(K-1,I)), 0.5) !limit should never get hit
          detturb  = 0.00008
          oow      = -0.060/MAX(1.0,(0.5*(Wn+UPW(K-1,I))))   !coef for dynamical detrainment rate
          detrate  = MIN(MAX(oow*(Wn-UPW(K-1,I))/dz(k), detturb), .0002) ! dynamical detrainment rate (m^-1)
          detrateUV= MIN(MAX(oow*(Wn-UPW(K-1,I))/dz(k), detturb), .0001) ! dynamical detrainment rate (m^-1) 
          envm_thl(k)=envm_thl(k) + (0.5*(thl_ent + UPTHL(K-1,I)) - thl(k))*detrate*aratio*MIN(dzp,dzpmax)
          qv_ent = 0.5*(MAX(qt_ent-qc_ent,0.) + MAX(UPQT(K-1,I)-UPQC(K-1,I),0.))
          envm_sqv(k)=envm_sqv(k) + (qv_ent-QV(K))*detrate*aratio*MIN(dzp,dzpmax)
          IF (UPQC(K-1,I) > 1E-8) THEN
             IF (QC(K) > 1E-6) THEN
                qc_grid = QC(K)
             ELSE
                qc_grid = cldfra_bl1d(k)*qc_bl1d(K)
             ENDIF
             envm_sqc(k)=envm_sqc(k) + MAX(UPA(K-1,I)*0.5*(QCn + UPQC(K-1,I)) - qc_grid, 0.0)*detrate*aratio*MIN(dzp,dzpmax)
          ENDIF
          envm_u(k)  =envm_u(k)   + (0.5*(Un + UPU(K-1,I)) - U(K))*detrateUV*aratio*MIN(dzp,dzpmax)
          envm_v(k)  =envm_v(k)   + (0.5*(Vn + UPV(K-1,I)) - V(K))*detrateUV*aratio*MIN(dzp,dzpmax)

          IF (Wn > 0.) THEN
             !Update plume variables at current k index
             UPW(K,I)=Wn  !sqrt(Wn2)
             UPTHV(K,I)=THVn
             UPTHL(K,I)=THLn
             UPQT(K,I)=QTn
             UPQC(K,I)=QCn
             UPU(K,I)=Un
             UPV(K,I)=Vn
             UPQKE(K,I)=QKEn
             UPQNC(K,I)=QNCn
             UPQNI(K,I)=QNIn
             UPQNWFA(K,I)=QNWFAn
             UPQNIFA(K,I)=QNIFAn
             UPQNBCA(K,I)=QNBCAn
             UPA(K,I)=UPA(K-1,I)
             IF ( mix_chem ) THEN
               do ic = 1,nchem
                 UPCHEM(k,I,ic) = chemn(ic)
               enddo
             ENDIF
             ktop = MAX(ktop,k)
          ELSE
             exit  !exit k-loop
          END IF
       ENDDO
       IF (debug_mf == 1) THEN
          IF (MAXVAL(UPW(:,I)) > 10.0 .OR. MINVAL(UPA(:,I)) < 0.0 .OR. &
              MAXVAL(UPA(:,I)) > Atot .OR. NUP2 > 10) THEN
             ! surface values
             print *,'flq:',flq,' fltv:',fltv2,' Nup2=',Nup2
             print *,'pblh:',pblh,' wstar:',wstar,' ktop=',ktop
             print *,'sigmaW=',sigmaW,' sigmaTH=',sigmaTH,' sigmaQT=',sigmaQT
             ! means
             print *,'u:',u
             print *,'v:',v
             print *,'thl:',thl
             print *,'UPA:',UPA(:,I)
             print *,'UPW:',UPW(:,I)
             print *,'UPTHL:',UPTHL(:,I)
             print *,'UPQT:',UPQT(:,I)
             print *,'ENT:',ENT(:,I)
          ENDIF
       ENDIF
    ENDDO
  ELSE
    !At least one of the conditions was not met for activating the MF scheme.
    NUP2=0.
  END IF !end criteria for mass-flux scheme

  ktop=MIN(ktop,KTE-1)  !  Just to be safe...
  IF (ktop == 0) THEN
     ztop = 0.0
  ELSE
     ztop=zw(ktop)
  ENDIF

  IF(nup2 > 0) THEN

    !Calculate the fluxes for each variable
    !All s_aw* variable are == 0 at k=1
    DO i=1,NUP !NUP2
      IF(I > NUP2) exit
      DO k=KTS,KTE-1
        IF(k > ktop) exit
        rho_int     = (rho(k)*DZ(k+1)+rho(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
        s_aw(k+1)   = s_aw(k+1)    + rho_int*UPA(K,i)*UPW(K,i)*Psig_w
        s_awthl(k+1)= s_awthl(k+1) + rho_int*UPA(K,i)*UPW(K,i)*UPTHL(K,i)*Psig_w
        s_awqt(k+1) = s_awqt(k+1)  + rho_int*UPA(K,i)*UPW(K,i)*UPQT(K,i)*Psig_w
        !to conform to grid mean properties, move qc to qv in grid mean
        !saturated layers, so total water fluxes are preserved but 
        !negative qc fluxes in unsaturated layers is reduced.
        IF (qc(k) > 1e-12 .OR. qc(k+1) > 1e-12) then
          qc_plume = UPQC(K,i)
        ELSE
          qc_plume = 0.0
        ENDIF
        s_awqc(k+1) = s_awqc(k+1)  + rho_int*UPA(K,i)*UPW(K,i)*qc_plume*Psig_w
        IF (momentum_opt > 0) THEN
          s_awu(k+1)  = s_awu(k+1)   + rho_int*UPA(K,i)*UPW(K,i)*UPU(K,i)*Psig_w
          s_awv(k+1)  = s_awv(k+1)   + rho_int*UPA(K,i)*UPW(K,i)*UPV(K,i)*Psig_w
        ENDIF
        IF (tke_opt > 0) THEN
          s_awqke(k+1)= s_awqke(k+1) + rho_int*UPA(K,i)*UPW(K,i)*UPQKE(K,i)*Psig_w
        ENDIF
        s_awqv(k+1) = s_awqt(k+1)  - s_awqc(k+1)
      ENDDO
    ENDDO

    IF ( mix_chem ) THEN
      DO k=KTS,KTE
        IF(k > KTOP) exit
        rho_int     = (rho(k)*DZ(k+1)+rho(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
        DO i=1,NUP !NUP2
          IF(I > NUP2) exit
          do ic = 1,nchem
            s_awchem(k+1,ic) = s_awchem(k+1,ic) + rho_int*UPA(K,i)*UPW(K,i)*UPCHEM(K,i,ic)*Psig_w
          enddo
        ENDDO
      ENDDO
    ENDIF

    IF (scalar_opt > 0) THEN
      DO k=KTS,KTE
        IF(k > KTOP) exit
        rho_int     = (rho(k)*DZ(k+1)+rho(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
        DO I=1,NUP !NUP2
          IF (I > NUP2) exit
          s_awqnc(k+1)= s_awqnc(K+1) + rho_int*UPA(K,i)*UPW(K,i)*UPQNC(K,i)*Psig_w
          s_awqni(k+1)= s_awqni(K+1) + rho_int*UPA(K,i)*UPW(K,i)*UPQNI(K,i)*Psig_w
          s_awqnwfa(k+1)= s_awqnwfa(K+1) + rho_int*UPA(K,i)*UPW(K,i)*UPQNWFA(K,i)*Psig_w
          s_awqnifa(k+1)= s_awqnifa(K+1) + rho_int*UPA(K,i)*UPW(K,i)*UPQNIFA(K,i)*Psig_w
          s_awqnbca(k+1)= s_awqnbca(K+1) + rho_int*UPA(K,i)*UPW(K,i)*UPQNBCA(K,i)*Psig_w
        ENDDO
      ENDDO
    ENDIF

    !Flux limiter: Check ratio of heat flux at top of first model layer
    !and at the surface. Make sure estimated flux out of the top of the
    !layer is < fluxportion*surface_heat_flux
    IF (s_aw(kts+1) /= 0.) THEN
       dzi(kts) = 0.5*(DZ(kts)+DZ(kts+1)) !dz centered at model interface
       flx1   = MAX(s_aw(kts+1)*(TH(kts)-TH(kts+1))/dzi(kts),1.0e-5)
    ELSE
       flx1 = 0.0
       !print*,"ERROR: s_aw(kts+1) == 0, NUP=",NUP," NUP2=",NUP2,&
       !       " superadiabatic=",superadiabatic," KTOP=",KTOP
    ENDIF
    adjustment=1.0
    !Print*,"Flux limiter in MYNN-EDMF, adjustment=",fluxportion*flt/dz(kts)/flx1
    !Print*,"flt/dz=",flt/dz(kts)," flx1=",flx1," s_aw(kts+1)=",s_aw(kts+1)
    IF (flx1 > fluxportion*flt/dz(kts) .AND. flx1>0.0) THEN
       adjustment= fluxportion*flt/dz(kts)/flx1
       s_aw   = s_aw*adjustment
       s_awthl= s_awthl*adjustment
       s_awqt = s_awqt*adjustment
       s_awqc = s_awqc*adjustment
       s_awqv = s_awqv*adjustment
       s_awqnc= s_awqnc*adjustment
       s_awqni= s_awqni*adjustment
       s_awqnwfa= s_awqnwfa*adjustment
       s_awqnifa= s_awqnifa*adjustment
       s_awqnbca= s_awqnbca*adjustment
       IF (momentum_opt > 0) THEN
          s_awu  = s_awu*adjustment
          s_awv  = s_awv*adjustment
       ENDIF
       IF (tke_opt > 0) THEN
          s_awqke= s_awqke*adjustment
       ENDIF
       IF ( mix_chem ) THEN
          s_awchem = s_awchem*adjustment
       ENDIF
       UPA = UPA*adjustment
    ENDIF
    !Print*,"adjustment=",adjustment," fluxportion=",fluxportion," flt=",flt

    !Calculate mean updraft properties for output:
    !all edmf_* variables at k=1 correspond to the interface at top of first model layer
    DO k=KTS,KTE-1
      IF(k > KTOP) exit
      rho_int     = (rho(k)*DZ(k+1)+rho(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
      DO I=1,NUP !NUP2
        IF(I > NUP2) exit
        edmf_a(K)  =edmf_a(K)  +UPA(K,i)
        edmf_w(K)  =edmf_w(K)  +rho_int*UPA(K,i)*UPW(K,i)
        edmf_qt(K) =edmf_qt(K) +rho_int*UPA(K,i)*UPQT(K,i)
        edmf_thl(K)=edmf_thl(K)+rho_int*UPA(K,i)*UPTHL(K,i)
        edmf_ent(K)=edmf_ent(K)+rho_int*UPA(K,i)*ENT(K,i)
        edmf_qc(K) =edmf_qc(K) +rho_int*UPA(K,i)*UPQC(K,i)
      ENDDO

      !Note that only edmf_a is multiplied by Psig_w. This takes care of the
      !scale-awareness of the subsidence below:
      IF (edmf_a(k)>0.) THEN
        edmf_w(k)=edmf_w(k)/edmf_a(k)
        edmf_qt(k)=edmf_qt(k)/edmf_a(k)
        edmf_thl(k)=edmf_thl(k)/edmf_a(k)
        edmf_ent(k)=edmf_ent(k)/edmf_a(k)
        edmf_qc(k)=edmf_qc(k)/edmf_a(k)
        edmf_a(k)=edmf_a(k)*Psig_w

        !FIND MAXIMUM MASS-FLUX IN THE COLUMN:
        IF(edmf_a(k)*edmf_w(k) > maxmf) maxmf = edmf_a(k)*edmf_w(k)
      ENDIF
    ENDDO ! end k

    !smoke/chem
    IF ( mix_chem ) THEN
      DO k=kts,kte-1
        IF(k > KTOP) exit
        rho_int     = (rho(k)*dz(k+1)+rho(k+1)*dz(k))/(dz(k+1)+dz(k))
        DO I=1,NUP !NUP2
          IF(I > NUP2) exit
          do ic = 1,nchem
            edmf_chem(k,ic) = edmf_chem(k,ic) + rho_int*UPA(K,I)*UPCHEM(k,i,ic)
          enddo
        ENDDO

        IF (edmf_a(k)>0.) THEN
          do ic = 1,nchem
            edmf_chem(k,ic) = edmf_chem(k,ic)/edmf_a(k)
          enddo
        ENDIF
      ENDDO ! end k
    ENDIF

    !Calculate the effects environmental subsidence.
    !All envi_*variables are valid at the interfaces, like the edmf_* variables
    IF (env_subs) THEN
       DO k=kts+1,kte-1
          !First, smooth the profiles of w & a, since sharp vertical gradients
          !in plume variables are not likely extended to env variables
          !Note1: w is treated as negative further below
          !Note2: both w & a will be transformed into env variables further below
          envi_w(k) = onethird*(edmf_w(k-1)+edmf_w(k)+edmf_w(k+1))
          envi_a(k) = onethird*(edmf_a(k-1)+edmf_a(k)+edmf_a(k+1))*adjustment
       ENDDO
       !define env variables at k=1 (top of first model layer)
       envi_w(kts) = edmf_w(kts)
       envi_a(kts) = edmf_a(kts)
       !define env variables at k=kte
       envi_w(kte) = 0.0
       envi_a(kte) = edmf_a(kte)
       !define env variables at k=kte+1
       envi_w(kte+1) = 0.0
       envi_a(kte+1) = edmf_a(kte)
       !Add limiter for very long time steps (i.e. dt > 300 s)
       !Note that this is not a robust check - only for violations in
       !   the first model level.
       IF (envi_w(kts) > 0.9*DZ(kts)/dt) THEN
          sublim = 0.9*DZ(kts)/dt/envi_w(kts)
       ELSE
          sublim = 1.0
       ENDIF
       !Transform w & a into env variables
       DO k=kts,kte
          temp=envi_a(k)
          envi_a(k)=1.0-temp
          envi_w(k)=csub*sublim*envi_w(k)*temp/(1.-temp)
       ENDDO
       !calculate tendencies from subsidence and detrainment valid at the middle of
       !each model layer. The lowest model layer uses an assumes w=0 at the surface.
       dzi(kts)    = 0.5*(dz(kts)+dz(kts+1))
       rho_int     = (rho(kts)*dz(kts+1)+rho(kts+1)*dz(kts))/(dz(kts+1)+dz(kts))
       sub_thl(kts)= 0.5*envi_w(kts)*envi_a(kts)*                               &
                     (rho(kts+1)*thl(kts+1)-rho(kts)*thl(kts))/dzi(kts)/rho_int
       sub_sqv(kts)= 0.5*envi_w(kts)*envi_a(kts)*                               &
                     (rho(kts+1)*qv(kts+1)-rho(kts)*qv(kts))/dzi(kts)/rho_int
       DO k=kts+1,kte-1
          dzi(k)    = 0.5*(dz(k)+dz(k+1))
          rho_int   = (rho(k)*dz(k+1)+rho(k+1)*dz(k))/(dz(k+1)+dz(k))
          sub_thl(k)= 0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                      (rho(k+1)*thl(k+1)-rho(k)*thl(k))/dzi(k)/rho_int
          sub_sqv(k)= 0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                      (rho(k+1)*qv(k+1)-rho(k)*qv(k))/dzi(k)/rho_int
       ENDDO

       DO k=KTS,KTE-1
          det_thl(k)=Cdet*(envm_thl(k)-thl(k))*envi_a(k)*Psig_w
          det_sqv(k)=Cdet*(envm_sqv(k)-qv(k))*envi_a(k)*Psig_w
          det_sqc(k)=Cdet*(envm_sqc(k)-qc(k))*envi_a(k)*Psig_w
       ENDDO

       IF (momentum_opt > 0) THEN
         rho_int     = (rho(kts)*dz(kts+1)+rho(kts+1)*dz(kts))/(dz(kts+1)+dz(kts))
         sub_u(kts)=0.5*envi_w(kts)*envi_a(kts)*                               &
                    (rho(kts+1)*u(kts+1)-rho(kts)*u(kts))/dzi(kts)/rho_int
         sub_v(kts)=0.5*envi_w(kts)*envi_a(kts)*                               &
                    (rho(kts+1)*v(kts+1)-rho(kts)*v(kts))/dzi(kts)/rho_int
         DO k=kts+1,kte-1
            rho_int   = (rho(k)*dz(k+1)+rho(k+1)*dz(k))/(dz(k+1)+dz(k))
            sub_u(k)=0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                      (rho(k+1)*u(k+1)-rho(k)*u(k))/dzi(k)/rho_int
            sub_v(k)=0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                      (rho(k+1)*v(k+1)-rho(k)*v(k))/dzi(k)/rho_int
         ENDDO

         DO k=KTS,KTE-1
           det_u(k) = Cdet*(envm_u(k)-u(k))*envi_a(k)*Psig_w
           det_v(k) = Cdet*(envm_v(k)-v(k))*envi_a(k)*Psig_w
         ENDDO
       ENDIF
    ENDIF !end subsidence/env detranment

    !First, compute exner, plume theta, and dz centered at interface
    !Here, k=1 is the top of the first model layer. These values do not 
    !need to be defined at k=kte (unused level).
    DO K=KTS,KTE-1
       exneri(k) = (exner(k)*DZ(k+1)+exner(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
       edmf_th(k)= edmf_thl(k) + xlvcp/exneri(k)*edmf_qc(K)
       dzi(k)    = 0.5*(DZ(k)+DZ(k+1))
    ENDDO

!JOE: ADD CLDFRA_bl1d, qc_bl1d. Note that they have already been defined in
!     mym_condensation. Here, a shallow-cu component is added, but no cumulus
!     clouds can be added at k=1 (start loop at k=2).
    DO K=KTS+1,KTE-2
       IF(k > KTOP) exit
         IF(0.5*(edmf_qc(k)+edmf_qc(k-1))>0.0)THEN
            !interpolate plume quantities to mass levels
            Aup = (edmf_a(k)*dzi(k-1)+edmf_a(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            THp = (edmf_th(k)*dzi(k-1)+edmf_th(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            QTp = (edmf_qt(k)*dzi(k-1)+edmf_qt(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            !convert TH to T
!            t = THp*exner(k)
            !SATURATED VAPOR PRESSURE
            esat = esat_blend(tk(k))
            !SATURATED SPECIFIC HUMIDITY
            qsl=ep_2*esat/max(1.e-7,(p(k)-ep_3*esat)) 

            !condensed liquid in the plume on mass levels
            IF (edmf_qc(k)>0.0 .AND. edmf_qc(k-1)>0.0)THEN
              QCp = (edmf_qc(k)*dzi(k-1)+edmf_qc(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            ELSE
              QCp = MAX(edmf_qc(k),edmf_qc(k-1))
            ENDIF

            !COMPUTE CLDFRA & QC_BL FROM MASS-FLUX SCHEME and recompute vt & vq
            xl = xl_blend(tk(k))                ! obtain blended heat capacity
            qsat_tk = qsat_blend(tk(k),p(k))    ! get saturation water vapor mixing ratio
                                                !   at t and p
            rsl = xl*qsat_tk / (r_v*tk(k)**2)   ! slope of C-C curve at t (abs temp)
                                                ! CB02, Eqn. 4
            cpm = cp + qt(k)*cpv                ! CB02, sec. 2, para. 1
            a   = 1./(1. + xl*rsl/cpm)          ! CB02 variable "a"
            b9  = a*rsl                         ! CB02 variable "b" 

            q2p  = xlvcp/exner(k)
            pt = thl(k) +q2p*QCp*Aup ! potential temp (env + plume)
            bb = b9*tk(k)/pt ! bb is "b9" in BCMT95.  Their "b9" differs from
                           ! "b9" in CB02 by a factor
                           ! of T/theta.  Strictly, b9 above is formulated in
                           ! terms of sat. mixing ratio, but bb in BCMT95 is
                           ! cast in terms of sat. specific humidity.  The
                           ! conversion is neglected here.
            qww   = 1.+0.61*qt(k)
            alpha = 0.61*pt
            beta  = pt*xl/(tk(k)*cp) - 1.61*pt
            !Buoyancy flux terms have been moved to the end of this section...

            !Now calculate convective component of the cloud fraction:
            if (a > 0.0) then
               f = MIN(1.0/a, 4.0)              ! f is vertical profile scaling function (CB2005)
            else
               f = 1.0
            endif

            !CB form:
            !sigq = 3.5E-3 * Aup * 0.5*(edmf_w(k)+edmf_w(k-1)) * f  ! convective component of sigma (CB2005)
            !sigq = SQRT(sigq**2 + sgm(k)**2)    ! combined conv + stratus components
            !Per S.DeRoode 2009?
            !sigq = 4. * Aup * (QTp - qt(k))
            sigq = 10. * Aup * (QTp - qt(k)) 
            !constrain sigq wrt saturation:
            sigq = max(sigq, qsat_tk*0.02 )
            sigq = min(sigq, qsat_tk*0.25 )

            qmq = a * (qt(k) - qsat_tk)           ! saturation deficit/excess;
            Q1  = qmq/sigq                        !   the numerator of Q1

            if ((landsea-1.5).GE.0) then      ! WATER
               !modified form from LES
               !mf_cf = min(max(0.5 + 0.36 * atan(1.20*(Q1+0.2)),0.01),0.6)
               !Original CB
               mf_cf = min(max(0.5 + 0.36 * atan(1.55*Q1),0.01),0.6)
               mf_cf = max(mf_cf, 1.2 * Aup)
               mf_cf = min(mf_cf, 5.0 * Aup)
            else                              ! LAND
               !LES form
               !mf_cf = min(max(0.5 + 0.36 * atan(1.20*(Q1+0.4)),0.01),0.6)
               !Original CB
               mf_cf = min(max(0.5 + 0.36 * atan(1.55*Q1),0.01),0.6)
               mf_cf = max(mf_cf, 1.75 * Aup)
               mf_cf = min(mf_cf, 5.0  * Aup)
            endif

            ! WA TEST 4/15/22 use fit to Aup rather than CB
            !IF (Aup > 0.1) THEN
            !   mf_cf = 2.5 * Aup
            !ELSE
            !   mf_cf = 1.8 * Aup
            !ENDIF

            !IF ( debug_code ) THEN
            !   print*,"In MYNN, StEM edmf"
            !   print*,"  CB: env qt=",qt(k)," qsat=",qsat_tk
            !   print*,"  k=",k," satdef=",QTp - qsat_tk," sgm=",sgm(k)
            !   print*,"  CB: sigq=",sigq," qmq=",qmq," tk=",tk(k)
            !   print*,"  CB: mf_cf=",mf_cf," cldfra_bl=",cldfra_bl1d(k)," edmf_a=",edmf_a(k)
            !ENDIF

            ! Update cloud fractions and specific humidities in grid cells
            ! where the mass-flux scheme is active. The specific humidities
            ! are converted to grid means (not in-cloud quantities).

            if ((landsea-1.5).GE.0) then     ! water
               !don't overwrite stratus CF & qc_bl - degrades marine stratus
               if (cldfra_bl1d(k) < cf_thresh) then
                  if (QCp * Aup > 5e-5) then
                     qc_bl1d(k) = 1.86 * (QCp * Aup) - 2.2e-5
                  else
                     qc_bl1d(k) = 1.18 * (QCp * Aup)
                  endif
                  if (mf_cf .ge. Aup) then
                    qc_bl1d(k) = qc_bl1d(k) / mf_cf
                  endif
                  cldfra_bl1d(k) = mf_cf
                  Ac_mf          = mf_cf
               endif
            else                             ! land
               if (QCp * Aup > 5e-5) then
                  qc_bl1d(k) = 1.86 * (QCp * Aup) - 2.2e-5
               else
                  qc_bl1d(k) = 1.18 * (QCp * Aup)
               endif
               if (mf_cf .ge. Aup) then
                  qc_bl1d(k) = qc_bl1d(k) / mf_cf
               endif
               cldfra_bl1d(k) = mf_cf
               Ac_mf          = mf_cf
            endif

            !Now recalculate the terms for the buoyancy flux for mass-flux clouds:
            !See mym_condensation for details on these formulations.
            !Use Bechtold and Siebesma (1998) piecewise estimation of Fng with 
            !limits ,since they really should be recalculated after all the other changes...:
            !Only overwrite vt & vq in non-stratus condition
            if (cldfra_bl1d(k) < cf_thresh) then
               !if ((landsea-1.5).GE.0) then      ! WATER
                  Q1=max(Q1,-2.25)
               !else
               !   Q1=max(Q1,-2.0)
               !endif

               if (Q1 .ge. 1.0) then
                  Fng = 1.0
               elseif (Q1 .ge. -1.7 .and. Q1 .lt. 1.0) then
                  Fng = EXP(-0.4*(Q1-1.0))
               elseif (Q1 .ge. -2.5 .and. Q1 .lt. -1.7) then
                  Fng = 3.0 + EXP(-3.8*(Q1+1.7))
               else
                  Fng = min(23.9 + EXP(-1.6*(Q1+2.5)), 60.)
               endif

               !link the buoyancy flux function to active clouds only (c*Aup):
               vt(k) = qww   - (1.5*Aup)*beta*bb*Fng - 1.
               vq(k) = alpha + (1.5*Aup)*beta*a*Fng  - tv0
            endif
         endif
      enddo !k-loop

    ENDIF  !end nup2 > 0

    !modify output (negative: dry plume, positive: moist plume)
    IF (ktop > 0) THEN
      maxqc = maxval(edmf_qc(1:ktop)) 
      IF ( maxqc < 1.E-8) maxmf = -1.0*maxmf
    ENDIF

!
! debugging   
!
IF (edmf_w(1) > 4.0) THEN 
! surface values
    print *,'flq:',flq,' fltv:',fltv2
    print *,'pblh:',pblh,' wstar:',wstar
    print *,'sigmaW=',sigmaW,' sigmaTH=',sigmaTH,' sigmaQT=',sigmaQT
! means
!   print *,'u:',u
!   print *,'v:',v  
!   print *,'thl:',thl
!   print *,'thv:',thv
!   print *,'qt:',qt
!   print *,'p:',p
 
! updrafts
! DO I=1,NUP2
!   print *,'up:A',i
!   print *,UPA(:,i)
!   print *,'up:W',i
!   print*,UPW(:,i)
!   print *,'up:thv',i
!   print *,UPTHV(:,i)
!   print *,'up:thl',i 
!   print *,UPTHL(:,i)
!   print *,'up:qt',i
!   print *,UPQT(:,i)
!   print *,'up:tQC',i
!   print *,UPQC(:,i)
!   print *,'up:ent',i
!   print *,ENT(:,i)   
! ENDDO
 
! mean updrafts
   print *,' edmf_a',edmf_a(1:14)
   print *,' edmf_w',edmf_w(1:14)
   print *,' edmf_qt:',edmf_qt(1:14)
   print *,' edmf_thl:',edmf_thl(1:14)
 
ENDIF !END Debugging


#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

END SUBROUTINE DMP_MF

END MODULE
