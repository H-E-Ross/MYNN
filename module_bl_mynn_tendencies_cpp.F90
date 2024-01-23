MODULE mynn_tend 

use module_bl_mynn_common, only: p608,r_d,xlscp,xlvcp
use mynn_functions

use iso_c_binding

IMPLICIT NONE

interface
subroutine moisture_check(kte,delt,dp,exner,qv,qc,qi,th,dqv,dqc,dqi,dth,xlscp,xlvcp)
import c_int, c_double
integer, intent(in)    :: kte
real, intent(in) :: delt,xlscp,xlvcp
real, dimension(kte), intent(in)     :: dp, exner
real, dimension(kte), intent(inout)  :: qv, qc, qi, th
real, dimension(kte), intent(inout)  :: dqv, dqc, dqi, dth
end subroutine moisture_check
end interface

CONTAINS
! ==================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine solves for tendencies of U, V, \f$\theta\f$, qv,
!! qc, and qi
  SUBROUTINE mynn_tendencies(kts,kte,i,    &
       &delt,dz,rho,                       &
       &u,v,th,tk,qv,qc,qi,qnc,qni,        &
       &psfc,p,exner,                      &
       &thl,sqv,sqc,sqi,sqw,               &
       &qnwfa,qnifa,qnbca,ozone,           &
       &ust,flt,flq,flqv,flqc,wspd,        &
       &uoce,voce,                         &
       &tsq,qsq,cov,                       &
       &tcd,qcd,                           &
       &dfm,dfh,dfq,                       &
       &Du,Dv,Dth,Dqv,Dqc,Dqi,Dqnc,Dqni,   &
       &Dqnwfa,Dqnifa,Dqnbca,Dozone,       &
       &diss_heat,                         &
       &s_aw,s_awthl,s_awqt,s_awqv,s_awqc, &
       &s_awu,s_awv,                       &
       &s_awqnc,s_awqni,                   &
       &s_awqnwfa,s_awqnifa,s_awqnbca,     &
       &sd_aw,sd_awthl,sd_awqt,sd_awqv,    &
       &sd_awqc,sd_awu,sd_awv,             &
       &sub_thl,sub_sqv,                   &
       &sub_u,sub_v,                       &
       &det_thl,det_sqv,det_sqc,           &
       &det_u,det_v,                       &
       &FLAG_QC,FLAG_QI,FLAG_QNC,FLAG_QNI, &
       &FLAG_QNWFA,FLAG_QNIFA,FLAG_QNBCA,  &
       &cldfra_bl1d,                       &
       &bl_mynn_cloudmix,                  &
       &bl_mynn_mixqt,                     &
       &bl_mynn_edmf,                      &
       &bl_mynn_edmf_mom,                  &
       &bl_mynn_mixscalars,                &
       &debug_code)

!-------------------------------------------------------------------
    INTEGER, INTENT(in) :: kts,kte,i

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    INTEGER, INTENT(in) :: bl_mynn_cloudmix,bl_mynn_mixqt,&
                           bl_mynn_edmf,bl_mynn_edmf_mom, &
                           bl_mynn_mixscalars
    LOGICAL, INTENT(IN) :: FLAG_QI,FLAG_QNI,FLAG_QC,FLAG_QNC,&
                           FLAG_QNWFA,FLAG_QNIFA,FLAG_QNBCA,&
                           debug_code

! thl - liquid water potential temperature
! qw - total water
! dfm,dfh,dfq - diffusivities i.e., dfh(k) = elq*sh(k) / dzk
! flt - surface flux of thl
! flq - surface flux of qw

! mass-flux plumes
    REAL, DIMENSION(kts:kte+1), INTENT(in) :: s_aw,s_awthl,s_awqt,&
         &s_awqnc,s_awqni,s_awqv,s_awqc,s_awu,s_awv,              &
         &s_awqnwfa,s_awqnifa,s_awqnbca,                          &
         &sd_aw,sd_awthl,sd_awqt,sd_awqv,sd_awqc,sd_awu,sd_awv
! tendencies from mass-flux environmental subsidence and detrainment
    REAL, DIMENSION(kts:kte), INTENT(in) :: sub_thl,sub_sqv,  &
         &sub_u,sub_v,det_thl,det_sqv,det_sqc,det_u,det_v
    REAL, DIMENSION(kts:kte), INTENT(in) :: u,v,th,tk,qv,qc,qi,qni,qnc,&
         &rho,p,exner,dfq,dz,tsq,qsq,cov,tcd,qcd,cldfra_bl1d,diss_heat
    REAL, DIMENSION(kts:kte), INTENT(inout) :: thl,sqw,sqv,sqc,sqi,&
         &qnwfa,qnifa,qnbca,ozone,dfm,dfh
    REAL, DIMENSION(kts:kte), INTENT(inout) :: du,dv,dth,dqv,dqc,dqi,&
         &dqni,dqnc,dqnwfa,dqnifa,dqnbca,dozone
    REAL, INTENT(IN) :: delt,ust,flt,flq,flqv,flqc,wspd,uoce,voce,&
         &psfc
    !debugging
    REAL ::wsp,wsp2,tk2,th2
    LOGICAL :: problem
    integer :: kproblem

!    REAL, INTENT(IN) :: gradu_top,gradv_top,gradth_top,gradqv_top

!local vars

    REAL, DIMENSION(kts:kte) :: dtz,dfhc,dfmc,delp
    REAL, DIMENSION(kts:kte) :: sqv2,sqc2,sqi2,sqw2,qni2,qnc2, & !AFTER MIXING
                                qnwfa2,qnifa2,qnbca2,ozone2
    REAL, DIMENSION(kts:kte) :: zfac,plumeKh,rhoinv
    REAL, DIMENSION(kts:kte) :: a,b,c,d,x
    REAL, DIMENSION(kts:kte+1) :: rhoz, & !rho on model interface
          &         khdz, kmdz
    REAL :: rhs,gfluxm,gfluxp,dztop,maxdfh,mindfh,maxcf,maxKh,zw
    REAL :: t,esat,qsl,onoff,kh,km,dzk,rhosfc
    REAL :: ustdrag,ustdiff,qvflux
    REAL :: th_new,portion_qc,portion_qi,condensate,qsat
    INTEGER :: k,kk

    !Activate nonlocal mixing from the mass-flux scheme for
    !number concentrations and aerosols (0.0 = no; 1.0 = yes)
    REAL, PARAMETER :: nonloc = 1.0

    dztop=.5*(dz(kte)+dz(kte-1))

    ! REGULATE THE MOMENTUM MIXING FROM THE MASS-FLUX SCHEME (on or off)
    ! Note that s_awu and s_awv already come in as 0.0 if bl_mynn_edmf_mom == 0, so
    ! we only need to zero-out the MF term
    IF (bl_mynn_edmf_mom == 0) THEN
       onoff=0.0
    ELSE
       onoff=1.0
    ENDIF

    !Prepare "constants" for diffusion equation.
    !khdz = rho*Kh/dz = rho*dfh
    rhosfc     = psfc/(R_d*(tk(kts)+p608*qv(kts)))
    dtz(kts)   =delt/dz(kts)
    rhoz(kts)  =rho(kts)
    rhoinv(kts)=1./rho(kts)
    khdz(kts)  =rhoz(kts)*dfh(kts)
    kmdz(kts)  =rhoz(kts)*dfm(kts)
    delp(kts)  = psfc - (p(kts+1)*dz(kts) + p(kts)*dz(kts+1))/(dz(kts)+dz(kts+1))
    DO k=kts+1,kte
       dtz(k)   =delt/dz(k)
       rhoz(k)  =(rho(k)*dz(k-1) + rho(k-1)*dz(k))/(dz(k-1)+dz(k))
       rhoz(k)  =  MAX(rhoz(k),1E-4)
       rhoinv(k)=1./MAX(rho(k),1E-4)
       dzk      = 0.5  *( dz(k)+dz(k-1) )
       khdz(k)  = rhoz(k)*dfh(k)
       kmdz(k)  = rhoz(k)*dfm(k)
    ENDDO
    DO k=kts+1,kte-1
       delp(k)  = (p(k)*dz(k-1) + p(k-1)*dz(k))/(dz(k)+dz(k-1)) - &
                  (p(k+1)*dz(k) + p(k)*dz(k+1))/(dz(k)+dz(k+1))
    ENDDO
    delp(kte)  =delp(kte-1)
    rhoz(kte+1)=rhoz(kte)
    khdz(kte+1)=rhoz(kte+1)*dfh(kte)
    kmdz(kte+1)=rhoz(kte+1)*dfm(kte)

    !stability criteria for mf
    DO k=kts+1,kte-1
       khdz(k) = MAX(khdz(k),  0.5*s_aw(k))
       khdz(k) = MAX(khdz(k), -0.5*(s_aw(k)-s_aw(k+1)))
       kmdz(k) = MAX(kmdz(k),  0.5*s_aw(k))
       kmdz(k) = MAX(kmdz(k), -0.5*(s_aw(k)-s_aw(k+1)))
    ENDDO

    ustdrag = MIN(ust*ust,0.99)/wspd  ! limit at ~ 20 m/s
    ustdiff = MIN(ust*ust,0.01)/wspd  ! limit at ~ 2 m/s
    dth(kts:kte) = 0.0  ! must initialize for moisture_check routine

!!============================================
!! u
!!============================================

    k=kts

!original approach (drag in b-vector):
!    a(1)=0.
!    b(1)=1. + dtz(k)*(dfm(k+1)+ust**2/wspd) - 0.5*dtz(k)*s_aw(k+1)*onoff
!    c(1)=-dtz(k)*dfm(k+1) - 0.5*dtz(k)*s_aw(k+1)*onoff
!    d(1)=u(k) + dtz(k)*uoce*ust**2/wspd - dtz(k)*s_awu(k+1)*onoff + &
!         sub_u(k)*delt + det_u(k)*delt

!rho-weighted (drag in b-vector):
    a(k)=  -dtz(k)*kmdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(kmdz(k+1)+rhosfc*ust**2/wspd)*rhoinv(k) &
           & - 0.5*dtz(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)*onoff
    c(k)=  -dtz(k)*kmdz(k+1)*rhoinv(k) &
           & - 0.5*dtz(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)*onoff
    d(k)=u(k)  + dtz(k)*uoce*ust**2/wspd - dtz(k)*s_awu(k+1)*onoff - &
       & dtz(k)*rhoinv(k)*sd_awu(k+1)*onoff + sub_u(k)*delt + det_u(k)*delt

!rho-weighted with drag term moved out of b-array
!    a(k)=  -dtz(k)*kmdz(k)*rhoinv(k)
!    b(k)=1.+dtz(k)*(kmdz(k+1))*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)*onoff
!    c(k)=  -dtz(k)*kmdz(k+1)*rhoinv(k)   - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)*onoff
!    d(k)=u(k)*(1.-ust**2/wspd*dtz(k)*rhosfc/rho(k)) + dtz(k)*uoce*ust**2/wspd - &
!    !!!d(k)=u(k)*(1.-ust**2/wspd*dtz(k)) + dtz(k)*uoce*ust**2/wspd - &
!      &  dtz(k)*rhoinv(k)*s_awu(k+1)*onoff - dtz(k)*rhoinv(k)*sd_awu(k+1)*onoff + sub_u(k)*delt + det_u(k)*delt

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*kmdz(k)*rhoinv(k)     + 0.5*dtz(k)*rhoinv(k)*s_aw(k)*onoff + 0.5*dtz(k)*rhoinv(k)*sd_aw(k)*onoff 
       b(k)=1.+dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k) + &
           &    0.5*dtz(k)*rhoinv(k)*(s_aw(k)-s_aw(k+1))*onoff + 0.5*dtz(k)*rhoinv(k)*(sd_aw(k)-sd_aw(k+1))*onoff
       c(k)=  -dtz(k)*kmdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)*onoff
       d(k)=u(k) + dtz(k)*rhoinv(k)*(s_awu(k)-s_awu(k+1))*onoff + dtz(k)*rhoinv(k)*(sd_awu(k)-sd_awu(k+1))*onoff + &
           &    sub_u(k)*delt + det_u(k)*delt
    ENDDO

!! no flux at the top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=0.

!! specified gradient at the top 
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=gradu_top*dztop

!! prescribed value
    a(kte)=0
    b(kte)=1.
    c(kte)=0.
    d(kte)=u(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
!       du(k)=(d(k-kts+1)-u(k))/delt
       du(k)=(x(k)-u(k))/delt
    ENDDO

!!============================================
!! v
!!============================================

    k=kts

!original approach (drag in b-vector):
!    a(1)=0.
!    b(1)=1. + dtz(k)*(dfm(k+1)+ust**2/wspd) - 0.5*dtz(k)*s_aw(k+1)*onoff
!    c(1)=   - dtz(k)*dfm(k+1)               - 0.5*dtz(k)*s_aw(k+1)*onoff
!    d(1)=v(k) + dtz(k)*voce*ust**2/wspd - dtz(k)*s_awv(k+1)*onoff + &
!          sub_v(k)*delt + det_v(k)*delt

!rho-weighted (drag in b-vector):
    a(k)=  -dtz(k)*kmdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(kmdz(k+1) + rhosfc*ust**2/wspd)*rhoinv(k) &
        &  - 0.5*dtz(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)*onoff
    c(k)=  -dtz(k)*kmdz(k+1)*rhoinv(k) - 0.5*dtz(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)*onoff
    d(k)=v(k)  + dtz(k)*voce*ust**2/wspd - dtz(k)*s_awv(k+1)*onoff - dtz(k)*rhoinv(k)*sd_awv(k+1)*onoff + &
       & sub_v(k)*delt + det_v(k)*delt

!rho-weighted with drag	term moved out of b-array
!    a(k)=  -dtz(k)*kmdz(k)*rhoinv(k)
!    b(k)=1.+dtz(k)*(kmdz(k+1))*rhoinv(k)  - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)*onoff
!    c(k)=  -dtz(k)*kmdz(k+1)*rhoinv(k)    - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)*onoff
!    d(k)=v(k)*(1.-ust**2/wspd*dtz(k)*rhosfc/rho(k)) + dtz(k)*voce*ust**2/wspd - &
!    !!!d(k)=v(k)*(1.-ust**2/wspd*dtz(k)) + dtz(k)*voce*ust**2/wspd - &
!      &  dtz(k)*rhoinv(k)*s_awv(k+1)*onoff - dtz(k)*rhoinv(k)*sd_awv(k+1)*onoff + sub_v(k)*delt + det_v(k)*delt

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*kmdz(k)*rhoinv(k)   + 0.5*dtz(k)*rhoinv(k)*s_aw(k)*onoff + 0.5*dtz(k)*rhoinv(k)*sd_aw(k)*onoff
       b(k)=1.+dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k) + &
           &    0.5*dtz(k)*rhoinv(k)*(s_aw(k)-s_aw(k+1))*onoff + 0.5*dtz(k)*rhoinv(k)*(sd_aw(k)-sd_aw(k+1))*onoff
       c(k)=  -dtz(k)*kmdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*onoff - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)*onoff 
       d(k)=v(k) + dtz(k)*rhoinv(k)*(s_awv(k)-s_awv(k+1))*onoff + dtz(k)*rhoinv(k)*(sd_awv(k)-sd_awv(k+1))*onoff + &
           &    sub_v(k)*delt + det_v(k)*delt
    ENDDO

!! no flux at the top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=0.

!! specified gradient at the top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=gradv_top*dztop

!! prescribed value
    a(kte)=0
    b(kte)=1.
    c(kte)=0.
    d(kte)=v(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
!       dv(k)=(d(k-kts+1)-v(k))/delt
       dv(k)=(x(k)-v(k))/delt
    ENDDO

!!============================================
!! thl tendency
!!============================================
    k=kts

!    a(k)=0.
!    b(k)=1.+dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1)
!    c(k)=  -dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1)
!    d(k)=thl(k) + dtz(k)*flt + tcd(k)*delt &
!        & -dtz(k)*s_awthl(kts+1) + diss_heat(k)*delt + &
!        & sub_thl(k)*delt + det_thl(k)*delt
!
!    DO k=kts+1,kte-1
!       a(k)=  -dtz(k)*dfh(k)            + 0.5*dtz(k)*s_aw(k)
!       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) + 0.5*dtz(k)*(s_aw(k)-s_aw(k+1))
!       c(k)=  -dtz(k)*dfh(k+1)          - 0.5*dtz(k)*s_aw(k+1)
!       d(k)=thl(k) + tcd(k)*delt + dtz(k)*(s_awthl(k)-s_awthl(k+1)) &
!           &       + diss_heat(k)*delt + &
!           &         sub_thl(k)*delt + det_thl(k)*delt
!    ENDDO

!rho-weighted: rhosfc*X*rhoinv(k)
    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1) - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)           - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1) - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)
    d(k)=thl(k)  + dtz(k)*rhosfc*flt*rhoinv(k) + tcd(k)*delt &
       & - dtz(k)*rhoinv(k)*s_awthl(k+1) -dtz(k)*rhoinv(k)*sd_awthl(k+1) + &
       & diss_heat(k)*delt + sub_thl(k)*delt + det_thl(k)*delt

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)     + 0.5*dtz(k)*rhoinv(k)*s_aw(k) + 0.5*dtz(k)*rhoinv(k)*sd_aw(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k) + &
           &   0.5*dtz(k)*rhoinv(k)*(s_aw(k)-s_aw(k+1)) + 0.5*dtz(k)*rhoinv(k)*(sd_aw(k)-sd_aw(k+1))
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1) - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)
       d(k)=thl(k) + tcd(k)*delt + &
          & dtz(k)*rhoinv(k)*(s_awthl(k)-s_awthl(k+1)) + dtz(k)*rhoinv(k)*(sd_awthl(k)-sd_awthl(k+1)) + &
          &       diss_heat(k)*delt + &
          &       sub_thl(k)*delt + det_thl(k)*delt
    ENDDO

!! no flux at the top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=0.

!! specified gradient at the top
!assume gradthl_top=gradth_top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=gradth_top*dztop

!! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=thl(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !thl(k)=d(k-kts+1)
       thl(k)=x(k)
    ENDDO

IF (bl_mynn_mixqt > 0) THEN
 !============================================
 ! MIX total water (sqw = sqc + sqv + sqi)
 ! NOTE: no total water tendency is output; instead, we must calculate
 !       the saturation specific humidity and then 
 !       subtract out the moisture excess (sqc & sqi)
 !============================================

    k=kts

!    a(k)=0.
!    b(k)=1.+dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1)
!    c(k)=  -dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1)
!    !rhs= qcd(k) !+ (gfluxp - gfluxm)/dz(k)&
!    d(k)=sqw(k) + dtz(k)*flq + qcd(k)*delt - dtz(k)*s_awqt(k+1)
!
!    DO k=kts+1,kte-1
!       a(k)=  -dtz(k)*dfh(k)            + 0.5*dtz(k)*s_aw(k)
!       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) + 0.5*dtz(k)*(s_aw(k)-s_aw(k+1))
!       c(k)=  -dtz(k)*dfh(k+1)          - 0.5*dtz(k)*s_aw(k+1)
!       d(k)=sqw(k) + qcd(k)*delt + dtz(k)*(s_awqt(k)-s_awqt(k+1))
!    ENDDO

!rho-weighted:
    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1) - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)           - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1) - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)
    d(k)=sqw(k)  + dtz(k)*rhosfc*flq*rhoinv(k) + qcd(k)*delt - dtz(k)*rhoinv(k)*s_awqt(k+1) - dtz(k)*rhoinv(k)*sd_awqt(k+1)

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)     + 0.5*dtz(k)*rhoinv(k)*s_aw(k) + 0.5*dtz(k)*rhoinv(k)*sd_aw(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k) + &
           &    0.5*dtz(k)*rhoinv(k)*(s_aw(k)-s_aw(k+1)) + 0.5*dtz(k)*rhoinv(k)*(sd_aw(k)-sd_aw(k+1))
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1) - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)
       d(k)=sqw(k) + qcd(k)*delt + dtz(k)*rhoinv(k)*(s_awqt(k)-s_awqt(k+1)) + dtz(k)*rhoinv(k)*(sd_awqt(k)-sd_awqt(k+1))
    ENDDO

!! no flux at the top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=0.
!! specified gradient at the top
!assume gradqw_top=gradqv_top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=gradqv_top*dztop
!! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=sqw(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,sqw2)
!    CALL tridiag3(kte,a,b,c,d,sqw2)

!    DO k=kts,kte
!       sqw2(k)=d(k-kts+1)
!    ENDDO
ELSE
    sqw2=sqw
ENDIF

IF (bl_mynn_mixqt == 0) THEN
!============================================
! cloud water ( sqc ). If mixing total water (bl_mynn_mixqt > 0),
! then sqc will be backed out of saturation check (below).
!============================================
  IF (bl_mynn_cloudmix > 0 .AND. FLAG_QC) THEN

    k=kts

!    a(k)=0.
!    b(k)=1.+dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1)
!    c(k)=  -dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1)
!    d(k)=sqc(k) + dtz(k)*flqc + qcd(k)*delt - &
!         dtz(k)*s_awqc(k+1)  + det_sqc(k)*delt
!
!    DO k=kts+1,kte-1
!       a(k)=  -dtz(k)*dfh(k)            + 0.5*dtz(k)*s_aw(k)
!       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) + 0.5*dtz(k)*(s_aw(k)-s_aw(k+1))
!       c(k)=  -dtz(k)*dfh(k+1)          - 0.5*dtz(k)*s_aw(k+1)
!       d(k)=sqc(k) + qcd(k)*delt + dtz(k)*(s_awqc(k)-s_awqc(k+1)) + &
!            det_sqc(k)*delt
!    ENDDO

!rho-weighted:
    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1) - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)           - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1) - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)
    d(k)=sqc(k)  + dtz(k)*rhosfc*flqc*rhoinv(k) + qcd(k)*delt &
       &  - dtz(k)*rhoinv(k)*s_awqc(k+1) - dtz(k)*rhoinv(k)*sd_awqc(k+1) + &
       &  det_sqc(k)*delt

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)     + 0.5*dtz(k)*rhoinv(k)*s_aw(k) + 0.5*dtz(k)*rhoinv(k)*sd_aw(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k) + &
           &    0.5*dtz(k)*rhoinv(k)*(s_aw(k)-s_aw(k+1)) + 0.5*dtz(k)*rhoinv(k)*(sd_aw(k)-sd_aw(k+1))
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1) - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)
       d(k)=sqc(k) + qcd(k)*delt + dtz(k)*rhoinv(k)*(s_awqc(k)-s_awqc(k+1)) + dtz(k)*rhoinv(k)*(sd_awqc(k)-sd_awqc(k+1)) + &
          & det_sqc(k)*delt
    ENDDO

! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=sqc(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,sqc2)
!    CALL tridiag3(kte,a,b,c,d,sqc2)

!    DO k=kts,kte
!       sqc2(k)=d(k-kts+1)
!    ENDDO
  ELSE
    !If not mixing clouds, set "updated" array equal to original array
    sqc2=sqc
  ENDIF
ENDIF

IF (bl_mynn_mixqt == 0) THEN
  !============================================
  ! MIX WATER VAPOR ONLY ( sqv ). If mixing total water (bl_mynn_mixqt > 0),
  ! then sqv will be backed out of saturation check (below).
  !============================================

    k=kts

!    a(k)=0.
!    b(k)=1.+dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1)
!    c(k)=  -dtz(k)*dfh(k+1) - 0.5*dtz(k)*s_aw(k+1)
!    d(k)=sqv(k) + dtz(k)*flqv + qcd(k)*delt - dtz(k)*s_awqv(k+1) + &
!       & sub_sqv(k)*delt + det_sqv(k)*delt
!
!    DO k=kts+1,kte-1
!       a(k)=  -dtz(k)*dfh(k)            + 0.5*dtz(k)*s_aw(k)
!       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1)) + 0.5*dtz(k)*(s_aw(k)-s_aw(k+1))
!       c(k)=  -dtz(k)*dfh(k+1)          - 0.5*dtz(k)*s_aw(k+1)
!       d(k)=sqv(k) + qcd(k)*delt + dtz(k)*(s_awqv(k)-s_awqv(k+1)) + &
!          & sub_sqv(k)*delt + det_sqv(k)*delt
!    ENDDO

    !limit unreasonably large negative fluxes:
    qvflux = flqv
    if (qvflux < 0.0) then
       !do not allow specified surface flux to reduce qv below 1e-8 kg/kg
       qvflux = max(qvflux, (min(0.9*sqv(kts) - 1e-8, 0.0)/dtz(kts)))
    endif

!rho-weighted:  rhosfc*X*rhoinv(k)
    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1) - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)           - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1) - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)
    d(k)=sqv(k)  + dtz(k)*rhosfc*qvflux*rhoinv(k) + qcd(k)*delt &
       &  - dtz(k)*rhoinv(k)*s_awqv(k+1) - dtz(k)*rhoinv(k)*sd_awqv(k+1) + &
       & sub_sqv(k)*delt + det_sqv(k)*delt

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)     + 0.5*dtz(k)*rhoinv(k)*s_aw(k) + 0.5*dtz(k)*rhoinv(k)*sd_aw(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k) + &
           &    0.5*dtz(k)*rhoinv(k)*(s_aw(k)-s_aw(k+1)) + 0.5*dtz(k)*rhoinv(k)*(sd_aw(k)-sd_aw(k+1))
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1) - 0.5*dtz(k)*rhoinv(k)*sd_aw(k+1)
       d(k)=sqv(k) + qcd(k)*delt + dtz(k)*rhoinv(k)*(s_awqv(k)-s_awqv(k+1)) + dtz(k)*rhoinv(k)*(sd_awqv(k)-sd_awqv(k+1)) + &
          & sub_sqv(k)*delt + det_sqv(k)*delt
    ENDDO

! no flux at the top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=0.

! specified gradient at the top
! assume gradqw_top=gradqv_top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=gradqv_top*dztop

! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=sqv(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,sqv2)
!    CALL tridiag3(kte,a,b,c,d,sqv2)

!    DO k=kts,kte
!       sqv2(k)=d(k-kts+1)
!    ENDDO
ELSE
    sqv2=sqv
ENDIF

!============================================
! MIX CLOUD ICE ( sqi )                      
!============================================
IF (bl_mynn_cloudmix > 0 .AND. FLAG_QI) THEN

    k=kts

!    a(k)=0.
!    b(k)=1.+dtz(k)*dfh(k+1)
!    c(k)=  -dtz(k)*dfh(k+1)
!    d(k)=sqi(k) !+ qcd(k)*delt !should we have qcd for ice?
!
!    DO k=kts+1,kte-1
!       a(k)=  -dtz(k)*dfh(k)
!       b(k)=1.+dtz(k)*(dfh(k)+dfh(k+1))
!       c(k)=  -dtz(k)*dfh(k+1)
!       d(k)=sqi(k) !+ qcd(k)*delt
!    ENDDO

!rho-weighted:
    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)
    d(k)=sqi(k)

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)
       d(k)=sqi(k)
    ENDDO

!! no flux at the top
!    a(kte)=-1.       
!    b(kte)=1.        
!    c(kte)=0.        
!    d(kte)=0.        

!! specified gradient at the top
!assume gradqw_top=gradqv_top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=gradqv_top*dztop

!! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=sqi(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,sqi2)
!    CALL tridiag3(kte,a,b,c,d,sqi2)

!    DO k=kts,kte
!       sqi2(k)=d(k-kts+1)
!    ENDDO
ELSE
   sqi2=sqi
ENDIF

!!============================================
!! cloud ice number concentration (qni)
!!============================================
IF (bl_mynn_cloudmix > 0 .AND. FLAG_QNI .AND. &
      bl_mynn_mixscalars > 0) THEN

    k=kts

    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*nonloc
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*nonloc
    d(k)=qni(k)  - dtz(k)*rhoinv(k)*s_awqni(k+1)*nonloc

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)     + 0.5*dtz(k)*rhoinv(k)*s_aw(k)*nonloc
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k) + &
           &    0.5*dtz(k)*rhoinv(k)*(s_aw(k)-s_aw(k+1))*nonloc
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*nonloc
       d(k)=qni(k) + dtz(k)*rhoinv(k)*(s_awqni(k)-s_awqni(k+1))*nonloc
    ENDDO

!! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=qni(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !qni2(k)=d(k-kts+1)
       qni2(k)=x(k)
    ENDDO

ELSE
    qni2=qni
ENDIF

!!============================================
!! cloud water number concentration (qnc)     
!! include non-local transport                
!!============================================
  IF (bl_mynn_cloudmix > 0 .AND. FLAG_QNC .AND. &
      bl_mynn_mixscalars > 0) THEN

    k=kts

    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*nonloc
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)           - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*nonloc
    d(k)=qnc(k)  - dtz(k)*rhoinv(k)*s_awqnc(k+1)*nonloc

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)     + 0.5*dtz(k)*rhoinv(k)*s_aw(k)*nonloc
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k) + &
           &    0.5*dtz(k)*rhoinv(k)*(s_aw(k)-s_aw(k+1))*nonloc
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*nonloc
       d(k)=qnc(k) + dtz(k)*rhoinv(k)*(s_awqnc(k)-s_awqnc(k+1))*nonloc
    ENDDO

!! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=qnc(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !qnc2(k)=d(k-kts+1)
       qnc2(k)=x(k)
    ENDDO

ELSE
    qnc2=qnc
ENDIF

!============================================
! Water-friendly aerosols ( qnwfa ).
!============================================
IF (bl_mynn_cloudmix > 0 .AND. FLAG_QNWFA .AND. &
      bl_mynn_mixscalars > 0) THEN

    k=kts

    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k) + khdz(k+1))*rhoinv(k) - &
           &    0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*nonloc
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*nonloc
    d(k)=qnwfa(k)  - dtz(k)*rhoinv(k)*s_awqnwfa(k+1)*nonloc

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)     + 0.5*dtz(k)*rhoinv(k)*s_aw(k)*nonloc
       b(k)=1.+dtz(k)*(khdz(k) + khdz(k+1))*rhoinv(k) + &
           &    0.5*dtz(k)*rhoinv(k)*(s_aw(k)-s_aw(k+1))*nonloc
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*nonloc
       d(k)=qnwfa(k) + dtz(k)*rhoinv(k)*(s_awqnwfa(k)-s_awqnwfa(k+1))*nonloc
    ENDDO

! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=qnwfa(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !qnwfa2(k)=d(k)
       qnwfa2(k)=x(k)
    ENDDO

ELSE
    !If not mixing aerosols, set "updated" array equal to original array
    qnwfa2=qnwfa
ENDIF

!============================================
! Ice-friendly aerosols ( qnifa ).
!============================================
IF (bl_mynn_cloudmix > 0 .AND. FLAG_QNIFA .AND. &
      bl_mynn_mixscalars > 0) THEN

   k=kts

    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k) + khdz(k+1))*rhoinv(k) - &
           &    0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*nonloc
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*nonloc
    d(k)=qnifa(k)  - dtz(k)*rhoinv(k)*s_awqnifa(k+1)*nonloc

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)     + 0.5*dtz(k)*rhoinv(k)*s_aw(k)*nonloc
       b(k)=1.+dtz(k)*(khdz(k) + khdz(k+1))*rhoinv(k) + &
           &    0.5*dtz(k)*rhoinv(k)*(s_aw(k)-s_aw(k+1))*nonloc
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*nonloc
       d(k)=qnifa(k) + dtz(k)*rhoinv(k)*(s_awqnifa(k)-s_awqnifa(k+1))*nonloc
    ENDDO

! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=qnifa(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !qnifa2(k)=d(k-kts+1)
       qnifa2(k)=x(k)
    ENDDO

ELSE
    !If not mixing aerosols, set "updated" array equal to original array
    qnifa2=qnifa
ENDIF

!============================================
! Black-carbon aerosols ( qnbca ).           
!============================================
IF (bl_mynn_cloudmix > 0 .AND. FLAG_QNBCA .AND. &
      bl_mynn_mixscalars > 0) THEN

   k=kts

    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k) + khdz(k+1))*rhoinv(k) - &
           &    0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*nonloc
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*nonloc
    d(k)=qnbca(k)  - dtz(k)*rhoinv(k)*s_awqnbca(k+1)*nonloc

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)     + 0.5*dtz(k)*rhoinv(k)*s_aw(k)*nonloc
       b(k)=1.+dtz(k)*(khdz(k) + khdz(k+1))*rhoinv(k) + &
           &    0.5*dtz(k)*rhoinv(k)*(s_aw(k)-s_aw(k+1))*nonloc
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*nonloc
       d(k)=qnbca(k) + dtz(k)*rhoinv(k)*(s_awqnbca(k)-s_awqnbca(k+1))*nonloc
    ENDDO

! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=qnbca(kte)

!    CALL tridiag(kte,a,b,c,d)
!    CALL tridiag2(kte,a,b,c,d,x)
    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !qnbca2(k)=d(k-kts+1)
       qnbca2(k)=x(k)
    ENDDO

ELSE
    !If not mixing aerosols, set "updated" array equal to original array
    qnbca2=qnbca
ENDIF

!============================================
! Ozone - local mixing only
!============================================

    k=kts

!rho-weighted:
    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)
    d(k)=ozone(k)

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)
       d(k)=ozone(k)
    ENDDO

! prescribed value                                                                                                           
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=ozone(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !ozone2(k)=d(k-kts+1)
       dozone(k)=(x(k)-ozone(k))/delt
    ENDDO

!!============================================
!! Compute tendencies and convert to mixing ratios for WRF.
!! Note that the momentum tendencies are calculated above.
!!============================================

   IF (bl_mynn_mixqt > 0) THEN 
      DO k=kts,kte
         !compute updated theta using updated thl and old condensate
         th_new = thl(k) + xlvcp/exner(k)*sqc(k) &
           &             + xlscp/exner(k)*sqi(k)

         t  = th_new*exner(k)
         qsat = qsat_blend(t,p(k)) 
         !SATURATED VAPOR PRESSURE
         !esat=esat_blend(t)
         !SATURATED SPECIFIC HUMIDITY
         !qsl=ep_2*esat/(p(k)-ep_3*esat)
         !qsl=ep_2*esat/max(1.e-4,(p(k)-ep_3*esat))

         IF (sqc(k) > 0.0 .or. sqi(k) > 0.0) THEN !initially saturated
            sqv2(k) = MIN(sqw2(k),qsat)
            portion_qc = sqc(k)/(sqc(k) + sqi(k))
            portion_qi = sqi(k)/(sqc(k) + sqi(k))
            condensate = MAX(sqw2(k) - qsat, 0.0)
            sqc2(k) = condensate*portion_qc
            sqi2(k) = condensate*portion_qi
         ELSE                     ! initially unsaturated -----
            sqv2(k) = sqw2(k)     ! let microphys decide what to do
            sqi2(k) = 0.0         ! if sqw2 > qsat 
            sqc2(k) = 0.0
         ENDIF
         !dqv(k) = (sqv2(k) - sqv(k))/delt
         !dqc(k) = (sqc2(k) - sqc(k))/delt
         !dqi(k) = (sqi2(k) - sqi(k))/delt
      ENDDO
   ENDIF


    !=====================
    ! WATER VAPOR TENDENCY
    !=====================
    DO k=kts,kte
       Dqv(k)=(sqv2(k)/(1.-sqv2(k)) - qv(k))/delt
       !if (sqv2(k) < 0.0)print*,"neg qv:",sqv2(k),k
    ENDDO

    IF (bl_mynn_cloudmix > 0) THEN
      !=====================
      ! CLOUD WATER TENDENCY
      !=====================
      !print*,"FLAG_QC:",FLAG_QC
      IF (FLAG_QC) THEN
         DO k=kts,kte
            Dqc(k)=(sqc2(k)/(1.-sqv2(k)) - qc(k))/delt
            !if (sqc2(k) < 0.0)print*,"neg qc:",sqc2(k),k
         ENDDO
      ELSE
         DO k=kts,kte
           Dqc(k) = 0.
         ENDDO
      ENDIF

      !===================
      ! CLOUD WATER NUM CONC TENDENCY
      !===================
      IF (FLAG_QNC .AND. bl_mynn_mixscalars > 0) THEN
         DO k=kts,kte
           Dqnc(k) = (qnc2(k)-qnc(k))/delt
           !IF(Dqnc(k)*delt + qnc(k) < 0.)Dqnc(k)=-qnc(k)/delt
         ENDDO 
      ELSE
         DO k=kts,kte
           Dqnc(k) = 0.
         ENDDO
      ENDIF

      !===================
      ! CLOUD ICE TENDENCY
      !===================
      IF (FLAG_QI) THEN
         DO k=kts,kte
           Dqi(k)=(sqi2(k)/(1.-sqv2(k)) - qi(k))/delt
           !if (sqi2(k) < 0.0)print*,"neg qi:",sqi2(k),k
         ENDDO
      ELSE
         DO k=kts,kte
           Dqi(k) = 0.
         ENDDO
      ENDIF

      !===================
      ! CLOUD ICE NUM CONC TENDENCY
      !===================
      IF (FLAG_QNI .AND. bl_mynn_mixscalars > 0) THEN
         DO k=kts,kte
           Dqni(k)=(qni2(k)-qni(k))/delt
           !IF(Dqni(k)*delt + qni(k) < 0.)Dqni(k)=-qni(k)/delt
         ENDDO
      ELSE
         DO k=kts,kte
           Dqni(k)=0.
         ENDDO
      ENDIF
    ELSE !-MIX CLOUD SPECIES?
      !CLOUDS ARE NOT NIXED (when bl_mynn_cloudmix == 0)
      DO k=kts,kte
         Dqc(k)=0.
         Dqnc(k)=0.
         Dqi(k)=0.
         Dqni(k)=0.
      ENDDO
    ENDIF

    !ensure non-negative moist species
    CALL moisture_check(kte, delt, delp, exner,  &
                        sqv2, sqc2, sqi2, thl,   &
                        dqv, dqc, dqi, dth,      &
                        xlscp,xlvcp )

    !=====================
    ! OZONE TENDENCY CHECK
    !=====================
    DO k=kts,kte
       IF(Dozone(k)*delt + ozone(k) < 0.) THEN
         Dozone(k)=-ozone(k)*0.99/delt
       ENDIF
    ENDDO

    !===================
    ! THETA TENDENCY
    !===================
    IF (FLAG_QI) THEN
      DO k=kts,kte
         Dth(k)=(thl(k) + xlvcp/exner(k)*sqc2(k) &
           &            + xlscp/exner(k)*sqi2(k) &
           &            - th(k))/delt
         !Use form from Tripoli and Cotton (1981) with their
         !suggested min temperature to improve accuracy:
         !Dth(k)=(thl(k)*(1.+ xlvcp/MAX(tk(k),TKmin)*sqc(k)  &
         !  &               + xlscp/MAX(tk(k),TKmin)*sqi(k)) &
         !  &               - th(k))/delt
      ENDDO
    ELSE
      DO k=kts,kte
         Dth(k)=(thl(k)+xlvcp/exner(k)*sqc2(k) - th(k))/delt
         !Use form from Tripoli and Cotton (1981) with their
         !suggested min temperature to improve accuracy.
         !Dth(k)=(thl(k)*(1.+ xlvcp/MAX(tk(k),TKmin)*sqc(k))  &
         !&               - th(k))/delt
      ENDDO
    ENDIF

    !===================
    ! AEROSOL TENDENCIES
    !===================
    IF (FLAG_QNWFA .AND. FLAG_QNIFA .AND. &
        bl_mynn_mixscalars > 0) THEN
       DO k=kts,kte
          !=====================
          ! WATER-friendly aerosols
          !=====================
          Dqnwfa(k)=(qnwfa2(k) - qnwfa(k))/delt
          !=====================
          ! Ice-friendly aerosols
          !=====================
          Dqnifa(k)=(qnifa2(k) - qnifa(k))/delt
          !=====================
          ! Black-carbon aerosols
          !=====================
          Dqnbca(k)=(qnbca2(k) - qnbca(k))/delt
       ENDDO
    ELSE
       DO k=kts,kte
          Dqnwfa(k)=0.
          Dqnifa(k)=0.
          Dqnbca(k)=0.
       ENDDO
    ENDIF

    !ensure non-negative moist species
    !note: if called down here, dth needs to be updated, but
    !      if called before the theta-tendency calculation, do not compute dth
    !CALL moisture_check(kte, delt, delp, exner,     &
    !                    sqv, sqc, sqi, thl,         &
    !                    dqv, dqc, dqi, dth )

    if (debug_code) then
       problem = .false.
       do k=kts,kte
          wsp  = sqrt(u(k)**2 + v(k)**2)
          wsp2 = sqrt((u(k)+du(k)*delt)**2 + (v(k)+du(k)*delt)**2)
          th2  = th(k) + Dth(k)*delt
          tk2  = th2*exner(k)
          if (wsp2 > 200. .or. tk2 > 360. .or. tk2 < 160.) then
             problem = .true.
             print*,"Outgoing problem at: i=",i," k=",k
             print*," incoming wsp=",wsp," outgoing wsp=",wsp2
             print*," incoming T=",th(k)*exner(k),"outgoing T:",tk2
             print*," du=",du(k)*delt," dv=",dv(k)*delt," dth=",dth(k)*delt
             print*," km=",kmdz(k)*dz(k)," kh=",khdz(k)*dz(k)
             print*," u*=",ust," wspd=",wspd,"rhosfc=",rhosfc
             print*," LH=",flq*rhosfc*1004.," HFX=",flt*rhosfc*1004.
             print*," drag term=",ust**2/wspd*dtz(k)*rhosfc/rho(kts)
             kproblem = k
          endif
       enddo
       if (problem) then
          print*,"==thl:",thl(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"===qv:",sqv2(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"===qc:",sqc2(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"===qi:",sqi2(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"====u:",u(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"====v:",v(max(kproblem-3,1):min(kproblem+3,kte))
       endif
    endif

#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mynn_tendencies

END MODULE
