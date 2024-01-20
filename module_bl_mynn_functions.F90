MODULE mynn_functions 

    use module_bl_mynn_common,only: tice,t0c,cice,cliq,cpv,xls,xlv 


! JAYMES-
!> Constants used for empirical calculations of saturation
!! vapor pressures (in function "esat") and saturation mixing ratios
!! (in function "qsat"), reproduced from module_mp_thompson.F,
!! v3.6
  REAL, PARAMETER:: J0= .611583699E03
  REAL, PARAMETER:: J1= .444606896E02
  REAL, PARAMETER:: J2= .143177157E01
  REAL, PARAMETER:: J3= .264224321E-1
  REAL, PARAMETER:: J4= .299291081E-3
  REAL, PARAMETER:: J5= .203154182E-5
  REAL, PARAMETER:: J6= .702620698E-8
  REAL, PARAMETER:: J7= .379534310E-11
  REAL, PARAMETER:: J8=-.321582393E-13

  REAL, PARAMETER:: K0= .609868993E03
  REAL, PARAMETER:: K1= .499320233E02
  REAL, PARAMETER:: K2= .184672631E01
  REAL, PARAMETER:: K3= .402737184E-1
  REAL, PARAMETER:: K4= .565392987E-3
  REAL, PARAMETER:: K5= .521693933E-5
  REAL, PARAMETER:: K6= .307839583E-7
  REAL, PARAMETER:: K7= .105785160E-9
  REAL, PARAMETER:: K8= .161444444E-12


CONTAINS
! ====================================================================

!>\ingroup gsd_mynn_edmf
!! This function extends function "esat" and returns a "blended"
!! saturation mixing ratio.
!!\author JAYMES
  FUNCTION qsat_blend(t, P, waterice)

      IMPLICIT NONE

      REAL, INTENT(IN):: t, P
      CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: waterice
      CHARACTER(LEN=1) :: wrt
      REAL :: qsat_blend,XC,ESL,ESI,RSLF,RSIF,chi

      IF ( .NOT. PRESENT(waterice) ) THEN 
          wrt = 'b'
      ELSE
          wrt = waterice
      ENDIF

      XC=MAX(-80.,t - t0c)

      IF ((t .GE. t0c) .OR. (wrt .EQ. 'w')) THEN
          ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8))))))) 
          qsat_blend = 0.622*ESL/max(P-ESL, 1e-5) 
!      ELSE IF (t .LE. 253.) THEN
      ELSE IF (t .LE. tice) THEN
          ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          qsat_blend = 0.622*ESI/max(P-ESI, 1e-5)
      ELSE
          ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
          ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          RSLF = 0.622*ESL/max(P-ESL, 1e-5)
          RSIF = 0.622*ESI/max(P-ESI, 1e-5)
!          chi  = (273.16-t)/20.16
          chi  = (t0c - t)/(t0c - tice) 
         qsat_blend = (1.-chi)*RSLF + chi*RSIF
      END IF

  END FUNCTION qsat_blend

! ===================================================================


! =====================================================================
!>\ingroup gsd_mynn_edmf
!! \author JAYMES- added 22 Apr 2015
!! This function calculates saturation vapor pressure.  Separate ice and liquid functions
!! are used (identical to those in module_mp_thompson.F, v3.6). Then, the
!! final returned value is a temperature-dependant "blend". Because the final
!! value is "phase-aware", this formulation may be preferred for use throughout
!! the module (replacing "svp").
  FUNCTION esat_blend(t)

      IMPLICIT NONE

      REAL, INTENT(IN):: t
      REAL :: esat_blend,XC,ESL,ESI,chi

      XC=MAX(-80.,t - t0c) !note t0c = 273.15, tice is set in module mynn_common

! For 253 < t < 273.16 K, the vapor pressures are "blended" as a function of temperature, 
! using the approach of Chaboureau and Bechtold (2002), JAS, p. 2363.  The resulting 
! values are returned from the function.
      IF (t .GE. t0c) THEN
          esat_blend = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
      ELSE IF (t .LE. tice) THEN
          esat_blend = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
      ELSE
          ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
          ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          chi  = (t0c - t)/(t0c - tice)
          esat_blend = (1.-chi)*ESL  + chi*ESI
      END IF

  END FUNCTION esat_blend

! ====================================================================

!>\ingroup gsd_mynn_edmf
!! This function interpolates the latent heats of vaporization and sublimation into
!! a single, temperature-dependent, "blended" value, following
!! Chaboureau and Bechtold (2002) \cite Chaboureau_2002, Appendix.
!!\author JAYMES
  FUNCTION xl_blend(t)

      IMPLICIT NONE

      REAL, INTENT(IN):: t
      REAL :: xl_blend,xlvt,xlst,chi
      !note: t0c = 273.15, tice is set in mynn_common

      IF (t .GE. t0c) THEN
          xl_blend = xlv + (cpv-cliq)*(t-t0c)  !vaporization/condensation
      ELSE IF (t .LE. tice) THEN
          xl_blend = xls + (cpv-cice)*(t-t0c)  !sublimation/deposition
      ELSE
          xlvt = xlv + (cpv-cliq)*(t-t0c)  !vaporization/condensation
          xlst = xls + (cpv-cice)*(t-t0c)  !sublimation/deposition
!          chi  = (273.16-t)/20.16
          chi  = (t0c - t)/(t0c - tice)
          xl_blend = (1.-chi)*xlvt + chi*xlst     !blended
      END IF

  END FUNCTION xl_blend

! ===================================================================

FUNCTION calc_phim(zet, cphh_unstm, cphm_unst)
     ! New stability function parameters for momentum (Puhales, 2020, WRF 4.2.1)
     ! The forms in unstable conditions (z/L < 0) use Grachev et al. (2000), which are a blend of
     ! the classical (Kansas) forms (i.e., Paulson 1970, Dyer and Hicks 1970), valid for weakly
     ! unstable conditions (-1 < z/L < 0). The stability functions for stable conditions use an
     ! updated form taken from Cheng and Brutsaert (2005), which extends the validity into very
     ! stable conditions [z/L ~ O(10)].
      IMPLICIT NONE

      REAL, INTENT(IN):: zet, cphh_unstm, cphm_unst
      REAL :: dummy_0,dummy_1,dummy_11,dummy_2,dummy_22,dummy_3,dummy_33,dummy_4,dummy_44,dummy_psi
      REAL, PARAMETER :: am_st=6.1, bm_st=2.5, rbm_st=1./bm_st
      REAL, PARAMETER :: ah_st=5.3, bh_st=1.1, rbh_st=1./bh_st
      REAL, PARAMETER :: am_unst=10., ah_unst=34.
      REAL :: phi_m,calc_phim

      if ( zet >= 0.0 ) then
         dummy_0=1+zet**bm_st
         dummy_1=zet+dummy_0**(rbm_st)
         dummy_11=1+dummy_0**(rbm_st-1)*zet**(bm_st-1)
         dummy_2=(-am_st/dummy_1)*dummy_11
         phi_m = 1-zet*dummy_2
      else
         dummy_0 = (1.0-cphm_unst*zet)**0.25
         phi_m = 1./dummy_0
         dummy_psi = 2.*log(0.5*(1.+dummy_0))+log(0.5*(1.+dummy_0**2))-2.*atan(dummy_0)+1.570796

         dummy_0=(1.-am_unst*zet)          ! parentesis arg
         dummy_1=dummy_0**0.333333         ! y
         dummy_11=-0.33333*am_unst*dummy_0**(-0.6666667) ! dy/dzet
         dummy_2 = 0.33333*(dummy_1**2.+dummy_1+1.)      ! f
         dummy_22 = 0.3333*dummy_11*(2.*dummy_1+1.)      ! df/dzet
         dummy_3 = 0.57735*(2.*dummy_1+1.) ! g
         dummy_33 = 1.1547*dummy_11        ! dg/dzet
         dummy_4 = 1.5*log(dummy_2)-1.73205*atan(dummy_3)+1.813799364 !psic
         dummy_44 = (1.5/dummy_2)*dummy_22-1.73205*dummy_33/(1.+dummy_3**2)! dpsic/dzet

         dummy_0 = zet**2
         dummy_1 = 1./(1.+dummy_0) ! denon
         dummy_11 = 2.*zet         ! denon/dzet
         dummy_2 = ((1-phi_m)/zet+dummy_11*dummy_4+dummy_0*dummy_44)*dummy_1
         dummy_22 = -dummy_11*(dummy_psi+dummy_0*dummy_4)*dummy_1**2

         phi_m = 1.-zet*(dummy_2+dummy_22)
      end if

      !phim = phi_m - zet
      calc_phim = phi_m

  END FUNCTION calc_phim
! ===================================================================

FUNCTION calc_phih(zet,cphh_unst)
    ! New stability function parameters for heat (Puhales, 2020, WRF 4.2.1)
    ! The forms in unstable conditions (z/L < 0) use Grachev et al. (2000), which are a blend of
    ! the classical (Kansas) forms (i.e., Paulson 1970, Dyer and Hicks 1970), valid for weakly
    ! unstable conditions (-1 < z/L < 0). The stability functions for stable conditions use an
    ! updated form taken from Cheng and Brutsaert (2005), which extends the validity into very
    ! stable conditions [z/L ~ O(10)].
      IMPLICIT NONE

      REAL, INTENT(IN):: zet, cphh_unst
      REAL :: dummy_0,dummy_1,dummy_11,dummy_2,dummy_22,dummy_3,dummy_33,dummy_4,dummy_44,dummy_psi
      REAL, PARAMETER :: am_st=6.1, bm_st=2.5, rbm_st=1./bm_st
      REAL, PARAMETER :: ah_st=5.3, bh_st=1.1, rbh_st=1./bh_st
      REAL, PARAMETER :: am_unst=10., ah_unst=34.
      REAL :: phh,calc_phih

      if ( zet >= 0.0 ) then
         dummy_0=1+zet**bh_st
         dummy_1=zet+dummy_0**(rbh_st)
         dummy_11=1+dummy_0**(rbh_st-1)*zet**(bh_st-1)
         dummy_2=(-ah_st/dummy_1)*dummy_11
         calc_phih = 1-zet*dummy_2
      else
         dummy_0 = (1.0-cphh_unst*zet)**0.5
         phh = 1./dummy_0
         dummy_psi = 2.*log(0.5*(1.+dummy_0))

         dummy_0=(1.-ah_unst*zet)          ! parentesis arg
         dummy_1=dummy_0**0.333333         ! y
         dummy_11=-0.33333*ah_unst*dummy_0**(-0.6666667) ! dy/dzet
         dummy_2 = 0.33333*(dummy_1**2.+dummy_1+1.)      ! f
         dummy_22 = 0.3333*dummy_11*(2.*dummy_1+1.)      ! df/dzet
         dummy_3 = 0.57735*(2.*dummy_1+1.) ! g
         dummy_33 = 1.1547*dummy_11        ! dg/dzet
         dummy_4 = 1.5*log(dummy_2)-1.73205*atan(dummy_3)+1.813799364 !psic
         dummy_44 = (1.5/dummy_2)*dummy_22-1.73205*dummy_33/(1.+dummy_3**2)! dpsic/dzet

         dummy_0 = zet**2
         dummy_1 = 1./(1.+dummy_0)         ! denon
         dummy_11 = 2.*zet                 ! ddenon/dzet
         dummy_2 = ((1-phh)/zet+dummy_11*dummy_4+dummy_0*dummy_44)*dummy_1
         dummy_22 = -dummy_11*(dummy_psi+dummy_0*dummy_4)*dummy_1**2

         calc_phih = 1.-zet*(dummy_2+dummy_22)
      end if

END FUNCTION calc_phih


! ===================================================================


END MODULE
