MODULE myn_tendencies 
IMPLICIT NONE

! ==================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine calculates hybrid diagnotic boundary-layer height (PBLH).
!!
!! NOTES ON THE PBLH FORMULATION: The 1.5-theta-increase method defines
!!PBL heights as the level at.
!!which the potential temperature first exceeds the minimum potential.
!!temperature within the boundary layer by 1.5 K. When applied to.
!!observed temperatures, this method has been shown to produce PBL-
!!height estimates that are unbiased relative to profiler-based.
!!estimates (Nielsen-Gammon et al. 2008 \cite Nielsen_Gammon_2008). 
!! However, their study did not
!!include LLJs. Banta and Pichugina (2008) \cite Pichugina_2008  show that a TKE-based.
!!threshold is a good estimate of the PBL height in LLJs. Therefore,
!!a hybrid definition is implemented that uses both methods, weighting
!!the TKE-method more during stable conditions (PBLH < 400 m).
!!A variable tke threshold (TKEeps) is used since no hard-wired
!!value could be found to work best in all conditions.
!>\section gen_get_pblh  GSD get_pblh General Algorithm
!> @{
  SUBROUTINE GET_PBLH(KTS,KTE,zi,thetav1D,qke1D,zw1D,dz1D,landsea,kzi)

    !---------------------------------------------------------------
    !             NOTES ON THE PBLH FORMULATION
    !
    !The 1.5-theta-increase method defines PBL heights as the level at 
    !which the potential temperature first exceeds the minimum potential 
    !temperature within the boundary layer by 1.5 K. When applied to 
    !observed temperatures, this method has been shown to produce PBL-
    !height estimates that are unbiased relative to profiler-based 
    !estimates (Nielsen-Gammon et al. 2008). However, their study did not
    !include LLJs. Banta and Pichugina (2008) show that a TKE-based 
    !threshold is a good estimate of the PBL height in LLJs. Therefore,
    !a hybrid definition is implemented that uses both methods, weighting
    !the TKE-method more during stable conditions (PBLH < 400 m).
    !A variable tke threshold (TKEeps) is used since no hard-wired
    !value could be found to work best in all conditions.
    !---------------------------------------------------------------

    INTEGER,INTENT(IN) :: KTS,KTE

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    REAL, INTENT(OUT) :: zi
    REAL, INTENT(IN) :: landsea
    REAL, DIMENSION(KTS:KTE), INTENT(IN) :: thetav1D, qke1D, dz1D
    REAL, DIMENSION(KTS:KTE+1), INTENT(IN) :: zw1D
    !LOCAL VARS
    REAL ::  PBLH_TKE,qtke,qtkem1,wt,maxqke,TKEeps,minthv
    REAL :: delt_thv   !delta theta-v; dependent on land/sea point
    REAL, PARAMETER :: sbl_lim  = 200. !upper limit of stable BL height (m).
    REAL, PARAMETER :: sbl_damp = 400. !transition length for blending (m).
    INTEGER :: I,J,K,kthv,ktke,kzi

    !Initialize KPBL (kzi)
    kzi = 2

    !> - FIND MIN THETAV IN THE LOWEST 200 M AGL
    k = kts+1
    kthv = 1
    minthv = 9.E9
    DO WHILE (zw1D(k) .LE. 200.)
    !DO k=kts+1,kte-1
       IF (minthv > thetav1D(k)) then
           minthv = thetav1D(k)
           kthv = k
       ENDIF
       k = k+1
       !IF (zw1D(k) .GT. sbl_lim) exit
    ENDDO

    !> - FIND THETAV-BASED PBLH (BEST FOR DAYTIME).
    zi=0.
    k = kthv+1
    IF((landsea-1.5).GE.0)THEN
        ! WATER
        delt_thv = 1.0
    ELSE
        ! LAND
        delt_thv = 1.25
    ENDIF

    zi=0.
    k = kthv+1
!    DO WHILE (zi .EQ. 0.) 
    DO k=kts+1,kte-1
       IF (thetav1D(k) .GE. (minthv + delt_thv))THEN
          zi = zw1D(k) - dz1D(k-1)* &
             & MIN((thetav1D(k)-(minthv + delt_thv))/ &
             & MAX(thetav1D(k)-thetav1D(k-1),1E-6),1.0)
       ENDIF
       !k = k+1
       IF (k .EQ. kte-1) zi = zw1D(kts+1) !EXIT SAFEGUARD
       IF (zi .NE. 0.0) exit
    ENDDO
    !print*,"IN GET_PBLH:",thsfc,zi

    !> - FOR STABLE BOUNDARY LAYERS, USE TKE METHOD TO COMPLEMENT THE
    !! THETAV-BASED DEFINITION (WHEN THE THETA-V BASED PBLH IS BELOW ~0.5 KM).
    !!THE TANH WEIGHTING FUNCTION WILL MAKE THE TKE-BASED DEFINITION NEGLIGIBLE 
    !!WHEN THE THETA-V-BASED DEFINITION IS ABOVE ~1 KM.
    ktke = 1
    maxqke = MAX(Qke1D(kts),0.)
    !Use 5% of tke max (Kosovic and Curry, 2000; JAS)
    !TKEeps = maxtke/20. = maxqke/40.
    TKEeps = maxqke/40.
    TKEeps = MAX(TKEeps,0.02) !0.025) 
    PBLH_TKE=0.

    k = ktke+1
!    DO WHILE (PBLH_TKE .EQ. 0.) 
    DO k=kts+1,kte-1
       !QKE CAN BE NEGATIVE (IF CKmod == 0)... MAKE TKE NON-NEGATIVE.
       qtke  =MAX(Qke1D(k)/2.,0.)      ! maximum TKE
       qtkem1=MAX(Qke1D(k-1)/2.,0.)
       IF (qtke .LE. TKEeps) THEN
           PBLH_TKE = zw1D(k) - dz1D(k-1)* &
             & MIN((TKEeps-qtke)/MAX(qtkem1-qtke, 1E-6), 1.0)
           !IN CASE OF NEAR ZERO TKE, SET PBLH = LOWEST LEVEL.
           PBLH_TKE = MAX(PBLH_TKE,zw1D(kts+1))
           !print *,"PBLH_TKE:",i,PBLH_TKE, Qke1D(k)/2., zw1D(kts+1)
       ENDIF
       !k = k+1
       IF (k .EQ. kte-1) PBLH_TKE = zw1D(kts+1) !EXIT SAFEGUARD
       IF (PBLH_TKE .NE. 0.) exit
    ENDDO

    !> - With TKE advection turned on, the TKE-based PBLH can be very large 
    !! in grid points with convective precipitation (> 8 km!),
    !! so an artificial limit is imposed to not let PBLH_TKE exceed the
    !!theta_v-based PBL height +/- 350 m.
    !!This has no impact on 98-99% of the domain, but is the simplest patch
    !!that adequately addresses these extremely large PBLHs.
    PBLH_TKE = MIN(PBLH_TKE,zi+350.)
    PBLH_TKE = MAX(PBLH_TKE,MAX(zi-350.,10.))

    wt=.5*TANH((zi - sbl_lim)/sbl_damp) + .5
    IF (maxqke <= 0.05) THEN
       !Cold pool situation - default to theta_v-based def
    ELSE
       !BLEND THE TWO PBLH TYPES HERE: 
       zi=PBLH_TKE*(1.-wt) + zi*wt
    ENDIF

    !Compute KPBL (kzi)
    DO k=kts+1,kte-1
       IF ( zw1D(k) >= zi) THEN
          kzi = k-1
          exit
       ENDIF
    ENDDO

#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE GET_PBLH
!> @}

END MODULE
