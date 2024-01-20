MODULE mynn_condensation
use module_bl_mynn_common,only:&
        cp,cpv,ep_2,ep_3,p608,r_d,r_v,tice,tv0,xlv,xlvcp
use mynn_functions
IMPLICIT NONE

CONTAINS
! ==================================================================
!     SUBROUTINE  mym_condensation:
!
!     Input variables:    see subroutine mym_initialize and turbulence
!       exner(nz)    : Perturbation of the Exner function    (J/kg K)
!                         defined on the walls of the grid boxes
!                         This is usually computed by integrating
!                         d(pi)/dz = h*g*tv/tref**2
!                         from the upper boundary, where tv is the
!                         virtual potential temperature minus tref.
!
!     Output variables:   see subroutine mym_initialize
!       cld(nx,nz,ny)   : Cloud fraction
!
!     Work arrays/variables:
!       qmq             : Q_w-Q_{sl}, where Q_{sl} is the saturation
!                         specific humidity at T=Tl
!       alp(nx,nz,ny)   : Functions in the condensation process
!       bet(nx,nz,ny)   : ditto
!       sgm(nx,nz,ny)   : Combined standard deviation sigma_s
!                         multiplied by 2/alp
!
!     # qmq, alp, bet and sgm are allowed to share storage units with
!       any four of other work arrays for saving memory.
!
!     # Results are sensitive particularly to values of cp and r_d.
!       Set these values to those adopted by you.
!
!-------------------------------------------------------------------
!>\ingroup gsd_mynn_edmf 
!! This subroutine calculates the nonconvective component of the 
!! subgrid cloud fraction and mixing ratio as well as the functions used to 
!! calculate the buoyancy flux. Different cloud PDFs can be selected by
!! use of the namelist parameter \p bl_mynn_cloudpdf .
  SUBROUTINE  mym_condensation (kts,kte,   &
    &            dx, dz, zw, xland,        &
    &            thl, qw, qv, qc, qi,      &
    &            p,exner,                  &
    &            tsq, qsq, cov,            &
    &            Sh, el, bl_mynn_cloudpdf, &
    &            qc_bl1D, qi_bl1D,         &
    &            cldfra_bl1D,              &
    &            PBLH1,HFX1,               &
    &            Vt, Vq, th, sgm, rmo,     &
    &            spp_pbl,rstoch_col,       &
    &            b2,rrp,rr2,tliq           )

!-------------------------------------------------------------------

    INTEGER, INTENT(IN)   :: kts,kte, bl_mynn_cloudpdf

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    REAL, INTENT(IN)      :: dx,PBLH1,HFX1,rmo,xland,b2,rrp,rr2,tliq
    REAL, DIMENSION(kts:kte), INTENT(IN) :: dz
    REAL, DIMENSION(kts:kte+1), INTENT(IN) :: zw
    REAL, DIMENSION(kts:kte), INTENT(IN) :: p,exner,thl,qw,qv,qc,qi, &
         &tsq, qsq, cov, th

    REAL, DIMENSION(kts:kte), INTENT(INOUT) :: vt,vq,sgm

    REAL, DIMENSION(kts:kte) :: alp,a,bet,b,ql,q1,RH
    REAL, DIMENSION(kts:kte), INTENT(OUT) :: qc_bl1D,qi_bl1D, &
                                             cldfra_bl1D
    DOUBLE PRECISION :: t3sq, r3sq, c3sq

    REAL :: qsl,esat,qsat,dqsl,cld0,q1k,qlk,eq1,qll,&
         &q2p,pt,rac,qt,t,xl,rsl,cpm,Fng,qww,alpha,beta,bb,&
         &ls,wt,cld_factor,fac_damp,liq_frac,ql_ice,ql_water,&
         &qmq,qsat_tk
    INTEGER :: i,j,k

    REAL :: erf

    !VARIABLES FOR ALTERNATIVE SIGMA
    REAL::dth,dtl,dqw,dzk,els
    REAL, DIMENSION(kts:kte), INTENT(IN) :: Sh,el

    !variables for SGS BL clouds
    REAL            :: zagl,damp,PBLH2
    REAL            :: cfmax

    !JAYMES:  variables for tropopause-height estimation
    REAL            :: theta1, theta2, ht1, ht2
    INTEGER         :: k_tropo

!   Stochastic
    INTEGER,  INTENT(IN)                          ::    spp_pbl
    REAL, DIMENSION(KTS:KTE)                      ::    rstoch_col
    REAL :: qw_pert

! First, obtain an estimate for the tropopause height (k), using the method employed in the
! Thompson subgrid-cloud scheme.  This height will be a consideration later when determining 
! the "final" subgrid-cloud properties.
! JAYMES:  added 3 Nov 2016, adapted from G. Thompson

    DO k = kte-3, kts, -1
       theta1 = th(k)
       theta2 = th(k+2)
       ht1 = 44307.692 * (1.0 - (p(k)/101325.)**0.190)
       ht2 = 44307.692 * (1.0 - (p(k+2)/101325.)**0.190)
       if ( (((theta2-theta1)/(ht2-ht1)) .lt. 10./1500. ) .AND.       &
     &                       (ht1.lt.19000.) .and. (ht1.gt.4000.) ) then 
          goto 86
       endif
    ENDDO
 86   continue
    k_tropo = MAX(kts+2, k+2)

    zagl = 0.

    SELECT CASE(bl_mynn_cloudpdf)

      CASE (0) ! ORIGINAL MYNN PARTIAL-CONDENSATION SCHEME

        DO k = kts,kte-1
           t  = th(k)*exner(k)

!x      if ( ct .gt. 0.0 ) then
!       a  =  17.27
!       b  = 237.3
!x      else
!x        a  =  21.87
!x        b  = 265.5
!x      end if
!
!   **  3.8 = 0.622*6.11 (hPa)  **

           !SATURATED VAPOR PRESSURE
           esat = esat_blend(t)
           !SATURATED SPECIFIC HUMIDITY
           !qsl=ep_2*esat/(p(k)-ep_3*esat)
           qsl=ep_2*esat/max(1.e-4,(p(k)-ep_3*esat))
           !dqw/dT: Clausius-Clapeyron
           dqsl = qsl*ep_2*xlv/( r_d*t**2 )

           alp(k) = 1.0/( 1.0+dqsl*xlvcp )
           bet(k) = dqsl*exner(k)

           !Sommeria and Deardorff (1977) scheme, as implemented
           !in Nakanishi and Niino (2009), Appendix B
           t3sq = MAX( tsq(k), 0.0 )
           r3sq = MAX( qsq(k), 0.0 )
           c3sq =      cov(k)
           c3sq = SIGN( MIN( ABS(c3sq), SQRT(t3sq*r3sq) ), c3sq )
           r3sq = r3sq +bet(k)**2*t3sq -2.0*bet(k)*c3sq
           !DEFICIT/EXCESS WATER CONTENT
           qmq  = qw(k) -qsl
           !ORIGINAL STANDARD DEVIATION
           sgm(k) = SQRT( MAX( r3sq, 1.0d-10 ))
           !NORMALIZED DEPARTURE FROM SATURATION
           q1(k)   = qmq / sgm(k)
           !CLOUD FRACTION. rr2 = 1/SQRT(2) = 0.707
           cldfra_bl1D(k) = 0.5*( 1.0+erf( q1(k)*rr2 ) )

           q1k  = q1(k)
           eq1  = rrp*EXP( -0.5*q1k*q1k )
           qll  = MAX( cldfra_bl1D(k)*q1k + eq1, 0.0 )
           !ESTIMATED LIQUID WATER CONTENT (UNNORMALIZED)
           ql(k) = alp(k)*sgm(k)*qll
           !LIMIT SPECIES TO TEMPERATURE RANGES
           liq_frac = min(1.0, max(0.0,(t-240.0)/29.0))
           qc_bl1D(k) = liq_frac*ql(k)
           qi_bl1D(k) = (1.0 - liq_frac)*ql(k)

           if(cldfra_bl1D(k)>0.01 .and. qc_bl1D(k)<1.E-6)qc_bl1D(k)=1.E-6
           if(cldfra_bl1D(k)>0.01 .and. qi_bl1D(k)<1.E-8)qi_bl1D(k)=1.E-8

           !Now estimate the buoyancy flux functions
           q2p = xlvcp/exner(k)
           pt = thl(k) +q2p*ql(k) ! potential temp

           !qt is a THETA-V CONVERSION FOR TOTAL WATER (i.e., THETA-V = qt*THETA)
           qt   = 1.0 +p608*qw(k) -(1.+p608)*(qc_bl1D(k)+qi_bl1D(k))*cldfra_bl1D(k)
           rac  = alp(k)*( cldfra_bl1D(K)-qll*eq1 )*( q2p*qt-(1.+p608)*pt )

           !BUOYANCY FACTORS: wherever vt and vq are used, there is a
           !"+1" and "+tv0", respectively, so these are subtracted out here.
           !vt is unitless and vq has units of K.
           vt(k) =      qt-1.0 -rac*bet(k)
           vq(k) = p608*pt-tv0 +rac

        END DO

      CASE (1, -1) !ALTERNATIVE FORM (Nakanishi & Niino 2004 BLM, eq. B6, and
                       !Kuwano-Yoshida et al. 2010 QJRMS, eq. 7):
        DO k = kts,kte-1
           t  = th(k)*exner(k)
           !SATURATED VAPOR PRESSURE
           esat = esat_blend(t)
           !SATURATED SPECIFIC HUMIDITY
           !qsl=ep_2*esat/(p(k)-ep_3*esat)
           qsl=ep_2*esat/max(1.e-4,(p(k)-ep_3*esat))
           !dqw/dT: Clausius-Clapeyron
           dqsl = qsl*ep_2*xlv/( r_d*t**2 )

           alp(k) = 1.0/( 1.0+dqsl*xlvcp )
           bet(k) = dqsl*exner(k)

           if (k .eq. kts) then 
             dzk = 0.5*dz(k)
           else
             dzk = dz(k)
           end if
           dth = 0.5*(thl(k+1)+thl(k)) - 0.5*(thl(k)+thl(MAX(k-1,kts)))
           dqw = 0.5*(qw(k+1) + qw(k)) - 0.5*(qw(k) + qw(MAX(k-1,kts)))
           sgm(k) = SQRT( MAX( (alp(k)**2 * MAX(el(k)**2,0.1) * &
                             b2 * MAX(Sh(k),0.03))/4. * &
                      (dqw/dzk - bet(k)*(dth/dzk ))**2 , 1.0e-10) )
           qmq   = qw(k) -qsl
           q1(k) = qmq / sgm(k)
           cldfra_bl1D(K) = 0.5*( 1.0+erf( q1(k)*rr2 ) )

           !now compute estimated lwc for PBL scheme's use 
           !qll IS THE NORMALIZED LIQUID WATER CONTENT (Sommeria and
           !Deardorff (1977, eq 29a). rrp = 1/(sqrt(2*pi)) = 0.3989
           q1k  = q1(k)
           eq1  = rrp*EXP( -0.5*q1k*q1k )
           qll  = MAX( cldfra_bl1D(K)*q1k + eq1, 0.0 )
           !ESTIMATED LIQUID WATER CONTENT (UNNORMALIZED)
           ql (k) = alp(k)*sgm(k)*qll
           liq_frac = min(1.0, max(0.0,(t-240.0)/29.0))
           qc_bl1D(k) = liq_frac*ql(k)
           qi_bl1D(k) = (1.0 - liq_frac)*ql(k)

           if(cldfra_bl1D(k)>0.01 .and. qc_bl1D(k)<1.E-6)qc_bl1D(k)=1.E-6
           if(cldfra_bl1D(k)>0.01 .and. qi_bl1D(k)<1.E-8)qi_bl1D(k)=1.E-8

           !Now estimate the buoyancy flux functions
           q2p = xlvcp/exner(k)
           pt = thl(k) +q2p*ql(k) ! potential temp

           !qt is a THETA-V CONVERSION FOR TOTAL WATER (i.e., THETA-V = qt*THETA)
           qt   = 1.0 +p608*qw(k) -(1.+p608)*(qc_bl1D(k)+qi_bl1D(k))*cldfra_bl1D(k)
           rac  = alp(k)*( cldfra_bl1D(K)-qll*eq1 )*( q2p*qt-(1.+p608)*pt )

           !BUOYANCY FACTORS: wherever vt and vq are used, there is a
           !"+1" and "+tv0", respectively, so these are subtracted out here.
           !vt is unitless and vq has units of K.
           vt(k) =      qt-1.0 -rac*bet(k)
           vq(k) = p608*pt-tv0 +rac

        END DO

      CASE (2, -2)

        !Diagnostic statistical scheme of Chaboureau and Bechtold (2002), JAS
        !but with use of higher-order moments to estimate sigma
        PBLH2=MAX(10.,PBLH1)
        zagl = 0.
        DO k = kts,kte-1
           zagl = zagl + dz(k)
           t  = th(k)*exner(k)

           xl = xl_blend(t)                  ! obtain latent heat
           qsat_tk = qsat_blend(t,  p(k))    ! saturation water vapor mixing ratio at tk and p
           rh(k)=MAX(MIN(1.0,qw(k)/MAX(1.E-8,qsat_tk)),0.001)

           !dqw/dT: Clausius-Clapeyron
           dqsl = qsat_tk*ep_2*xlv/( r_d*t**2 )
           alp(k) = 1.0/( 1.0+dqsl*xlvcp )
           bet(k) = dqsl*exner(k)
 
           rsl = xl*qsat_tk / (r_v*t**2)     ! slope of C-C curve at t (=abs temperature)
                                             ! CB02, Eqn. 4
           cpm = cp + qw(k)*cpv              ! CB02, sec. 2, para. 1
           a(k) = 1./(1. + xl*rsl/cpm)       ! CB02 variable "a"
           b(k) = a(k)*rsl                   ! CB02 variable "b"

           !SPP
           qw_pert = qw(k) + qw(k)*0.5*rstoch_col(k)*real(spp_pbl)

           !This form of qmq (the numerator of Q1) no longer uses the a(k) factor
           qmq    = qw_pert - qsat_tk          ! saturation deficit/excess;

           !Use the form of Eq. (6) in Chaboureau and Bechtold (2002)
           !except neglect all but the first term for sig_r
           r3sq = MAX( qsq(k), 0.0 )
           !Calculate sigma using higher-order moments:
           sgm(k) = SQRT( r3sq )
           !Set limits on sigma relative to saturation water vapor
           sgm(k) = MIN( sgm(k), qsat_tk*0.666 ) !500 )
           sgm(k) = MAX( sgm(k), qsat_tk*0.035 ) !Note: 0.02 results in SWDOWN similar
                                                 !to the first-order version of sigma
           q1(k) = qmq  / sgm(k)  ! Q1, the normalized saturation
           q1k   = q1(k)          ! backup Q1 for later modification

           ! Specify cloud fraction
           !Original C-B cloud fraction, allows cloud fractions out to q1 = -3.5
           !cldfra_bl1D(K) = max(0., min(1., 0.5+0.36*atan(1.55*q1(k)))) ! Eq. 7 in CB02
           !Waynes LES fit  - over-diffuse, when limits removed from vt & vq & fng
           !cldfra_bl1D(K) = max(0., min(1., 0.5+0.36*atan(1.2*(q1(k)+0.4))))
           !Best compromise: Improves marine stratus without adding much cold bias.
           cldfra_bl1D(K) = max(0., min(1., 0.5+0.36*atan(1.8*(q1(k)+0.2))))

           ! Specify hydrometeors
           ! JAYMES- this option added 8 May 2015
           ! The cloud water formulations are taken from CB02, Eq. 8.
           IF (q1k < 0.) THEN        !unsaturated
              ql_water = sgm(k)*EXP(1.2*q1k-1)
              ql_ice   = sgm(k)*EXP(1.2*q1k-1.)
           ELSE IF (q1k > 2.) THEN   !supersaturated
              ql_water = sgm(k)*q1k
              ql_ice   = sgm(k)*q1k
           ELSE                      !slightly saturated (0 > q1 < 2)
              ql_water = sgm(k)*(EXP(-1.) + 0.66*q1k + 0.086*q1k**2)
              ql_ice   = sgm(k)*(EXP(-1.) + 0.66*q1k + 0.086*q1k**2)
           ENDIF

           !In saturated grid cells, use average of SGS and resolved values
           if ( qc(k) > 1.e-6 ) ql_water = 0.5 * ( ql_water + qc(k) ) 
           !since ql_ice is actually the total frozen condensate (snow+ice),
           !do not average with grid-scale ice alone 
           !if ( qi(k) > 1.e-9 ) ql_ice = 0.5 * ( ql_ice + qi(k) )

           if (cldfra_bl1D(k) < 0.01) then
              ql_ice   = 0.0
              ql_water = 0.0
              cldfra_bl1D(k) = 0.0
           endif

           !PHASE PARTITIONING: currently commented out since we are moving towards prognostic sgs clouds
           !Make some inferences about the relative amounts of
           !subgrid cloud water vs. ice based on collocated explicit clouds.  Otherise, 
           !use a simple temperature-dependent partitioning.
           ! IF ( qc(k) + qi(k) > 0.0 ) THEN ! explicit condensate exists, retain its phase partitioning
           !    IF ( qi(k) == 0.0 ) THEN       ! explicit contains no ice; assume subgrid liquid
           !      liq_frac = 1.0
           !    ELSE IF ( qc(k) == 0.0 ) THEN  ! explicit contains no liquid; assume subgrid ice
           !      liq_frac = 0.0
           !    ELSE IF ( (qc(k) >= 1.E-10) .AND. (qi(k) >= 1.E-10) ) THEN  ! explicit contains mixed phase of workably 
           !                                                                ! large amounts; assume subgrid follows 
           !                                                               ! same partioning
           !      liq_frac = qc(k) / ( qc(k) + qi(k) )
           !    ELSE
           !      liq_frac = MIN(1.0, MAX(0.0, (t-tice)/(t0c-tice))) ! explicit contains mixed phase, but at least one 
           !                                                         ! species is very small, so make a temperature-
           !                                                         ! depedent guess
           !    ENDIF
           ! ELSE                          ! no explicit condensate, so make a temperature-dependent guess
             liq_frac = MIN(1.0, MAX(0.0, (t-tice)/(tliq-tice)))
           ! ENDIF

           qc_bl1D(k) = liq_frac*ql_water       ! apply liq_frac to ql_water and ql_ice
           qi_bl1D(k) = (1.0-liq_frac)*ql_ice

           !Above tropopause:  eliminate subgrid clouds from CB scheme
           if (k .ge. k_tropo-1) then
              cldfra_bl1D(K) = 0.
              qc_bl1D(k)     = 0.
              qi_bl1D(k)     = 0.
           endif

           !Buoyancy-flux-related calculations follow...
           !limiting Q1 to avoid too much diffusion in cloud layers
           q1k=max(Q1(k),-2.0)

           ! "Fng" represents the non-Gaussian transport factor
           ! (non-dimensional) from Bechtold et al. 1995 
           ! (hereafter BCMT95), section 3(c).  Their suggested 
           ! forms for Fng (from their Eq. 20) are:
           !IF (q1k < -2.) THEN
           !  Fng = 2.-q1k
           !ELSE IF (q1k > 0.) THEN
           !  Fng = 1.
           !ELSE
           !  Fng = 1.-1.5*q1k
           !ENDIF
           ! Use the form of "Fng" from Bechtold and Siebesma (1998, JAS)
           IF (q1k .GE. 1.0) THEN
              Fng = 1.0
           ELSEIF (q1k .GE. -1.7 .AND. q1k .LT. 1.0) THEN
              Fng = EXP(-0.4*(q1k-1.0))
           ELSEIF (q1k .GE. -2.5 .AND. q1k .LT. -1.7) THEN
              Fng = 3.0 + EXP(-3.8*(q1k+1.7))
           ELSE
              Fng = MIN(23.9 + EXP(-1.6*(q1k+2.5)), 60.)
           ENDIF

           cfmax= min(cldfra_bl1D(k), 0.5)
           bb = b(k)*t/th(k) ! bb is "b" in BCMT95.  Their "b" differs from 
                             ! "b" in CB02 (i.e., b(k) above) by a factor 
                             ! of T/theta.  Strictly, b(k) above is formulated in
                             ! terms of sat. mixing ratio, but bb in BCMT95 is
                             ! cast in terms of sat. specific humidity.  The
                             ! conversion is neglected here. 
           qww   = 1.+0.61*qw(k)
           alpha = 0.61*th(k)
           beta  = (th(k)/t)*(xl/cp) - 1.61*th(k)
           vt(k) = qww   - cfmax*beta*bb*Fng   - 1.
           vq(k) = alpha + cfmax*beta*a(k)*Fng - tv0
           ! vt and vq correspond to beta-theta and beta-q, respectively,  
           ! in NN09, Eq. B8.  They also correspond to the bracketed
           ! expressions in BCMT95, Eq. 15, since (s*ql/sigma^2) = cldfra*Fng
           ! The "-1" and "-tv0" terms are included for consistency with 
           ! the legacy vt and vq formulations (above).

           ! dampen amplification factor where need be
           fac_damp = min(zagl * 0.0025, 1.0)
           !cld_factor = 1.0 + fac_damp*MAX(0.0, ( RH(k) - 0.75 ) / 0.26 )**1.9 !HRRRv4
           !cld_factor = 1.0 + fac_damp*min((max(0.0, ( RH(k) - 0.92 )) / 0.25 )**2, 0.3)
           cld_factor = 1.0 + fac_damp*min((max(0.0, ( RH(k) - 0.92 )) / 0.145)**2, 0.4)
           cldfra_bl1D(K) = MIN( 1., cld_factor*cldfra_bl1D(K) )
        enddo

      END SELECT !end cloudPDF option

      !For testing purposes only, option for isolating on the mass-flux clouds.
      IF (bl_mynn_cloudpdf .LT. 0) THEN
         DO k = kts,kte-1
            cldfra_bl1D(k) = 0.0
            qc_bl1D(k) = 0.0
            qi_bl1D(k) = 0.0
         END DO
      ENDIF
!
      ql(kte) = ql(kte-1)
      vt(kte) = vt(kte-1)
      vq(kte) = vq(kte-1)
      qc_bl1D(kte)=0.
      qi_bl1D(kte)=0.
      cldfra_bl1D(kte)=0.
    RETURN

#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mym_condensation

! ==================================================================

END MODULE
 
