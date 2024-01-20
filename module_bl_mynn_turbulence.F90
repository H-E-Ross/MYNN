MODULE turbulence
use module_bl_mynn_common,only:tv0,gtr
use mynn_functions
IMPLICIT NONE
CONTAINS
! ==================================================================
!     SUBROUTINE  mym_turbulence:
!
!     Input variables:    see subroutine mym_initialize
!       closure        : closure level (2.5, 2.6, or 3.0)
!
!     # ql, vt, vq, qke, tsq, qsq and cov are changed to input variables.
!
!     Output variables:   see subroutine mym_initialize
!       dfm(nx,nz,ny) : Diffusivity coefficient for momentum,
!                         divided by dz (not dz*h(i,j))            (m/s)
!       dfh(nx,nz,ny) : Diffusivity coefficient for heat,
!                         divided by dz (not dz*h(i,j))            (m/s)
!       dfq(nx,nz,ny) : Diffusivity coefficient for q^2,
!                         divided by dz (not dz*h(i,j))            (m/s)
!       tcd(nx,nz,ny)   : Countergradient diffusion term for Theta_l
!                                                                  (K/s)
!       qcd(nx,nz,ny)   : Countergradient diffusion term for Q_w
!                                                              (kg/kg s)
!       pd?(nx,nz,ny) : Half of the production terms
!
!       Only tcd and qcd are defined at the center of the grid boxes
!
!     # DO NOT forget that tcd and qcd are added on the right-hand side
!       of the equations for Theta_l and Q_w, respectively.
!
!     Work arrays:        see subroutine mym_initialize and level2
!
!     # dtl, dqw, dtv, gm and gh are allowed to share storage units with
!       dfm, dfh, dfq, tcd and qcd, respectively, for saving memory.
!
!>\ingroup gsd_mynn_edmf
!! This subroutine calculates the vertical diffusivity coefficients and the 
!! production terms for the turbulent quantities.      
!>\section gen_mym_turbulence GSD mym_turbulence General Algorithm
!! Two subroutines mym_level2() and mym_length() are called within this
!!subrouine to collect variable to carry out successive calculations:
!! - mym_level2() calculates the level 2 nondimensional wind shear \f$G_M\f$
!! and vertical temperature gradient \f$G_H\f$ as well as the level 2 stability
!! functions \f$S_h\f$ and \f$S_m\f$.
!! - mym_length() calculates the mixing lengths.
!! - The stability criteria from Helfand and Labraga (1989) are applied.
!! - The stability functions for level 2.5 or level 3.0 are calculated.
!! - If level 3.0 is used, counter-gradient terms are calculated.
!! - Production terms of TKE,\f$\theta^{'2}\f$,\f$q^{'2}\f$, and \f$\theta^{'}q^{'}\f$
!! are calculated.
!! - Eddy diffusivity \f$K_h\f$ and eddy viscosity \f$K_m\f$ are calculated.
!! - TKE budget terms are calculated (if the namelist parameter \p tke_budget 
!! is set to True)
  SUBROUTINE  mym_turbulence (                                &
    &            kts,kte,                                     &
    &            xland,closure,                               &
    &            dz, dx, zw,                                  &
    &            u, v, thl, thetav, ql, qw,                   &
    &            thlsg, qwsg,                                 &
    &            qke, tsq, qsq, cov,                          &
    &            vt, vq,                                      &
    &            rmo, flt, flq,                               &
    &            zi,theta,                                    &
    &            sh, sm,                                      &
    &            El,                                          &
    &            Dfm, Dfh, Dfq, Tcd, Qcd, Pdk, Pdt, Pdq, Pdc, &
    &		 qWT1D,qSHEAR1D,qBUOY1D,qDISS1D,              &
    &            tke_budget,                                  &
    &            Psig_bl,Psig_shcu,cldfra_bl1D,               &
    &            bl_mynn_mixlength,                           &
    &            edmf_w1,edmf_a1,                             &
    &            TKEprodTD,                                   &
    &            spp_pbl,rstoch_col,                          &
    &            qmin, zmax, CKmod,                           &
    &            a1,a2,b1,b2,c1,cc3,                           &
    &            debug_code,                                  &
    &            e1c,e2c,e3c,e4c,e5c                          )

!-------------------------------------------------------------------
!
    INTEGER, INTENT(IN)   :: kts,kte

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    INTEGER, INTENT(IN)   :: bl_mynn_mixlength,tke_budget
    REAL, INTENT(IN)      :: closure, qmin,zmax,CKmod
    REAL, INTENT(IN)      :: a1,a2,b1,b2,c1,cc3,e1c,e2c,e3c,e4c,e5c
    LOGICAL,INTENT(IN)    :: debug_code
    REAL, DIMENSION(kts:kte), INTENT(in) :: dz
    REAL, DIMENSION(kts:kte+1), INTENT(in) :: zw
    REAL, INTENT(in) :: rmo,flt,flq,Psig_bl,Psig_shcu,dx,xland,zi
    REAL, DIMENSION(kts:kte), INTENT(in) :: u,v,thl,thetav,qw,& 
         &ql,vt,vq,qke,tsq,qsq,cov,cldfra_bl1D,edmf_w1,edmf_a1,&
         &TKEprodTD,thlsg,qwsg

    REAL, DIMENSION(kts:kte), INTENT(out) :: dfm,dfh,dfq,&
         &pdk,pdt,pdq,pdc,tcd,qcd,el

    REAL, DIMENSION(kts:kte), INTENT(inout) :: &
         qWT1D,qSHEAR1D,qBUOY1D,qDISS1D
    REAL :: q3sq_old,dlsq1,qWTP_old,qWTP_new
    REAL :: dudz,dvdz,dTdz,&
            upwp,vpwp,Tpwp

    REAL, DIMENSION(kts:kte) :: qkw,dtl,dqw,dtv,gm,gh,sm,sh

    INTEGER :: k
!    REAL :: cc2,cc3,e1c,e2c,e3c,e4c,e5c
    REAL :: e6c,dzk,afk,abk,vtt,vqq,&
         &cw25,clow,cupp,gamt,gamq,smd,gamv,elq,elh

    REAL :: cldavg
    REAL, DIMENSION(kts:kte), INTENT(in) :: theta

    REAL ::  a2fac, duz, ri !JOE-Canuto/Kitamura mod

    REAL:: auh,aum,adh,adm,aeh,aem,Req,Rsl,Rsl2,&
           gmelq,sm20,sh20,sm25max,sh25max,sm25min,sh25min,&
           sm_pbl,sh_pbl,zi2,wt,slht,wtpr

    DOUBLE PRECISION  q2sq, t2sq, r2sq, c2sq, elsq, gmel, ghel
    DOUBLE PRECISION  q3sq, t3sq, r3sq, c3sq, dlsq, qdiv
    DOUBLE PRECISION  e1, e2, e3, e4, enum, eden, wden

!   Stochastic
    INTEGER,  INTENT(IN)                          ::    spp_pbl
    REAL, DIMENSION(KTS:KTE)                      ::    rstoch_col
    REAL :: Prnum, Prlim
    REAL, PARAMETER :: Prlimit = 5.0


!
!    tv0 = 0.61*tref
!    gtr = 9.81/tref
!
!    cc2 =  1.0-c2
!    cc3 =  1.0-c3
!    e1c =  3.0*a2*b2*cc3
!    e2c =  9.0*a1*a2*cc2
!    e3c =  9.0*a2*a2*cc2*( 1.0-c5 )
!    e4c = 12.0*a1*a2*cc2
!    e5c =  6.0*a1*a1
!

    CALL mym_level2 (kts,kte,                   &
    &            dz,                            &
    &            u, v, thl, thetav, qw,         &
    &            thlsg, qwsg,                   &
    &            ql, vt, vq,                    &
    &            dtl, dqw, dtv, gm, gh, sm, sh, &
    &            CKmod                          )
!
    CALL mym_length (                           &
    &            kts,kte,xland,                 &
    &            dz, dx, zw,                    &
    &            rmo, flt, flq,                 &
    &            vt, vq,                        &
    &            u, v, qke,                     &
    &            dtv,                           &
    &            el,                            &
    &            zi,theta,                      &
    &            qkw,Psig_bl,cldfra_bl1D,       &
    &            bl_mynn_mixlength,             &
    &            edmf_w1,edmf_a1,               &
    &            qmin,zmax                      )
!

    DO k = kts+1,kte
       dzk = 0.5  *( dz(k)+dz(k-1) )
       afk = dz(k)/( dz(k)+dz(k-1) )
       abk = 1.0 -afk
       elsq = el (k)**2
       q3sq = qkw(k)**2
       q2sq = b1*elsq*( sm(k)*gm(k)+sh(k)*gh(k) )

       sh20 = MAX(sh(k), 1e-5)
       sm20 = MAX(sm(k), 1e-5)
       sh(k)= MAX(sh(k), 1e-5)

       !Canuto/Kitamura mod
       duz = ( u(k)-u(k-1) )**2 +( v(k)-v(k-1) )**2
       duz =   duz                    /dzk**2
       !   **  Gradient Richardson number  **
       ri = -gh(k)/MAX( duz, 1.0e-10 )
       IF (CKmod .eq. 1) THEN
          a2fac = 1./(1. + MAX(ri,0.0))
       ELSE
          a2fac = 1.
       ENDIF
       !end Canuto/Kitamura mod

       !level 2.0 Prandtl number
       !Prnum = MIN(sm20/sh20, 4.0)
       !The form of Zilitinkevich et al. (2006) but modified
       !half-way towards Esau and Grachev (2007, Wind Eng)
       !Prnum = MIN(0.76 + 3.0*MAX(ri,0.0), Prlimit)
       Prnum = MIN(0.76 + 4.0*MAX(ri,0.0), Prlimit)
       !Prnum = MIN(0.76 + 5.0*MAX(ri,0.0), Prlimit)
!
!  Modified: Dec/22/2005, from here, (dlsq -> elsq)
       gmel = gm (k)*elsq
       ghel = gh (k)*elsq
!  Modified: Dec/22/2005, up to here

       ! Level 2.0 debug prints
       IF ( debug_code ) THEN
         IF (sh(k)<0.0 .OR. sm(k)<0.0) THEN
           print*,"MYNN; mym_turbulence 2.0; sh=",sh(k)," k=",k
           print*," gm=",gm(k)," gh=",gh(k)," sm=",sm(k)
           print*," q2sq=",q2sq," q3sq=",q3sq," q3/q2=",q3sq/q2sq
           print*," qke=",qke(k)," el=",el(k)," ri=",ri
           print*," PBLH=",zi," u=",u(k)," v=",v(k)
         ENDIF
       ENDIF

!     **  Since qkw is set to more than 0.0, q3sq > 0.0.  **

!     new stability criteria in level 2.5 (as well as level 3) - little/no impact
!     **  Limitation on q, instead of L/q  **
       dlsq =  elsq
       IF ( q3sq/dlsq .LT. -gh(k) ) q3sq = -dlsq*gh(k)

       IF ( q3sq .LT. q2sq ) THEN
          !Apply Helfand & Labraga mod
          qdiv = SQRT( q3sq/q2sq )   !HL89: (1-alfa)
!
          !Use level 2.5 stability functions
          !e1   = q3sq - e1c*ghel*a2fac
          !e2   = q3sq - e2c*ghel*a2fac
          !e3   = e1   + e3c*ghel*a2fac**2
          !e4   = e1   - e4c*ghel*a2fac
          !eden = e2*e4 + e3*e5c*gmel
          !eden = MAX( eden, 1.0d-20 )
          !sm(k) = q3sq*a1*( e3-3.0*c1*e4       )/eden
          !!JOE-Canuto/Kitamura mod
          !!sh(k) = q3sq*a2*( e2+3.0*c1*e5c*gmel )/eden
          !sh(k) = q3sq*(a2*a2fac)*( e2+3.0*c1*e5c*gmel )/eden
          !sm(k) = Prnum*sh(k)
          !sm(k) = sm(k) * qdiv

          !Use level 2.0 functions as in original MYNN
          sh(k) = sh(k) * qdiv
          sm(k) = sm(k) * qdiv
        !  !sm_pbl = sm(k) * qdiv
        !
        !  !Or, use the simple Pr relationship
        !  sm(k) = Prnum*sh(k)
        !
        !  !or blend them:
        !  zi2   = MAX(zi, 300.)
        !  wt    =.5*TANH((zw(k) - zi2)/200.) + .5
        !  sm(k) = sm_pbl*(1.-wt) + sm(k)*wt

          !Recalculate terms for later use
          !JOE-Canuto/Kitamura mod
          !e1   = q3sq - e1c*ghel * qdiv**2
          !e2   = q3sq - e2c*ghel * qdiv**2
          !e3   = e1   + e3c*ghel * qdiv**2
          !e4   = e1   - e4c*ghel * qdiv**2
          e1   = q3sq - e1c*ghel*a2fac * qdiv**2
          e2   = q3sq - e2c*ghel*a2fac * qdiv**2
          e3   = e1   + e3c*ghel*a2fac**2 * qdiv**2
          e4   = e1   - e4c*ghel*a2fac * qdiv**2
          eden = e2*e4 + e3*e5c*gmel * qdiv**2
          eden = MAX( eden, 1.0d-20 )
          !!JOE-Canuto/Kitamura mod
          !!sh(k) = q3sq*a2*( e2+3.0*c1*e5c*gmel )/eden  - retro 5
          !sh(k) = q3sq*(a2*a2fac)*( e2+3.0*c1*e5c*gmel )/eden
          !sm(k) = Prnum*sh(k)
       ELSE
          !JOE-Canuto/Kitamura mod
          !e1   = q3sq - e1c*ghel
          !e2   = q3sq - e2c*ghel
          !e3   = e1   + e3c*ghel
          !e4   = e1   - e4c*ghel
          e1   = q3sq - e1c*ghel*a2fac
          e2   = q3sq - e2c*ghel*a2fac
          e3   = e1   + e3c*ghel*a2fac**2
          e4   = e1   - e4c*ghel*a2fac
          eden = e2*e4 + e3*e5c*gmel
          eden = MAX( eden, 1.0d-20 )

          qdiv = 1.0
          !Use level 2.5 stability functions
          sm(k) = q3sq*a1*( e3-3.0*c1*e4       )/eden
        !  sm_pbl = q3sq*a1*( e3-3.0*c1*e4       )/eden
          !!JOE-Canuto/Kitamura mod
          !!sh(k) = q3sq*a2*( e2+3.0*c1*e5c*gmel )/eden
          sh(k) = q3sq*(a2*a2fac)*( e2+3.0*c1*e5c*gmel )/eden
        !  sm(k) = Prnum*sh(k)

        !  !or blend them:
        !  zi2   = MAX(zi, 300.)
        !  wt    = .5*TANH((zw(k) - zi2)/200.) + .5
        !  sm(k) = sm_pbl*(1.-wt) + sm(k)*wt
       END IF !end Helfand & Labraga check

       !Impose broad limits on Sh and Sm:
       gmelq    = MAX(gmel/q3sq, 1e-8)
       sm25max  = 4.  !MIN(sm20*3.0, SQRT(.1936/gmelq))
       sh25max  = 4.  !MIN(sh20*3.0, 0.76*b2)
       sm25min  = 0.0 !MAX(sm20*0.1, 1e-6)
       sh25min  = 0.0 !MAX(sh20*0.1, 1e-6)

       !JOE: Level 2.5 debug prints
       ! HL88 , lev2.5 criteria from eqs. 3.17, 3.19, & 3.20
       IF ( debug_code ) THEN
         IF ((sh(k)<sh25min .OR. sm(k)<sm25min .OR. &
              sh(k)>sh25max .OR. sm(k)>sm25max) ) THEN
           print*,"In mym_turbulence 2.5: k=",k
           print*," sm=",sm(k)," sh=",sh(k)
           print*," ri=",ri," Pr=",sm(k)/MAX(sh(k),1e-8)
           print*," gm=",gm(k)," gh=",gh(k)
           print*," q2sq=",q2sq," q3sq=",q3sq, q3sq/q2sq
           print*," qke=",qke(k)," el=",el(k)
           print*," PBLH=",zi," u=",u(k)," v=",v(k)
           print*," SMnum=",q3sq*a1*( e3-3.0*c1*e4)," SMdenom=",eden
           print*," SHnum=",q3sq*(a2*a2fac)*( e2+3.0*c1*e5c*gmel ),&
                  " SHdenom=",eden
         ENDIF
       ENDIF

       !Enforce constraints for level 2.5 functions
       IF ( sh(k) > sh25max ) sh(k) = sh25max
       IF ( sh(k) < sh25min ) sh(k) = sh25min
       !IF ( sm(k) > sm25max ) sm(k) = sm25max
       !IF ( sm(k) < sm25min ) sm(k) = sm25min
       !sm(k) = Prnum*sh(k)

       !surface layer PR
       !slht  = zi*0.1
       !wtpr  = min( max( (slht - zw(k))/slht, 0.0), 1.0) ! 1 at z=0, 0 above sfc layer
       !Prlim = 1.0*wtpr + (1.0 - wtpr)*Prlimit
       !Prlim = 2.0*wtpr + (1.0 - wtpr)*Prlimit
       !sm(k) = MIN(sm(k), Prlim*Sh(k))
       !Pending more testing, keep same Pr limit in sfc layer
       sm(k) = MIN(sm(k), Prlimit*Sh(k))

!   **  Level 3 : start  **
       IF ( closure .GE. 3.0 ) THEN
          t2sq = qdiv*b2*elsq*sh(k)*dtl(k)**2
          r2sq = qdiv*b2*elsq*sh(k)*dqw(k)**2
          c2sq = qdiv*b2*elsq*sh(k)*dtl(k)*dqw(k)
          t3sq = MAX( tsq(k)*abk+tsq(k-1)*afk, 0.0 )
          r3sq = MAX( qsq(k)*abk+qsq(k-1)*afk, 0.0 )
          c3sq =      cov(k)*abk+cov(k-1)*afk

!  Modified: Dec/22/2005, from here
          c3sq = SIGN( MIN( ABS(c3sq), SQRT(t3sq*r3sq) ), c3sq )
!
          vtt  = 1.0 +vt(k)*abk +vt(k-1)*afk
          vqq  = tv0 +vq(k)*abk +vq(k-1)*afk

          t2sq = vtt*t2sq +vqq*c2sq
          r2sq = vtt*c2sq +vqq*r2sq
          c2sq = MAX( vtt*t2sq+vqq*r2sq, 0.0d0 )
          t3sq = vtt*t3sq +vqq*c3sq
          r3sq = vtt*c3sq +vqq*r3sq
          c3sq = MAX( vtt*t3sq+vqq*r3sq, 0.0d0 )
!
          cw25 = e1*( e2 + 3.0*c1*e5c*gmel*qdiv**2 )/( 3.0*eden )
!
!     **  Limitation on q, instead of L/q  **
          dlsq =  elsq
          IF ( q3sq/dlsq .LT. -gh(k) ) q3sq = -dlsq*gh(k)
!
!     **  Limitation on c3sq (0.12 =< cw =< 0.76) **
          ! Use Janjic's (2001; p 13-17) methodology (eqs 4.11-414 and 5.7-5.10)
          ! to calculate an exact limit for c3sq:
          auh = 27.*a1*((a2*a2fac)**2)*b2*(gtr)**2
          aum = 54.*(a1**2)*(a2*a2fac)*b2*c1*(gtr)
          adh = 9.*a1*((a2*a2fac)**2)*(12.*a1 + 3.*b2)*(gtr)**2
          adm = 18.*(a1**2)*(a2*a2fac)*(b2 - 3.*(a2*a2fac))*(gtr)

          aeh = (9.*a1*((a2*a2fac)**2)*b1 +9.*a1*((a2*a2fac)**2)* &
                (12.*a1 + 3.*b2))*(gtr)
          aem = 3.*a1*(a2*a2fac)*b1*(3.*(a2*a2fac) + 3.*b2*c1 + &
                (18.*a1*c1 - b2)) + &
                (18.)*(a1**2)*(a2*a2fac)*(b2 - 3.*(a2*a2fac))

          Req = -aeh/aem
          Rsl = (auh + aum*Req)/(3.*adh + 3.*adm*Req)
          !For now, use default values, since tests showed little/no sensitivity
          Rsl = .12             !lower limit
          Rsl2= 1.0 - 2.*Rsl    !upper limit
          !IF (k==2)print*,"Dynamic limit RSL=",Rsl
          !IF (Rsl < 0.10 .OR. Rsl > 0.18) THEN
          !   print*,'--- ERROR: MYNN: Dynamic Cw '// &
          !        'limit exceeds reasonable limits'
          !   print*," MYNN: Dynamic Cw limit needs attention=",Rsl
          !ENDIF

          !JOE-Canuto/Kitamura mod
          !e2   = q3sq - e2c*ghel * qdiv**2
          !e3   = q3sq + e3c*ghel * qdiv**2
          !e4   = q3sq - e4c*ghel * qdiv**2
          e2   = q3sq - e2c*ghel*a2fac * qdiv**2
          e3   = q3sq + e3c*ghel*a2fac**2 * qdiv**2
          e4   = q3sq - e4c*ghel*a2fac * qdiv**2
          eden = e2*e4  + e3 *e5c*gmel * qdiv**2

          !JOE-Canuto/Kitamura mod
          !wden = cc3*gtr**2 * dlsq**2/elsq * qdiv**2 &
          !     &        *( e2*e4c - e3c*e5c*gmel * qdiv**2 )
          wden = cc3*gtr**2 * dlsq**2/elsq * qdiv**2 &
               &        *( e2*e4c*a2fac - e3c*e5c*gmel*a2fac**2 * qdiv**2 )

          IF ( wden .NE. 0.0 ) THEN
             !JOE: test dynamic limits
             clow = q3sq*( 0.12-cw25 )*eden/wden
             cupp = q3sq*( 0.76-cw25 )*eden/wden
             !clow = q3sq*( Rsl -cw25 )*eden/wden
             !cupp = q3sq*( Rsl2-cw25 )*eden/wden
!
             IF ( wden .GT. 0.0 ) THEN
                c3sq  = MIN( MAX( c3sq, c2sq+clow ), c2sq+cupp )
             ELSE
                c3sq  = MAX( MIN( c3sq, c2sq+clow ), c2sq+cupp )
             END IF
          END IF
!
          e1   = e2 + e5c*gmel * qdiv**2
          eden = MAX( eden, 1.0d-20 )
!  Modified: Dec/22/2005, up to here

          !JOE-Canuto/Kitamura mod
          !e6c  = 3.0*a2*cc3*gtr * dlsq/elsq
          e6c  = 3.0*(a2*a2fac)*cc3*gtr * dlsq/elsq

          !============================
          !     **  for Gamma_theta  **
          !!          enum = qdiv*e6c*( t3sq-t2sq )
          IF ( t2sq .GE. 0.0 ) THEN
             enum = MAX( qdiv*e6c*( t3sq-t2sq ), 0.0d0 )
          ELSE
             enum = MIN( qdiv*e6c*( t3sq-t2sq ), 0.0d0 )
          ENDIF
          gamt =-e1  *enum    /eden

          !============================
          !     **  for Gamma_q  **
          !!          enum = qdiv*e6c*( r3sq-r2sq )
          IF ( r2sq .GE. 0.0 ) THEN
             enum = MAX( qdiv*e6c*( r3sq-r2sq ), 0.0d0 )
          ELSE
             enum = MIN( qdiv*e6c*( r3sq-r2sq ), 0.0d0 )
          ENDIF
          gamq =-e1  *enum    /eden

          !============================
          !     **  for Sm' and Sh'd(Theta_V)/dz  **
          !!          enum = qdiv*e6c*( c3sq-c2sq )
          enum = MAX( qdiv*e6c*( c3sq-c2sq ), 0.0d0)

          !JOE-Canuto/Kitamura mod
          !smd  = dlsq*enum*gtr/eden * qdiv**2 * (e3c+e4c)*a1/a2
          smd  = dlsq*enum*gtr/eden * qdiv**2 * (e3c*a2fac**2 + &
               & e4c*a2fac)*a1/(a2*a2fac)

          gamv = e1  *enum*gtr/eden
          sm(k) = sm(k) +smd

          !============================
          !     **  For elh (see below), qdiv at Level 3 is reset to 1.0.  **
          qdiv = 1.0

          ! Level 3 debug prints
          IF ( debug_code ) THEN
            IF (sh(k)<-0.3 .OR. sm(k)<-0.3 .OR. &
              qke(k) < -0.1 .or. ABS(smd) .gt. 2.0) THEN
              print*," MYNN; mym_turbulence3.0; sh=",sh(k)," k=",k
              print*," gm=",gm(k)," gh=",gh(k)," sm=",sm(k)
              print*," q2sq=",q2sq," q3sq=",q3sq," q3/q2=",q3sq/q2sq
              print*," qke=",qke(k)," el=",el(k)," ri=",ri
              print*," PBLH=",zi," u=",u(k)," v=",v(k)
            ENDIF
          ENDIF

!   **  Level 3 : end  **

       ELSE
!     **  At Level 2.5, qdiv is not reset.  **
          gamt = 0.0
          gamq = 0.0
          gamv = 0.0
       END IF
!
!      Add min background stability function (diffusivity) within model levels
!      with active plumes and clouds.
       cldavg = 0.5*(cldfra_bl1D(k-1) + cldfra_bl1D(k))
       IF (edmf_a1(k) > 0.001 .OR. cldavg > 0.02) THEN
           ! for mass-flux columns
           sm(k) = MAX(sm(k), 0.03*MIN(10.*edmf_a1(k)*edmf_w1(k),1.0) )
           sh(k) = MAX(sh(k), 0.03*MIN(10.*edmf_a1(k)*edmf_w1(k),1.0) )
           ! for clouds
           sm(k) = MAX(sm(k), 0.05*MIN(cldavg,1.0) )
           sh(k) = MAX(sh(k), 0.05*MIN(cldavg,1.0) )
       ENDIF
!
       elq = el(k)*qkw(k)
       elh = elq*qdiv

       ! Production of TKE (pdk), T-variance (pdt),
       ! q-variance (pdq), and covariance (pdc)
       pdk(k) = elq*( sm(k)*gm(k)                &
            &        +sh(k)*gh(k)+gamv ) +       &
            &   TKEprodTD(k)
       pdt(k) = elh*( sh(k)*dtl(k)+gamt )*dtl(k)
       pdq(k) = elh*( sh(k)*dqw(k)+gamq )*dqw(k)
       pdc(k) = elh*( sh(k)*dtl(k)+gamt )        &
            &   *dqw(k)*0.5                      &
            & + elh*( sh(k)*dqw(k)+gamq )*dtl(k)*0.5

       ! Contergradient terms
       tcd(k) = elq*gamt
       qcd(k) = elq*gamq

       ! Eddy Diffusivity/Viscosity divided by dz
       dfm(k) = elq*sm(k) / dzk
       dfh(k) = elq*sh(k) / dzk
!  Modified: Dec/22/2005, from here
!   **  In sub.mym_predict, dfq for the TKE and scalar variance **
!   **  are set to 3.0*dfm and 1.0*dfm, respectively. (Sqfac)   **
       dfq(k) =     dfm(k)
!  Modified: Dec/22/2005, up to here

   IF (tke_budget .eq. 1) THEN
       !TKE BUDGET
!       dudz = ( u(k)-u(k-1) )/dzk
!       dvdz = ( v(k)-v(k-1) )/dzk
!       dTdz = ( thl(k)-thl(k-1) )/dzk

!       upwp = -elq*sm(k)*dudz
!       vpwp = -elq*sm(k)*dvdz
!       Tpwp = -elq*sh(k)*dTdz
!       Tpwp = SIGN(MAX(ABS(Tpwp),1.E-6),Tpwp)

       
!!  TKE budget  (Puhales, 2020, WRF 4.2.1)  << EOB   

       !!!Shear Term
       !!!qSHEAR1D(k)=-(upwp*dudz + vpwp*dvdz)
       qSHEAR1D(k) = elq*sm(k)*gm(k) !staggered

       !!!Buoyancy Term    
       !!!qBUOY1D(k)=grav*Tpwp/thl(k)
       !qBUOY1D(k)= elq*(sh(k)*gh(k) + gamv)
       !qBUOY1D(k) = elq*(sh(k)*(-dTdz*grav/thl(k)) + gamv) !! ORIGINAL CODE
       
       !! Buoyncy term takes the TKEprodTD(k) production now
       qBUOY1D(k) = elq*(sh(k)*gh(k)+gamv)+TKEprodTD(k) !staggered

       !!!Dissipation Term (now it evaluated on mym_predict)
       !qDISS1D(k) = (q3sq**(3./2.))/(b1*MAX(el(k),1.)) !! ORIGINAL CODE
       
       !! >> EOB
    ENDIF

    END DO
!

    dfm(kts) = 0.0
    dfh(kts) = 0.0
    dfq(kts) = 0.0
    tcd(kts) = 0.0
    qcd(kts) = 0.0

    tcd(kte) = 0.0
    qcd(kte) = 0.0

!
    DO k = kts,kte-1
       dzk = dz(k)
       tcd(k) = ( tcd(k+1)-tcd(k) )/( dzk )
       qcd(k) = ( qcd(k+1)-qcd(k) )/( dzk )
    END DO
!


    if (spp_pbl==1) then
       DO k = kts,kte
          dfm(k)= dfm(k) + dfm(k)* rstoch_col(k) * 1.5 * MAX(exp(-MAX(zw(k)-8000.,0.0)/2000.),0.001)
          dfh(k)= dfh(k) + dfh(k)* rstoch_col(k) * 1.5 * MAX(exp(-MAX(zw(k)-8000.,0.0)/2000.),0.001)
       END DO
    endif

!    RETURN
#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mym_turbulence
      

END MODULE
 
