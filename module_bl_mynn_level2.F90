MODULE mynn_level2

use module_bl_mynn_common,only: gtr,tv0
IMPLICIT NONE


! Closure constants
  REAL, PARAMETER ::  &
       &pr  =  0.74,  &
       &g1  =  0.235, &  ! NN2009 = 0.235
       &b1  = 24.0,   &
       &b2  = 15.0,   &  ! CKmod     NN2009
       &c2  =  0.729, &  ! 0.729, & !0.75, &
       &c3  =  0.340, &  ! 0.340, & !0.352, &
       &c4  =  0.0,   &
       &c5  =  0.2,   &
       &a1  = b1*( 1.0-3.0*g1 )/6.0, &
!       &c1  = g1 -1.0/( 3.0*a1*b1**(1.0/3.0) ), &
       &c1  = g1 -1.0/( 3.0*a1*2.88449914061481660), &
       &a2  = a1*( g1-c1 )/( g1*pr ), &
       &g2  = b2/b1*( 1.0-c3 ) +2.0*a1/b1*( 3.0-2.0*c2 )
       
  REAL, PARAMETER :: &
       &cc2 =  1.0-c2, &
       &cc3 =  1.0-c3, &
       &e1c =  3.0*a2*b2*cc3, &
       &e2c =  9.0*a1*a2*cc2, &
       &e3c =  9.0*a2*a2*cc2*( 1.0-c5 ), &
       &e4c = 12.0*a1*a2*cc2, &
       &e5c =  6.0*a1*a1 


CONTAINS
!
! ==================================================================
!     SUBROUTINE  mym_level2:
!
!     Input variables:    see subroutine mym_initialize
!
!     Output variables:
!       dtl(nx,nz,ny) : Vertical gradient of Theta_l             (K/m)
!       dqw(nx,nz,ny) : Vertical gradient of Q_w
!       dtv(nx,nz,ny) : Vertical gradient of Theta_V             (K/m)
!       gm (nx,nz,ny) : G_M divided by L^2/q^2                (s^(-2))
!       gh (nx,nz,ny) : G_H divided by L^2/q^2                (s^(-2))
!       sm (nx,nz,ny) : Stability function for momentum, at Level 2
!       sh (nx,nz,ny) : Stability function for heat, at Level 2
!
!       These are defined on the walls of the grid boxes.
!

!>\ingroup gsd_mynn_edmf
!! This subroutine calculates the level 2, non-dimensional wind shear
!! \f$G_M\f$ and vertical temperature gradient \f$G_H\f$ as well as 
!! the level 2 stability funcitons \f$S_h\f$ and \f$S_m\f$.
!!\param kts    horizontal dimension
!!\param kte    vertical dimension
!!\param dz     vertical grid spacings (\f$m\f$)
!!\param u      west-east component of the horizontal wind (\f$m s^{-1}\f$)
!!\param v      south-north component of the horizontal wind (\f$m s^{-1}\f$)
!!\param thl    liquid water potential temperature
!!\param qw     total water content \f$Q_w\f$
!!\param ql     liquid water content (\f$kg kg^{-1}\f$)
!!\param vt
!!\param vq
!!\param dtl     vertical gradient of \f$\theta_l\f$ (\f$K m^{-1}\f$)
!!\param dqw     vertical gradient of \f$Q_w\f$
!!\param dtv     vertical gradient of \f$\theta_V\f$ (\f$K m^{-1}\f$)
!!\param gm      \f$G_M\f$ divided by \f$L^{2}/q^{2}\f$ (\f$s^{-2}\f$)
!!\param gh      \f$G_H\f$ divided by \f$L^{2}/q^{2}\f$ (\f$s^{-2}\f$)
!!\param sm      stability function for momentum, at Level 2
!!\param sh      stability function for heat, at Level 2
!!\section gen_mym_level2 GSD MYNN-EDMF mym_level2 General Algorithm
!! @ {
  SUBROUTINE  mym_level2 (kts,kte,                &
       &            dz,                           &
       &            u, v, thl, thetav, qw,        &
       &            thlsg, qwsg,                  &
       &            ql, vt, vq,                   &
       &            dtl, dqw, dtv, gm, gh, sm, sh,&
       &            cKmod                         )
!
!-------------------------------------------------------------------

    INTEGER, INTENT(IN)   :: kts,kte
    REAL, INTENT(IN) :: CKmod
#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    REAL, DIMENSION(kts:kte), INTENT(in) :: dz
    REAL, DIMENSION(kts:kte), INTENT(in) :: u,v,thl,qw,ql,vt,vq,&
                                            thetav,thlsg,qwsg
    REAL, DIMENSION(kts:kte), INTENT(out) :: &
         &dtl,dqw,dtv,gm,gh,sm,sh

    INTEGER :: k

    REAL :: rfc,f1,f2,rf1,rf2,smc,shc,&
         &ri1,ri2,ri3,ri4,duz,dtz,dqz,vtt,vqq,dtq,dzk,afk,abk,ri,rf

    REAL ::   a2fac

!    ev  = 2.5e6
!    tv0 = 0.61*tref
!    tv1 = 1.61*tref
!    gtr = 9.81/tref
!
    rfc = g1/( g1+g2 )
    f1  = b1*( g1-c1 ) +3.0*a2*( 1.0    -c2 )*( 1.0-c5 ) &
    &                   +2.0*a1*( 3.0-2.0*c2 )
    f2  = b1*( g1+g2 ) -3.0*a1*( 1.0    -c2 )
    rf1 = b1*( g1-c1 )/f1
    rf2 = b1*  g1     /f2
    smc = a1 /a2*  f1/f2
    shc = 3.0*a2*( g1+g2 )
!
    ri1 = 0.5/smc
    ri2 = rf1*smc
    ri3 = 4.0*rf2*smc -2.0*ri2
    ri4 = ri2**2
!
    DO k = kts+1,kte
       dzk = 0.5  *( dz(k)+dz(k-1) )
       afk = dz(k)/( dz(k)+dz(k-1) )
       abk = 1.0 -afk
       duz = ( u(k)-u(k-1) )**2 +( v(k)-v(k-1) )**2
       duz =   duz                    /dzk**2
       dtz = ( thl(k)-thl(k-1) )/( dzk )
       !Alternatively, use SGS clouds for thl
       !dtz = ( thlsg(k)-thlsg(k-1) )/( dzk )
       dqz = ( qw(k)-qw(k-1) )/( dzk )
       !Alternatively, use SGS clouds for qw
       !dqz = ( qwsg(k)-qwsg(k-1) )/( dzk )
!
       vtt =  1.0 +vt(k)*abk +vt(k-1)*afk  ! Beta-theta in NN09, Eq. 39
       vqq =  tv0 +vq(k)*abk +vq(k-1)*afk  ! Beta-q
       dtq =  vtt*dtz +vqq*dqz
       !Alternatively, use theta-v without the SGS clouds
       !dtq = ( thetav(k)-thetav(k-1) )/( dzk )
!
       dtl(k) =  dtz
       dqw(k) =  dqz
       dtv(k) =  dtq
!?      dtv(i,j,k) =  dtz +tv0*dqz
!?   :              +( xlv/pi0(i,j,k)-tv1 )
!?   :              *( ql(i,j,k)-ql(i,j,k-1) )/( dzk*h(i,j) )
!
       gm (k) =  duz
       gh (k) = -dtq*gtr
!
!   **  Gradient Richardson number  **
       ri = -gh(k)/MAX( duz, 1.0e-10 )

    !a2fac is needed for the Canuto/Kitamura mod
    IF (CKmod .eq. 1) THEN
       a2fac = 1./(1. + MAX(ri,0.0))
    ELSE
       a2fac = 1.
    ENDIF

       rfc = g1/( g1+g2 )
       f1  = b1*( g1-c1 ) +3.0*a2*a2fac *( 1.0    -c2 )*( 1.0-c5 ) &
    &                     +2.0*a1*( 3.0-2.0*c2 )
       f2  = b1*( g1+g2 ) -3.0*a1*( 1.0    -c2 )
       rf1 = b1*( g1-c1 )/f1
       rf2 = b1*  g1     /f2
       smc = a1 /(a2*a2fac)*  f1/f2
       shc = 3.0*(a2*a2fac)*( g1+g2 )

       ri1 = 0.5/smc
       ri2 = rf1*smc
       ri3 = 4.0*rf2*smc -2.0*ri2
       ri4 = ri2**2

!   **  Flux Richardson number  **
       rf = MIN( ri1*( ri + ri2-SQRT(ri**2 - ri3*ri + ri4) ), rfc )
!
       sh (k) = shc*( rfc-rf )/( 1.0-rf )
       sm (k) = smc*( rf1-rf )/( rf2-rf ) * sh(k)
    END DO
!
!    RETURN

#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mym_level2


END MODULE
 
