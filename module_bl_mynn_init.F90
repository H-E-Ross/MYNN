MODULE mym_init

use module_bl_mynn_common,only: karman, twothirds
        
IMPLICIT NONE


CONTAINS
!=======================================================================
!     SUBROUTINE  mym_initialize:
!
!     Input variables:
!       iniflag         : <>0; turbulent quantities will be initialized
!                         = 0; turbulent quantities have been already
!                              given, i.e., they will not be initialized
!       nx, nz          : Dimension sizes of the
!                         x and z directions, respectively
!       tref            : Reference temperature                      (K)
!       dz(nz)          : Vertical grid spacings                     (m)
!                         # dz(nz)=dz(nz-1)
!       zw(nz+1)        : Heights of the walls of the grid boxes     (m)
!                         # zw(1)=0.0 and zw(k)=zw(k-1)+dz(k-1)
!       exner(nx,nz)    : Exner function at zw*h+zg             (J/kg K)
!                         defined by c_p*( p_basic/1000hPa )^kappa
!                         This is usually computed by integrating
!                         d(pi0)/dz = -h*g/tref.
!       rmo(nx)         : Inverse of the Obukhov length         (m^(-1))
!       flt, flq(nx)    : Turbulent fluxes of potential temperature and
!                         total water, respectively:
!                                    flt=-u_*Theta_*             (K m/s)
!                                    flq=-u_*qw_*            (kg/kg m/s)
!       ust(nx)         : Friction velocity                        (m/s)
!       pmz(nx)         : phi_m-zeta at z1*h+z0, where z1 (=0.5*dz(1))
!                         is the first grid point above the surafce, z0
!                         the roughness length and zeta=(z1*h+z0)*rmo
!       phh(nx)         : phi_h at z1*h+z0
!       u, v(nx,nz)     : Components of the horizontal wind        (m/s)
!       thl(nx,nz)      : Liquid water potential temperature
!                                                                    (K)
!       qw(nx,nz)       : Total water content Q_w                (kg/kg)
!
!     Output variables:
!       ql(nx,nz)       : Liquid water content                   (kg/kg)
!       vt, vq(nx,nz)   : Functions for computing the buoyancy flux
!       qke(nx,nz)      : Twice the turbulent kinetic energy q^2
!                                                              (m^2/s^2)
!       tsq(nx,nz)      : Variance of Theta_l                      (K^2)
!       qsq(nx,nz)      : Variance of Q_w
!       cov(nx,nz)      : Covariance of Theta_l and Q_w              (K)
!       el(nx,nz)       : Master length scale L                      (m)
!                         defined on the walls of the grid boxes
!
!     Work arrays:        see subroutine mym_level2
!       pd?(nx,nz,ny) : Half of the production terms at Level 2
!                         defined on the walls of the grid boxes
!       qkw(nx,nz,ny) : q on the walls of the grid boxes         (m/s)
!
!     # As to dtl, ...gh, see subroutine mym_turbulence.
!
!-------------------------------------------------------------------

!>\ingroup gsd_mynn_edmf
!! This subroutine initializes the mixing length, TKE, \f$\theta^{'2}\f$,
!! \f$q^{'2}\f$, and \f$\theta^{'}q^{'}\f$.
!!\section gen_mym_ini GSD MYNN-EDMF mym_initialize General Algorithm 
!> @{
  SUBROUTINE  mym_initialize (                                & 
       &            kts,kte,xland,                            &
       &            dz, dx, zw,                               &
       &            u, v, thl, qw,                            &
       &            thlsg, qwsg,                              &
!       &            ust, rmo, pmz, phh, flt, flq,             &
       &            zi, theta, thetav, sh, sm,                &
       &            ust, rmo, el,                             &
       &            Qke, Tsq, Qsq, Cov, Psig_bl, cldfra_bl1D, &
       &            bl_mynn_mixlength,                        &
       &            edmf_w1,edmf_a1,                          &
       &            INITIALIZE_QKE,                           &
       &            spp_pbl,rstoch_col,                       &
       &            b1, b2, qkemin, CKmod                     )
!
!-------------------------------------------------------------------
    
    INTEGER, INTENT(IN)   :: kts,kte
    INTEGER, INTENT(IN)   :: bl_mynn_mixlength
    LOGICAL, INTENT(IN)   :: INITIALIZE_QKE
!    REAL, INTENT(IN)   :: ust, rmo, pmz, phh, flt, flq
    REAL, INTENT(IN)   :: ust, rmo, Psig_bl, dx, xland, b1, b2, qkemin, CKmod
    REAL, DIMENSION(kts:kte), INTENT(in) :: dz
    REAL, DIMENSION(kts:kte+1), INTENT(in) :: zw
    REAL, DIMENSION(kts:kte), INTENT(in) :: u,v,thl,qw,cldfra_bl1D,&
                                          edmf_w1,edmf_a1
    REAL, DIMENSION(kts:kte), INTENT(out) :: tsq,qsq,cov
    REAL, DIMENSION(kts:kte), INTENT(inout) :: el,qke
    REAL, DIMENSION(kts:kte) :: &
         &ql,pdk,pdt,pdq,pdc,dtl,dqw,dtv,&
         &gm,gh,sm,sh,qkw,vt,vq
    INTEGER :: k,l,lmax
    REAL :: phm,vkz,elq,elv,b1l,b2l,pmz=1.,phh=1.,flt=0.,flq=0.,tmpq
    REAL :: zi
    REAL, DIMENSION(kts:kte) :: theta,thetav,thlsg,qwsg

    REAL, DIMENSION(kts:kte) :: rstoch_col
    INTEGER ::spp_pbl

!> - At first ql, vt and vq are set to zero.
    DO k = kts,kte
       ql(k) = 0.0
       vt(k) = 0.0
       vq(k) = 0.0
    END DO
!
!> - Call mym_level2() to calculate the stability functions at level 2.
    CALL mym_level2 ( kts,kte,                      &
         &            dz,                           &
         &            u, v, thl, thetav, qw,        &
         &            thlsg, qwsg,                  &
         &            ql, vt, vq,                   &
         &            dtl, dqw, dtv, gm, gh, sm, sh,&
         &            CKmod                         )
!
!   **  Preliminary setting  **

    el (kts) = 0.0
    IF (INITIALIZE_QKE) THEN
       !qke(kts) = ust**2 * ( b1*pmz )**(2.0/3.0)
       qke(kts) = 1.5 * ust**2 * ( b1*pmz )**(2.0/3.0)
       DO k = kts+1,kte
          !qke(k) = 0.0
          !linearly taper off towards top of pbl
          qke(k)=qke(kts)*MAX((ust*700. - zw(k))/(MAX(ust,0.01)*700.), 0.01)
       ENDDO
    ENDIF
!
    phm      = phh*b2 / ( b1*pmz )**(1.0/3.0)
    tsq(kts) = phm*( flt/ust )**2
    qsq(kts) = phm*( flq/ust )**2
    cov(kts) = phm*( flt/ust )*( flq/ust )
!
    DO k = kts+1,kte
       vkz = karman*zw(k)
       el (k) = vkz/( 1.0 + vkz/100.0 )
!       qke(k) = 0.0
!
       tsq(k) = 0.0
       qsq(k) = 0.0
       cov(k) = 0.0
    END DO
!
!   **  Initialization with an iterative manner          **
!   **  lmax is the iteration count. This is arbitrary.  **
    lmax = 5
!
    DO l = 1,lmax
!
!> - call mym_length() to calculate the master length scale.
       CALL mym_length (                          &
            &            kts,kte,xland,           &
            &            dz, dx, zw,              &
            &            rmo, flt, flq,           &
            &            vt, vq,                  &
            &            u, v, qke,               &
            &            dtv,                     &
            &            el,                      &
            &            zi,theta,                &
            &            qkw,Psig_bl,cldfra_bl1D, &
            &            bl_mynn_mixlength,       &
            &            edmf_w1,edmf_a1          )
!
       DO k = kts+1,kte
          elq = el(k)*qkw(k)
          pdk(k) = elq*( sm(k)*gm(k) + &
               &         sh(k)*gh(k) )
          pdt(k) = elq*  sh(k)*dtl(k)**2
          pdq(k) = elq*  sh(k)*dqw(k)**2
          pdc(k) = elq*  sh(k)*dtl(k)*dqw(k)
       END DO
!
!   **  Strictly, vkz*h(i,j) -> karman*( 0.5*dz(1)*h(i,j)+z0 )  **
       vkz = karman*0.5*dz(kts)
       elv = 0.5*( el(kts+1)+el(kts) ) /  vkz
       IF (INITIALIZE_QKE)THEN 
          !qke(kts) = ust**2 * ( b1*pmz*elv    )**(2.0/3.0)
          qke(kts) = 1.0 * MAX(ust,0.02)**2 * ( b1*pmz*elv    )**(2.0/3.0) 
       ENDIF

       phm      = phh*b2 / ( b1*pmz/elv**2 )**(1.0/3.0)
       tsq(kts) = phm*( flt/ust )**2
       qsq(kts) = phm*( flq/ust )**2
       cov(kts) = phm*( flt/ust )*( flq/ust )

       DO k = kts+1,kte-1
          b1l = b1*0.25*( el(k+1)+el(k) )
          !tmpq=MAX(b1l*( pdk(k+1)+pdk(k) ),qkemin)
          !add MIN to limit unreasonable QKE
          tmpq=MIN(MAX(b1l*( pdk(k+1)+pdk(k) ),qkemin),125.)
!          PRINT *,'tmpqqqqq',tmpq,pdk(k+1),pdk(k)
          IF (INITIALIZE_QKE)THEN
             qke(k) = tmpq**twothirds
          ENDIF

          IF ( qke(k) .LE. 0.0 ) THEN
             b2l = 0.0
          ELSE
             b2l = b2*( b1l/b1 ) / SQRT( qke(k) )
          END IF

          tsq(k) = b2l*( pdt(k+1)+pdt(k) )
          qsq(k) = b2l*( pdq(k+1)+pdq(k) )
          cov(k) = b2l*( pdc(k+1)+pdc(k) )
       END DO

    END DO

!!    qke(kts)=qke(kts+1)
!!    tsq(kts)=tsq(kts+1)
!!    qsq(kts)=qsq(kts+1)
!!    cov(kts)=cov(kts+1)

    IF (INITIALIZE_QKE)THEN
       qke(kts)=0.5*(qke(kts)+qke(kts+1))
       qke(kte)=qke(kte-1)
    ENDIF
    tsq(kte)=tsq(kte-1)
    qsq(kte)=qsq(kte-1)
    cov(kte)=cov(kte-1)

!
!    RETURN

  END SUBROUTINE mym_initialize


END MODULE

