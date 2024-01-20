MODULE module_mym_length
use module_bl_mynn_common,only: grav, gtr,karman,onethird,tv0        
IMPLICIT NONE
CONTAINS


! ==================================================================
!     SUBROUTINE  mym_length:
!
!     Input variables:    see subroutine mym_initialize
!
!     Output variables:   see subroutine mym_initialize
!
!     Work arrays:
!       elt(nx,ny)      : Length scale depending on the PBL depth    (m)
!       vsc(nx,ny)      : Velocity scale q_c                       (m/s)
!                         at first, used for computing elt
!
!     NOTE: the mixing lengths are meant to be calculated at the full-
!           sigmal levels (or interfaces beween the model layers).
!
!>\ingroup gsd_mynn_edmf
!! This subroutine calculates the mixing lengths.
  SUBROUTINE  mym_length (                     & 
    &            kts,kte,xland,                &
    &            dz, dx, zw,                   &
    &            rmo, flt, flq,                &
    &            vt, vq,                       &
    &            u1, v1, qke,                  &
    &            dtv,                          &
    &            el,                           &
    &            zi, theta, qkw,               &
    &            Psig_bl, cldfra_bl1D,         &
    &            bl_mynn_mixlength,            &
    &            edmf_w1,edmf_a1,              &
    &            qmin, zmax    )
    
!-------------------------------------------------------------------

    INTEGER, INTENT(IN)   :: kts,kte

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    INTEGER, INTENT(IN)   :: bl_mynn_mixlength
    REAL, DIMENSION(kts:kte), INTENT(in)   :: dz
    REAL, DIMENSION(kts:kte+1), INTENT(in) :: zw
    REAL, INTENT(in) :: rmo,flt,flq,Psig_bl,dx,xland,qmin,zmax
    REAL, DIMENSION(kts:kte), INTENT(IN)   :: u1,v1,qke,vt,vq,cldfra_bl1D,&
                                          edmf_w1,edmf_a1
    REAL, DIMENSION(kts:kte), INTENT(out)  :: qkw, el
    REAL, DIMENSION(kts:kte), INTENT(in)   :: dtv

    REAL :: elt,vsc

    REAL, DIMENSION(kts:kte), INTENT(IN) :: theta
    REAL, DIMENSION(kts:kte) :: qtke,elBLmin,elBLavg,thetaw
    REAL :: wt,wt2,zi,zi2,h1,h2,hs,elBLmin0,elBLavg0,cldavg

    ! THE FOLLOWING CONSTANTS ARE IMPORTANT FOR REGULATING THE
    ! MIXING LENGTHS:
    REAL :: cns,   &   !< for surface layer (els) in stable conditions
            alp1,  &   !< for turbulent length scale (elt)
            alp2,  &   !< for buoyancy length scale (elb)
            alp3,  &   !< for buoyancy enhancement factor of elb
            alp4,  &   !< for surface layer (els) in unstable conditions
            alp5,  &   !< for BouLac mixing length or above PBLH
            alp6       !< for mass-flux/

    !THE FOLLOWING LIMITS DO NOT DIRECTLY AFFECT THE ACTUAL PBLH.
    !THEY ONLY IMPOSE LIMITS ON THE CALCULATION OF THE MIXING LENGTH 
    !SCALES SO THAT THE BOULAC MIXING LENGTH (IN FREE ATMOS) DOES
    !NOT ENCROACH UPON THE BOUNDARY LAYER MIXING LENGTH (els, elb & elt).
    REAL, PARAMETER :: minzi = 300.  !< min mixed-layer height
    REAL, PARAMETER :: maxdz = 750.  !< max (half) transition layer depth
                                     !! =0.3*2500 m PBLH, so the transition
                                     !! layer stops growing for PBLHs > 2.5 km.
    REAL, PARAMETER :: mindz = 300.  !< 300  !min (half) transition layer depth

    !SURFACE LAYER LENGTH SCALE MODS TO REDUCE IMPACT IN UPPER BOUNDARY LAYER
    REAL, PARAMETER :: ZSLH = 100. !< Max height correlated to surface conditions (m)
    REAL, PARAMETER :: CSL = 2.    !< CSL = constant of proportionality to L O(1)


    INTEGER :: i,j,k
    REAL :: afk,abk,zwk,zwk1,dzk,qdz,vflx,bv,tau_cloud,wstar,elb,els, &
            & elf,el_stab,el_mf,el_stab_mf,elb_mf,                    &
            & PBLH_PLUS_ENT,Uonset,Ugrid,wt_u,el_les
    REAL, PARAMETER :: ctau = 1000. !constant for tau_cloud

!    tv0 = 0.61*tref
!    gtr = 9.81/tref

    SELECT CASE(bl_mynn_mixlength)

      CASE (0) ! ORIGINAL MYNN MIXING LENGTH + BouLac

        cns  = 2.7
        alp1 = 0.23
        alp2 = 1.0
        alp3 = 5.0
        alp4 = 100.
        alp5 = 0.3

        ! Impose limits on the height integration for elt and the transition layer depth
        zi2  = MIN(10000.,zw(kte-2))  !originally integrated to model top, not just 10 km.
        h1=MAX(0.3*zi2,mindz)
        h1=MIN(h1,maxdz)         ! 1/2 transition layer depth
        h2=h1/2.0                ! 1/4 transition layer depth

        qkw(kts) = SQRT(MAX(qke(kts),1.0e-10))
        DO k = kts+1,kte
           afk = dz(k)/( dz(k)+dz(k-1) )
           abk = 1.0 -afk
           qkw(k) = SQRT(MAX(qke(k)*abk+qke(k-1)*afk,1.0e-3))
        END DO

        elt = 1.0e-5
        vsc = 1.0e-5        

        !   **  Strictly, zwk*h(i,j) -> ( zwk*h(i,j)+z0 )  **
        k = kts+1
        zwk = zw(k)
        DO WHILE (zwk .LE. zi2+h1)
           dzk = 0.5*( dz(k)+dz(k-1) )
           qdz = MAX( qkw(k)-qmin, 0.03 )*dzk
           elt = elt +qdz*zwk
           vsc = vsc +qdz
           k   = k+1
           zwk = zw(k)
        END DO

        elt =  alp1*elt/vsc
        vflx = ( vt(kts)+1.0 )*flt +( vq(kts)+tv0 )*flq
        vsc = ( gtr*elt*MAX( vflx, 0.0 ) )**(1.0/3.0)

        !   **  Strictly, el(i,k=1) is not zero.  **
        el(kts) = 0.0
        zwk1    = zw(kts+1)

        DO k = kts+1,kte
           zwk = zw(k)              !full-sigma levels

           !   **  Length scale limited by the buoyancy effect  **
           IF ( dtv(k) .GT. 0.0 ) THEN
              bv  = SQRT( gtr*dtv(k) )
              elb = alp2*qkw(k) / bv &
                  &       *( 1.0 + alp3/alp2*&
                  &SQRT( vsc/( bv*elt ) ) )
              elf = alp2 * qkw(k)/bv

           ELSE
              elb = 1.0e10
              elf = elb
           ENDIF

           !   **  Length scale in the surface layer  **
           IF ( rmo .GT. 0.0 ) THEN
              els  = karman*zwk/(1.0+cns*MIN( zwk*rmo, zmax ))
           ELSE
              els  =  karman*zwk*( 1.0 - alp4* zwk*rmo )**0.2
           END IF

           !   ** HARMONC AVERGING OF MIXING LENGTH SCALES:
           !       el(k) =      MIN(elb/( elb/elt+elb/els+1.0 ),elf)
           !       el(k) =      elb/( elb/elt+elb/els+1.0 )

           wt=.5*TANH((zwk - (zi2+h1))/h2) + .5

           el(k) = MIN(elb/( elb/elt+elb/els+1.0 ),elf)

        END DO

      CASE (1) !NONLOCAL (using BouLac) FORM OF MIXING LENGTH

        ugrid = sqrt(u1(kts)**2 + v1(kts)**2)
        uonset= 15. 
        wt_u = (1.0 - min(max(ugrid - uonset, 0.0)/30.0, 0.5)) 
        cns  = 3.5
        alp1 = 0.22 !was 0.21
        alp2 = 0.25 !was 0.3
        alp3 = 2.0 * wt_u !taper off bouyancy enhancement in shear-driven pbls
        alp4 = 5.0
        alp5 = 0.3
        alp6 = 50.

        ! Impose limits on the height integration for elt and the transition layer depth
        zi2=MAX(zi,200.) !minzi)
        h1=MAX(0.3*zi2,200.)
        h1=MIN(h1,500.)          ! 1/2 transition layer depth
        h2=h1/2.0                ! 1/4 transition layer depth

        qtke(kts)=MAX(0.5*qke(kts), 0.01) !tke at full sigma levels
        thetaw(kts)=theta(kts)            !theta at full-sigma levels
        qkw(kts) = SQRT(MAX(qke(kts),1.0e-10))

        DO k = kts+1,kte
           afk = dz(k)/( dz(k)+dz(k-1) )
           abk = 1.0 -afk
           qkw(k) = SQRT(MAX(qke(k)*abk+qke(k-1)*afk,1.0e-3))
           qtke(k) = 0.5*(qkw(k)**2)     ! q -> TKE
           thetaw(k)= theta(k)*abk + theta(k-1)*afk
        END DO

        elt = 1.0e-5
        vsc = 1.0e-5

        !   **  Strictly, zwk*h(i,j) -> ( zwk*h(i,j)+z0 )  **
        k = kts+1
        zwk = zw(k)
        DO WHILE (zwk .LE. zi2+h1)
           dzk = 0.5*( dz(k)+dz(k-1) )
           qdz = min(max( qkw(k)-qmin, 0.03 ), 30.0)*dzk
           elt = elt +qdz*zwk
           vsc = vsc +qdz
           k   = k+1
           zwk = zw(k)
        END DO

        elt = MIN( MAX( alp1*elt/vsc, 10.), 400.)
        !avoid use of buoyancy flux functions which are ill-defined at the surface
        !vflx = ( vt(kts)+1.0 )*flt + ( vq(kts)+tv0 )*flq
        vflx = flt
        vsc = ( gtr*elt*MAX( vflx, 0.0 ) )**onethird

        !   **  Strictly, el(i,j,1) is not zero.  **
        el(kts) = 0.0
        zwk1    = zw(kts+1)              !full-sigma levels

        ! COMPUTE BouLac mixing length
        CALL boulac_length(kts,kte,zw,dz,qtke,thetaw,elBLmin,elBLavg)

        DO k = kts+1,kte
           zwk = zw(k)              !full-sigma levels

           !   **  Length scale limited by the buoyancy effect  **
           IF ( dtv(k) .GT. 0.0 ) THEN
              bv  = max( sqrt( gtr*dtv(k) ), 0.001)
              elb = MAX(alp2*qkw(k),                          &
                  &    alp6*edmf_a1(k-1)*edmf_w1(k-1)) / bv   &
                  &  *( 1.0 + alp3*SQRT( vsc/(bv*elt) ) )
              elb = MIN(elb, zwk)
              elf = 0.65 * qkw(k)/bv
              elBLavg(k) = MAX(elBLavg(k), alp6*edmf_a1(k-1)*edmf_w1(k-1)/bv)
           ELSE
              elb = 1.0e10
              elf = elb
           ENDIF

           !   **  Length scale in the surface layer  **
           IF ( rmo .GT. 0.0 ) THEN
              els  = karman*zwk/(1.0+cns*MIN( zwk*rmo, zmax ))
           ELSE
              els  =  karman*zwk*( 1.0 - alp4* zwk*rmo )**0.2
           END IF

           !   ** NOW BLEND THE MIXING LENGTH SCALES:
           wt=.5*TANH((zwk - (zi2+h1))/h2) + .5

           !add blending to use BouLac mixing length in free atmos;
           !defined relative to the PBLH (zi) + transition layer (h1)
           !el(k) = MIN(elb/( elb/elt+elb/els+1.0 ),elf)
           !try squared-blending
           el(k) = SQRT( els**2/(1. + (els**2/elt**2) +(els**2/elb**2)))
           el(k) = MIN (el(k), elf)
           el(k) = el(k)*(1.-wt) + alp5*elBLavg(k)*wt

           ! include scale-awareness, except for original MYNN
           el(k) = el(k)*Psig_bl

         END DO

      CASE (2) !Local (mostly) mixing length formulation

        Uonset = 3.5 + dz(kts)*0.1
        Ugrid  = sqrt(u1(kts)**2 + v1(kts)**2)
        cns  = 3.5 !JOE-test  * (1.0 - MIN(MAX(Ugrid - Uonset, 0.0)/10.0, 1.0))
        alp1 = 0.22 !0.21
        alp2 = 0.25 !0.30
        alp3 = 2.0  !1.5
        alp4 = 5.0
        alp5 = alp2 !like alp2, but for free atmosphere
        alp6 = 50.0 !used for MF mixing length

        ! Impose limits on the height integration for elt and the transition layer depth
        !zi2=MAX(zi,minzi)
        zi2=MAX(zi,    200.)
        !h1=MAX(0.3*zi2,mindz)
        !h1=MIN(h1,maxdz)         ! 1/2 transition layer depth
        h1=MAX(0.3*zi2,200.)
        h1=MIN(h1,500.)
        h2=h1*0.5                ! 1/4 transition layer depth

        qtke(kts)=MAX(0.5*qke(kts),0.01) !tke at full sigma levels
        qkw(kts) = SQRT(MAX(qke(kts),1.0e-4))

        DO k = kts+1,kte
           afk = dz(k)/( dz(k)+dz(k-1) )
           abk = 1.0 -afk
           qkw(k) = SQRT(MAX(qke(k)*abk+qke(k-1)*afk,1.0e-3))
           qtke(k) = 0.5*qkw(k)**2  ! qkw -> TKE
        END DO

        elt = 1.0e-5
        vsc = 1.0e-5

        !   **  Strictly, zwk*h(i,j) -> ( zwk*h(i,j)+z0 )  **
        PBLH_PLUS_ENT = MAX(zi+h1, 100.)
        k = kts+1
        zwk = zw(k)
        DO WHILE (zwk .LE. PBLH_PLUS_ENT)
           dzk = 0.5*( dz(k)+dz(k-1) )
           qdz = min(max( qkw(k)-qmin, 0.03 ), 30.0)*dzk
           elt = elt +qdz*zwk
           vsc = vsc +qdz
           k   = k+1
           zwk = zw(k)
        END DO

        elt = MIN( MAX(alp1*elt/vsc, 10.), 400.)
        !avoid use of buoyancy flux functions which are ill-defined at the surface
        !vflx = ( vt(kts)+1.0 )*flt +( vq(kts)+tv0 )*flq
        vflx = flt
        vsc = ( gtr*elt*MAX( vflx, 0.0 ) )**onethird

        !   **  Strictly, el(i,j,1) is not zero.  **
        el(kts) = 0.0
        zwk1    = zw(kts+1)

        DO k = kts+1,kte
           zwk = zw(k)              !full-sigma levels
           dzk = 0.5*( dz(k)+dz(k-1) )
           cldavg = 0.5*(cldfra_bl1D(k-1)+cldfra_bl1D(k))

           !   **  Length scale limited by the buoyancy effect  **
           IF ( dtv(k) .GT. 0.0 ) THEN
              !impose min value on bv
              bv  = MAX( SQRT( gtr*dtv(k) ), 0.001)  
              !elb_mf = alp2*qkw(k) / bv  &
              elb_mf = MAX(alp2*qkw(k),                    &
                  &    alp6*edmf_a1(k-1)*edmf_w1(k-1)) / bv    &
                  &  *( 1.0 + alp3*SQRT( vsc/( bv*elt ) ) )
              elb = MIN(MAX(alp5*qkw(k), alp6*edmf_a1(k)*edmf_w1(k))/bv, zwk)

              !tau_cloud = MIN(MAX(0.5*zi/((gtr*zi*MAX(vflx,1.0e-4))**onethird),30.),150.)
              wstar = 1.25*(gtr*zi*MAX(vflx,1.0e-4))**onethird
              tau_cloud = MIN(MAX(ctau * wstar/grav, 30.), 150.)
              !minimize influence of surface heat flux on tau far away from the PBLH.
              wt=.5*TANH((zwk - (zi2+h1))/h2) + .5
              tau_cloud = tau_cloud*(1.-wt) + 50.*wt
              elf = MIN(MAX(tau_cloud*SQRT(MIN(qtke(k),40.)), &
                  &         alp6*edmf_a1(k)*edmf_w1(k)/bv), zwk)

              !IF (zwk > zi .AND. elf > 400.) THEN
              !   ! COMPUTE BouLac mixing length
              !   !CALL boulac_length0(k,kts,kte,zw,dz,qtke,thetaw,elBLmin0,elBLavg0)
              !   !elf = alp5*elBLavg0
              !   elf = MIN(MAX(50.*SQRT(qtke(k)), 400.), zwk)
              !ENDIF

           ELSE
              ! use version in development for RAP/HRRR 2016
              ! JAYMES-
              ! tau_cloud is an eddy turnover timescale;
              ! see Teixeira and Cheinet (2004), Eq. 1, and
              ! Cheinet and Teixeira (2003), Eq. 7.  The
              ! coefficient 0.5 is tuneable. Expression in
              ! denominator is identical to vsc (a convective
              ! velocity scale), except that elt is relpaced
              ! by zi, and zero is replaced by 1.0e-4 to
              ! prevent division by zero.
              !tau_cloud = MIN(MAX(0.5*zi/((gtr*zi*MAX(vflx,1.0e-4))**onethird),50.),150.)
              wstar = 1.25*(gtr*zi*MAX(vflx,1.0e-4))**onethird
              tau_cloud = MIN(MAX(ctau * wstar/grav, 50.), 200.)
              !minimize influence of surface heat flux on tau far away from the PBLH.
              wt=.5*TANH((zwk - (zi2+h1))/h2) + .5
              !tau_cloud = tau_cloud*(1.-wt) + 50.*wt
              tau_cloud = tau_cloud*(1.-wt) + MAX(100.,dzk*0.25)*wt

              elb = MIN(tau_cloud*SQRT(MIN(qtke(k),40.)), zwk)
              !elf = elb
              elf = elb !/(1. + (elb/800.))  !bound free-atmos mixing length to < 800 m.
              elb_mf = elb
         END IF
         elf    = elf/(1. + (elf/800.))  !bound free-atmos mixing length to < 800 m.
         elb_mf = MAX(elb_mf, 0.01) !to avoid divide-by-zero below

         !   **  Length scale in the surface layer  **
         IF ( rmo .GT. 0.0 ) THEN
            els  = karman*zwk/(1.0+cns*MIN( zwk*rmo, zmax ))
         ELSE
            els  =  karman*zwk*( 1.0 - alp4* zwk*rmo )**0.2
         END IF

         !   ** NOW BLEND THE MIXING LENGTH SCALES:
         wt=.5*TANH((zwk - (zi2+h1))/h2) + .5

         !try squared-blending
         el(k) = SQRT( els**2/(1. + (els**2/elt**2) +(els**2/elb_mf**2)))
         el(k) = el(k)*(1.-wt) + elf*wt

         ! include scale-awareness. For now, use simple asymptotic kz -> 12 m (should be ~dz).
         el_les= MIN(els/(1. + (els/12.)), elb_mf)
         el(k) = el(k)*Psig_bl + (1.-Psig_bl)*el_les

       END DO

    END SELECT


#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mym_length

END MODULE
 
