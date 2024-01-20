MODULE mynn_predict
use module_bl_mynn_common,only:karman
IMPLICIT NONE
CONTAINS
! ==================================================================
!     SUBROUTINE  mym_predict:
!
!     Input variables:    see subroutine mym_initialize and turbulence
!       qke(nx,nz,ny) : qke at (n)th time level
!       tsq, ...cov     : ditto
!
!     Output variables:
!       qke(nx,nz,ny) : qke at (n+1)th time level
!       tsq, ...cov     : ditto
!
!     Work arrays:
!       qkw(nx,nz,ny)   : q at the center of the grid boxes        (m/s)
!       bp (nx,nz,ny)   : = 1/2*F,     see below
!       rp (nx,nz,ny)   : = P-1/2*F*Q, see below
!
!     # The equation for a turbulent quantity Q can be expressed as
!          dQ/dt + Ah + Av = Dh + Dv + P - F*Q,                      (1)
!       where A is the advection, D the diffusion, P the production,
!       F*Q the dissipation and h and v denote horizontal and vertical,
!       respectively. If Q is q^2, F is 2q/B_1L.
!       Using the Crank-Nicholson scheme for Av, Dv and F*Q, a finite
!       difference equation is written as
!          Q{n+1} - Q{n} = dt  *( Dh{n}   - Ah{n}   + P{n} )
!                        + dt/2*( Dv{n}   - Av{n}   - F*Q{n}   )
!                        + dt/2*( Dv{n+1} - Av{n+1} - F*Q{n+1} ),    (2)
!       where n denotes the time level.
!       When the advection and diffusion terms are discretized as
!          dt/2*( Dv - Av ) = a(k)Q(k+1) - b(k)Q(k) + c(k)Q(k-1),    (3)
!       Eq.(2) can be rewritten as
!          - a(k)Q(k+1) + [ 1 + b(k) + dt/2*F ]Q(k) - c(k)Q(k-1)
!                 = Q{n} + dt  *( Dh{n}   - Ah{n}   + P{n} )
!                        + dt/2*( Dv{n}   - Av{n}   - F*Q{n}   ),    (4)
!       where Q on the left-hand side is at (n+1)th time level.
!
!       In this subroutine, a(k), b(k) and c(k) are obtained from
!       subprogram coefvu and are passed to subprogram tinteg via
!       common. 1/2*F and P-1/2*F*Q are stored in bp and rp,
!       respectively. Subprogram tinteg solves Eq.(4).
!
!       Modify this subroutine according to your numerical integration
!       scheme (program).
!
!-------------------------------------------------------------------
!>\ingroup gsd_mynn_edmf
!! This subroutine predicts the turbulent quantities at the next step.
  SUBROUTINE  mym_predict (kts,kte,                                     &
       &            closure,                                            &
       &            delt,                                               &
       &            dz,                                                 &
       &            ust, flt, flq, pmz, phh,                            &
       &            el, dfq, rho,                                       &
       &            pdk, pdt, pdq, pdc,                                 &
       &            qke, tsq, qsq, cov,                                 &
       &            s_aw,s_awqke,bl_mynn_edmf_tke,                      &
       &            qWT1D, qDISS1D,tke_budget,                           &
       &            b1,b2,sqfac)  !! TKE budget  (Puhales, 2020)
       
!-------------------------------------------------------------------
    INTEGER, INTENT(IN) :: kts,kte    

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    REAL, INTENT(IN)    :: closure,b1,b2,sqfac
    INTEGER, INTENT(IN) :: bl_mynn_edmf_tke, tke_budget
    REAL, INTENT(IN)    :: delt
    REAL, DIMENSION(kts:kte), INTENT(IN) :: dz, dfq, el, rho
    REAL, DIMENSION(kts:kte), INTENT(INOUT) :: pdk, pdt, pdq, pdc
    REAL, INTENT(IN)    ::  flt, flq, ust, pmz, phh
    REAL, DIMENSION(kts:kte), INTENT(INOUT) :: qke,tsq, qsq, cov
! WA 8/3/15
    REAL, DIMENSION(kts:kte+1), INTENT(INOUT) :: s_awqke,s_aw
    
    !!  TKE budget  (Puhales, 2020, WRF 4.2.1)  << EOB 
    REAL, DIMENSION(kts:kte), INTENT(OUT) :: qWT1D, qDISS1D  
    REAL, DIMENSION(kts:kte) :: tke_up,dzinv  
    !! >> EOB
    
    INTEGER :: k
    REAL, DIMENSION(kts:kte) :: qkw, bp, rp, df3q
    REAL :: vkz,pdk1,phm,pdt1,pdq1,pdc1,b1l,b2l,onoff
    REAL, DIMENSION(kts:kte) :: dtz
    REAL, DIMENSION(kts:kte) :: a,b,c,d,x

    REAL, DIMENSION(kts:kte) :: rhoinv
    REAL, DIMENSION(kts:kte+1) :: rhoz,kqdz,kmdz

    ! REGULATE THE MOMENTUM MIXING FROM THE MASS-FLUX SCHEME (on or off)
    IF (bl_mynn_edmf_tke == 0) THEN
       onoff=0.0
    ELSE
       onoff=1.0
    ENDIF

!   **  Strictly, vkz*h(i,j) -> karman*( 0.5*dz(1)*h(i,j)+z0 )  **
    vkz = karman*0.5*dz(kts)
!
!   **  dfq for the TKE is 3.0*dfm.  **
!
    DO k = kts,kte
!!       qke(k) = MAX(qke(k), 0.0)
       qkw(k) = SQRT( MAX( qke(k), 0.0 ) )
       df3q(k)=Sqfac*dfq(k)
       dtz(k)=delt/dz(k)
    END DO
!
!JOE-add conservation + stability criteria
    !Prepare "constants" for diffusion equation.
    !khdz = rho*Kh/dz = rho*dfh
    rhoz(kts)  =rho(kts)
    rhoinv(kts)=1./rho(kts)
    kqdz(kts)  =rhoz(kts)*df3q(kts)
    kmdz(kts)  =rhoz(kts)*dfq(kts)
    DO k=kts+1,kte
       rhoz(k)  =(rho(k)*dz(k-1) + rho(k-1)*dz(k))/(dz(k-1)+dz(k))
       rhoz(k)  =  MAX(rhoz(k),1E-4)
       rhoinv(k)=1./MAX(rho(k),1E-4)
       kqdz(k)  = rhoz(k)*df3q(k) ! for TKE
       kmdz(k)  = rhoz(k)*dfq(k)  ! for T'2, q'2, and T'q'
    ENDDO
    rhoz(kte+1)=rhoz(kte)
    kqdz(kte+1)=rhoz(kte+1)*df3q(kte)
    kmdz(kte+1)=rhoz(kte+1)*dfq(kte)

    !stability criteria for mf
    DO k=kts+1,kte-1
       kqdz(k) = MAX(kqdz(k),  0.5* s_aw(k))
       kqdz(k) = MAX(kqdz(k), -0.5*(s_aw(k)-s_aw(k+1)))
       kmdz(k) = MAX(kmdz(k),  0.5* s_aw(k))
       kmdz(k) = MAX(kmdz(k), -0.5*(s_aw(k)-s_aw(k+1)))
    ENDDO
!JOE-end conservation mods

    pdk1 = 2.0*ust**3*pmz/( vkz )
    phm  = 2.0/ust   *phh/( vkz )
    pdt1 = phm*flt**2
    pdq1 = phm*flq**2
    pdc1 = phm*flt*flq
!
!   **  pdk(i,j,1)+pdk(i,j,2) corresponds to pdk1.  **
    pdk(kts) = pdk1 -pdk(kts+1)

!!    pdt(kts) = pdt1 -pdt(kts+1)
!!    pdq(kts) = pdq1 -pdq(kts+1)
!!    pdc(kts) = pdc1 -pdc(kts+1)
    pdt(kts) = pdt(kts+1)
    pdq(kts) = pdq(kts+1)
    pdc(kts) = pdc(kts+1)
!
!   **  Prediction of twice the turbulent kinetic energy  **
!!    DO k = kts+1,kte-1
    DO k = kts,kte-1
       b1l = b1*0.5*( el(k+1)+el(k) )
       bp(k) = 2.*qkw(k) / b1l
       rp(k) = pdk(k+1) + pdk(k)
    END DO

!!    a(1)=0.
!!    b(1)=1.
!!    c(1)=-1.
!!    d(1)=0.

! Since df3q(kts)=0.0, a(1)=0.0 and b(1)=1.+dtz(k)*df3q(k+1)+bp(k)*delt.
    DO k=kts,kte-1
!       a(k-kts+1)=-dtz(k)*df3q(k)
!       b(k-kts+1)=1.+dtz(k)*(df3q(k)+df3q(k+1))+bp(k)*delt
!       c(k-kts+1)=-dtz(k)*df3q(k+1)
!       d(k-kts+1)=rp(k)*delt + qke(k)
! WA 8/3/15 add EDMF contribution
!       a(k)=   - dtz(k)*df3q(k) + 0.5*dtz(k)*s_aw(k)*onoff
!       b(k)=1. + dtz(k)*(df3q(k)+df3q(k+1)) &
!               + 0.5*dtz(k)*(s_aw(k)-s_aw(k+1))*onoff + bp(k)*delt
!       c(k)=   - dtz(k)*df3q(k+1) - 0.5*dtz(k)*s_aw(k+1)*onoff
!       d(k)=rp(k)*delt + qke(k) + dtz(k)*(s_awqke(k)-s_awqke(k+1))*onoff
!JOE 8/22/20 improve conservation
       a(k)=   - dtz(k)*kqdz(k)*rhoinv(k)                       &
           &   + 0.5*dtz(k)*rhoinv(k)*s_aw(k)*onoff
       b(k)=1. + dtz(k)*(kqdz(k)+kqdz(k+1))*rhoinv(k)           &
           &   + 0.5*dtz(k)*rhoinv(k)*(s_aw(k)-s_aw(k+1))*onoff &
           &   + bp(k)*delt
       c(k)=   - dtz(k)*kqdz(k+1)*rhoinv(k)                     &
           &   - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)*onoff
       d(k)=rp(k)*delt + qke(k)                                 &
           &   + dtz(k)*rhoinv(k)*(s_awqke(k)-s_awqke(k+1))*onoff
    ENDDO

!!    DO k=kts+1,kte-1
!!       a(k-kts+1)=-dtz(k)*df3q(k)
!!       b(k-kts+1)=1.+dtz(k)*(df3q(k)+df3q(k+1))
!!       c(k-kts+1)=-dtz(k)*df3q(k+1)
!!       d(k-kts+1)=rp(k)*delt + qke(k) - qke(k)*bp(k)*delt
!!    ENDDO

!! "no flux at top"
!    a(kte)=-1. !0.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=0.
!! "prescribed value"
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=qke(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)

    DO k=kts,kte
!       qke(k)=max(d(k-kts+1), 1.e-4)
       qke(k)=max(x(k), 1.e-4)
       qke(k)=min(qke(k), 150.)
    ENDDO
      
   
!!  TKE budget  (Puhales, 2020, WRF 4.2.1)  << EOB 
    IF (tke_budget .eq. 1) THEN
       !! TKE Vertical transport << EOBvt
        tke_up=0.5*qke
        dzinv=1./dz
        k=kts
        qWT1D(k)=dzinv(k)*(                                    &
            &  (kqdz(k+1)*(tke_up(k+1)-tke_up(k))-kqdz(k)*tke_up(k)) &
            &  + 0.5*rhoinv(k)*(s_aw(k+1)*tke_up(k+1)          &
            &  +      (s_aw(k+1)-s_aw(k))*tke_up(k)            &
            &  +      (s_awqke(k)-s_awqke(k+1)))*onoff) !unstaggered
        DO k=kts+1,kte-1
            qWT1D(k)=dzinv(k)*(                                &
            & (kqdz(k+1)*(tke_up(k+1)-tke_up(k))-kqdz(k)*(tke_up(k)-tke_up(k-1))) &
            &  + 0.5*rhoinv(k)*(s_aw(k+1)*tke_up(k+1)          &
            &  +      (s_aw(k+1)-s_aw(k))*tke_up(k)            &
            &  -                  s_aw(k)*tke_up(k-1)          &
            &  +      (s_awqke(k)-s_awqke(k+1)))*onoff) !unstaggered
        ENDDO
        k=kte
        qWT1D(k)=dzinv(k)*(-kqdz(k)*(tke_up(k)-tke_up(k-1)) &
            &  + 0.5*rhoinv(k)*(-s_aw(k)*tke_up(k)-s_aw(k)*tke_up(k-1)+s_awqke(k))*onoff) !unstaggared
        !!  >> EOBvt
        qDISS1D=bp*tke_up !! TKE dissipation rate !unstaggered
    END IF
!! >> EOB 
   
    IF ( closure > 2.5 ) THEN

       !   **  Prediction of the moisture variance  **
       DO k = kts,kte-1
          b2l   = b2*0.5*( el(k+1)+el(k) )
          bp(k) = 2.*qkw(k) / b2l
          rp(k) = pdq(k+1) + pdq(k)
       END DO

       !zero gradient for qsq at bottom and top
       !a(1)=0.
       !b(1)=1.
       !c(1)=-1.
       !d(1)=0.

       ! Since dfq(kts)=0.0, a(1)=0.0 and b(1)=1.+dtz(k)*dfq(k+1)+bp(k)*delt.
       DO k=kts,kte-1
          a(k)=   - dtz(k)*kmdz(k)*rhoinv(k)
          b(k)=1. + dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k) + bp(k)*delt
          c(k)=   - dtz(k)*kmdz(k+1)*rhoinv(k)
          d(k)=rp(k)*delt + qsq(k)
       ENDDO

       a(kte)=-1. !0.
       b(kte)=1.
       c(kte)=0.
       d(kte)=0.

!       CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
       
       DO k=kts,kte
          !qsq(k)=d(k-kts+1)
          qsq(k)=MAX(x(k),1e-17)
       ENDDO
    ELSE
       !level 2.5 - use level 2 diagnostic
       DO k = kts,kte-1
          IF ( qkw(k) .LE. 0.0 ) THEN
             b2l = 0.0
          ELSE
             b2l = b2*0.25*( el(k+1)+el(k) )/qkw(k)
          END IF
          qsq(k) = b2l*( pdq(k+1)+pdq(k) )
       END DO
       qsq(kte)=qsq(kte-1)
    END IF
!!!!!!!!!!!!!!!!!!!!!!end level 2.6   

    IF ( closure .GE. 3.0 ) THEN
!
!   **  dfq for the scalar variance is 1.0*dfm.  **
!
!   **  Prediction of the temperature variance  **
!!       DO k = kts+1,kte-1
       DO k = kts,kte-1
          b2l = b2*0.5*( el(k+1)+el(k) )
          bp(k) = 2.*qkw(k) / b2l
          rp(k) = pdt(k+1) + pdt(k) 
       END DO
       
!zero gradient for tsq at bottom and top
       
!!       a(1)=0.
!!       b(1)=1.
!!       c(1)=-1.
!!       d(1)=0.

! Since dfq(kts)=0.0, a(1)=0.0 and b(1)=1.+dtz(k)*dfq(k+1)+bp(k)*delt.
       DO k=kts,kte-1
          !a(k-kts+1)=-dtz(k)*dfq(k)
          !b(k-kts+1)=1.+dtz(k)*(dfq(k)+dfq(k+1))+bp(k)*delt
          !c(k-kts+1)=-dtz(k)*dfq(k+1)
          !d(k-kts+1)=rp(k)*delt + tsq(k)
!JOE 8/22/20 improve conservation
          a(k)=   - dtz(k)*kmdz(k)*rhoinv(k)
          b(k)=1. + dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k) + bp(k)*delt
          c(k)=   - dtz(k)*kmdz(k+1)*rhoinv(k)
          d(k)=rp(k)*delt + tsq(k)
       ENDDO

!!       DO k=kts+1,kte-1
!!          a(k-kts+1)=-dtz(k)*dfq(k)
!!          b(k-kts+1)=1.+dtz(k)*(dfq(k)+dfq(k+1))
!!          c(k-kts+1)=-dtz(k)*dfq(k+1)
!!          d(k-kts+1)=rp(k)*delt + tsq(k) - tsq(k)*bp(k)*delt
!!       ENDDO

       a(kte)=-1. !0.
       b(kte)=1.
       c(kte)=0.
       d(kte)=0.
       
!       CALL tridiag(kte,a,b,c,d)
       CALL tridiag2(kte,a,b,c,d,x)

       DO k=kts,kte
!          tsq(k)=d(k-kts+1)
           tsq(k)=x(k)
       ENDDO

!   **  Prediction of the temperature-moisture covariance  **
!!       DO k = kts+1,kte-1
       DO k = kts,kte-1
          b2l = b2*0.5*( el(k+1)+el(k) )
          bp(k) = 2.*qkw(k) / b2l
          rp(k) = pdc(k+1) + pdc(k) 
       END DO
       
!zero gradient for tqcov at bottom and top
       
!!       a(1)=0.
!!       b(1)=1.
!!       c(1)=-1.
!!       d(1)=0.

! Since dfq(kts)=0.0, a(1)=0.0 and b(1)=1.+dtz(k)*dfq(k+1)+bp(k)*delt.
       DO k=kts,kte-1
          !a(k-kts+1)=-dtz(k)*dfq(k)
          !b(k-kts+1)=1.+dtz(k)*(dfq(k)+dfq(k+1))+bp(k)*delt
          !c(k-kts+1)=-dtz(k)*dfq(k+1)
          !d(k-kts+1)=rp(k)*delt + cov(k)
!JOE 8/22/20 improve conservation
          a(k)=   - dtz(k)*kmdz(k)*rhoinv(k)
          b(k)=1. + dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k) + bp(k)*delt
          c(k)=   - dtz(k)*kmdz(k+1)*rhoinv(k)
          d(k)=rp(k)*delt + cov(k)
       ENDDO

!!       DO k=kts+1,kte-1
!!          a(k-kts+1)=-dtz(k)*dfq(k)
!!          b(k-kts+1)=1.+dtz(k)*(dfq(k)+dfq(k+1))
!!          c(k-kts+1)=-dtz(k)*dfq(k+1)
!!          d(k-kts+1)=rp(k)*delt + cov(k) - cov(k)*bp(k)*delt
!!       ENDDO

       a(kte)=-1. !0.
       b(kte)=1.
       c(kte)=0.
       d(kte)=0.

!       CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
       
       DO k=kts,kte
!          cov(k)=d(k-kts+1)
          cov(k)=x(k)
       ENDDO
       
    ELSE

       !Not level 3 - default to level 2 diagnostic
       DO k = kts,kte-1
          IF ( qkw(k) .LE. 0.0 ) THEN
             b2l = 0.0
          ELSE
             b2l = b2*0.25*( el(k+1)+el(k) )/qkw(k)
          END IF
!
          tsq(k) = b2l*( pdt(k+1)+pdt(k) )
          cov(k) = b2l*( pdc(k+1)+pdc(k) )
       END DO
       
       tsq(kte)=tsq(kte-1)
       cov(kte)=cov(kte-1)
      
    END IF

#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mym_predict

END MODULE
 
