MODULE myn_coeffs 

        IMPLICIT NONE
CONTAINS
! ==================================================================
!>\ingroup gsd_mynn_edmf
  SUBROUTINE retrieve_exchange_coeffs(kts,kte,&
       &dfm,dfh,dz,K_m,K_h)

!-------------------------------------------------------------------

    INTEGER , INTENT(in) :: kts,kte

    REAL, DIMENSION(KtS:KtE), INTENT(in) :: dz,dfm,dfh

    REAL, DIMENSION(KtS:KtE), INTENT(out) :: K_m, K_h


    INTEGER :: k
    REAL :: dzk

    K_m(kts)=0.
    K_h(kts)=0.

    DO k=kts+1,kte
       dzk = 0.5  *( dz(k)+dz(k-1) )
       K_m(k)=dfm(k)*dzk
       K_h(k)=dfh(k)*dzk
    ENDDO

  END SUBROUTINE retrieve_exchange_coeffs

! ==================================================================

END MODULE
