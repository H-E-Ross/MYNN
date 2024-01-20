MODULE myn_mix_chem 
IMPLICIT NONE
CONTAINS
! ==================================================================

  SUBROUTINE mynn_mix_chem(kts,kte,i,     &
       delt,dz,pblh,                      &
       nchem, kdvel, ndvel,               &
       chem1, vd1,                        &
       rho,                               &
       flt, tcd, qcd,                     &
       dfh,                               &
       s_aw, s_awchem,                    &
       emis_ant_no, frp, rrfs_sd,         &
       enh_mix, smoke_dbg                 )

!-------------------------------------------------------------------
    INTEGER, INTENT(in) :: kts,kte,i
    REAL, DIMENSION(kts:kte), INTENT(IN)    :: dfh,dz,tcd,qcd
    REAL, DIMENSION(kts:kte), INTENT(INOUT) :: rho
    REAL, INTENT(IN)    :: delt,flt,pblh
    INTEGER, INTENT(IN) :: nchem, kdvel, ndvel
    REAL, DIMENSION( kts:kte+1), INTENT(IN) :: s_aw
    REAL, DIMENSION( kts:kte, nchem ), INTENT(INOUT) :: chem1
    REAL, DIMENSION( kts:kte+1,nchem), INTENT(IN) :: s_awchem
    REAL, DIMENSION( ndvel ), INTENT(IN) :: vd1
    REAL, INTENT(IN) :: emis_ant_no,frp
    LOGICAL, INTENT(IN) :: rrfs_sd,enh_mix,smoke_dbg
!local vars

    REAL, DIMENSION(kts:kte)     :: dtz
    REAL, DIMENSION(kts:kte) :: a,b,c,d,x
    REAL :: rhs,dztop
    REAL :: t,dzk
    REAL :: hght 
    REAL :: khdz_old, khdz_back
    INTEGER :: k,kk,kmaxfire                         ! JLS 12/21/21
    INTEGER :: ic  ! Chemical array loop index
    
    INTEGER, SAVE :: icall

    REAL, DIMENSION(kts:kte) :: rhoinv
    REAL, DIMENSION(kts:kte+1) :: rhoz,khdz
    REAL, PARAMETER :: NO_threshold    = 10.0     ! For anthropogenic sources
    REAL, PARAMETER :: frp_threshold   = 10.0     ! RAR 02/11/22: I increased the frp threshold to enhance mixing over big fires
    REAL, PARAMETER :: pblh_threshold  = 100.0

    dztop=.5*(dz(kte)+dz(kte-1))

    DO k=kts,kte
       dtz(k)=delt/dz(k)
    ENDDO

    !Prepare "constants" for diffusion equation.
    !khdz = rho*Kh/dz = rho*dfh
    rhoz(kts)  =rho(kts)
    rhoinv(kts)=1./rho(kts)
    khdz(kts)  =rhoz(kts)*dfh(kts)

    DO k=kts+1,kte
       rhoz(k)  =(rho(k)*dz(k-1) + rho(k-1)*dz(k))/(dz(k-1)+dz(k))
       rhoz(k)  =  MAX(rhoz(k),1E-4)
       rhoinv(k)=1./MAX(rho(k),1E-4)
       dzk      = 0.5  *( dz(k)+dz(k-1) )
       khdz(k)  = rhoz(k)*dfh(k)
    ENDDO
    rhoz(kte+1)=rhoz(kte)
    khdz(kte+1)=rhoz(kte+1)*dfh(kte)

    !stability criteria for mf
    DO k=kts+1,kte-1
       khdz(k) = MAX(khdz(k),  0.5*s_aw(k))
       khdz(k) = MAX(khdz(k), -0.5*(s_aw(k)-s_aw(k+1)))
    ENDDO

    !Enhanced mixing over fires
    IF ( rrfs_sd .and. enh_mix ) THEN
       DO k=kts+1,kte-1
          khdz_old  = khdz(k)
          khdz_back = pblh * 0.15 / dz(k)
          !Modify based on anthropogenic emissions of NO and FRP
          IF ( pblh < pblh_threshold ) THEN
             IF ( emis_ant_no > NO_threshold ) THEN
                khdz(k) = MAX(1.1*khdz(k),sqrt((emis_ant_no / NO_threshold)) / dz(k) * rhoz(k)) ! JLS 12/21/21
!                khdz(k) = MAX(khdz(k),khdz_back)
             ENDIF
             IF ( frp > frp_threshold ) THEN
                kmaxfire = ceiling(log(frp))
                khdz(k) = MAX(1.1*khdz(k), (1. - k/(kmaxfire*2.)) * ((log(frp))**2.- 2.*log(frp)) / dz(k)*rhoz(k)) ! JLS 12/21/21
!                khdz(k) = MAX(khdz(k),khdz_back)
             ENDIF
          ENDIF
       ENDDO
    ENDIF

  !============================================
  ! Patterned after mixing of water vapor in mynn_tendencies.
  !============================================

    DO ic = 1,nchem
       k=kts

       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
       b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)           - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)
       d(k)=chem1(k,ic) & !dtz(k)*flt  !neglecting surface sources 
            & - dtz(k)*vd1(ic)*chem1(k,ic) &
            & - dtz(k)*rhoinv(k)*s_awchem(k+1,ic)

       DO k=kts+1,kte-1
          a(k)=  -dtz(k)*khdz(k)*rhoinv(k)     + 0.5*dtz(k)*rhoinv(k)*s_aw(k)
          b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k) + &
             &    0.5*dtz(k)*rhoinv(k)*(s_aw(k)-s_aw(k+1))
          c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw(k+1)
          d(k)=chem1(k,ic) + dtz(k)*rhoinv(k)*(s_awchem(k,ic)-s_awchem(k+1,ic))
       ENDDO

      ! prescribed value at top
       a(kte)=0.
       b(kte)=1.
       c(kte)=0.
       d(kte)=chem1(kte,ic)

       CALL tridiag3(kte,a,b,c,d,x)

       DO k=kts,kte
          chem1(k,ic)=x(k)
       ENDDO
    ENDDO

  END SUBROUTINE mynn_mix_chem

END MODULE
