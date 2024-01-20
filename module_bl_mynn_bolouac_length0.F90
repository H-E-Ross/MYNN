MODULE mynn_bolouac_length0
use module_bl_mynn_common,only:gtr
IMPLICIT NONE
CONTAINS
! ==================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine was taken from the BouLac scheme in WRF-ARW and modified for
!! integration into the MYNN PBL scheme. WHILE loops were added to reduce the
!! computational expense. This subroutine computes the length scales up and down
!! and then computes the min, average of the up/down length scales, and also
!! considers the distance to the surface.
!\param dlu  the distance a parcel can be lifted upwards give a finite
!  amount of TKE.
!\param dld  the distance a parcel can be displaced downwards given a
!  finite amount of TKE.
!\param lb1  the minimum of the length up and length down
!\param lb2  the average of the length up and length down
  SUBROUTINE boulac_length0(k,kts,kte,zw,dz,qtke,theta,lb1,lb2)
!
!    NOTE: This subroutine was taken from the BouLac scheme in WRF-ARW
!          and modified for integration into the MYNN PBL scheme.
!          WHILE loops were added to reduce the computational expense.
!          This subroutine computes the length scales up and down
!          and then computes the min, average of the up/down
!          length scales, and also considers the distance to the
!          surface.
!
!      dlu = the distance a parcel can be lifted upwards give a finite
!            amount of TKE.
!      dld = the distance a parcel can be displaced downwards given a
!            finite amount of TKE.
!      lb1 = the minimum of the length up and length down
!      lb2 = the average of the length up and length down
!-------------------------------------------------------------------

     INTEGER, INTENT(IN) :: k,kts,kte
     REAL, DIMENSION(kts:kte), INTENT(IN) :: qtke,dz,theta
     REAL, INTENT(OUT) :: lb1,lb2
     REAL, DIMENSION(kts:kte+1), INTENT(IN) :: zw

     !LOCAL VARS
     INTEGER :: izz, found
     REAL :: dlu,dld
     REAL :: dzt, zup, beta, zup_inf, bbb, tl, zdo, zdo_sup, zzz


     !----------------------------------
     ! FIND DISTANCE UPWARD             
     !----------------------------------
     zup=0.
     dlu=zw(kte+1)-zw(k)-dz(k)*0.5
     zzz=0.
     zup_inf=0.
     beta=gtr           !Buoyancy coefficient (g/tref)

     !print*,"FINDING Dup, k=",k," zw=",zw(k)

     if (k .lt. kte) then      !cant integrate upwards from highest level
        found = 0
        izz=k
        DO WHILE (found .EQ. 0)

           if (izz .lt. kte) then
              dzt=dz(izz)                   ! layer depth above
              zup=zup-beta*theta(k)*dzt     ! initial PE the parcel has at k
              !print*,"  ",k,izz,theta(izz),dz(izz)
              zup=zup+beta*(theta(izz+1)+theta(izz))*dzt*0.5 ! PE gained by lifting a parcel to izz+1
              zzz=zzz+dzt                   ! depth of layer k to izz+1
              !print*,"  PE=",zup," TKE=",qtke(k)," z=",zw(izz)
              if (qtke(k).lt.zup .and. qtke(k).ge.zup_inf) then
                 bbb=(theta(izz+1)-theta(izz))/dzt
                 if (bbb .ne. 0.) then
                    !fractional distance up into the layer where TKE becomes < PE
                    tl=(-beta*(theta(izz)-theta(k)) + &
                      & sqrt( max(0.,(beta*(theta(izz)-theta(k)))**2 + &
                      &       2.*bbb*beta*(qtke(k)-zup_inf))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(k))then
                       tl=(qtke(k)-zup_inf)/(beta*(theta(izz)-theta(k)))
                    else
                       tl=0.
                    endif
                 endif
                 dlu=zzz-dzt+tl
                 !print*,"  FOUND Dup:",dlu," z=",zw(izz)," tl=",tl
                 found =1
              endif
              zup_inf=zup
              izz=izz+1
           ELSE
              found = 1
           ENDIF

        ENDDO

     endif

     !----------------------------------
     ! FIND DISTANCE DOWN               
     !----------------------------------
     zdo=0.
     zdo_sup=0.
     dld=zw(k)
     zzz=0.

     !print*,"FINDING Ddown, k=",k," zwk=",zw(k)
     if (k .gt. kts) then  !cant integrate downwards from lowest level

        found = 0
        izz=k
        DO WHILE (found .EQ. 0)

           if (izz .gt. kts) then
              dzt=dz(izz-1)
              zdo=zdo+beta*theta(k)*dzt
              !print*,"  ",k,izz,theta(izz),dz(izz-1)
              zdo=zdo-beta*(theta(izz-1)+theta(izz))*dzt*0.5
              zzz=zzz+dzt
              !print*,"  PE=",zdo," TKE=",qtke(k)," z=",zw(izz)
              if (qtke(k).lt.zdo .and. qtke(k).ge.zdo_sup) then
                 bbb=(theta(izz)-theta(izz-1))/dzt
                 if (bbb .ne. 0.) then
                    tl=(beta*(theta(izz)-theta(k))+ &
                      & sqrt( max(0.,(beta*(theta(izz)-theta(k)))**2 + &
                      &       2.*bbb*beta*(qtke(k)-zdo_sup))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(k)) then
                       tl=(qtke(k)-zdo_sup)/(beta*(theta(izz)-theta(k)))
                    else
                       tl=0.
                    endif
                 endif
                 dld=zzz-dzt+tl
                 !print*,"  FOUND Ddown:",dld," z=",zw(izz)," tl=",tl
                 found = 1
              endif
              zdo_sup=zdo
              izz=izz-1
           ELSE
              found = 1
           ENDIF
        ENDDO

     endif

     !----------------------------------
     ! GET MINIMUM (OR AVERAGE)         
     !----------------------------------
     !The surface layer length scale can exceed z for large z/L,
     !so keep maximum distance down > z.
     dld = min(dld,zw(k+1))!not used in PBL anyway, only free atmos
     lb1 = min(dlu,dld)     !minimum
     !JOE-fight floating point errors
     dlu=MAX(0.1,MIN(dlu,1000.))
     dld=MAX(0.1,MIN(dld,1000.))
     lb2 = sqrt(dlu*dld)    !average - biased towards smallest
     !lb2 = 0.5*(dlu+dld)   !average

     if (k .eq. kte) then
        lb1 = 0.
        lb2 = 0.
     endif
     !print*,"IN MYNN-BouLac",k,lb1
     !print*,"IN MYNN-BouLac",k,dld,dlu

  END SUBROUTINE boulac_length0

END MODULE
 
