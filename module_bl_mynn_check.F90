MODULE myn_check

use module_bl_mynn_common,only:xlscp,xlvcp 
IMPLICIT NONE

CONTAINS
! ==================================================================
  SUBROUTINE moisture_check(kte, delt, dp, exner, &
                            qv, qc, qi, th,       &
                            dqv, dqc, dqi, dth )

  ! This subroutine was adopted from the CAM-UW ShCu scheme and 
  ! adapted for use here.
  !
  ! If qc < qcmin, qi < qimin, or qv < qvmin happens in any layer,
  ! force them to be larger than minimum value by (1) condensating 
  ! water vapor into liquid or ice, and (2) by transporting water vapor 
  ! from the very lower layer.
  ! 
  ! We then update the final state variables and tendencies associated
  ! with this correction. If any condensation happens, update theta too.
  ! Note that (qv,qc,qi,th) are the final state variables after
  ! applying corresponding input tendencies and corrective tendencies.

    implicit none
    integer,  intent(in)     :: kte
    real, intent(in)         :: delt
    real, dimension(kte), intent(in)     :: dp, exner
    real, dimension(kte), intent(inout)  :: qv, qc, qi, th
    real, dimension(kte), intent(inout)  :: dqv, dqc, dqi, dth
    integer   k
    real ::  dqc2, dqi2, dqv2, sum, aa, dum
    real, parameter :: qvmin = 1e-20,   &
                       qcmin = 0.0,     &
                       qimin = 0.0

    do k = kte, 1, -1  ! From the top to the surface
       dqc2 = max(0.0, qcmin-qc(k)) !qc deficit (>=0)
       dqi2 = max(0.0, qimin-qi(k)) !qi deficit (>=0)

       !fix tendencies
       dqc(k) = dqc(k) +  dqc2/delt
       dqi(k) = dqi(k) +  dqi2/delt
       dqv(k) = dqv(k) - (dqc2+dqi2)/delt
       dth(k) = dth(k) + xlvcp/exner(k)*(dqc2/delt) + &
                         xlscp/exner(k)*(dqi2/delt)
       !update species
       qc(k)  = qc(k)  +  dqc2
       qi(k)  = qi(k)  +  dqi2
       qv(k)  = qv(k)  -  dqc2 - dqi2
       th(k)  = th(k)  +  xlvcp/exner(k)*dqc2 + &
                          xlscp/exner(k)*dqi2

       !then fix qv
       dqv2   = max(0.0, qvmin-qv(k)) !qv deficit (>=0)
       dqv(k) = dqv(k) + dqv2/delt
       qv(k)  = qv(k)  + dqv2
       if( k .ne. 1 ) then
           qv(k-1)   = qv(k-1)  - dqv2*dp(k)/dp(k-1)
           dqv(k-1)  = dqv(k-1) - dqv2*dp(k)/dp(k-1)/delt
       endif
       qv(k) = max(qv(k),qvmin)
       qc(k) = max(qc(k),qcmin)
       qi(k) = max(qi(k),qimin)
    end do
    ! Extra moisture used to satisfy 'qv(1)>=qvmin' is proportionally
    ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
    ! preserves column moisture.
    if( dqv2 .gt. 1.e-20 ) then
        sum = 0.0
        do k = 1, kte
           if( qv(k) .gt. 2.0*qvmin ) sum = sum + qv(k)*dp(k)
        enddo
        aa = dqv2*dp(1)/max(1.e-20,sum)
        if( aa .lt. 0.5 ) then
            do k = 1, kte
               if( qv(k) .gt. 2.0*qvmin ) then
                   dum    = aa*qv(k)
                   qv(k)  = qv(k) - dum
                   dqv(k) = dqv(k) - dum/delt
               endif
            enddo
        else
        ! For testing purposes only (not yet found in any output):
        !    write(*,*) 'Full moisture conservation is impossible'
        endif
    endif

    return

  END SUBROUTINE moisture_check

END MODULE
