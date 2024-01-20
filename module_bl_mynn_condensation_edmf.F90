MODULE myn_condensation_edmf

use module_bl_mynn_common,only:&
        cp, p1000mb, rcp, rvovrd, xlv, p1000mb, rcp, rvovrd,xlvcp 
use mynn_functions
IMPLICIT NONE

CONTAINS
!=================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine 
subroutine condensation_edmf(QT,THL,P,zagl,THV,QC)
!
! zero or one condensation for edmf: calculates THV and QC
!
real,intent(in)   :: QT,THL,P,zagl
real,intent(out)  :: THV
real,intent(inout):: QC

integer :: niter,i
real :: diff,exn,t,th,qs,qcold

! constants used from module_model_constants.F
! p1000mb
! rcp ... Rd/cp
! xlv ... latent heat for water (2.5e6)
! cp
! rvord .. r_v/r_d  (1.6) 

! number of iterations
  niter=50
! minimum difference (usually converges in < 8 iterations with diff = 2e-5)
  diff=1.e-6

  EXN=(P/p1000mb)**rcp
  !QC=0.  !better first guess QC is incoming from lower level, do not set to zero
  do i=1,NITER
     T=EXN*THL + xlvcp*QC
     QS=qsat_blend(T,P)
     QCOLD=QC
     QC=0.5*QC + 0.5*MAX((QT-QS),0.)
     if (abs(QC-QCOLD)<Diff) exit
  enddo

  T=EXN*THL + xlvcp*QC
  QS=qsat_blend(T,P)
  QC=max(QT-QS,0.)

  !Do not allow saturation below 100 m
  if(zagl < 100.)QC=0.

  !THV=(THL+xlv/cp*QC).*(1+(1-rvovrd)*(QT-QC)-QC);
  THV=(THL+xlvcp*QC)*(1.+QT*(rvovrd-1.)-rvovrd*QC)

!  IF (QC > 0.0) THEN
!    PRINT*,"EDMF SAT, p:",p," iterations:",i
!    PRINT*," T=",T," THL=",THL," THV=",THV
!    PRINT*," QS=",QS," QT=",QT," QC=",QC,"ratio=",qc/qs
!  ENDIF

  !THIS BASICALLY GIVE THE SAME RESULT AS THE PREVIOUS LINE
  !TH = THL + xlv/cp/EXN*QC
  !THV= TH*(1. + 0.608*QT)

  !print *,'t,p,qt,qs,qc'
  !print *,t,p,qt,qs,qc 


end subroutine condensation_edmf

!===============================================================

subroutine condensation_edmf_r(QT,THL,P,zagl,THV,QC)
!                                                                                                
! zero or one condensation for edmf: calculates THL and QC                                       
! similar to condensation_edmf but with different inputs                                         
!                                                                                                
real,intent(in)   :: QT,THV,P,zagl
real,intent(out)  :: THL, QC

integer :: niter,i
real :: diff,exn,t,th,qs,qcold

! number of iterations                                                                           
  niter=50
! minimum difference                                                                             
  diff=2.e-5

  EXN=(P/p1000mb)**rcp
  ! assume first that th = thv                                                                   
  T = THV*EXN
  !QS = qsat_blend(T,P)                                                                          
  !QC = QS - QT                                                                                  

  QC=0.

  do i=1,NITER
     QCOLD = QC
     T = EXN*THV/(1.+QT*(rvovrd-1.)-rvovrd*QC)
     QS=qsat_blend(T,P)
     QC= MAX((QT-QS),0.)
     if (abs(QC-QCOLD)<Diff) exit
  enddo
  THL = (T - xlv/cp*QC)/EXN

end subroutine condensation_edmf_r


END MODULE
