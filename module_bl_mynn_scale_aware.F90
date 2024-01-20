MODULE myn_scale_aware
IMPLICIT NONE
CONTAINS 
!===============================================================


SUBROUTINE SCALE_AWARE(dx,PBL1,Psig_bl,Psig_shcu)

    !---------------------------------------------------------------
    !             NOTES ON SCALE-AWARE FORMULATION
    !
    !JOE: add scale-aware factor (Psig) here, taken from Honnert et al. (2011,
    !     JAS) and/or from Hyeyum Hailey Shin and Song-You Hong (2013, JAS)
    !
    ! Psig_bl tapers local mixing
    ! Psig_shcu tapers nonlocal mixing

    REAL,INTENT(IN) :: dx,PBL1
    REAL, INTENT(OUT) :: Psig_bl,Psig_shcu
    REAL :: dxdh

    Psig_bl=1.0
    Psig_shcu=1.0
    dxdh=MAX(2.5*dx,10.)/MIN(PBL1,3000.)
    ! Honnert et al. 2011, TKE in PBL  *** original form used until 201605
    !Psig_bl= ((dxdh**2) + 0.07*(dxdh**0.667))/((dxdh**2) + &
    !         (3./21.)*(dxdh**0.67) + (3./42.))
    ! Honnert et al. 2011, TKE in entrainment layer
    !Psig_bl= ((dxdh**2) + (4./21.)*(dxdh**0.667))/((dxdh**2) + &
     !        (3./20.)*(dxdh**0.67) + (7./21.))
    ! New form to preseve parameterized mixing - only down 5% at dx = 750 m
     Psig_bl= ((dxdh**2) + 0.106*(dxdh**0.667))/((dxdh**2) +0.066*(dxdh**0.667) + 0.071)

    !assume a 500 m cloud depth for shallow-cu clods
    dxdh=MAX(2.5*dx,10.)/MIN(PBL1+500.,3500.)
    ! Honnert et al. 2011, TKE in entrainment layer *** original form used until 201605
    !Psig_shcu= ((dxdh**2) + (4./21.)*(dxdh**0.667))/((dxdh**2) + &
    !         (3./20.)*(dxdh**0.67) + (7./21.))

    ! Honnert et al. 2011, TKE in cumulus
    !Psig(i)= ((dxdh**2) + 1.67*(dxdh**1.4))/((dxdh**2) +1.66*(dxdh**1.4) +
    !0.2)

    ! Honnert et al. 2011, w'q' in PBL
    !Psig(i)= 0.5 + 0.5*((dxdh**2) + 0.03*(dxdh**1.4) -
    !(4./13.))/((dxdh**2) + 0.03*(dxdh**1.4) + (4./13.))
    ! Honnert et al. 2011, w'q' in cumulus
    !Psig(i)= ((dxdh**2) - 0.07*(dxdh**1.4))/((dxdh**2) -0.07*(dxdh**1.4) +
    !0.02)

    ! Honnert et al. 2011, q'q' in PBL
    !Psig(i)= 0.5 + 0.5*((dxdh**2) + 0.25*(dxdh**0.667) -0.73)/((dxdh**2)
    !-0.03*(dxdh**0.667) + 0.73)
    ! Honnert et al. 2011, q'q' in cumulus
    !Psig(i)= ((dxdh**2) - 0.34*(dxdh**1.4))/((dxdh**2) - 0.35*(dxdh**1.4)
    !+ 0.37)

    ! Hyeyum Hailey Shin and Song-You Hong 2013, TKE in PBL (same as Honnert's above)
    !Psig_shcu= ((dxdh**2) + 0.070*(dxdh**0.667))/((dxdh**2)
    !+0.142*(dxdh**0.667) + 0.071)
    ! Hyeyum Hailey Shin and Song-You Hong 2013, TKE in entrainment zone  *** switch to this form 201605
    Psig_shcu= ((dxdh**2) + 0.145*(dxdh**0.667))/((dxdh**2) +0.172*(dxdh**0.667) + 0.170)

    ! Hyeyum Hailey Shin and Song-You Hong 2013, w'theta' in PBL
    !Psig(i)= 0.5 + 0.5*((dxdh**2) -0.098)/((dxdh**2) + 0.106) 
    ! Hyeyum Hailey Shin and Song-You Hong 2013, w'theta' in entrainment zone
    !Psig(i)= 0.5 + 0.5*((dxdh**2) - 0.112*(dxdh**0.25) -0.071)/((dxdh**2)
    !+ 0.054*(dxdh**0.25) + 0.10)

    !print*,"in scale_aware; dx, dxdh, Psig(i)=",dx,dxdh,Psig(i)
    !If(Psig_bl(i) < 0.0 .OR. Psig(i) > 1.)print*,"dx, dxdh, Psig(i)=",dx,dxdh,Psig_bl(i) 
    If(Psig_bl > 1.0) Psig_bl=1.0
    If(Psig_bl < 0.0) Psig_bl=0.0

    If(Psig_shcu > 1.0) Psig_shcu=1.0
    If(Psig_shcu < 0.0) Psig_shcu=0.0

  END SUBROUTINE SCALE_AWARE

! =====================================================================

END MODULE
 
