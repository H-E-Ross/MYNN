MODULE myn_topdown_cloudrad

use module_bl_mynn_common,only:&
        cp, ep_2, grav, karman, p608, r_d, xlv, xlvcp
IMPLICIT NONE
CONTAINS
! ==================================================================
 SUBROUTINE topdown_cloudrad(kts,kte,dz1,zw,xland,kpbl,PBLH,  &
               &sqc,sqi,sqw,thl,th1,ex1,p1,rho1,thetav,       &
               &cldfra_bl1D,rthraten,                         &
               &maxKHtopdown,KHtopdown,TKEprodTD              )

    !input
    integer, intent(in) :: kte,kts
    real, dimension(kts:kte), intent(in) :: dz1,sqc,sqi,sqw,&
          thl,th1,ex1,p1,rho1,thetav,cldfra_bl1D,rthraten
    real, dimension(kts:kte+1), intent(in) :: zw
    real, intent(in) :: pblh,xland
    integer,intent(in) :: kpbl
    !output
    real, intent(out) :: maxKHtopdown
    real, dimension(kts:kte), intent(out) :: KHtopdown,TKEprodTD
    !local
    real, dimension(kts:kte) :: zfac,wscalek2,zfacent
    real :: bfx0,sflux,wm2,wm3,h1,h2,bfxpbl,dthvx,tmp1
    real :: temps,templ,zl1,wstar3_2
    real :: ent_eff,radsum,radflux,we,rcldb,rvls,minrad,zminrad
    real, parameter :: pfac =2.0, zfmin = 0.01, phifac=8.0
    integer :: k,kk,kminrad
    logical :: cloudflg

    cloudflg=.false.
    minrad=100.
    kminrad=kpbl
    zminrad=PBLH
    KHtopdown(kts:kte)=0.0
    TKEprodTD(kts:kte)=0.0
    maxKHtopdown=0.0

    !CHECK FOR STRATOCUMULUS-TOPPED BOUNDARY LAYERS
    DO kk = MAX(1,kpbl-2),kpbl+3
       if (sqc(kk).gt. 1.e-6 .OR. sqi(kk).gt. 1.e-6 .OR. &
           cldfra_bl1D(kk).gt.0.5) then
          cloudflg=.true.
       endif
       if (rthraten(kk) < minrad)then
          minrad=rthraten(kk)
          kminrad=kk
          zminrad=zw(kk) + 0.5*dz1(kk)
       endif
    ENDDO

    IF (MAX(kminrad,kpbl) < 2)cloudflg = .false.
    IF (cloudflg) THEN
       zl1 = dz1(kts)
       k = MAX(kpbl-1, kminrad-1)
       !Best estimate of height of TKE source (top of downdrafts):
       !zminrad = 0.5*pblh(i) + 0.5*zminrad

       templ=thl(k)*ex1(k)
       !rvls is ws at full level
       rvls=100.*6.112*EXP(17.67*(templ-273.16)/(templ-29.65))*(ep_2/p1(k+1))
       temps=templ + (sqw(k)-rvls)/(cp/xlv  +  ep_2*xlv*rvls/(r_d*templ**2))
       rvls=100.*6.112*EXP(17.67*(temps-273.15)/(temps-29.65))*(ep_2/p1(k+1))
       rcldb=max(sqw(k)-rvls,0.)

       !entrainment efficiency
       dthvx     = (thl(k+2) + th1(k+2)*p608*sqw(k+2)) &
                 - (thl(k)   + th1(k)  *p608*sqw(k))
       dthvx     = max(dthvx,0.1)
       tmp1      = xlvcp * rcldb/(ex1(k)*dthvx)
       !Originally from Nichols and Turton (1986), where a2 = 60, but lowered
       !here to 8, as in Grenier and Bretherton (2001).
       ent_eff   = 0.2 + 0.2*8.*tmp1

       radsum=0.
       DO kk = MAX(1,kpbl-3),kpbl+3
          radflux=rthraten(kk)*ex1(kk)         !converts theta/s to temp/s
          radflux=radflux*cp/grav*(p1(kk)-p1(kk+1)) ! converts temp/s to W/m^2
          if (radflux < 0.0 ) radsum=abs(radflux)+radsum
       ENDDO

       !More strict limits over land to reduce stable-layer mixouts
       if ((xland-1.5).GE.0)THEN      ! WATER
          radsum=MIN(radsum,90.0)
          bfx0 = max(radsum/rho1(k)/cp,0.)
       else                           ! LAND
          radsum=MIN(0.25*radsum,30.0)!practically turn off over land
          bfx0 = max(radsum/rho1(k)/cp - max(sflux,0.0),0.)
       endif

       !entrainment from PBL top thermals
       wm3    = grav/thetav(k)*bfx0*MIN(pblh,1500.) ! this is wstar3(i)
       wm2    = wm2 + wm3**h2
       bfxpbl = - ent_eff * bfx0
       dthvx  = max(thetav(k+1)-thetav(k),0.1)
       we     = max(bfxpbl/dthvx,-sqrt(wm3**h2))

       DO kk = kts,kpbl+3
          !Analytic vertical profile
          zfac(kk) = min(max((1.-(zw(kk+1)-zl1)/(zminrad-zl1)),zfmin),1.)
          zfacent(kk) = 10.*MAX((zminrad-zw(kk+1))/zminrad,0.0)*(1.-zfac(kk))**3

          !Calculate an eddy diffusivity profile (not used at the moment)
          wscalek2(kk) = (phifac*karman*wm3*(zfac(kk)))**h1
          !Modify shape of Kh to be similar to Lock et al (2000): use pfac = 3.0
          KHtopdown(kk) = wscalek2(kk)*karman*(zminrad-zw(kk+1))*(1.-zfac(kk))**3 !pfac
          KHtopdown(kk) = MAX(KHtopdown(kk),0.0)

          !Calculate TKE production = 2(g/TH)(w'TH'), where w'TH' = A(TH/g)wstar^3/PBLH,
          !A = ent_eff, and wstar is associated with the radiative cooling at top of PBL.
          !An analytic profile controls the magnitude of this TKE prod in the vertical.
          TKEprodTD(kk)=2.*ent_eff*wm3/MAX(pblh,100.)*zfacent(kk)
          TKEprodTD(kk)= MAX(TKEprodTD(kk),0.0)
       ENDDO
    ENDIF !end cloud check
    maxKHtopdown=MAXVAL(KHtopdown(:))

 END SUBROUTINE topdown_cloudrad
! ==================================================================
! ===================================================================
! ===================================================================

END MODULE
