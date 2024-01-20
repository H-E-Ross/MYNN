MODULE myn_ddmf_jpl
use module_bl_mynn_common,only:cp,grav,onethird,svp1
IMPLICIT NONE
CONTAINS
! ===================================================================
! This is the downdraft mass flux scheme - analogus to edmf_JPL but  
! flipped updraft to downdraft. This scheme is currently only tested 
! for Stratocumulus cloud conditions. For a detailed desctiption of the
! model, see paper.

SUBROUTINE DDMF_JPL(kts,kte,dt,zw,dz,p,              &
              &u,v,th,thl,thv,tk,qt,qv,qc,           &
              &rho,exner,                            &
              &ust,wthl,wqt,pblh,kpbl,               &
              &edmf_a_dd,edmf_w_dd, edmf_qt_dd,      &
              &edmf_thl_dd,edmf_ent_dd,edmf_qc_dd,   &
              &sd_aw,sd_awthl,sd_awqt,               &
              &sd_awqv,sd_awqc,sd_awu,sd_awv,        &
              &sd_awqke,                             &
              &qc_bl1d,cldfra_bl1d,                  &
              &rthraten                              )

        INTEGER, INTENT(IN) :: KTS,KTE,KPBL
        REAL,DIMENSION(KTS:KTE), INTENT(IN) :: U,V,TH,THL,TK,QT,QV,QC,&
            THV,P,rho,exner,rthraten,dz
        ! zw .. heights of the downdraft levels (edges of boxes)
        REAL,DIMENSION(KTS:KTE+1), INTENT(IN) :: ZW
        REAL, INTENT(IN) :: DT,UST,WTHL,WQT,PBLH

  ! outputs - downdraft properties
        REAL,DIMENSION(KTS:KTE), INTENT(OUT) :: edmf_a_dd,edmf_w_dd,   &
                      & edmf_qt_dd,edmf_thl_dd, edmf_ent_dd,edmf_qc_dd

  ! outputs - variables needed for solver (sd_aw - sum ai*wi, sd_awphi - sum ai*wi*phii)
        REAL,DIMENSION(KTS:KTE+1) :: sd_aw, sd_awthl, sd_awqt, sd_awu, &
                            sd_awv, sd_awqc, sd_awqv, sd_awqke, sd_aw2

        REAL,DIMENSION(KTS:KTE), INTENT(IN) :: qc_bl1d, cldfra_bl1d

        INTEGER, PARAMETER :: NDOWN=5, debug_mf=0 !fixing number of plumes to 5
  ! draw downdraft starting height randomly between cloud base and cloud top
        INTEGER, DIMENSION(1:NDOWN) :: DD_initK
        REAL   , DIMENSION(1:NDOWN) :: randNum
  ! downdraft properties
        REAL,DIMENSION(KTS:KTE+1,1:NDOWN) :: DOWNW,DOWNTHL,DOWNQT,&
                    DOWNQC,DOWNA,DOWNU,DOWNV,DOWNTHV

  ! entrainment variables
        REAl,DIMENSION(KTS+1:KTE+1,1:NDOWN) :: ENT,ENTf
        INTEGER,DIMENSION(KTS+1:KTE+1,1:NDOWN) :: ENTi

  ! internal variables
        INTEGER :: K,I,ki, kminrad, qlTop, p700_ind, qlBase
        REAL :: wthv,wstar,qstar,thstar,sigmaW,sigmaQT,sigmaTH,z0, &
            pwmin,pwmax,wmin,wmax,wlv,wtv,went,mindownw
        REAL :: B,QTn,THLn,THVn,QCn,Un,Vn,QKEn,Wn2,Wn,THVk,Pk, &
                EntEXP,EntW, Beta_dm, EntExp_M, rho_int
        REAL :: jump_thetav, jump_qt, jump_thetal, &
                refTHL, refTHV, refQT
  ! DD specific internal variables
        REAL :: minrad,zminrad, radflux, F0, wst_rad, wst_dd
        logical :: cloudflg

        REAL :: sigq,xl,rsl,cpm,a,mf_cf,diffqt,&
               Fng,qww,alpha,beta,bb,f,pt,t,q2p,b9,satvp,rhgrid

  ! w parameters
        REAL,PARAMETER :: &
            &Wa=1., &
            &Wb=1.5,&
            &Z00=100.,&
            &BCOEFF=0.2
  ! entrainment parameters
        REAL,PARAMETER :: &
        & L0=80,&
        & ENT0=0.2

   pwmin=-3. ! drawing from the negative tail -3sigma to -1sigma
   pwmax=-1.

  ! initialize downdraft properties
   DOWNW=0.
   DOWNTHL=0.
   DOWNTHV=0.
   DOWNQT=0.
   DOWNA=0.
   DOWNU=0.
   DOWNV=0.
   DOWNQC=0.
   ENT=0.
   DD_initK=0

   edmf_a_dd  =0.
   edmf_w_dd  =0.
   edmf_qt_dd =0.
   edmf_thl_dd=0.
   edmf_ent_dd=0.
   edmf_qc_dd =0.

   sd_aw=0.
   sd_awthl=0.
   sd_awqt=0.
   sd_awqv=0.
   sd_awqc=0.
   sd_awu=0.
   sd_awv=0.
   sd_awqke=0.

  ! FIRST, CHECK FOR STRATOCUMULUS-TOPPED BOUNDARY LAYERS
   cloudflg=.false.
   minrad=100.
   kminrad=kpbl
   zminrad=PBLH
   qlTop = 1 !initialize at 0
   qlBase = 1
   wthv=wthl+svp1*wqt
   do k = MAX(3,kpbl-2),kpbl+3
      if (qc(k).gt. 1.e-6 .AND. cldfra_bl1D(k).gt.0.5) then
          cloudflg=.true. ! found Sc cloud
          qlTop = k       ! index for Sc cloud top
      endif
   enddo

   do k = qlTop, kts, -1
      if (qc(k) .gt. 1E-6) then
         qlBase = k ! index for Sc cloud base
      endif
   enddo
   qlBase = (qlTop+qlBase)/2 ! changed base to half way through the cloud

!   call init_random_seed_1()
!   call RANDOM_NUMBER(randNum)
   do i=1,NDOWN
      ! downdraft starts somewhere between cloud base to cloud top
      ! the probability is equally distributed
      DD_initK(i) = qlTop ! nint(randNum(i)*REAL(qlTop-qlBase)) + qlBase
   enddo

   ! LOOP RADFLUX
   F0 = 0.
   do k = 1, qlTop ! Snippet from YSU, YSU loops until qlTop - 1
      radflux = rthraten(k) * exner(k) ! Converts theta/s to temperature/s
      radflux = radflux * cp / grav * ( p(k) - p(k+1) ) ! Converts K/s to W/m^2
      if ( radflux < 0.0 ) F0 = abs(radflux) + F0
   enddo
   F0 = max(F0, 1.0)
   !found Sc cloud and cloud not at surface, trigger downdraft
   if (cloudflg) then

!      !get entrainent coefficient
!      do i=1,NDOWN
!         do k=kts+1,kte
!            ENTf(k,i)=(ZW(k+1)-ZW(k))/L0
!         enddo
!      enddo
!
!      ! get Poisson P(dz/L0)
!      call Poisson(1,NDOWN,kts+1,kte,ENTf,ENTi)


      ! entrainent: Ent=Ent0/dz*P(dz/L0)
      do i=1,NDOWN
         do k=kts+1,kte
!            ENT(k,i)=real(ENTi(k,i))*Ent0/(ZW(k+1)-ZW(k))
            ENT(k,i) = 0.002
            ENT(k,i) = min(ENT(k,i),0.9/(ZW(k+1)-ZW(k)))
         enddo
      enddo

      !!![EW: INVJUMP] find 700mb height then subtract trpospheric lapse rate!!!
      p700_ind = MINLOC(ABS(p-70000),1)!p1D is 70000
      jump_thetav = thv(p700_ind) - thv(1) - (thv(p700_ind)-thv(qlTop+3))/(ZW(p700_ind)-ZW(qlTop+3))*(ZW(p700_ind)-ZW(qlTop))
      jump_qt = qc(p700_ind) + qv(p700_ind) - qc(1) - qv(1)
      jump_thetal = thl(p700_ind) - thl(1) - (thl(p700_ind)-thl(qlTop+3))/(ZW(p700_ind)-ZW(qlTop+3))*(ZW(p700_ind)-ZW(qlTop))

      refTHL = thl(qlTop) !sum(thl(1:qlTop)) / (qlTop) ! avg over BL for now or just at qlTop
      refTHV = thv(qlTop) !sum(thv(1:qlTop)) / (qlTop)
      refQT  = qt(qlTop)  !sum(qt(1:qlTop))  / (qlTop)

      ! wstar_rad, following Lock and MacVean (1999a)
      wst_rad = ( grav * zw(qlTop) * F0 / (refTHL * rho(qlTop) * cp) ) ** (0.333)
      wst_rad = max(wst_rad, 0.1)
      wstar   = max(0.,(grav/thv(1)*wthv*pblh)**(onethird))
      went    = thv(1) / ( grav * jump_thetav * zw(qlTop) ) * &
                (0.15 * (wstar**3 + 5*ust**3) + 0.35 * wst_rad**3 )
      qstar  = abs(went*jump_qt/wst_rad)
      thstar = F0/rho(qlTop)/cp/wst_rad - went*jump_thetav/wst_rad
      !wstar_dd = mixrad + surface wst
      wst_dd = (0.15 * (wstar**3 + 5*ust**3) + 0.35 * wst_rad**3 ) ** (0.333)

      print*,"qstar=",qstar," thstar=",thstar," wst_dd=",wst_dd
      print*,"F0=",F0," wst_rad=",wst_rad," jump_thv=",jump_thetav
      print*,"entrainment velocity=",went

      sigmaW  = 0.2*wst_dd  ! 0.8*wst_dd !wst_rad tuning parameter ! 0.5 was good
      sigmaQT = 40  * qstar ! 50 was good
      sigmaTH = 1.0 * thstar! 0.5 was good

      wmin=sigmaW*pwmin
      wmax=sigmaW*pwmax
      !print*,"sigw=",sigmaW," wmin=",wmin," wmax=",wmax

      do I=1,NDOWN !downdraft now starts at different height
         ki = DD_initK(I)

         wlv=wmin+(wmax-wmin)/REAL(NDOWN)*(i-1)
         wtv=wmin+(wmax-wmin)/REAL(NDOWN)*i

         !DOWNW(ki,I)=0.5*(wlv+wtv)
         DOWNW(ki,I)=wlv
         !DOWNA(ki,I)=0.5*ERF(wtv/(sqrt(2.)*sigmaW))-0.5*ERF(wlv/(sqrt(2.)*sigmaW))
         DOWNA(ki,I)=.1/REAL(NDOWN)
         DOWNU(ki,I)=(u(ki-1)*DZ(ki) + u(ki)*DZ(ki-1)) /(DZ(ki)+DZ(ki-1))
         DOWNV(ki,I)=(v(ki-1)*DZ(ki) + v(ki)*DZ(ki-1)) /(DZ(ki)+DZ(ki-1))

         !reference now depends on where dd starts
!         refTHL = 0.5 * (thl(ki) + thl(ki-1))
!         refTHV = 0.5 * (thv(ki) + thv(ki-1))
!         refQT  = 0.5 * (qt(ki)  + qt(ki-1) )

         refTHL = (thl(ki-1)*DZ(ki) + thl(ki)*DZ(ki-1)) /(DZ(ki)+DZ(ki-1))
         refTHV = (thv(ki-1)*DZ(ki) + thv(ki)*DZ(ki-1)) /(DZ(ki)+DZ(ki-1))
         refQT  = (qt(ki-1)*DZ(ki)  + qt(ki)*DZ(ki-1))  /(DZ(ki)+DZ(ki-1))

         !DOWNQC(ki,I) = 0.0
         DOWNQC(ki,I) = (qc(ki-1)*DZ(ki) + qc(ki)*DZ(ki-1)) /(DZ(ki)+DZ(ki-1))
         DOWNQT(ki,I) = refQT  !+ 0.5  *DOWNW(ki,I)*sigmaQT/sigmaW
         DOWNTHV(ki,I)= refTHV + 0.01 *DOWNW(ki,I)*sigmaTH/sigmaW
         DOWNTHL(ki,I)= refTHL + 0.01 *DOWNW(ki,I)*sigmaTH/sigmaW

         !input :: QT,THV,P,zagl,  output :: THL, QC
!         Pk  =(P(ki-1)*DZ(ki)+P(ki)*DZ(ki-1))/(DZ(ki)+DZ(ki-1))
!         call condensation_edmf_r(DOWNQT(ki,I),   &
!              &        DOWNTHL(ki,I),Pk,ZW(ki),   &
!              &     DOWNTHV(ki,I),DOWNQC(ki,I)    )

      enddo


      !print*, " Begin integration of downdrafts:"
      DO I=1,NDOWN
         !print *, "Plume # =", I,"======================="
         DO k=DD_initK(I)-1,KTS+1,-1
            !starting at the first interface level below cloud top
            !EntExp=exp(-ENT(K,I)*dz(k))
            !EntExp_M=exp(-ENT(K,I)/3.*dz(k))
            EntExp  =ENT(K,I)*dz(k)
            EntExp_M=ENT(K,I)*0.333*dz(k)

            QTn =DOWNQT(k+1,I) *(1.-EntExp) + QT(k)*EntExp
            THLn=DOWNTHL(k+1,I)*(1.-EntExp) + THL(k)*EntExp
            Un  =DOWNU(k+1,I)  *(1.-EntExp) + U(k)*EntExp_M
            Vn  =DOWNV(k+1,I)  *(1.-EntExp) + V(k)*EntExp_M
            !QKEn=DOWNQKE(k-1,I)*(1.-EntExp) + QKE(k)*EntExp

!            QTn =DOWNQT(K+1,I) +(QT(K) -DOWNQT(K+1,I)) *(1.-EntExp)
!            THLn=DOWNTHL(K+1,I)+(THL(K)-DOWNTHL(K+1,I))*(1.-EntExp)
!            Un  =DOWNU(K+1,I)  +(U(K)  -DOWNU(K+1,I))*(1.-EntExp_M)
!            Vn  =DOWNV(K+1,I)  +(V(K)  -DOWNV(K+1,I))*(1.-EntExp_M)

            ! given new p & z, solve for thvn & qcn
            Pk  =(P(k-1)*DZ(k)+P(k)*DZ(k-1))/(DZ(k)+DZ(k-1))
            call condensation_edmf(QTn,THLn,Pk,ZW(k),THVn,QCn)
!            B=grav*(0.5*(THVn+DOWNTHV(k+1,I))/THV(k)-1.)
            THVk  =(THV(k-1)*DZ(k)+THV(k)*DZ(k-1))/(DZ(k)+DZ(k-1))
            B=grav*(THVn/THVk - 1.0)
!            Beta_dm = 2*Wb*ENT(K,I) + 0.5/(ZW(k)-dz(k)) * &
!                 &    max(1. - exp((ZW(k) -dz(k))/Z00 - 1. ) , 0.)
!            EntW=exp(-Beta_dm * dz(k))
            EntW=EntExp
!            if (Beta_dm >0) then
!               Wn2=DOWNW(K+1,I)**2*EntW - Wa*B/Beta_dm * (1. - EntW)
!            else
!               Wn2=DOWNW(K+1,I)**2      - 2.*Wa*B*dz(k)
!            end if

            mindownw = MIN(DOWNW(K+1,I),-0.2)
            Wn = DOWNW(K+1,I) + (-2.*ENT(K,I)*DOWNW(K+1,I) - &
                    BCOEFF*B/mindownw)*MIN(dz(k), 250.)

            !Do not allow a parcel to accelerate more than 1.25 m/s over 200 m.
            !Add max increase of 2.0 m/s for coarse vertical resolution.
            IF (Wn < DOWNW(K+1,I) - MIN(1.25*dz(k)/200., 2.0))THEN
                Wn = DOWNW(K+1,I) - MIN(1.25*dz(k)/200., 2.0)
            ENDIF
            !Add symmetrical max decrease in w
            IF (Wn > DOWNW(K+1,I) + MIN(1.25*dz(k)/200., 2.0))THEN
                Wn = DOWNW(K+1,I) + MIN(1.25*dz(k)/200., 2.0)
            ENDIF
            Wn = MAX(MIN(Wn,0.0), -3.0)

            !print *, "  k       =",      k,      " z    =", ZW(k)
            !print *, "  entw    =",ENT(K,I),     " Bouy =", B
            !print *, "  downthv =",   THVn,      " thvk =", thvk
            !print *, "  downthl =",   THLn,      " thl  =", thl(k)
            !print *, "  downqt  =",   QTn ,      " qt   =", qt(k)
            !print *, "  downw+1 =",DOWNW(K+1,I), " Wn2  =", Wn

            IF (Wn .lt. 0.) THEN !terminate when velocity is too small
               DOWNW(K,I)  = Wn !-sqrt(Wn2)
               DOWNTHV(K,I)= THVn
               DOWNTHL(K,I)= THLn
               DOWNQT(K,I) = QTn
               DOWNQC(K,I) = QCn
               DOWNU(K,I)  = Un
               DOWNV(K,I)  = Vn
               DOWNA(K,I)  = DOWNA(K+1,I)
            ELSE
               !plumes must go at least 2 levels
               if (DD_initK(I) - K .lt. 2) then
                  DOWNW(:,I)  = 0.0
                  DOWNTHV(:,I)= 0.0
                  DOWNTHL(:,I)= 0.0
                  DOWNQT(:,I) = 0.0
                  DOWNQC(:,I) = 0.0
                  DOWNU(:,I)  = 0.0
                  DOWNV(:,I)  = 0.0
               endif
               exit
            ENDIF
         ENDDO
      ENDDO
   endif ! end cloud flag

   DOWNW(1,:) = 0. !make sure downdraft does not go to the surface
   DOWNA(1,:) = 0.

   ! Combine both moist and dry plume, write as one averaged plume
   ! Even though downdraft starts at different height, average all up to qlTop
   DO k=qlTop,KTS,-1
      DO I=1,NDOWN
         IF (I > NDOWN) exit
         edmf_a_dd(K)  =edmf_a_dd(K)  +DOWNA(K-1,I)
         edmf_w_dd(K)  =edmf_w_dd(K)  +DOWNA(K-1,I)*DOWNW(K-1,I)
         edmf_qt_dd(K) =edmf_qt_dd(K) +DOWNA(K-1,I)*DOWNQT(K-1,I)
         edmf_thl_dd(K)=edmf_thl_dd(K)+DOWNA(K-1,I)*DOWNTHL(K-1,I)
         edmf_ent_dd(K)=edmf_ent_dd(K)+DOWNA(K-1,I)*ENT(K-1,I)
         edmf_qc_dd(K) =edmf_qc_dd(K) +DOWNA(K-1,I)*DOWNQC(K-1,I)
      ENDDO

      IF (edmf_a_dd(k) >0.) THEN
          edmf_w_dd(k)  =edmf_w_dd(k)  /edmf_a_dd(k)
          edmf_qt_dd(k) =edmf_qt_dd(k) /edmf_a_dd(k)
          edmf_thl_dd(k)=edmf_thl_dd(k)/edmf_a_dd(k)
          edmf_ent_dd(k)=edmf_ent_dd(k)/edmf_a_dd(k)
          edmf_qc_dd(k) =edmf_qc_dd(k) /edmf_a_dd(k)
      ENDIF
   ENDDO

   !
   ! computing variables needed for solver
   !

   DO k=KTS,qlTop
      rho_int = (rho(k)*DZ(k+1)+rho(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
      DO I=1,NDOWN
         sd_aw(k)   =sd_aw(k)   +rho_int*DOWNA(k,i)*DOWNW(k,i)
         sd_awthl(k)=sd_awthl(k)+rho_int*DOWNA(k,i)*DOWNW(k,i)*DOWNTHL(k,i)
         sd_awqt(k) =sd_awqt(k) +rho_int*DOWNA(k,i)*DOWNW(k,i)*DOWNQT(k,i)
         sd_awqc(k) =sd_awqc(k) +rho_int*DOWNA(k,i)*DOWNW(k,i)*DOWNQC(k,i)
         sd_awu(k)  =sd_awu(k)  +rho_int*DOWNA(k,i)*DOWNW(k,i)*DOWNU(k,i)
         sd_awv(k)  =sd_awv(k)  +rho_int*DOWNA(k,i)*DOWNW(k,i)*DOWNV(k,i)
      ENDDO
      sd_awqv(k) = sd_awqt(k)  - sd_awqc(k)
   ENDDO

END SUBROUTINE DDMF_JPL
!===============================================================

END MODULE
