!>\file module_bl_mynn.F90
!! This file contains the entity of MYNN-EDMF PBL scheme.
! **********************************************************************
! *   An improved Mellor-Yamada turbulence closure model               *
! *                                                                    *
! *      Original author: M. Nakanishi (N.D.A), naka@nda.ac.jp         *
! *      Translated into F90 and implemented in WRF-ARW by:            *
! *                       Mariusz Pagowski (NOAA-GSL)                  *
! *      Subsequently developed by:                                    *
! *                 Joseph Olson, Jaymes Kenyon (NOAA/GSL),            *
! *                 Wayne Angevine (NOAA/CSL), Kay Suselj (NASA/JPL),  *
! *                 Franciano Puhales (UFSM), Laura Fowler (NCAR),     *
! *                 Elynn Wu (UCSD), and Jordan Schnell (NOAA/GSL)     *
! *                                                                    *
! *   Contents:                                                        *
! *                                                                    *
! *   mynn_bl_driver - main subroutine which calls all other routines  *
! *   --------------                                                   *
! *     1. mym_initialize  (to be called once initially)               *
! *        gives the closure constants and initializes the turbulent   *
! *        quantities.                                                 *
! *     2. get_pblh                                                    *
! *        Calculates the boundary layer height                        *
! *     3. scale_aware                                                 *
! *        Calculates scale-adaptive tapering functions                *
! *     4. mym_condensation                                            *
! *        determines the liquid water content and the cloud fraction  *
! *        diagnostically.                                             *
! *     5. dmp_mf                                                      *
! *        Calls the (nonlocal) mass-flux component                    *
! *     6. ddmf_jpl                                                    *
! *        Calls the downdraft mass-flux component                     *
! *    (-) mym_level2      (called in the other subroutines)           *
! *        calculates the stability functions at Level 2.              *
! *    (-) mym_length      (called in the other subroutines)           *
! *        calculates the master length scale.                         *
! *     7. mym_turbulence                                              *
! *        calculates the vertical diffusivity coefficients and the    *
! *        production terms for the turbulent quantities.              *
! *     8. mym_predict                                                 *
! *        predicts the turbulent quantities at the next step.         *
! *                                                                    *
! *             call mym_initialize                                    *
! *                  |                                                 *
! *                  |<----------------+                               *
! *                  |                 |                               *
! *             call get_pblh          |                               *
! *             call scale_aware       |                               *
! *             call mym_condensation  |                               *
! *             call dmp_mf            |                               *
! *             call ddmf_jpl          |                               *
! *             call mym_turbulence    |                               *
! *             call mym_predict       |                               *
! *                  |                 |                               *
! *                  |-----------------+                               *
! *                  |                                                 *
! *                 end                                                *
! *                                                                    *
! *   Variables worthy of special mention:                             *
! *     tref   : Reference temperature                                 *
! *     thl    : Liquid water potential temperature                    *
! *     qw     : Total water (water vapor+liquid water) content        *
! *     ql     : Liquid water content                                  *
! *     vt, vq : Functions for computing the buoyancy flux             *
! *     qke    : 2 * TKE                                               *
! *     el     : mixing length                                         *
! *                                                                    *
! *     If the water contents are unnecessary, e.g., in the case of    *
! *     ocean models, thl is the potential temperature and qw, ql, vt  *
! *     and vq are all zero.                                           *
! *                                                                    *
! *   Grid arrangement:                                                *
! *             k+1 +---------+                                        *
! *                 |         |     i = 1 - nx                         *
! *             (k) |    *    |     k = 1 - nz                         *
! *                 |         |                                        *
! *              k  +---------+                                        *
! *                 i   (i)  i+1                                       *
! *                                                                    *
! *     All the predicted variables are defined at the center (*) of   *
! *     the grid boxes. The diffusivity coefficients and two of their  *
! *     components (el and stability functions sh & sm) are, however,  *
! *     defined on the walls of the grid boxes.                        *
! *     # Upper boundary values are given at k=nz.                     *
! *                                                                    *
! *   References:                                                      *
! *     1. Nakanishi, M., 2001:                                        *
! *        Boundary-Layer Meteor., 99, 349-378.                        *
! *     2. Nakanishi, M. and H. Niino, 2004:                           *
! *        Boundary-Layer Meteor., 112, 1-31.                          *
! *     3. Nakanishi, M. and H. Niino, 2006:                           *
! *        Boundary-Layer Meteor., 119, 397-407.                       *
! *     4. Nakanishi, M. and H. Niino, 2009:                           *
! *        Jour. Meteor. Soc. Japan, 87, 895-912.                      *
! *     5. Olson J. and coauthors, 2019: A description of the          *
! *        MYNN-EDMF scheme and coupling to other components in        *
! *        WRF-ARW. NOAA Tech. Memo. OAR GSD, 61, 37 pp.,              *
! *        https://doi.org/10.25923/n9wm-be49.                         * 
! *     6. Puhales, Franciano S. and coauthors, 2020: Turbulent        *
! *        Kinetic Energy Budget for MYNN-EDMF PBL Scheme in WRF model.*
! *        Universidade Federal de Santa Maria Technical Note. 9 pp.   *
! **********************************************************************
! ==================================================================
! Notes on original implementation into WRF-ARW
! changes to original code:
! 1. code is 1D (in z)
! 2. option to advect TKE, but not the covariances and variances
! 3. Cranck-Nicholson replaced with the implicit scheme
! 4. removed terrain-dependent grid since input in WRF in actual
!    distances in z[m]
! 5. cosmetic changes to adhere to WRF standard (remove common blocks,
!            intent etc)
!-------------------------------------------------------------------
! Further modifications post-implementation
!
! 1. Addition of BouLac mixing length in the free atmosphere.
! 2. Changed the turbulent mixing length to be integrated from the
!    surface to the top of the BL + a transition layer depth.
! v3.4.1:    Option to use Kitamura/Canuto modification which removes 
!            the critical Richardson number and negative TKE (default).
!            Hybrid PBL height diagnostic, which blends a theta-v-based
!            definition in neutral/convective BL and a TKE-based definition
!            in stable conditions.
!            TKE budget output option (bl_mynn_tkebudget)
! v3.5.0:    TKE advection option (bl_mynn_tkeadvect)
! v3.5.1:    Fog deposition related changes.
! v3.6.0:    Removed fog deposition from the calculation of tendencies
!            Added mixing of qc, qi, qni
!            Added output for wstar, delta, TKE_PBL, & KPBL for correct 
!                   coupling to shcu schemes  
! v3.8.0:    Added subgrid scale cloud output for coupling to radiation
!            schemes (activated by setting icloud_bl =1 in phys namelist).
!            Added WRF_DEBUG prints (at level 3000)
!            Added Tripoli and Cotton (1981) correction.
!            Added namelist option bl_mynn_cloudmix to test effect of mixing
!                cloud species (default = 1: on). 
!            Added mass-flux option (bl_mynn_edmf, = 1 for DMP mass-flux, 0: off).
!                Related options: 
!                 bl_mynn_edmf_mom = 1 : activate momentum transport in MF scheme
!                 bl_mynn_edmf_tke = 1 : activate TKE transport in MF scheme
!            Added mixing length option (bl_mynn_mixlength, see notes below)
!            Added more sophisticated saturation checks, following Thompson scheme
!            Added new cloud PDF option (bl_mynn_cloudpdf = 2) from Chaboureau
!                and Bechtold (2002, JAS, with mods) 
!            Added capability to mix chemical species when env variable
!                WRF_CHEM = 1, thanks to Wayne Angevine.
!            Added scale-aware mixing length, following Junshi Ito's work
!                Ito et al. (2015, BLM).
! v3.9.0    Improvement to the mass-flux scheme (dynamic number of plumes,
!                better plume/cloud depth, significant speed up, better cloud
!                fraction). 
!            Added Stochastic Parameter Perturbation (SPP) implementation.
!            Many miscellaneous tweaks to the mixing lengths and stratus
!                component of the subgrid clouds.
! v.4.0      Removed or added alternatives to WRF-specific functions/modules
!                for the sake of portability to other models.
!                the sake of portability to other models.
!            Further refinement of mass-flux scheme from SCM experiments with
!                Wayne Angevine: switch to linear entrainment and back to
!                Simpson and Wiggert-type w-equation.
!            Addition of TKE production due to radiation cooling at top of 
!                clouds (proto-version); not activated by default.
!            Some code rewrites to move if-thens out of loops in an attempt to
!                improve computational efficiency.
!            New tridiagonal solver, which is supposedly 14% faster and more
!                conservative. Impact seems very small.
!            Many miscellaneous tweaks to the mixing lengths and stratus
!                component of the subgrid-scale (SGS) clouds.
! v4.1       Big improvements in downward SW radiation due to revision of subgrid clouds
!                - better cloud fraction and subgrid scale mixing ratios.
!                - may experience a small cool bias during the daytime now that high 
!                  SW-down bias is greatly reduced...
!            Some tweaks to increase the turbulent mixing during the daytime for
!                bl_mynn_mixlength option 2 to alleviate cool bias (very small impact).
!            Improved ensemble spread from changes to SPP in MYNN
!                - now perturbing eddy diffusivity and eddy viscosity directly
!                - now perturbing background rh (in SGS cloud calc only)
!                - now perturbing entrainment rates in mass-flux scheme
!            Added IF checks (within IFDEFS) to protect mixchem code from being used
!                when HRRR smoke is used (no impact on regular non-wrf chem use)
!            Important bug fix for wrf chem when transporting chemical species in MF scheme
!            Removed 2nd mass-flux scheme (no only bl_mynn_edmf = 1, no option 2)
!            Removed unused stochastic code for mass-flux scheme
!            Changed mass-flux scheme to be integrated on interface levels instead of
!                mass levels - impact is small
!            Added option to mix 2nd moments in MYNN as opposed to the scalar_pblmix option.
!                - activated with bl_mynn_mixscalars = 1; this sets scalar_pblmix = 0
!                - added tridagonal solver used in scalar_pblmix option to duplicate tendencies
!                - this alone changes the interface call considerably from v4.0.
!            Slight revision to TKE production due to radiation cooling at top of clouds
!            Added the non-Guassian buoyancy flux function of Bechtold and Siebesma (1998, JAS).
!                - improves TKE in SGS clouds
!            Added heating due to dissipation of TKE (small impact, maybe + 0.1 C daytime PBL temp)
!            Misc changes made for FV3/MPAS compatibility
! v4.2       A series of small tweaks to help reduce a cold bias in the PBL:
!                - slight increase in diffusion in convective conditions
!                - relaxed criteria for mass-flux activation/strength
!                - added capability to cycle TKE for continuity in hourly updating HRRR
!                - added effects of compensational environmental subsidence in mass-flux scheme,
!                  which resulted in tweaks to detrainment rates.
!            Bug fix for diagnostic-decay of SGS clouds - noticed by Greg Thompson. This has
!                a very small, but primarily  positive, impact on SW-down biases.
!            Tweak to calculation of KPBL - urged by Laura Fowler - to make more intuitive.
!            Tweak to temperature range of blending for saturation check (water to ice). This
!                slightly reduces excessive SGS clouds in polar region. No impact warm clouds. 
!            Added namelist option bl_mynn_output (0 or 1) to suppress or activate the
!                allocation and output of 10 3D variables. Most people will want this
!                set to 0 (default) to save memory and disk space.
!            Added new array qi_bl as opposed to using qc_bl for both SGS qc and qi. This
!                gives us more control of the magnitudes which can be confounded by using
!                a single array. As a results, many subroutines needed to be modified,
!                especially mym_condensation.
!            Added the blending of the stratus component of the SGS clouds to the mass-flux
!                clouds to account for situations where stratus and cumulus may exist in the
!                grid cell.
!            Misc small-impact bugfixes:
!                1) dz was incorrectly indexed in mym_condensation
!                2) configurations with icloud_bl = 0 were using uninitialized arrays
! v4.5 / CCPP
!            This version includes many modifications that proved valuable in the global
!            framework and removes some key lingering bugs in the mixing of chemical species.
!            TKE Budget output fixed (Puhales, 2020-12)
!            New option for stability function: (Puhales, 2020-12)
!                bl_mynn_stfunc = 0 (original, Kansas-type function, Paulson, 1970 )
!                bl_mynn_stfunc = 1 (expanded range, same as used for Jimenez et al (MWR)
!                see the Technical Note for this implementation (small impact).
!            Improved conservation of momentum and higher-order moments.
!            Important bug fixes for mixing of chemical species.
!            Addition of pressure-gradient effects on updraft momentum transport.
!            Addition of bl_mynn_closure option = 2.5, 2.6, or 3.0
!            Addition of higher-order moments for sigma when using 
!                bl_mynn_cloudpdf = 2 (Chab-Becht).
!            Removed WRF_CHEM dependencies.
!            Many miscellaneous tweaks.
!
! Many of these changes are now documented in references listed above.
!====================================================================

MODULE module_bl_mynn

  use module_bl_mynn_common,only: &
        cp        , cpv       , cliq       , cice      , &
        p608      , ep_2      , ep_3       , gtr       , &
        grav      , g_inv     , karman     , p1000mb   , &
        rcp       , r_d       , r_v        , rk        , &
        rvovrd    , svp1      , svp2       , svp3      , &
        xlf       , xlv       , xls        , xlscp     , &
        xlvcp     , tv0       , tv1        , tref      , &
        zero      , half      , one        , two       , &
        onethird  , twothirds , tkmin      , t0c       , &
        tice      , kind_phys

  use mynn_functions
  use mym_init
  use mynn_tend
  use myn_ddmf_jpl
  use mynn_dmp_mf
  use mynn_condensation
  use myn_condensation_edmf
  use myn_tridiag
  use myn_mix_chem
  use mynn_bolouac_length
  use mynn_level2
  use myn_topdown_cloudrad
  use mynn_bolouac_length0
  use myn_coeffs
!  use myn_scale_aware
  use turbulence
  !use myn_check
  use mynn_predict

  use iso_c_binding
  IMPLICIT NONE


!===================================================================
!INTERFACES
interface
subroutine SCALE_AWARE(dx, PBL1, Psig_bl, Psig_shcu)
import c_float
real(kind=c_float), intent(in) :: dx, PBL1
real(kind=c_float), intent(inout) ::Psig_bl,Psig_shcu 
end subroutine SCALE_AWARE
end interface

!===================================================================
! From here on, these are MYNN-specific parameters:
! The parameters below depend on stability functions of module_sf_mynn.
  REAL, PARAMETER :: cphm_st=5.0, cphm_unst=16.0, &
                     cphh_st=5.0, cphh_unst=16.0

! Constants for min tke in elt integration (qmin), max z/L in els (zmax), 
! and factor for eddy viscosity for TKE (Kq = Sqfac*Km):
  REAL, PARAMETER :: qmin=0.0, zmax=1.0, Sqfac=3.0
! Note that the following mixing-length constants are now specified in mym_length
!      &cns=3.5, alp1=0.23, alp2=0.3, alp3=3.0, alp4=10.0, alp5=0.2

  REAL, PARAMETER :: gpw=5./3., qcgmin=1.e-8, qkemin=1.e-12
  REAL, PARAMETER :: tliq = 269. !all hydrometeors are liquid when T > tliq

! Constants for cloud PDF (mym_condensation)
  REAL, PARAMETER :: rr2=0.7071068, rrp=0.3989423

  !>Use Canuto/Kitamura mod (remove Ric and negative TKE) (1:yes, 0:no)
  !!For more info, see Canuto et al. (2008 JAS) and Kitamura (Journal of the 
  !!Meteorological Society of Japan, Vol. 88, No. 5, pp. 857-864, 2010).
  !!Note that this change required further modification of other parameters
  !!above (c2, c3). If you want to remove this option, set c2 and c3 constants 
  !!(above) back to NN2009 values (see commented out lines next to the
  !!parameters above). This only removes the negative TKE problem
  !!but does not necessarily improve performance - neutral impact.
  REAL, PARAMETER :: CKmod=1.

  !>Use Ito et al. (2015, BLM) scale-aware (0: no, 1: yes). Note that this also has impacts
  !!on the cloud PDF and mass-flux scheme, using Honnert et al. (2011) similarity function
  !!for TKE in the upper PBL/cloud layer.
  REAL, PARAMETER :: scaleaware=1.

  !>Of the following the options, use one OR the other, not both.
  !>Adding top-down diffusion driven by cloud-top radiative cooling
  INTEGER, PARAMETER :: bl_mynn_topdown = 0
  !>Option to activate downdrafts, from Elynn Wu (0: deactive, 1: active)
  INTEGER, PARAMETER :: bl_mynn_edmf_dd = 0

  !>Option to activate heating due to dissipation of TKE (to activate, set to 1.0)
  INTEGER, PARAMETER :: dheat_opt = 1

  !Option to activate environmental subsidence in mass-flux scheme
  LOGICAL, PARAMETER :: env_subs = .false.

  !Option to switch flux-profile relationship for surface (from Puhales et al. 2020)
  !0: use original Dyer-Hicks, 1: use Cheng-Brustaert and Blended COARE
  INTEGER, PARAMETER :: bl_mynn_stfunc = 1

  !option to print out more stuff for debugging purposes
  LOGICAL, PARAMETER :: debug_code = .false.
  INTEGER, PARAMETER :: idbg = 23 !specific i-point to write out
  
  ! Used in WRF-ARW module_physics_init.F
  INTEGER :: mynn_level


CONTAINS

! ==================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine is the GSD MYNN-EDNF PBL driver routine,which
!! encompassed the majority of the subroutines that comprise the 
!! procedures that ultimately solve for tendencies of 
!! \f$U, V, \theta, q_v, q_c, and q_i\f$.
!!\section gen_mynn_bl_driver GSD mynn_bl_driver General Algorithm
!> @{
  SUBROUTINE mynn_bl_driver(            &
       &initflag,restart,cycling,       &
       &delt,dz,dx,znt,                 &
       &u,v,w,th,sqv3d,sqc3d,sqi3d,     &
       &qnc,qni,                        &
       &qnwfa,qnifa,qnbca,ozone,        &
       &p,exner,rho,t3d,                &
       &xland,ts,qsfc,ps,               &
       &ust,ch,hfx,qfx,rmol,wspd,       &
       &uoce,voce,                      & !ocean current
       &qke,qke_adv,                    &
       &sh3d,sm3d,                      &
       &nchem,kdvel,ndvel,              & !smoke/chem variables
       &chem3d,vdep,                    &
       &frp,emis_ant_no,                &
       &mix_chem,enh_mix,               & !note: these arrays/flags are still under development
       &rrfs_sd,smoke_dbg,              & !end smoke/chem variables
       &tsq,qsq,cov,                    &
       &rublten,rvblten,rthblten,       &
       &rqvblten,rqcblten,rqiblten,     &
       &rqncblten,rqniblten,            &
       &rqnwfablten,rqnifablten,        &
       &rqnbcablten,dozone,             &
       &exch_h,exch_m,                  &
       &pblh,kpbl,                      & 
       &el_pbl,                         &
       &dqke,qwt,qshear,qbuoy,qdiss,    &
       &qc_bl,qi_bl,cldfra_bl,          &
       &bl_mynn_tkeadvect,              &
       &tke_budget,                     &
       &bl_mynn_cloudpdf,               &
       &bl_mynn_mixlength,              &
       &icloud_bl,                      &
       &closure,                        &
       &bl_mynn_edmf,                   &
       &bl_mynn_edmf_mom,               &
       &bl_mynn_edmf_tke,               &
       &bl_mynn_mixscalars,             &
       &bl_mynn_output,                 &
       &bl_mynn_cloudmix,bl_mynn_mixqt, &
       &edmf_a,edmf_w,edmf_qt,          &
       &edmf_thl,edmf_ent,edmf_qc,      &
       &sub_thl3D,sub_sqv3D,            &
       &det_thl3D,det_sqv3D,            &
       &nupdraft,maxMF,ktop_plume,      &
       &spp_pbl,pattern_spp_pbl,        &
       &rthraten,                       &
       &FLAG_QC,FLAG_QI,FLAG_QNC,       &
       &FLAG_QNI,FLAG_QNWFA,FLAG_QNIFA, &
       &FLAG_QNBCA,                     &
       &IDS,IDE,JDS,JDE,KDS,KDE,        &
       &IMS,IME,JMS,JME,KMS,KME,        &
       &ITS,ITE,JTS,JTE,KTS,KTE         )
    
!-------------------------------------------------------------------

    INTEGER, INTENT(in) :: initflag
    !INPUT NAMELIST OPTIONS:
    LOGICAL, INTENT(in) :: restart,cycling
    INTEGER, INTENT(in) :: tke_budget
    INTEGER, INTENT(in) :: bl_mynn_cloudpdf
    INTEGER, INTENT(in) :: bl_mynn_mixlength
    INTEGER, INTENT(in) :: bl_mynn_edmf
    LOGICAL, INTENT(in) :: bl_mynn_tkeadvect
    INTEGER, INTENT(in) :: bl_mynn_edmf_mom
    INTEGER, INTENT(in) :: bl_mynn_edmf_tke
    INTEGER, INTENT(in) :: bl_mynn_mixscalars
    INTEGER, INTENT(in) :: bl_mynn_output
    INTEGER, INTENT(in) :: bl_mynn_cloudmix
    INTEGER, INTENT(in) :: bl_mynn_mixqt
    INTEGER, INTENT(in) :: icloud_bl
    REAL,    INTENT(in) :: closure

    LOGICAL, INTENT(in) :: FLAG_QI,FLAG_QNI,FLAG_QC,FLAG_QNC,&
                           FLAG_QNWFA,FLAG_QNIFA,FLAG_QNBCA

    LOGICAL, INTENT(IN) :: mix_chem,enh_mix,rrfs_sd,smoke_dbg

    INTEGER, INTENT(in) :: &
         & IDS,IDE,JDS,JDE,KDS,KDE &
         &,IMS,IME,JMS,JME,KMS,KME &
         &,ITS,ITE,JTS,JTE,KTS,KTE

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

! initflag > 0  for TRUE
! else        for FALSE
!       closure       : <= 2.5;  Level 2.5
!                  2.5< and <3;  Level 2.6
!                        =   3;  Level 3
    
    REAL, INTENT(in) :: delt
    REAL, DIMENSION(IMS:IME), INTENT(in) :: dx
    REAL, DIMENSION(IMS:IME,KMS:KME), INTENT(in) :: dz,      &
         &u,v,w,th,sqv3D,p,exner,rho,t3d
    REAL, DIMENSION(IMS:IME,KMS:KME), OPTIONAL, INTENT(in):: &
         &sqc3D,sqi3D,qni,qnc,qnwfa,qnifa,qnbca
    REAL, DIMENSION(IMS:IME,KMS:KME), OPTIONAL, INTENT(in):: ozone
    REAL, DIMENSION(IMS:IME), INTENT(in) :: xland,ust,       &
         &ch,ts,qsfc,ps,hfx,qfx,wspd,uoce,voce,znt

    REAL, DIMENSION(IMS:IME,KMS:KME), INTENT(inout) ::       &
         &qke,tsq,qsq,cov,qke_adv

    REAL, DIMENSION(IMS:IME,KMS:KME), INTENT(inout) ::       &
         &rublten,rvblten,rthblten,rqvblten,rqcblten,        &
         &rqiblten,rqniblten,rqncblten,                      &
         &rqnwfablten,rqnifablten,rqnbcablten
    REAL, DIMENSION(IMS:IME,KMS:KME), INTENT(inout) :: dozone

    REAL, DIMENSION(IMS:IME,KMS:KME), INTENT(in)    :: rthraten

    REAL, DIMENSION(IMS:IME,KMS:KME), INTENT(out)   ::       &
         &exch_h,exch_m

   !These 10 arrays are only allocated when bl_mynn_output > 0
   REAL, DIMENSION(IMS:IME,KMS:KME), OPTIONAL, INTENT(inout) :: &
         & edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc,  &
         & sub_thl3D,sub_sqv3D,det_thl3D,det_sqv3D

!   REAL, DIMENSION(IMS:IME,KMS:KME)   :: &
!         & edmf_a_dd,edmf_w_dd,edmf_qt_dd,edmf_thl_dd,edmf_ent_dd,edmf_qc_dd

    REAL, DIMENSION(IMS:IME), INTENT(inout) :: pblh,rmol

    REAL, DIMENSION(IMS:IME) :: psig_bl,psig_shcu

    INTEGER,DIMENSION(IMS:IME),INTENT(INOUT) ::             &
         &kpbl,nupdraft,ktop_plume

    REAL, DIMENSION(IMS:IME), INTENT(OUT) ::                &
         &maxmf

    REAL, DIMENSION(IMS:IME,KMS:KME), INTENT(inout) ::      &
         &el_pbl

    REAL, DIMENSION(IMS:IME,KMS:KME), optional, INTENT(out) :: &
         &qwt,qshear,qbuoy,qdiss,dqke
    ! 3D budget arrays are not allocated when tke_budget == 0
    ! 1D (local) budget arrays are used for passing between subroutines.
    REAL, DIMENSION(kts:kte) :: qwt1,qshear1,qbuoy1,qdiss1, &
         &dqke1,diss_heat

    REAL, DIMENSION(IMS:IME,KMS:KME), intent(out) :: Sh3D,Sm3D

    REAL, DIMENSION(IMS:IME,KMS:KME), INTENT(inout) ::      &
         &qc_bl,qi_bl,cldfra_bl
    REAL, DIMENSION(KTS:KTE) :: qc_bl1D,qi_bl1D,cldfra_bl1D,&
                    qc_bl1D_old,qi_bl1D_old,cldfra_bl1D_old

! smoke/chemical arrays
    INTEGER, INTENT(IN   ) ::   nchem, kdvel, ndvel
    REAL,    DIMENSION(ims:ime, kms:kme, nchem), INTENT(INOUT), optional :: chem3d
    REAL,    DIMENSION(ims:ime, ndvel),   INTENT(IN),   optional :: vdep
    REAL,    DIMENSION(ims:ime),     INTENT(IN),    optional :: frp,EMIS_ANT_NO
    !local
    REAL,    DIMENSION(kts:kte  ,nchem) :: chem1
    REAL,    DIMENSION(kts:kte+1,nchem) :: s_awchem1
    REAL,    DIMENSION(ndvel)           :: vd1
    INTEGER :: ic

!local vars
    INTEGER :: ITF,JTF,KTF, IMD,JMD
    INTEGER :: i,j,k,kproblem
    REAL, DIMENSION(KTS:KTE) :: thl,thvl,tl,qv1,qc1,qi1,sqw,&
         &el, dfm, dfh, dfq, tcd, qcd, pdk, pdt, pdq, pdc,  &
         &vt, vq, sgm, thlsg, sqwsg
    REAL, DIMENSION(KTS:KTE) :: thetav,sh,sm,u1,v1,w1,p1,   &
         &ex1,dz1,th1,tk1,rho1,qke1,tsq1,qsq1,cov1,         &
         &sqv,sqi,sqc,du1,dv1,dth1,dqv1,dqc1,dqi1,ozone1,   &
         &k_m1,k_h1,qni1,dqni1,qnc1,dqnc1,qnwfa1,qnifa1,    &
         &qnbca1,dqnwfa1,dqnifa1,dqnbca1,dozone1

    !mass-flux variables
    REAL, DIMENSION(KTS:KTE) :: dth1mf,dqv1mf,dqc1mf,du1mf,dv1mf
    REAL, DIMENSION(KTS:KTE) :: edmf_a1,edmf_w1,edmf_qt1,   &
         &edmf_thl1,edmf_ent1,edmf_qc1
    REAL, DIMENSION(KTS:KTE) :: edmf_a_dd1,edmf_w_dd1,      &
         &edmf_qt_dd1,edmf_thl_dd1,                         &
         &edmf_ent_dd1,edmf_qc_dd1
    REAL, DIMENSION(KTS:KTE) :: sub_thl,sub_sqv,sub_u,sub_v,&
                        det_thl,det_sqv,det_sqc,det_u,det_v
    REAL,DIMENSION(KTS:KTE+1) :: s_aw1,s_awthl1,s_awqt1,    &
                  s_awqv1,s_awqc1,s_awu1,s_awv1,s_awqke1,   &
                  s_awqnc1,s_awqni1,s_awqnwfa1,s_awqnifa1,  &
                  s_awqnbca1
    REAL,DIMENSION(KTS:KTE+1) :: sd_aw1,sd_awthl1,sd_awqt1, &
                  sd_awqv1,sd_awqc1,sd_awu1,sd_awv1,sd_awqke1

    REAL, DIMENSION(KTS:KTE+1) :: zw
    REAL :: cpm,sqcg,flt,fltv,flq,flqv,flqc,pmz,phh,exnerg,zet,phi_m,&
          & afk,abk,ts_decay, qc_bl2, qi_bl2,                        &
          & th_sfc,ztop_plume,sqc9,sqi9,wsp

    !top-down diffusion
    REAL, DIMENSION(ITS:ITE) :: maxKHtopdown
    REAL, DIMENSION(KTS:KTE) :: KHtopdown,TKEprodTD

    LOGICAL :: INITIALIZE_QKE,problem

    ! Stochastic fields 
    INTEGER,  INTENT(IN)                                     ::spp_pbl
    REAL, DIMENSION( ims:ime, kms:kme), INTENT(IN),OPTIONAL  ::pattern_spp_pbl
    REAL, DIMENSION(KTS:KTE)                                 ::rstoch_col

    ! Substepping TKE
    INTEGER :: nsub
    real    :: delt2


    if (debug_code) then !check incoming values
      do i=its,ite
        problem = .false.
        do k=kts,kte
          wsp  = sqrt(u(i,k)**2 + v(i,k)**2)
          if (abs(hfx(i)) > 1200. .or. abs(qfx(i)) > 0.001 .or.         &
              wsp > 200. .or. t3d(i,k) > 360. .or. t3d(i,k) < 160. .or. &
              sqv3d(i,k)< 0.0 .or. sqc3d(i,k)< 0.0 ) then
             kproblem = k
             problem = .true.
             print*,"Incoming problem at: i=",i," k=1"
             print*," QFX=",qfx(i)," HFX=",hfx(i)
             print*," wsp=",wsp," T=",t3d(i,k)
             print*," qv=",sqv3d(i,k)," qc=",sqc3d(i,k)
             print*," u*=",ust(i)," wspd=",wspd(i)
             print*," xland=",xland(i)," ts=",ts(i)
             print*," z/L=",0.5*dz(i,1)*rmol(i)," ps=",ps(i)
             print*," znt=",znt(i)," dx=",dx(i)
          endif
        enddo
        if (problem) then
          print*,"===tk:",t3d(i,max(kproblem-3,1):min(kproblem+3,kte))
          print*,"===qv:",sqv3d(i,max(kproblem-3,1):min(kproblem+3,kte))
          print*,"===qc:",sqc3d(i,max(kproblem-3,1):min(kproblem+3,kte))
          print*,"===qi:",sqi3d(i,max(kproblem-3,1):min(kproblem+3,kte))
          print*,"====u:",u(i,max(kproblem-3,1):min(kproblem+3,kte))
          print*,"====v:",v(i,max(kproblem-3,1):min(kproblem+3,kte))
        endif
      enddo
    endif

!***  Begin debugging
    IMD=(IMS+IME)/2
    JMD=(JMS+JME)/2
!***  End debugging 

    JTF=JTE
    ITF=ITE
    KTF=KTE

    IF (bl_mynn_output > 0) THEN !research mode
       edmf_a(its:ite,kts:kte)=0.
       edmf_w(its:ite,kts:kte)=0.
       edmf_qt(its:ite,kts:kte)=0.
       edmf_thl(its:ite,kts:kte)=0.
       edmf_ent(its:ite,kts:kte)=0.
       edmf_qc(its:ite,kts:kte)=0.
       sub_thl3D(its:ite,kts:kte)=0.
       sub_sqv3D(its:ite,kts:kte)=0.
       det_thl3D(its:ite,kts:kte)=0.
       det_sqv3D(its:ite,kts:kte)=0.

       !edmf_a_dd(its:ite,kts:kte)=0.
       !edmf_w_dd(its:ite,kts:kte)=0.
       !edmf_qt_dd(its:ite,kts:kte)=0.
       !edmf_thl_dd(its:ite,kts:kte)=0.
       !edmf_ent_dd(its:ite,kts:kte)=0.
       !edmf_qc_dd(its:ite,kts:kte)=0.
    ENDIF
    ktop_plume(its:ite)=0   !int
    nupdraft(its:ite)=0     !int
    maxmf(its:ite)=0.
    maxKHtopdown(its:ite)=0.

    ! DH* CHECK HOW MUCH OF THIS INIT IF-BLOCK IS ACTUALLY NEEDED FOR RESTARTS
!> - Within the MYNN-EDMF, there is a dependecy check for the first time step,
!! If true, a three-dimensional initialization loop is entered. Within this loop,
!! several arrays are initialized and k-oriented (vertical) subroutines are called 
!! at every i and j point, corresponding to the x- and y- directions, respectively.  
    IF (initflag > 0 .and. .not.restart) THEN

       !Test to see if we want to initialize qke
       IF ( (restart .or. cycling)) THEN
          IF (MAXVAL(QKE(its:ite,kts)) < 0.0002) THEN
             INITIALIZE_QKE = .TRUE.
             !print*,"QKE is too small, must initialize"
          ELSE
             INITIALIZE_QKE = .FALSE.
             !print*,"Using background QKE, will not initialize"
          ENDIF
       ELSE ! not cycling or restarting:
          INITIALIZE_QKE = .TRUE.
          !print*,"not restart nor cycling, must initialize QKE"
       ENDIF
 
       if (.not.restart .or. .not.cycling) THEN
         Sh3D(its:ite,kts:kte)=0.
         Sm3D(its:ite,kts:kte)=0.
         el_pbl(its:ite,kts:kte)=0.
         tsq(its:ite,kts:kte)=0.
         qsq(its:ite,kts:kte)=0.
         cov(its:ite,kts:kte)=0.
         cldfra_bl(its:ite,kts:kte)=0.
         qc_bl(its:ite,kts:kte)=0.
         qke(its:ite,kts:kte)=0.
       else
         qc_bl1D(kts:kte)=0.0
         qi_bl1D(kts:kte)=0.0
         cldfra_bl1D(kts:kte)=0.0
       end if
       dqc1(kts:kte)=0.0
       dqi1(kts:kte)=0.0
       dqni1(kts:kte)=0.0
       dqnc1(kts:kte)=0.0
       dqnwfa1(kts:kte)=0.0
       dqnifa1(kts:kte)=0.0
       dqnbca1(kts:kte)=0.0
       dozone1(kts:kte)=0.0
       qc_bl1D_old(kts:kte)=0.0
       cldfra_bl1D_old(kts:kte)=0.0
       edmf_a1(kts:kte)=0.0
       edmf_w1(kts:kte)=0.0
       edmf_qc1(kts:kte)=0.0
       edmf_a_dd1(kts:kte)=0.0
       edmf_w_dd1(kts:kte)=0.0
       edmf_qc_dd1(kts:kte)=0.0
       sgm(kts:kte)=0.0
       vt(kts:kte)=0.0
       vq(kts:kte)=0.0

       DO k=KTS,KTE
          DO i=ITS,ITF
             exch_m(i,k)=0.
             exch_h(i,k)=0.
          ENDDO
       ENDDO

       IF (tke_budget .eq. 1) THEN
          DO k=KTS,KTE
             DO i=ITS,ITF
                qWT(i,k)=0.
                qSHEAR(i,k)=0.
                qBUOY(i,k)=0.
                qDISS(i,k)=0.
                dqke(i,k)=0.
             ENDDO
          ENDDO
       ENDIF

       DO i=ITS,ITF
          DO k=KTS,KTE !KTF
                dz1(k)=dz(i,k)
                u1(k) = u(i,k)
                v1(k) = v(i,k)
                w1(k) = w(i,k)
                th1(k)=th(i,k)
                tk1(k)=T3D(i,k)
                ex1(k)=exner(i,k)
                rho1(k)=rho(i,k)
                sqc(k)=sqc3D(i,k) !/(1.+qv(i,k))
                sqv(k)=sqv3D(i,k) !/(1.+qv(i,k))
                thetav(k)=th(i,k)*(1.+0.608*sqv(k))
                IF (icloud_bl > 0) THEN
                   CLDFRA_BL1D(k)=CLDFRA_BL(i,k)
                   QC_BL1D(k)=QC_BL(i,k)
                   QI_BL1D(k)=QI_BL(i,k)
                ENDIF
                IF (FLAG_QI ) THEN
                   sqi(k)=sqi3D(i,k) !/(1.+qv(i,k))
                   sqw(k)=sqv(k)+sqc(k)+sqi(k)
                   thl(k)=th1(k) - xlvcp/ex1(k)*sqc(k) &
                       &         - xlscp/ex1(k)*sqi(k)
                   !Use form from Tripoli and Cotton (1981) with their
                   !suggested min temperature to improve accuracy.
                   !thl(k)=th(i,k)*(1.- xlvcp/MAX(tk1(k),TKmin)*sqc(k) &
                   !    &               - xlscp/MAX(tk1(k),TKmin)*sqi(k))
                   !COMPUTE THL USING SGS CLOUDS FOR PBLH DIAG
                   IF(sqc(k)<1e-6 .and. sqi(k)<1e-8 .and. CLDFRA_BL1D(k)>0.001)THEN
                      sqc9=QC_BL1D(k)*CLDFRA_BL1D(k)
                      sqi9=QI_BL1D(k)*CLDFRA_BL1D(k)
                   ELSE
                      sqc9=sqc(k)
                      sqi9=sqi(k)
                   ENDIF
                   thlsg(k)=th1(k) - xlvcp/ex1(k)*sqc9 &
                         &         - xlscp/ex1(k)*sqi9
                   sqwsg(k)=sqv(k)+sqc9+sqi9
                ELSE
                   sqi(k)=0.0
                   sqw(k)=sqv(k)+sqc(k)
                   thl(k)=th1(k)-xlvcp/ex1(k)*sqc(k)
                   !Use form from Tripoli and Cotton (1981) with their 
                   !suggested min temperature to improve accuracy.      
                   !thl(k)=th(i,k)*(1.- xlvcp/MAX(tk1(k),TKmin)*sqc(k))
                   !COMPUTE THL USING SGS CLOUDS FOR PBLH DIAG
                   IF(sqc(k)<1e-6 .and. CLDFRA_BL1D(k)>0.001)THEN
                            sqc9=QC_BL1D(k)*CLDFRA_BL1D(k)
                      sqi9=0.0
                   ELSE
                      sqc9=sqc(k)
                      sqi9=0.0
                   ENDIF
                   thlsg(k)=th1(k) - xlvcp/ex1(k)*sqc9 &
                         &         - xlscp/ex1(k)*sqi9
                   sqwsg(k)=sqv(k)+sqc9+sqi9
                ENDIF
                thvl(k)=thlsg(k)*(1.+0.61*sqv(k))

                IF (k==kts) THEN
                   zw(k)=0.
                ELSE
                   zw(k)=zw(k-1)+dz(i,k-1)
                ENDIF
                IF (INITIALIZE_QKE) THEN
                   !Initialize tke for initial PBLH calc only - using 
                   !simple PBLH form of Koracin and Berkowicz (1988, BLM)
                   !to linearly taper off tke towards top of PBL.
                   qke1(k)=5.*ust(i) * MAX((ust(i)*700. - zw(k))/(MAX(ust(i),0.01)*700.), 0.01)
                ELSE
                   qke1(k)=qke(i,k)
                ENDIF
                el(k)=el_pbl(i,k)
                sh(k)=Sh3D(i,k)
                sm(k)=Sm3D(i,k)
                tsq1(k)=tsq(i,k)
                qsq1(k)=qsq(i,k)
                cov1(k)=cov(i,k)
                if (spp_pbl==1) then
                    rstoch_col(k)=pattern_spp_pbl(i,k)
                else
                    rstoch_col(k)=0.0
                endif

             ENDDO

             zw(kte+1)=zw(kte)+dz(i,kte)

!>  - Call get_pblh() to calculate hybrid (\f$\theta_{vli}-TKE\f$) PBL height.
             CALL GET_PBLH(KTS,KTE,PBLH(i),thetav,&
               &  Qke1,zw,dz1,xland(i),KPBL(i))
             
!>  - Call scale_aware() to calculate similarity functions for scale-adaptive control
!! (\f$P_{\sigma-PBL}\f$ and \f$P_{\sigma-shcu}\f$).
             IF (scaleaware > 0.) THEN
                CALL SCALE_AWARE(dx(i),PBLH(i),Psig_bl(i),Psig_shcu(i))
             ELSE
                Psig_bl(i)=1.0
                Psig_shcu(i)=1.0
             ENDIF

             ! DH* CHECK IF WE CAN DO WITHOUT CALLING THIS ROUTINE FOR RESTARTS
!>  - Call mym_initialize() to initializes the mixing length, TKE, \f$\theta^{'2}\f$,
!! \f$q^{'2}\f$, and \f$\theta^{'}q^{'}\f$. These variables are calculated after 
!! obtaining prerequisite variables by calling the following subroutines from 
!! within mym_initialize(): mym_level2() and mym_length().
             CALL mym_initialize (                & 
                  &kts,kte,xland(i),              &
                  &dz1, dx(i), zw,                &
                  &u1, v1, thl, sqv,              &
                  &thlsg, sqwsg,                  &
                  &PBLH(i), th1, thetav, sh, sm,  &
                  &ust(i), rmol(i),               &
                  &el, Qke1, Tsq1, Qsq1, Cov1,    &
                  &Psig_bl(i), cldfra_bl1D,       &
                  &bl_mynn_mixlength,             &
                  &edmf_w1,edmf_a1,               &
                  &INITIALIZE_QKE,                &
                  &spp_pbl,rstoch_col,            &
                  &b1,b2,qkemin,CKmod             )

             IF (.not.restart) THEN
                !UPDATE 3D VARIABLES
                DO k=KTS,KTE !KTF
                   el_pbl(i,k)=el(k)
                   sh3d(i,k)=sh(k)
                   sm3d(i,k)=sm(k)
                   qke(i,k)=qke1(k)
                   tsq(i,k)=tsq1(k)
                   qsq(i,k)=qsq1(k)
                   cov(i,k)=cov1(k)
                ENDDO
                !initialize qke_adv array if using advection
                IF (bl_mynn_tkeadvect) THEN
                   DO k=KTS,KTE
                      qke_adv(i,k)=qke1(k)
                   ENDDO
                ENDIF
             ENDIF

!***  Begin debugging
!             IF(I==IMD .AND. J==JMD)THEN
!               PRINT*,"MYNN DRIVER INIT: k=",1," sh=",sh(k)
!               PRINT*," sqw=",sqw(k)," thl=",thl(k)," k_m=",exch_m(i,k)
!               PRINT*," xland=",xland(i)," rmol=",rmol(i)," ust=",ust(i)
!               PRINT*," qke=",qke(i,k)," el=",el_pbl(i,k)," tsq=",Tsq(i,k)
!               PRINT*," PBLH=",PBLH(i)," u=",u(i,k)," v=",v(i,k)
!             ENDIF
!***  End debugging

       ENDDO !end i-loop

    ENDIF ! end initflag

!> - After initializing all required variables, the regular procedures 
!! performed at every time step are ready for execution.
    !ACF- copy qke_adv array into qke if using advection
    IF (bl_mynn_tkeadvect) THEN
       qke=qke_adv
    ENDIF

    DO i=ITS,ITF
       DO k=KTS,KTE !KTF
            !JOE-TKE BUDGET
             IF (tke_budget .eq. 1) THEN
                dqke(i,k)=qke(i,k)
             END IF
             IF (icloud_bl > 0) THEN
                CLDFRA_BL1D(k)=CLDFRA_BL(i,k)
                QC_BL1D(k)=QC_BL(i,k)
                QI_BL1D(k)=QI_BL(i,k)
                cldfra_bl1D_old(k)=cldfra_bl(i,k)
                qc_bl1D_old(k)=qc_bl(i,k)
                qi_bl1D_old(k)=qi_bl(i,k)
             else
                CLDFRA_BL1D(k)=0.0
                QC_BL1D(k)=0.0
                QI_BL1D(k)=0.0
                cldfra_bl1D_old(k)=0.0
                qc_bl1D_old(k)=0.0
                qi_bl1D_old(k)=0.0
             ENDIF
             dz1(k)= dz(i,k)
             u1(k) = u(i,k)
             v1(k) = v(i,k)
             w1(k) = w(i,k)
             th1(k)= th(i,k)
             tk1(k)=T3D(i,k)
             p1(k) = p(i,k)
             ex1(k)= exner(i,k)
             rho1(k)=rho(i,k)
             sqv(k)= sqv3D(i,k) !/(1.+qv(i,k))
             sqc(k)= sqc3D(i,k) !/(1.+qv(i,k))
             qv1(k)= sqv(k)/(1.-sqv(k))
             qc1(k)= sqc(k)/(1.-sqv(k))
             dqc1(k)=0.0
             dqi1(k)=0.0
             dqni1(k)=0.0
             dqnc1(k)=0.0
             dqnwfa1(k)=0.0
             dqnifa1(k)=0.0
             dqnbca1(k)=0.0
             dozone1(k)=0.0
             IF(FLAG_QI)THEN
                sqi(k)= sqi3D(i,k) !/(1.+qv(i,k))
                qi1(k)= sqi(k)/(1.-sqv(k))
                sqw(k)= sqv(k)+sqc(k)+sqi(k)
                thl(k)= th1(k) - xlvcp/ex1(k)*sqc(k) &
                     &         - xlscp/ex1(k)*sqi(k)
                !Use form from Tripoli and Cotton (1981) with their
                !suggested min temperature to improve accuracy.    
                !thl(k)=th(i,k)*(1.- xlvcp/MAX(tk1(k),TKmin)*sqc(k) &
                !    &               - xlscp/MAX(tk1(k),TKmin)*sqi(k))
                !COMPUTE THL USING SGS CLOUDS FOR PBLH DIAG
                IF(sqc(k)<1e-6 .and. sqi(k)<1e-8 .and. CLDFRA_BL1D(k)>0.001)THEN
                   sqc9=QC_BL1D(k)*CLDFRA_BL1D(k)
                   sqi9=QI_BL1D(k)*CLDFRA_BL1D(k)
                ELSE
                   sqc9=sqc(k)
                   sqi9=sqi(k)
                ENDIF
                thlsg(k)=th1(k) - xlvcp/ex1(k)*sqc9 &
                      &         - xlscp/ex1(k)*sqi9
                sqwsg(k)=sqv(k)+sqc9+sqi9
             ELSE
                qi1(k)=0.0
                sqi(k)=0.0
                sqw(k)= sqv(k)+sqc(k)
                thl(k)= th1(k)-xlvcp/ex1(k)*sqc(k)
                !Use form from Tripoli and Cotton (1981) with their
                !suggested min temperature to improve accuracy.    
                !thl(k)=th(i,k)*(1.- xlvcp/MAX(tk1(k),TKmin)*sqc(k))
                !COMPUTE THL USING SGS CLOUDS FOR PBLH DIAG
                IF(sqc(k)<1e-6 .and. CLDFRA_BL1D(k)>0.001)THEN
                   sqc9=QC_BL1D(k)*CLDFRA_BL1D(k)
                   sqi9=QI_BL1D(k)*CLDFRA_BL1D(k)
                ELSE
                   sqc9=sqc(k)
                   sqi9=0.0
                ENDIF
                thlsg(k)=th1(k) - xlvcp/ex1(k)*sqc9 &
                      &         - xlscp/ex1(k)*sqi9 
            ENDIF
            thetav(k)=th1(k)*(1.+0.608*sqv(k))
            thvl(k)  =thlsg(k) *(1.+0.608*sqv(k))

             IF (FLAG_QNI ) THEN
                qni1(k)=qni(i,k)
             ELSE
                qni1(k)=0.0
             ENDIF
             IF (FLAG_QNC ) THEN
                qnc1(k)=qnc(i,k)
             ELSE
                qnc1(k)=0.0
             ENDIF
             IF (FLAG_QNWFA ) THEN
                qnwfa1(k)=qnwfa(i,k)
             ELSE
                qnwfa1(k)=0.0
             ENDIF
             IF (FLAG_QNIFA ) THEN
                qnifa1(k)=qnifa(i,k)
             ELSE
                qnifa1(k)=0.0
             ENDIF
             IF (FLAG_QNBCA .and. PRESENT(qnbca)) THEN
                qnbca1(k)=qnbca(i,k)
             ELSE
                qnbca1(k)=0.0
             ENDIF
             IF (PRESENT(ozone)) THEN
                ozone1(k)=ozone(i,k)
             ELSE
                ozone1(k)=0.0
             ENDIF
             el(k) = el_pbl(i,k)
             qke1(k)=qke(i,k)
             sh(k)  =sh3d(i,k)
             sm(k)  =sm3d(i,k)
             tsq1(k)=tsq(i,k)
             qsq1(k)=qsq(i,k)
             cov1(k)=cov(i,k)
             if (spp_pbl==1) then
                rstoch_col(k)=pattern_spp_pbl(i,k)
             else
                rstoch_col(k)=0.0
             endif

             !edmf
             edmf_a1(k)=0.0
             edmf_w1(k)=0.0
             edmf_qc1(k)=0.0
             s_aw1(k)=0.
             s_awthl1(k)=0.
             s_awqt1(k)=0.
             s_awqv1(k)=0.
             s_awqc1(k)=0.
             s_awu1(k)=0.
             s_awv1(k)=0.
             s_awqke1(k)=0.
             s_awqnc1(k)=0.
             s_awqni1(k)=0.
             s_awqnwfa1(k)=0.
             s_awqnifa1(k)=0.
             s_awqnbca1(k)=0.
             ![EWDD]
             edmf_a_dd1(k)=0.0
             edmf_w_dd1(k)=0.0
             edmf_qc_dd1(k)=0.0
             sd_aw1(k)=0.
             sd_awthl1(k)=0.
             sd_awqt1(k)=0.
             sd_awqv1(k)=0.
             sd_awqc1(k)=0.
             sd_awu1(k)=0.
             sd_awv1(k)=0.
             sd_awqke1(k)=0.
             sub_thl(k)=0.
             sub_sqv(k)=0.
             sub_u(k)=0.
             sub_v(k)=0.
             det_thl(k)=0.
             det_sqv(k)=0.
             det_sqc(k)=0.
             det_u(k)=0.
             det_v(k)=0.

             IF (k==kts) THEN
                zw(k)=0.
             ELSE
                zw(k)=zw(k-1)+dz(i,k-1)
             ENDIF
          ENDDO ! end k

          !initialize smoke/chem arrays (if used):
             if  ( mix_chem ) then
                do ic = 1,ndvel
                   vd1(ic) = vdep(i,ic) ! dry deposition velocity
                enddo
                do k = kts,kte
                   do ic = 1,nchem
                      chem1(k,ic) = chem3d(i,k,ic)
                      s_awchem1(k,ic)=0.
                   enddo
                enddo
             else
                do ic = 1,ndvel
                   vd1(ic) = 0. ! dry deposition velocity
                enddo
                do k = kts,kte
                   do ic = 1,nchem
                      chem1(k,ic) = 0.
                      s_awchem1(k,ic)=0.
                   enddo
                enddo
             endif

          zw(kte+1)=zw(kte)+dz(i,kte)
          !EDMF
          s_aw1(kte+1)=0.
          s_awthl1(kte+1)=0.
          s_awqt1(kte+1)=0.
          s_awqv1(kte+1)=0.
          s_awqc1(kte+1)=0.
          s_awu1(kte+1)=0.
          s_awv1(kte+1)=0.
          s_awqke1(kte+1)=0.
          s_awqnc1(kte+1)=0.
          s_awqni1(kte+1)=0.
          s_awqnwfa1(kte+1)=0.
          s_awqnifa1(kte+1)=0.
          s_awqnbca1(kte+1)=0.
          sd_aw1(kte+1)=0.
          sd_awthl1(kte+1)=0.
          sd_awqt1(kte+1)=0.
          sd_awqv1(kte+1)=0.
          sd_awqc1(kte+1)=0.
          sd_awu1(kte+1)=0.
          sd_awv1(kte+1)=0.
          sd_awqke1(kte+1)=0.
          IF ( mix_chem ) THEN
             DO ic = 1,nchem
                s_awchem1(kte+1,ic)=0.
             ENDDO
          ENDIF

!>  - Call get_pblh() to calculate the hybrid \f$\theta_{vli}-TKE\f$
!! PBL height diagnostic.
          CALL GET_PBLH(KTS,KTE,PBLH(i),thetav,&
          & Qke1,zw,dz1,xland(i),KPBL(i))

!>  - Call scale_aware() to calculate the similarity functions,
!! \f$P_{\sigma-PBL}\f$ and \f$P_{\sigma-shcu}\f$, to control 
!! the scale-adaptive behaviour for the local and nonlocal 
!! components, respectively.
          IF (scaleaware > 0.) THEN
             CALL SCALE_AWARE(dx(i),PBLH(i),Psig_bl(i),Psig_shcu(i))
          ELSE
             Psig_bl(i)=1.0
             Psig_shcu(i)=1.0
          ENDIF

          sqcg= 0.0   !ill-defined variable; qcg has been removed
          cpm=cp*(1.+0.84*qv1(kts))
          exnerg=(ps(i)/p1000mb)**rcp

          !-----------------------------------------------------
          !ORIGINAL CODE
          !flt = hfx(i)/( rho(i,kts)*cpm ) &
          ! +xlvcp*ch(i)*(sqc(kts)/exner(i,kts) -sqcg/exnerg)
          !flq = qfx(i)/  rho(i,kts)       &
          !    -ch(i)*(sqc(kts)   -sqcg )
          !-----------------------------------------------------
          flqv   = qfx(i)/rho1(kts)
          flqc   = 0.0 !currently no sea-spray fluxes, fog settling hangled elsewhere
          th_sfc = ts(i)/ex1(kts)

          ! TURBULENT FLUX FOR TKE BOUNDARY CONDITIONS
          flq =flqv+flqc                  !! LATENT
          flt =hfx(i)/(rho1(kts)*cpm )-xlvcp*flqc/ex1(kts)  !! Temperature flux
          fltv=flt + flqv*p608*th_sfc     !! Virtual temperature flux

          ! Update 1/L using updated sfc heat flux and friction velocity
          rmol(i) = -karman*gtr*fltv/max(ust(i)**3,1.0e-6)
          zet = 0.5*dz(i,kts)*rmol(i)
          zet = MAX(zet, -20.)
          zet = MIN(zet,  20.)
          !if(i.eq.idbg)print*,"updated z/L=",zet
          if (bl_mynn_stfunc == 0) then
             !Original Kansas-type stability functions
             if ( zet >= 0.0 ) then
                pmz = 1.0 + (cphm_st-1.0) * zet
                phh = 1.0 +  cphh_st      * zet
             else
                pmz = 1.0/    (1.0-cphm_unst*zet)**0.25 - zet
                phh = 1.0/SQRT(1.0-cphh_unst*zet)
             end if
          else
             !Updated stability functions (Puhales, 2020)
             phi_m = calc_phim(zet,cphh_unst,cphm_unst)
             pmz   = phi_m - zet
             phh   = calc_phih(zet,cphh_unst)
          end if

!>  - Call mym_condensation() to calculate the nonconvective component
!! of the subgrid cloud fraction and mixing ratio as well as the functions
!! used to calculate the buoyancy flux. Different cloud PDFs can be
!! selected by use of the namelist parameter \p bl_mynn_cloudpdf.

          CALL  mym_condensation ( kts,kte,      &
               &dx(i),dz1,zw,xland(i),           &
               &thl,sqw,sqv,sqc,sqi,             &
               &p1,ex1,tsq1,qsq1,cov1,           &
               &Sh,el,bl_mynn_cloudpdf,          &
               &qc_bl1D,qi_bl1D,cldfra_bl1D,     &
               &PBLH(i),HFX(i),                  &
               &Vt, Vq, th1, sgm, rmol(i),       &
               &spp_pbl, rstoch_col,             &
               &b2,rrp,rr2,tliq                  )

!>  - Add TKE source driven by cloud top cooling
!!  Calculate the buoyancy production of TKE from cloud-top cooling when
!! \p bl_mynn_topdown =1.
          IF (bl_mynn_topdown.eq.1)then
             CALL topdown_cloudrad(kts,kte,dz1,zw,          &
                &xland(i),kpbl(i),PBLH(i),                  &
                &sqc,sqi,sqw,thl,th1,ex1,p1,rho1,thetav,    &
                &cldfra_bl1D,rthraten(i,:),                 &
                &maxKHtopdown(i),KHtopdown,TKEprodTD        )
          ELSE
             maxKHtopdown(i)  = 0.0
             KHtopdown(kts:kte) = 0.0
             TKEprodTD(kts:kte) = 0.0
          ENDIF

          IF (bl_mynn_edmf > 0) THEN
            !PRINT*,"Calling DMP Mass-Flux: i= ",i
            CALL DMP_mf(                          &
               &kts,kte,delt,zw,dz1,p1,rho1,      &
               &bl_mynn_edmf_mom,                 &
               &bl_mynn_edmf_tke,                 &
               &bl_mynn_mixscalars,               &
               &u1,v1,w1,th1,thl,thetav,tk1,      &
               &sqw,sqv,sqc,qke1,                 &
               &qnc1,qni1,qnwfa1,qnifa1,qnbca1,   &
               &ex1,Vt,Vq,sgm,                    &
               &ust(i),flt,fltv,flq,flqv,         &
               &PBLH(i),KPBL(i),DX(i),            &
               &xland(i),th_sfc,                  &
            ! now outputs - tendencies
            ! &,dth1mf,dqv1mf,dqc1mf,du1mf,dv1mf  &
            ! outputs - updraft properties
               & edmf_a1,edmf_w1,edmf_qt1,        &
               & edmf_thl1,edmf_ent1,edmf_qc1,    &
            ! for the solver
               & s_aw1,s_awthl1,s_awqt1,          &
               & s_awqv1,s_awqc1,                 &
               & s_awu1,s_awv1,s_awqke1,          &
               & s_awqnc1,s_awqni1,               &
               & s_awqnwfa1,s_awqnifa1,s_awqnbca1,&
               & sub_thl,sub_sqv,                 &
               & sub_u,sub_v,                     &
               & det_thl,det_sqv,det_sqc,         &
               & det_u,det_v,                     &
            ! chem/smoke mixing
               & nchem,chem1,s_awchem1,           &
               & mix_chem,                        &
               & qc_bl1D,cldfra_bl1D,             &
               & qc_bl1D_old,cldfra_bl1D_old,     &
               & FLAG_QC,FLAG_QI,                 &
               & FLAG_QNC,FLAG_QNI,               &
               & FLAG_QNWFA,FLAG_QNIFA,FLAG_QNBCA,&
               & Psig_shcu(i),                    &
               & nupdraft(i),ktop_plume(i),       &
               & maxmf(i),ztop_plume,             &
               & spp_pbl,rstoch_col,              &
               & env_subs                         )
          ENDIF

          IF (bl_mynn_edmf_dd == 1) THEN
            CALL DDMF_JPL(kts,kte,delt,zw,dz1,p1, &
              &u1,v1,th1,thl,thetav,tk1,          &
              sqw,sqv,sqc,rho1,ex1,               &
              &ust(i),flt,flq,                    &
              &PBLH(i),KPBL(i),                   &
              &edmf_a_dd1,edmf_w_dd1,edmf_qt_dd1, &
              &edmf_thl_dd1,edmf_ent_dd1,         &
              &edmf_qc_dd1,                       &
              &sd_aw1,sd_awthl1,sd_awqt1,         &
              &sd_awqv1,sd_awqc1,sd_awu1,sd_awv1, &
              &sd_awqke1,                         &
              &qc_bl1d,cldfra_bl1d,               &
              &rthraten(i,:)                      )
          ENDIF

          !Capability to substep the eddy-diffusivity portion
          !do nsub = 1,2
          delt2 = delt !*0.5    !only works if topdown=0

          CALL mym_turbulence (                  & 
               &kts,kte,xland(i),closure,        &
               &dz1, DX(i), zw,                  &
               &u1, v1, thl, thetav, sqc, sqw,   &
               &thlsg, sqwsg,                    &
               &qke1, tsq1, qsq1, cov1,          &
               &vt, vq,                          &
               &rmol(i), flt, flq,               &
               &PBLH(i),th1,                     &
               &Sh,Sm,el,                        &
               &Dfm,Dfh,Dfq,                     &
               &Tcd,Qcd,Pdk,                     &
               &Pdt,Pdq,Pdc,                     &
               &qWT1,qSHEAR1,qBUOY1,qDISS1,      &
               &tke_budget,                      &
               &Psig_bl(i),Psig_shcu(i),         &
               &cldfra_bl1D,bl_mynn_mixlength,   &
               &edmf_w1,edmf_a1,                 &
               &TKEprodTD,                       &
               &spp_pbl,rstoch_col,              &
               &qmin, zmax, CKmod,               &                            
               &a1,a2,b1,b2,c1,cc3,              &   
               &debug_code,                      &
               &e1c,e2c,e3c,e4c,e5c              )

!>  - Call mym_predict() to solve TKE and 
!! \f$\theta^{'2}, q^{'2}, and \theta^{'}q^{'}\f$
!! for the following time step.
          CALL mym_predict (kts,kte,closure,      &
               &delt2, dz1,                       &
               &ust(i), flt, flq, pmz, phh,       &
               &el, dfq, rho1, pdk, pdt, pdq, pdc,&
               &Qke1, Tsq1, Qsq1, Cov1,           &
               &s_aw1, s_awqke1, bl_mynn_edmf_tke,&
               &qWT1, qDISS1,tke_budget,          &
               &b1,b2,sqfac                       ) !! TKE budget  (Puhales, 2020)

          if (dheat_opt > 0) then
             DO k=kts,kte-1
                ! Set max dissipative heating rate to 7.2 K per hour
                diss_heat(k) = MIN(MAX(1.0*(qke1(k)**1.5)/(b1*MAX(0.5*(el(k)+el(k+1)),1.))/cp, 0.0),0.002)
                ! Limit heating above 100 mb:
                diss_heat(k) = diss_heat(k) * exp(-10000./MAX(p1(k),1.)) 
             ENDDO
             diss_heat(kte) = 0.
          else
             diss_heat(1:kte) = 0.
          endif

!>  - Call mynn_tendencies() to solve for tendencies of 
!! \f$U, V, \theta, q_{v}, q_{c}, and q_{i}\f$.
          CALL mynn_tendencies(kts,kte,i,        &
               &delt, dz1, rho1,                 &
               &u1, v1, th1, tk1, qv1,           &
               &qc1, qi1, qnc1, qni1,            &
               &ps(i), p1, ex1, thl,             &
               &sqv, sqc, sqi, sqw,              &
               &qnwfa1, qnifa1, qnbca1, ozone1,  &
               &ust(i),flt,flq,flqv,flqc,        &
               &wspd(i),uoce(i),voce(i),         &
               &tsq1, qsq1, cov1,                &
               &tcd, qcd,                        &
               &dfm, dfh, dfq,                   &
               &Du1, Dv1, Dth1, Dqv1,            &
               &Dqc1, Dqi1, Dqnc1, Dqni1,        &
               &Dqnwfa1, Dqnifa1, Dqnbca1,       &
               &Dozone1,                         &
               &diss_heat,                       &
               ! mass flux components
               &s_aw1,s_awthl1,s_awqt1,          &
               &s_awqv1,s_awqc1,s_awu1,s_awv1,   &
               &s_awqnc1,s_awqni1,               &
               &s_awqnwfa1,s_awqnifa1,s_awqnbca1,&
               &sd_aw1,sd_awthl1,sd_awqt1,       &
               &sd_awqv1,sd_awqc1,               &
               sd_awu1,sd_awv1,                  &
               &sub_thl,sub_sqv,                 &
               &sub_u,sub_v,                     &
               &det_thl,det_sqv,det_sqc,         &
               &det_u,det_v,                     &
               &FLAG_QC,FLAG_QI,FLAG_QNC,        &
               &FLAG_QNI,FLAG_QNWFA,FLAG_QNIFA,  &
               &FLAG_QNBCA,                      &
               &cldfra_bl1d,                     &
               &bl_mynn_cloudmix,                &
               &bl_mynn_mixqt,                   &
               &bl_mynn_edmf,                    &
               &bl_mynn_edmf_mom,                &
               &bl_mynn_mixscalars,              &
               &debug_code                       )


          IF ( mix_chem ) THEN
            IF ( rrfs_sd ) THEN 
             CALL mynn_mix_chem(kts,kte,i,       &
                  &delt, dz1, pblh(i),           &
                  &nchem, kdvel, ndvel,          &
                  &chem1, vd1,                   &
                  &rho1,flt,                     &
                  &tcd, qcd,                     &
                  &dfh,                          &
                  &s_aw1,s_awchem1,              &
                  &emis_ant_no(i),               &
                  &frp(i), rrfs_sd,              &
                  &enh_mix, smoke_dbg            )
             ELSE
              CALL mynn_mix_chem(kts,kte,i,       &
                   &delt, dz1, pblh(i),           &
                   &nchem, kdvel, ndvel,          &
                   &chem1, vd1,                   &
                   &rho1,flt,                     &
                   &tcd, qcd,                     &
                   &dfh,                          &
                   &s_aw1,s_awchem1,              &
                   &zero,                         &
                   &zero, rrfs_sd,                &
                   &enh_mix, smoke_dbg            )
             ENDIF
             DO ic = 1,nchem
                DO k = kts,kte
                   chem3d(i,k,ic) = max(1.e-12, chem1(k,ic))
                ENDDO
             ENDDO
          ENDIF
 
          CALL retrieve_exchange_coeffs(kts,kte,&
               &dfm, dfh, dz1, K_m1, K_h1)

          !UPDATE 3D ARRAYS
          do k=kts,kte
             exch_m(i,k)=K_m1(k)
             exch_h(i,k)=K_h1(k)
             rublten(i,k)=du1(k)
             rvblten(i,k)=dv1(k)
             rthblten(i,k)=dth1(k)
             rqvblten(i,k)=dqv1(k)
             if (bl_mynn_cloudmix > 0) then
               if (present(sqc3D) .and. flag_qc) rqcblten(i,k)=dqc1(k)
               if (present(sqi3D) .and. flag_qi) rqiblten(i,k)=dqi1(k)
             else
               if (present(sqc3D) .and. flag_qc) rqcblten(i,k)=0.
               if (present(sqi3D) .and. flag_qi) rqiblten(i,k)=0.
             endif
             if (bl_mynn_cloudmix > 0 .and. bl_mynn_mixscalars > 0) then
               if (present(qnc) .and. flag_qnc) rqncblten(i,k)=dqnc1(k)
               if (present(qni) .and. flag_qni) rqniblten(i,k)=dqni1(k)
               if (present(qnwfa) .and. flag_qnwfa) rqnwfablten(i,k)=dqnwfa1(k)
               if (present(qnifa) .and. flag_qnifa) rqnifablten(i,k)=dqnifa1(k)
               if (present(qnbca) .and. flag_qnbca) rqnbcablten(i,k)=dqnbca1(k)
             else
               if (present(qnc) .and. flag_qnc) rqncblten(i,k)=0.
               if (present(qni) .and. flag_qni) rqniblten(i,k)=0.
               if (present(qnwfa) .and. flag_qnwfa) rqnwfablten(i,k)=0.
               if (present(qnifa) .and. flag_qnifa) rqnifablten(i,k)=0.
               if (present(qnbca) .and. flag_qnbca) rqnbcablten(i,k)=0.
             endif
             dozone(i,k)=dozone1(k)

             if (icloud_bl > 0) then
                qc_bl(i,k)=qc_bl1D(k)
                qi_bl(i,k)=qi_bl1D(k)
                cldfra_bl(i,k)=cldfra_bl1D(k)
             endif

             el_pbl(i,k)=el(k)
             qke(i,k)=qke1(k)
             tsq(i,k)=tsq1(k)
             qsq(i,k)=qsq1(k)
             cov(i,k)=cov1(k)
             sh3d(i,k)=sh(k)
             sm3d(i,k)=sm(k)
          enddo !end-k

          if (tke_budget .eq. 1) then
             !! TKE budget is now given in m**2/s**-3 (Puhales, 2020)
             !! Lower boundary condtions (using similarity relationships such as the prognostic equation for Qke)
             k=kts
             qSHEAR1(k)=4.*(ust(i)**3*phi_m/(karman*dz(i,k)))-qSHEAR1(k+1) !! staggered
             qBUOY1(k)=4.*(-ust(i)**3*zet/(karman*dz(i,k)))-qBUOY1(k+1) !! staggered
             !! unstaggering SHEAR and BUOY and trasfering all TKE budget to 3D array               
             do k = kts,kte-1
                qSHEAR(i,k)=0.5*(qSHEAR1(k)+qSHEAR1(k+1)) !!! unstaggering in z
                qBUOY(i,k)=0.5*(qBUOY1(k)+qBUOY1(k+1)) !!! unstaggering in z
                qWT(i,k)=qWT1(k)
                qDISS(i,k)=qDISS1(k)
                dqke(i,k)=(qke1(k)-dqke(i,k))*0.5/delt
             enddo
             !! Upper boundary conditions               
             k=kte
             qSHEAR(i,k)=0.
             qBUOY(i,k)=0.
             qWT(i,k)=0.
             qDISS(i,k)=0.
             dqke(i,k)=0.
          endif

          !update updraft/downdraft properties
          if (bl_mynn_output > 0) THEN !research mode == 1
             if (bl_mynn_edmf > 0) THEN
                DO k = kts,kte
                   edmf_a(i,k)=edmf_a1(k)
                   edmf_w(i,k)=edmf_w1(k)
                   edmf_qt(i,k)=edmf_qt1(k)
                   edmf_thl(i,k)=edmf_thl1(k)
                   edmf_ent(i,k)=edmf_ent1(k)
                   edmf_qc(i,k)=edmf_qc1(k)
                   sub_thl3D(i,k)=sub_thl(k)
                   sub_sqv3D(i,k)=sub_sqv(k)
                   det_thl3D(i,k)=det_thl(k)
                   det_sqv3D(i,k)=det_sqv(k)
                ENDDO
             endif
!             if (bl_mynn_edmf_dd > 0) THEN
!                DO k = kts,kte
!                   edmf_a_dd(i,k)=edmf_a_dd1(k)
!                   edmf_w_dd(i,k)=edmf_w_dd1(k)
!                   edmf_qt_dd(i,k)=edmf_qt_dd1(k)
!                   edmf_thl_dd(i,k)=edmf_thl_dd1(k)
!                   edmf_ent_dd(i,k)=edmf_ent_dd1(k)
!                   edmf_qc_dd(i,k)=edmf_qc_dd1(k)
!                ENDDO
!             ENDIF
          ENDIF

          !***  Begin debug prints
          IF ( debug_code .and. (i .eq. idbg)) THEN
             IF ( ABS(QFX(i))>.001)print*,&
                "SUSPICIOUS VALUES AT: i=",i," QFX=",QFX(i)
             IF ( ABS(HFX(i))>1100.)print*,&
                "SUSPICIOUS VALUES AT: i=",i," HFX=",HFX(i)
             DO k = kts,kte
               IF ( sh(k) < 0. .OR. sh(k)> 200.)print*,&
                  "SUSPICIOUS VALUES AT: i,k=",i,k," sh=",sh(k)
               IF ( ABS(vt(k)) > 2.0 )print*,&
                  "SUSPICIOUS VALUES AT: i,k=",i,k," vt=",vt(k)
               IF ( ABS(vq(k)) > 7000.)print*,&
                  "SUSPICIOUS VALUES AT: i,k=",i,k," vq=",vq(k)
               IF ( qke(i,k) < -1. .OR. qke(i,k)> 200.)print*,&
                  "SUSPICIOUS VALUES AT: i,k=",i,k," qke=",qke(i,k)
               IF ( el_pbl(i,k) < 0. .OR. el_pbl(i,k)> 1500.)print*,&
                  "SUSPICIOUS VALUES AT: i,k=",i,k," el_pbl=",el_pbl(i,k)
               IF ( exch_m(i,k) < 0. .OR. exch_m(i,k)> 2000.)print*,&
                  "SUSPICIOUS VALUES AT: i,k=",i,k," exxch_m=",exch_m(i,k)
               IF (icloud_bl > 0) then
                  IF( cldfra_bl(i,k) < 0.0 .OR. cldfra_bl(i,k)> 1.)THEN
                  PRINT*,"SUSPICIOUS VALUES: CLDFRA_BL=",cldfra_bl(i,k)," qc_bl=",QC_BL(i,k)
                  ENDIF
               ENDIF

               !IF (I==IMD .AND. J==JMD) THEN
               !   PRINT*,"MYNN DRIVER END: k=",k," sh=",sh(k)
               !   PRINT*," sqw=",sqw(k)," thl=",thl(k)," exch_m=",exch_m(i,k)
               !   PRINT*," xland=",xland(i)," rmol=",rmol(i)," ust=",ust(i)
               !   PRINT*," qke=",qke(i,k)," el=",el_pbl(i,k)," tsq=",tsq(i,k)
               !   PRINT*," PBLH=",PBLH(i)," u=",u(i,k)," v=",v(i,k)
               !   PRINT*," vq=",vq(k)," vt=",vt(k)
               !ENDIF
             ENDDO !end-k
          ENDIF
          !***  End debug prints

          !JOE-add tke_pbl for coupling w/shallow-cu schemes (TKE_PBL = QKE/2.)
          !    TKE_PBL is defined on interfaces, while QKE is at middle of layer.
          !tke_pbl(i,kts) = 0.5*MAX(qke(i,kts),1.0e-10)
          !DO k = kts+1,kte
          !   afk = dz1(k)/( dz1(k)+dz1(k-1) )
          !   abk = 1.0 -afk
          !   tke_pbl(i,k) = 0.5*MAX(qke(i,k)*abk+qke(i,k-1)*afk,1.0e-3)
          !ENDDO

    ENDDO !end i-loop

!ACF copy qke into qke_adv if using advection
    IF (bl_mynn_tkeadvect) THEN
       qke_adv=qke
    ENDIF
!ACF-end

#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mynn_bl_driver
!> @}

END MODULE 
