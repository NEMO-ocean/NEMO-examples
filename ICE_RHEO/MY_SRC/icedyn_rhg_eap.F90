MODULE icedyn_rhg_eap
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg_eap  ***
   !!   Sea-Ice dynamics : rheology Elasto-Viscous-Plastic
   !!======================================================================
   !! History :   -   !  2007-03  (M.A. Morales Maqueda, S. Bouillon) Original code
   !!            3.0  !  2008-03  (M. Vancoppenolle) adaptation to new model
   !!             -   !  2008-11  (M. Vancoppenolle, S. Bouillon, Y. Aksenov) add surface tilt in ice rheolohy
   !!            3.3  !  2009-05  (G.Garric)    addition of the evp case
   !!            3.4  !  2011-01  (A. Porter)   dynamical allocation
   !!            3.5  !  2012-08  (R. Benshila) AGRIF
   !!            3.6  !  2016-06  (C. Rousset)  Rewriting + landfast ice + mEVP (Bouillon 2013)
   !!            3.7  !  2017     (C. Rousset)  add aEVP (Kimmritz 2016-2017)
   !!            4.0  !  2018     (many people) SI3 [aka Sea Ice cube]
   !!                 !  2019     (S. Rynders, Y. Aksenov, C. Rousset)  change into eap rheology from
   !!                                           CICE code (Tsamados, Heorton)
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_dyn_rhg_eap : computes ice velocities from EVP rheology
   !!   rhg_eap_rst     : read/write EVP fields in ice restart
   !!----------------------------------------------------------------------
   USE phycst         ! Physical constant
   USE dom_oce        ! Ocean domain
   USE sbc_oce , ONLY : ln_ice_embd, nn_fsbc, ssh_m
   USE sbc_ice , ONLY : utau_ice, vtau_ice, snwice_mass, snwice_mass_b
   USE ice            ! sea-ice: ice variables
   USE icevar         ! ice_var_sshdyn
   USE icedyn_rdgrft  ! sea-ice: ice strength
   USE bdy_oce , ONLY : ln_bdy
   USE bdyice
#if defined key_agrif
   USE agrif_ice_interp
#endif
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE prtctl         ! Print control

   USE netcdf         ! NetCDF library for convergence test
   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_rhg_eap   ! called by icedyn_rhg.F90
   PUBLIC   rhg_eap_rst       ! called by icedyn_rhg.F90

   REAL(wp), PARAMETER ::   pphi = 3.141592653589793_wp/12._wp    ! diamond shaped floe smaller angle (default phi = 30 deg)

   ! look-up table for calculating structure tensor
   INTEGER, PARAMETER ::   nx_yield = 41
   INTEGER, PARAMETER ::   ny_yield = 41
   INTEGER, PARAMETER ::   na_yield = 21

   REAL(wp), DIMENSION(nx_yield, ny_yield, na_yield) ::   s11r, s12r, s22r, s11s, s12s, s22s
   REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   fimask   ! mask at F points for the ice

   !! for convergence tests
   INTEGER ::   ncvgid   ! netcdf file id
   INTEGER ::   nvarid   ! netcdf variable id

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icedyn_rhg_eap.F90 11536 2019-09-11 13:54:18Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_rhg_eap( kt, Kmm, pstress1_i, pstress2_i, pstress12_i, pshear_i, pdivu_i, pdelta_i, paniso_11, paniso_12, prdg_conv )
      !!-------------------------------------------------------------------
      !!                 ***  SUBROUTINE ice_dyn_rhg_eap  ***
      !!                             EAP-C-grid
      !!
      !! ** purpose : determines sea ice drift from wind stress, ice-ocean
      !!  stress and sea-surface slope. Ice-ice interaction is described by
      !!  a non-linear elasto-anisotropic-plastic (EAP) law including shear
      !!  strength and a bulk rheology .
      !!
      !!  The points in the C-grid look like this, dear reader
      !!
      !!                              (ji,jj)
      !!                                 |
      !!                                 |
      !!                      (ji-1,jj)  |  (ji,jj)
      !!                             ---------
      !!                            |         |
      !!                            | (ji,jj) |------(ji,jj)
      !!                            |         |
      !!                             ---------
      !!                     (ji-1,jj-1)     (ji,jj-1)
      !!
      !! ** Inputs  : - wind forcing (stress), oceanic currents
      !!                ice total volume (vt_i) per unit area
      !!                snow total volume (vt_s) per unit area
      !!
      !! ** Action  : - compute u_ice, v_ice : the components of the
      !!                sea-ice velocity vector
      !!              - compute delta_i, shear_i, divu_i, which are inputs
      !!                of the ice thickness distribution
      !!
      !! ** Steps   : 0) compute mask at F point
      !!              1) Compute ice snow mass, ice strength
      !!              2) Compute wind, oceanic stresses, mass terms and
      !!                 coriolis terms of the momentum equation
      !!              3) Solve the momentum equation (iterative procedure)
      !!              4) Recompute delta, shear and divergence
      !!                 (which are inputs of the ITD) & store stress
      !!                 for the next time step
      !!              5) Diagnostics including charge ellipse
      !!
      !! ** Notes   : There is the possibility to use aEVP from the nice work of Kimmritz et al. (2016 & 2017)
      !!              by setting up ln_aEVP=T (i.e. changing alpha and beta parameters).
      !!              This is an upgraded version of mEVP from Bouillon et al. 2013
      !!              (i.e. more stable and better convergence)
      !!
      !! References : Hunke and Dukowicz, JPO97
      !!              Bouillon et al., Ocean Modelling 2009
      !!              Bouillon et al., Ocean Modelling 2013
      !!              Kimmritz et al., Ocean Modelling 2016 & 2017
      !!-------------------------------------------------------------------
      INTEGER                 , INTENT(in   ) ::   kt                                    ! time step
      INTEGER                 , INTENT(in   ) ::   Kmm                                   ! ocean time level index
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pstress1_i, pstress2_i, pstress12_i   !
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   pshear_i  , pdivu_i   , pdelta_i      !
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   paniso_11 , paniso_12                 ! structure tensor components
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   prdg_conv                             ! for ridging
      !!
      INTEGER ::   ji, jj       ! dummy loop indices
      INTEGER ::   jter         ! local integers
      !
      REAL(wp) ::   zrhoco                                              ! rau0 * rn_cio
      REAL(wp) ::   zdtevp, z1_dtevp                                    ! time step for subcycling
      REAL(wp) ::   ecc2, z1_ecc2                                       ! square of yield ellipse eccenticity
      REAL(wp) ::   zalph1, z1_alph1, zalph2, z1_alph2                  ! alpha coef from Bouillon 2009 or Kimmritz 2017
      REAl(wp) ::   zbetau, zbetav
      REAL(wp) ::   zm1, zm2, zm3, zmassU, zmassV, zvU, zvV             ! ice/snow mass and volume
      REAL(wp) ::   zds2, zdt, zdt2, zdiv, zdiv2, zdsT                  ! temporary scalars
      REAL(wp) ::   zTauO, zTauB, zRHS, zvel                            ! temporary scalars
      REAL(wp) ::   zkt                                                 ! isotropic tensile strength for landfast ice
      REAL(wp) ::   zvCr                                                ! critical ice volume above which ice is landfast
      !
      REAL(wp) ::   zintb, zintn                                        ! dummy argument
      REAL(wp) ::   zfac_x, zfac_y
      REAL(wp) ::   zshear, zdum1, zdum2
      REAL(wp) ::   zstressptmp, zstressmtmp, zstress12tmpF             ! anisotropic stress tensor components
      REAL(wp) ::   zalphar, zalphas                                    ! for mechanical redistribution
      REAL(wp) ::   zmresult11, zmresult12, z1dtevpkth, zp5kth, z1_dtevp_A  ! for structure tensor evolution
      REAL(wp) ::   zinvw                                                ! for test case

      !
      REAL(wp), DIMENSION(jpi,jpj) ::   zstress12tmp                    ! anisotropic stress tensor component for regridding
      REAL(wp), DIMENSION(jpi,jpj) ::   zyield11, zyield22, zyield12    ! yield surface tensor for history
      REAL(wp), DIMENSION(jpi,jpj) ::   zdelta, zp_delt                 ! delta and P/delta at T points
      REAL(wp), DIMENSION(jpi,jpj) ::   zten_i                          ! tension
      REAL(wp), DIMENSION(jpi,jpj) ::   zbeta                           ! beta coef from Kimmritz 2017
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   zdt_m                           ! (dt / ice-snow_mass) on T points
      REAL(wp), DIMENSION(jpi,jpj) ::   zaU  , zaV                      ! ice fraction on U/V points
      REAL(wp), DIMENSION(jpi,jpj) ::   zmU_t, zmV_t                    ! (ice-snow_mass / dt) on U/V points
      REAL(wp), DIMENSION(jpi,jpj) ::   zmf                             ! coriolis parameter at T points
      REAL(wp), DIMENSION(jpi,jpj) ::   v_oceU, u_oceV, v_iceU, u_iceV  ! ocean/ice u/v component on V/U points
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   zds                             ! shear
      REAL(wp), DIMENSION(jpi,jpj) ::   zs1, zs2, zs12                  ! stress tensor components
      REAL(wp), DIMENSION(jpi,jpj) ::   zsshdyn                         ! array used for the calculation of ice surface slope:
      !                                                                 !    ocean surface (ssh_m) if ice is not embedded
      !                                                                 !    ice bottom surface if ice is embedded
      REAL(wp), DIMENSION(jpi,jpj) ::   zfU  , zfV                      ! internal stresses
      REAL(wp), DIMENSION(jpi,jpj) ::   zspgU, zspgV                    ! surface pressure gradient at U/V points
      REAL(wp), DIMENSION(jpi,jpj) ::   zCorU, zCorV                    ! Coriolis stress array
      REAL(wp), DIMENSION(jpi,jpj) ::   ztaux_ai, ztauy_ai              ! ice-atm. stress at U-V points
      REAL(wp), DIMENSION(jpi,jpj) ::   ztaux_oi, ztauy_oi              ! ice-ocean stress at U-V points
      REAL(wp), DIMENSION(jpi,jpj) ::   ztaux_bi, ztauy_bi              ! ice-OceanBottom stress at U-V points (landfast)
      REAL(wp), DIMENSION(jpi,jpj) ::   ztaux_base, ztauy_base          ! ice-bottom stress at U-V points (landfast)
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   zmsk00, zmsk15
      REAL(wp), DIMENSION(jpi,jpj) ::   zmsk01x, zmsk01y                ! dummy arrays
      REAL(wp), DIMENSION(jpi,jpj) ::   zmsk00x, zmsk00y                ! mask for ice presence

      REAL(wp), PARAMETER          ::   zepsi  = 1.0e-20_wp             ! tolerance parameter
      REAL(wp), PARAMETER          ::   zmmin  = 1._wp                  ! ice mass (kg/m2)  below which ice velocity becomes very small
      REAL(wp), PARAMETER          ::   zamin  = 0.001_wp               ! ice concentration below which ice velocity becomes very small
      !! --- check convergence
      REAL(wp), DIMENSION(jpi,jpj) ::   zu_ice, zv_ice
      !! --- diags
      REAL(wp) ::   zsig1, zsig2, zsig12, zfac, z1_strength
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zsig_I, zsig_II, zsig1_p, zsig2_p
      !! --- SIMIP diags
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zdiag_xmtrp_ice ! X-component of ice mass transport (kg/s)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zdiag_ymtrp_ice ! Y-component of ice mass transport (kg/s)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zdiag_xmtrp_snw ! X-component of snow mass transport (kg/s)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zdiag_ymtrp_snw ! Y-component of snow mass transport (kg/s)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zdiag_xatrp     ! X-component of area transport (m2/s)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zdiag_yatrp     ! Y-component of area transport (m2/s)
      !!-------------------------------------------------------------------

      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- ice_dyn_rhg_eap: EAP sea-ice rheology'
      !
      ! for diagnostics and convergence tests
      DO_2D( 1, 1, 1, 1 )
         zmsk00(ji,jj) = MAX( 0._wp , SIGN( 1._wp , at_i(ji,jj) - epsi06  ) ) ! 1 if ice    , 0 if no ice
      END_2D
      IF( nn_rhg_chkcvg > 0 ) THEN
         DO_2D( 1, 1, 1, 1 )
            zmsk15(ji,jj) = MAX( 0._wp , SIGN( 1._wp , at_i(ji,jj) - 0.15_wp ) ) ! 1 if 15% ice, 0 if less
         END_2D
      ENDIF
      !
      !------------------------------------------------------------------------------!
      ! 0) mask at F points for the ice
      !------------------------------------------------------------------------------!
      IF( kt == nit000 ) THEN
         ! ocean/land mask
         ALLOCATE( fimask(jpi,jpj) )
         IF( rn_ishlat == 0._wp ) THEN
            DO_2D( 0, 0, 0, 0 )
               fimask(ji,jj) = tmask(ji,jj,1) * tmask(ji+1,jj,1) * tmask(ji,jj+1,1) * tmask(ji+1,jj+1,1)
            END_2D
         ELSE
            DO_2D( 0, 0, 0, 0 )
               fimask(ji,jj) = tmask(ji,jj,1) * tmask(ji+1,jj,1) * tmask(ji,jj+1,1) * tmask(ji+1,jj+1,1)
               ! Lateral boundary conditions on velocity (modify fimask)
               IF( fimask(ji,jj) == 0._wp ) THEN
                  fimask(ji,jj) = rn_ishlat * MIN( 1._wp , MAX( umask(ji,jj,1), umask(ji,jj+1,1), &
                     &                                          vmask(ji,jj,1), vmask(ji+1,jj,1) ) )
               ENDIF
            END_2D
         ENDIF
         CALL lbc_lnk( 'icedyn_rhg_eap', fimask, 'F', 1.0_wp )
      ENDIF

      !------------------------------------------------------------------------------!
      ! 1) define some variables and initialize arrays
      !------------------------------------------------------------------------------!
      zrhoco = rho0 * rn_cio
!extra code for test case boundary conditions
      zinvw=1._wp/(zrhoco*0.5_wp)

      ! ecc2: square of yield ellipse eccenticrity
      ecc2    = rn_ecc * rn_ecc
      z1_ecc2 = 1._wp / ecc2

      ! alpha parameters (Bouillon 2009)
      IF( .NOT. ln_aEVP ) THEN
         zdtevp   = rDt_ice / REAL( nn_nevp )
         zalph1 =   2._wp * rn_relast * REAL( nn_nevp )
         zalph2 = zalph1 * z1_ecc2

         z1_alph1 = 1._wp / ( zalph1 + 1._wp )
         z1_alph2 = 1._wp / ( zalph2 + 1._wp )
      ELSE
         zdtevp   = rdt_ice
         ! zalpha parameters set later on adaptatively
      ENDIF
      z1_dtevp = 1._wp / zdtevp

      ! Initialise stress tensor
      zs1 (:,:) = pstress1_i (:,:)
      zs2 (:,:) = pstress2_i (:,:)
      zs12(:,:) = pstress12_i(:,:)

      ! constants for structure tensor
      z1_dtevp_A = z1_dtevp/10.0_wp
      z1dtevpkth = 1._wp / (z1_dtevp_A + 0.00002_wp)
      zp5kth = 0.5_wp * 0.00002_wp

      ! Ice strength
      CALL ice_strength

      ! landfast param from Lemieux(2016): add isotropic tensile strength (following Konig Beatty and Holland, 2010)
      IF( ln_landfast_L16 ) THEN   ;   zkt = rn_lf_tensile
      ELSE                         ;   zkt = 0._wp
      ENDIF
      !
      !------------------------------------------------------------------------------!
      ! 2) Wind / ocean stress, mass terms, coriolis terms
      !------------------------------------------------------------------------------!
      ! sea surface height
      !    embedded sea ice: compute representative ice top surface
      !    non-embedded sea ice: use ocean surface for slope calculation
      zsshdyn(:,:) = ice_var_sshdyn( ssh_m, snwice_mass, snwice_mass_b)

      DO_2D( 0, 0, 0, 0 )

         ! ice fraction at U-V points
         zaU(ji,jj) = 0.5_wp * ( at_i(ji,jj) * e1e2t(ji,jj) + at_i(ji+1,jj) * e1e2t(ji+1,jj) ) * r1_e1e2u(ji,jj) * umask(ji,jj,1)
         zaV(ji,jj) = 0.5_wp * ( at_i(ji,jj) * e1e2t(ji,jj) + at_i(ji,jj+1) * e1e2t(ji,jj+1) ) * r1_e1e2v(ji,jj) * vmask(ji,jj,1)

         ! Ice/snow mass at U-V points
         zm1 = ( rhos * vt_s(ji  ,jj  ) + rhoi * vt_i(ji  ,jj  ) )
         zm2 = ( rhos * vt_s(ji+1,jj  ) + rhoi * vt_i(ji+1,jj  ) )
         zm3 = ( rhos * vt_s(ji  ,jj+1) + rhoi * vt_i(ji  ,jj+1) )
         zmassU = 0.5_wp * ( zm1 * e1e2t(ji,jj) + zm2 * e1e2t(ji+1,jj) ) * r1_e1e2u(ji,jj) * umask(ji,jj,1)
         zmassV = 0.5_wp * ( zm1 * e1e2t(ji,jj) + zm3 * e1e2t(ji,jj+1) ) * r1_e1e2v(ji,jj) * vmask(ji,jj,1)

         ! Ocean currents at U-V points
         v_oceU(ji,jj)   = 0.25_wp * ( v_oce(ji,jj) + v_oce(ji,jj-1) + v_oce(ji+1,jj) + v_oce(ji+1,jj-1) ) * umask(ji,jj,1)
         u_oceV(ji,jj)   = 0.25_wp * ( u_oce(ji,jj) + u_oce(ji-1,jj) + u_oce(ji,jj+1) + u_oce(ji-1,jj+1) ) * vmask(ji,jj,1)

         ! Coriolis at T points (m*f)
         zmf(ji,jj)      = zm1 * ff_t(ji,jj)

         ! dt/m at T points (for alpha and beta coefficients)
         zdt_m(ji,jj)    = zdtevp / MAX( zm1, zmmin )

         ! m/dt
         zmU_t(ji,jj)    = zmassU * z1_dtevp
         zmV_t(ji,jj)    = zmassV * z1_dtevp

         ! Drag ice-atm.
         ztaux_ai(ji,jj) = zaU(ji,jj) * utau_ice(ji,jj)
         ztauy_ai(ji,jj) = zaV(ji,jj) * vtau_ice(ji,jj)

         ! Surface pressure gradient (- m*g*GRAD(ssh)) at U-V points
         zspgU(ji,jj)    = - zmassU * grav * ( zsshdyn(ji+1,jj) - zsshdyn(ji,jj) ) * r1_e1u(ji,jj)
         zspgV(ji,jj)    = - zmassV * grav * ( zsshdyn(ji,jj+1) - zsshdyn(ji,jj) ) * r1_e2v(ji,jj)

         ! masks
         zmsk00x(ji,jj) = 1._wp - MAX( 0._wp, SIGN( 1._wp, -zmassU ) )  ! 0 if no ice
         zmsk00y(ji,jj) = 1._wp - MAX( 0._wp, SIGN( 1._wp, -zmassV ) )  ! 0 if no ice

         ! switches
         IF( zmassU <= zmmin .AND. zaU(ji,jj) <= zamin ) THEN   ;   zmsk01x(ji,jj) = 0._wp
         ELSE                                                   ;   zmsk01x(ji,jj) = 1._wp   ;   ENDIF
         IF( zmassV <= zmmin .AND. zaV(ji,jj) <= zamin ) THEN   ;   zmsk01y(ji,jj) = 0._wp
         ELSE                                                   ;   zmsk01y(ji,jj) = 1._wp   ;   ENDIF

      END_2D
      CALL lbc_lnk( 'icedyn_rhg_eap', zmf, 'T', 1.0_wp, zdt_m, 'T', 1.0_wp )
      !
      !                                  !== Landfast ice parameterization ==!
      !
      IF( ln_landfast_L16 ) THEN         !-- Lemieux 2016
         DO_2D( 0, 0, 0, 0 )
            ! ice thickness at U-V points
            zvU = 0.5_wp * ( vt_i(ji,jj) * e1e2t(ji,jj) + vt_i(ji+1,jj) * e1e2t(ji+1,jj) ) * r1_e1e2u(ji,jj) * umask(ji,jj,1)
            zvV = 0.5_wp * ( vt_i(ji,jj) * e1e2t(ji,jj) + vt_i(ji,jj+1) * e1e2t(ji,jj+1) ) * r1_e1e2v(ji,jj) * vmask(ji,jj,1)
            ! ice-bottom stress at U points
            zvCr = zaU(ji,jj) * rn_lf_depfra * hu(ji,jj,Kmm) * ( 1._wp - icb_mask(ji,jj) ) ! if grounded icebergs are read: ocean depth = 0
            ztaux_base(ji,jj) = - rn_lf_bfr * MAX( 0._wp, zvU - zvCr ) * EXP( -rn_crhg * ( 1._wp - zaU(ji,jj) ) )
            ! ice-bottom stress at V points
            zvCr = zaV(ji,jj) * rn_lf_depfra * hv(ji,jj,Kmm) * ( 1._wp - icb_mask(ji,jj) ) ! if grounded icebergs are read: ocean depth = 0
            ztauy_base(ji,jj) = - rn_lf_bfr * MAX( 0._wp, zvV - zvCr ) * EXP( -rn_crhg * ( 1._wp - zaV(ji,jj) ) )
            ! ice_bottom stress at T points
            zvCr = at_i(ji,jj) * rn_lf_depfra * ht(ji,jj) * ( 1._wp - icb_mask(ji,jj) )    ! if grounded icebergs are read: ocean depth = 0
            tau_icebfr(ji,jj) = - rn_lf_bfr * MAX( 0._wp, vt_i(ji,jj) - zvCr ) * EXP( -rn_crhg * ( 1._wp - at_i(ji,jj) ) )
         END_2D
         CALL lbc_lnk( 'icedyn_rhg_eap', tau_icebfr(:,:), 'T', 1.0_wp )
         !
      ELSE                               !-- no landfast
         DO_2D( 0, 0, 0, 0 )
            ztaux_base(ji,jj) = 0._wp
            ztauy_base(ji,jj) = 0._wp
         END_2D
      ENDIF

      !------------------------------------------------------------------------------!
      ! 3) Solution of the momentum equation, iterative procedure
      !------------------------------------------------------------------------------!
      !
      !                                               ! ==================== !
      DO jter = 1 , nn_nevp                           !    loop over jter    !
         !                                            ! ==================== !
         l_full_nf_update = jter == nn_nevp   ! false: disable full North fold update (performances) for iter = 1 to nn_nevp-1
         !
         ! convergence test
         IF( nn_rhg_chkcvg == 1 .OR. nn_rhg_chkcvg == 2  ) THEN
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               zu_ice(ji,jj) = u_ice(ji,jj) * umask(ji,jj,1) ! velocity at previous time step
               zv_ice(ji,jj) = v_ice(ji,jj) * vmask(ji,jj,1)
            END_2D
         ENDIF

         ! --- divergence, tension & shear (Appendix B of Hunke & Dukowicz, 2002) --- !
         DO_2D( 1, 0, 1, 0 )

            ! shear at F points
            zds(ji,jj) = ( ( u_ice(ji,jj+1) * r1_e1u(ji,jj+1) - u_ice(ji,jj) * r1_e1u(ji,jj) ) * e1f(ji,jj) * e1f(ji,jj)   &
               &         + ( v_ice(ji+1,jj) * r1_e2v(ji+1,jj) - v_ice(ji,jj) * r1_e2v(ji,jj) ) * e2f(ji,jj) * e2f(ji,jj)   &
               &         ) * r1_e1e2f(ji,jj) * fimask(ji,jj)

         END_2D

         DO_2D( 0, 0, 0, 0 )

            ! shear**2 at T points (doc eq. A16)
            zds2 = ( zds(ji,jj  ) * zds(ji,jj  ) * e1e2f(ji,jj  ) + zds(ji-1,jj  ) * zds(ji-1,jj  ) * e1e2f(ji-1,jj  )  &
               &   + zds(ji,jj-1) * zds(ji,jj-1) * e1e2f(ji,jj-1) + zds(ji-1,jj-1) * zds(ji-1,jj-1) * e1e2f(ji-1,jj-1)  &
               &   ) * 0.25_wp * r1_e1e2t(ji,jj)

            ! divergence at T points
            zdiv  = ( e2u(ji,jj) * u_ice(ji,jj) - e2u(ji-1,jj) * u_ice(ji-1,jj)   &
               &    + e1v(ji,jj) * v_ice(ji,jj) - e1v(ji,jj-1) * v_ice(ji,jj-1)   &
               &    ) * r1_e1e2t(ji,jj)
            zdiv2 = zdiv * zdiv

            ! tension at T points
            zdt  = ( ( u_ice(ji,jj) * r1_e2u(ji,jj) - u_ice(ji-1,jj) * r1_e2u(ji-1,jj) ) * e2t(ji,jj) * e2t(ji,jj)   &
               &   - ( v_ice(ji,jj) * r1_e1v(ji,jj) - v_ice(ji,jj-1) * r1_e1v(ji,jj-1) ) * e1t(ji,jj) * e1t(ji,jj)   &
               &   ) * r1_e1e2t(ji,jj)
            zdt2 = zdt * zdt

            ! delta at T points
            zdelta(ji,jj) = SQRT( zdiv2 + ( zdt2 + zds2 ) * z1_ecc2 )

         END_2D
         CALL lbc_lnk( 'icedyn_rhg_eap', zdelta, 'T', 1.0_wp )

         ! P/delta at T points
         DO_2D( 1, 1, 1, 1 )
            zp_delt(ji,jj) = strength(ji,jj) / ( zdelta(ji,jj) + rn_creepl )
         END_2D

         DO_2D( 0, 1, 0, 1 )   ! loop ends at jpi,jpj so that no lbc_lnk are needed for zs1 and zs2

             ! shear at T points
            zdsT = ( zds(ji,jj  ) * e1e2f(ji,jj  ) + zds(ji-1,jj  ) * e1e2f(ji-1,jj  )  &
               &   + zds(ji,jj-1) * e1e2f(ji,jj-1) + zds(ji-1,jj-1) * e1e2f(ji-1,jj-1)  &
               &   ) * 0.25_wp * r1_e1e2t(ji,jj)

           ! divergence at T points (duplication to avoid communications)
            zdiv  = ( e2u(ji,jj) * u_ice(ji,jj) - e2u(ji-1,jj) * u_ice(ji-1,jj)   &
               &    + e1v(ji,jj) * v_ice(ji,jj) - e1v(ji,jj-1) * v_ice(ji,jj-1)   &
               &    ) * r1_e1e2t(ji,jj)

            ! tension at T points (duplication to avoid communications)
            zdt  = ( ( u_ice(ji,jj) * r1_e2u(ji,jj) - u_ice(ji-1,jj) * r1_e2u(ji-1,jj) ) * e2t(ji,jj) * e2t(ji,jj)   &
               &   - ( v_ice(ji,jj) * r1_e1v(ji,jj) - v_ice(ji,jj-1) * r1_e1v(ji,jj-1) ) * e1t(ji,jj) * e1t(ji,jj)   &
               &   ) * r1_e1e2t(ji,jj)

            ! --- anisotropic stress calculation --- !
            CALL update_stress_rdg (jter, nn_nevp, zdiv, zdt, zdsT, paniso_11(ji,jj), paniso_12(ji,jj), &
                                    zstressptmp, zstressmtmp, zstress12tmp(ji,jj), strength(ji,jj), zalphar, zalphas)

            ! structure tensor update
               CALL calc_ffrac(zstressptmp, zstressmtmp, zstress12tmp(ji,jj), paniso_11(ji,jj), paniso_12(ji,jj), zmresult11,  zmresult12)

               paniso_11(ji,jj) = (paniso_11(ji,jj)  + 0.5*2.e-5*zdtevp + zmresult11*zdtevp) / (1. + 2.e-5*zdtevp) ! implicit
               paniso_12(ji,jj) = (paniso_12(ji,jj)                     + zmresult12*zdtevp) / (1. + 2.e-5*zdtevp) ! implicit

            IF (jter == nn_nevp) THEN
               zyield11(ji,jj) = 0.5_wp * (zstressptmp + zstressmtmp)
               zyield22(ji,jj) = 0.5_wp * (zstressptmp - zstressmtmp)
               zyield12(ji,jj) = zstress12tmp(ji,jj)
               prdg_conv(ji,jj) = -min( zalphar, 0._wp)
            ENDIF

            ! alpha for aEVP
            !   gamma = 0.5*P/(delta+creepl) * (c*pi)**2/Area * dt/m
            !   alpha = beta = sqrt(4*gamma)
            IF( ln_aEVP ) THEN
               zalph1   = MAX( 50._wp, rpi * SQRT( 0.5_wp * zp_delt(ji,jj) * r1_e1e2t(ji,jj) * zdt_m(ji,jj) ) )
               z1_alph1 = 1._wp / ( zalph1 + 1._wp )
               zalph2   = zalph1
               z1_alph2 = z1_alph1
               ! explicit:
               ! z1_alph1 = 1._wp / zalph1
               ! z1_alph2 = 1._wp / zalph1
               ! zalph1 = zalph1 - 1._wp
               ! zalph2 = zalph1
            ENDIF

            ! stress at T points (zkt/=0 if landfast)
            zs1(ji,jj) = ( zs1(ji,jj) * zalph1 + zstressptmp ) * z1_alph1
            zs2(ji,jj) = ( zs2(ji,jj) * zalph1 + zstressmtmp ) * z1_alph1
         END_2D
         CALL lbc_lnk( 'icedyn_rhg_eap', zstress12tmp, 'T', 1.0_wp , paniso_11, 'T', 1.0_wp , paniso_12, 'T', 1.0_wp)

        ! Save beta at T-points for further computations
         IF( ln_aEVP ) THEN
            DO_2D( 1, 1, 1, 1 )
               zbeta(ji,jj) = MAX( 50._wp, rpi * SQRT( 0.5_wp * zp_delt(ji,jj) * r1_e1e2t(ji,jj) * zdt_m(ji,jj) ) )
            END_2D
         ENDIF

         DO_2D( 1, 0, 1, 0 )
            ! stress12tmp at F points
            zstress12tmpF = ( zstress12tmp(ji,jj+1) * e1e2t(ji,jj+1) + zstress12tmp(ji+1,jj+1) * e1e2t(ji+1,jj+1)  &
               &            + zstress12tmp(ji,jj  ) * e1e2t(ji,jj  ) + zstress12tmp(ji+1,jj  ) * e1e2t(ji+1,jj  )  &
               &            ) * 0.25_wp * r1_e1e2f(ji,jj)

            ! alpha for aEVP
            IF( ln_aEVP ) THEN
               zalph2   = MAX( zbeta(ji,jj), zbeta(ji+1,jj), zbeta(ji,jj+1), zbeta(ji+1,jj+1) )
               z1_alph2 = 1._wp / ( zalph2 + 1._wp )
               ! explicit:
               ! z1_alph2 = 1._wp / zalph2
               ! zalph2 = zalph2 - 1._wp
            ENDIF

            ! stress at F points (zkt/=0 if landfast)
            zs12(ji,jj) = ( zs12(ji,jj) * zalph1 + zstress12tmpF ) * z1_alph1

         END_2D
         CALL lbc_lnk( 'icedyn_rhg_eap', zs1, 'T', 1.0_wp, zs2, 'T', 1.0_wp, zs12, 'F', 1.0_wp )

         ! --- Ice internal stresses (Appendix C of Hunke and Dukowicz, 2002) --- !
         DO_2D( 0, 0, 0, 0 )
            !                   !--- U points
            zfU(ji,jj) = 0.5_wp * ( ( zs1(ji+1,jj) - zs1(ji,jj) ) * e2u(ji,jj)                                             &
               &                  + ( zs2(ji+1,jj) * e2t(ji+1,jj) * e2t(ji+1,jj) - zs2(ji,jj) * e2t(ji,jj) * e2t(ji,jj)    &
               &                    ) * r1_e2u(ji,jj)                                                                      &
               &                  + ( zs12(ji,jj) * e1f(ji,jj) * e1f(ji,jj) - zs12(ji,jj-1) * e1f(ji,jj-1) * e1f(ji,jj-1)  &
               &                    ) * 2._wp * r1_e1u(ji,jj)                                                              &
               &                  ) * r1_e1e2u(ji,jj)
            !
            !                !--- V points
            zfV(ji,jj) = 0.5_wp * ( ( zs1(ji,jj+1) - zs1(ji,jj) ) * e1v(ji,jj)                                             &
               &                  - ( zs2(ji,jj+1) * e1t(ji,jj+1) * e1t(ji,jj+1) - zs2(ji,jj) * e1t(ji,jj) * e1t(ji,jj)    &
               &                    ) * r1_e1v(ji,jj)                                                                      &
               &                  + ( zs12(ji,jj) * e2f(ji,jj) * e2f(ji,jj) - zs12(ji-1,jj) * e2f(ji-1,jj) * e2f(ji-1,jj)  &
               &                    ) * 2._wp * r1_e2v(ji,jj)                                                              &
               &                  ) * r1_e1e2v(ji,jj)
            !
            !                !--- ice currents at U-V point
            v_iceU(ji,jj) = 0.25_wp * ( v_ice(ji,jj) + v_ice(ji,jj-1) + v_ice(ji+1,jj) + v_ice(ji+1,jj-1) ) * umask(ji,jj,1)
            u_iceV(ji,jj) = 0.25_wp * ( u_ice(ji,jj) + u_ice(ji-1,jj) + u_ice(ji,jj+1) + u_ice(ji-1,jj+1) ) * vmask(ji,jj,1)
            !
         END_2D
         !
         ! --- Computation of ice velocity --- !
         !  Bouillon et al. 2013 (eq 47-48) => unstable unless alpha, beta vary as in Kimmritz 2016 & 2017
         !  Bouillon et al. 2009 (eq 34-35) => stable
         IF( MOD(jter,2) == 0 ) THEN ! even iterations
            !
            DO_2D( 0, 0, 0, 0 )
               !                 !--- tau_io/(v_oce - v_ice)
               zTauO = zaV(ji,jj) * zrhoco * SQRT( ( v_ice (ji,jj) - v_oce (ji,jj) ) * ( v_ice (ji,jj) - v_oce (ji,jj) )  &
                  &                              + ( u_iceV(ji,jj) - u_oceV(ji,jj) ) * ( u_iceV(ji,jj) - u_oceV(ji,jj) ) )
               !                 !--- Ocean-to-Ice stress
               ztauy_oi(ji,jj) = zTauO * ( v_oce(ji,jj) - v_ice(ji,jj) )
               !
               !                 !--- tau_bottom/v_ice
               zvel  = 5.e-05_wp + SQRT( v_ice(ji,jj) * v_ice(ji,jj) + u_iceV(ji,jj) * u_iceV(ji,jj) )
               zTauB = ztauy_base(ji,jj) / zvel
               !                 !--- OceanBottom-to-Ice stress
               ztauy_bi(ji,jj) = zTauB * v_ice(ji,jj)
               !
               !                 !--- Coriolis at V-points (energy conserving formulation)
               zCorV(ji,jj)  = - 0.25_wp * r1_e2v(ji,jj) *  &
                  &    ( zmf(ji,jj  ) * ( e2u(ji,jj  ) * u_ice(ji,jj  ) + e2u(ji-1,jj  ) * u_ice(ji-1,jj  ) )  &
                  &    + zmf(ji,jj+1) * ( e2u(ji,jj+1) * u_ice(ji,jj+1) + e2u(ji-1,jj+1) * u_ice(ji-1,jj+1) ) )
               !
               !                 !--- Sum of external forces (explicit solution) = F + tau_ia + Coriolis + spg + tau_io
               zRHS = zfV(ji,jj) + ztauy_ai(ji,jj) + zCorV(ji,jj) + zspgV(ji,jj) + ztauy_oi(ji,jj)
               !
               !                 !--- landfast switch => 0 = static  friction : TauB > RHS & sign(TauB) /= sign(RHS)
               !                                         1 = sliding friction : TauB < RHS
               rswitch = 1._wp - MIN( 1._wp, ABS( SIGN( 1._wp, zRHS + ztauy_base(ji,jj) ) - SIGN( 1._wp, zRHS ) ) )
               !
               IF( ln_aEVP ) THEN !--- ice velocity using aEVP (Kimmritz et al 2016 & 2017)
                  zbetav = MAX( zbeta(ji,jj), zbeta(ji,jj+1) )
                  v_ice(ji,jj) = ( (          rswitch   * ( zmV_t(ji,jj) * ( zbetav * v_ice(ji,jj) + v_ice_b(ji,jj) )         & ! previous velocity
                     &                                    + zRHS + zTauO * v_ice(ji,jj)                                       & ! F + tau_ia + Coriolis + spg + tau_io(only ocean part)
                     &                                    ) / MAX( zepsi, zmV_t(ji,jj) * ( zbetav + 1._wp ) + zTauO - zTauB ) & ! m/dt + tau_io(only ice part) + landfast
                     &            + ( 1._wp - rswitch ) * (  v_ice_b(ji,jj)                                                   &
                     &                                     + v_ice  (ji,jj) * MAX( 0._wp, zbetav - zdtevp * rn_lf_relax )     & ! static friction => slow decrease to v=0
                     &                                    ) / ( zbetav + 1._wp )                                              &
                     &             ) * zmsk01y(ji,jj) + v_oce(ji,jj) * 0.01_wp * ( 1._wp - zmsk01y(ji,jj) )                   & ! v_ice = v_oce/100 if mass < zmmin & conc < zamin
                     &           )   * zmsk00y(ji,jj)
               ELSE               !--- ice velocity using EVP implicit formulation (cf Madec doc & Bouillon 2009)
                  v_ice(ji,jj) = ( (         rswitch   * ( zmV_t(ji,jj)  * v_ice(ji,jj)                                       & ! previous velocity
                     &                                    + zRHS + zTauO * v_ice(ji,jj)                                       & ! F + tau_ia + Coriolis + spg + tau_io(only ocean part)
                     &                                    ) / MAX( zepsi, zmV_t(ji,jj) + zTauO - zTauB )                      & ! m/dt + tau_io(only ice part) + landfast
                     &            + ( 1._wp - rswitch ) *   v_ice(ji,jj) * MAX( 0._wp, 1._wp - zdtevp * rn_lf_relax )         & ! static friction => slow decrease to v=0
                     &              ) * zmsk01y(ji,jj) + v_oce(ji,jj) * 0.01_wp * ( 1._wp - zmsk01y(ji,jj) )                  & ! v_ice = v_oce/100 if mass < zmmin & conc < zamin
                     &            )   * zmsk00y(ji,jj)
               ENDIF
!extra code for test case boundary conditions
               IF (mjg(jj)<25 .or. mjg(jj)>975 .or. mig(ji)<25 .or. mig(ji)>975) THEN
                  v_ice(ji,jj) = zinvw*(ztauy_ai(ji,jj) + zCorV(ji,jj) + zspgV(ji,jj) + ztauy_oi(ji,jj))
               END IF

            END_2D
            CALL lbc_lnk( 'icedyn_rhg_eap', v_ice, 'V', -1.0_wp )
            !
#if defined key_agrif
!!          CALL agrif_interp_ice( 'V', jter, nn_nevp )
            CALL agrif_interp_ice( 'V' )
#endif
            IF( ln_bdy )   CALL bdy_ice_dyn( 'V' )
            !
            DO_2D( 0, 0, 0, 0 )
               !                 !--- tau_io/(u_oce - u_ice)
               zTauO = zaU(ji,jj) * zrhoco * SQRT( ( u_ice (ji,jj) - u_oce (ji,jj) ) * ( u_ice (ji,jj) - u_oce (ji,jj) )  &
                  &                              + ( v_iceU(ji,jj) - v_oceU(ji,jj) ) * ( v_iceU(ji,jj) - v_oceU(ji,jj) ) )
               !                 !--- Ocean-to-Ice stress
               ztaux_oi(ji,jj) = zTauO * ( u_oce(ji,jj) - u_ice(ji,jj) )
               !
               !                 !--- tau_bottom/u_ice
               zvel  = 5.e-05_wp + SQRT( v_iceU(ji,jj) * v_iceU(ji,jj) + u_ice(ji,jj) * u_ice(ji,jj) )
               zTauB = ztaux_base(ji,jj) / zvel
               !                 !--- OceanBottom-to-Ice stress
               ztaux_bi(ji,jj) = zTauB * u_ice(ji,jj)
               !
               !                 !--- Coriolis at U-points (energy conserving formulation)
               zCorU(ji,jj)  =   0.25_wp * r1_e1u(ji,jj) *  &
                  &    ( zmf(ji  ,jj) * ( e1v(ji  ,jj) * v_ice(ji  ,jj) + e1v(ji  ,jj-1) * v_ice(ji  ,jj-1) )  &
                  &    + zmf(ji+1,jj) * ( e1v(ji+1,jj) * v_ice(ji+1,jj) + e1v(ji+1,jj-1) * v_ice(ji+1,jj-1) ) )
               !
               !                 !--- Sum of external forces (explicit solution) = F + tau_ia + Coriolis + spg + tau_io
               zRHS = zfU(ji,jj) + ztaux_ai(ji,jj) + zCorU(ji,jj) + zspgU(ji,jj) + ztaux_oi(ji,jj)
               !
               !                 !--- landfast switch => 0 = static  friction : TauB > RHS & sign(TauB) /= sign(RHS)
               !                                         1 = sliding friction : TauB < RHS
               rswitch = 1._wp - MIN( 1._wp, ABS( SIGN( 1._wp, zRHS + ztaux_base(ji,jj) ) - SIGN( 1._wp, zRHS ) ) )
               !
               IF( ln_aEVP ) THEN !--- ice velocity using aEVP (Kimmritz et al 2016 & 2017)
                  zbetau = MAX( zbeta(ji,jj), zbeta(ji+1,jj) )
                  u_ice(ji,jj) = ( (          rswitch   * ( zmU_t(ji,jj) * ( zbetau * u_ice(ji,jj) + u_ice_b(ji,jj) )         & ! previous velocity
                     &                                    + zRHS + zTauO * u_ice(ji,jj)                                       & ! F + tau_ia + Coriolis + spg + tau_io(only ocean part)
                     &                                    ) / MAX( zepsi, zmU_t(ji,jj) * ( zbetau + 1._wp ) + zTauO - zTauB ) & ! m/dt + tau_io(only ice part) + landfast
                     &            + ( 1._wp - rswitch ) * (  u_ice_b(ji,jj)                                                   &
                     &                                     + u_ice  (ji,jj) * MAX( 0._wp, zbetau - zdtevp * rn_lf_relax )     & ! static friction => slow decrease to v=0
                     &                                    ) / ( zbetau + 1._wp )                                              &
                     &             ) * zmsk01x(ji,jj) + u_oce(ji,jj) * 0.01_wp * ( 1._wp - zmsk01x(ji,jj) )                   & ! v_ice = v_oce/100 if mass < zmmin & conc < zamin
                     &           )   * zmsk00x(ji,jj)
               ELSE               !--- ice velocity using EVP implicit formulation (cf Madec doc & Bouillon 2009)
                  u_ice(ji,jj) = ( (         rswitch   * ( zmU_t(ji,jj)  * u_ice(ji,jj)                                       & ! previous velocity
                     &                                    + zRHS + zTauO * u_ice(ji,jj)                                       & ! F + tau_ia + Coriolis + spg + tau_io(only ocean part)
                     &                                    ) / MAX( zepsi, zmU_t(ji,jj) + zTauO - zTauB )                      & ! m/dt + tau_io(only ice part) + landfast
                     &            + ( 1._wp - rswitch ) *   u_ice(ji,jj) * MAX( 0._wp, 1._wp - zdtevp * rn_lf_relax )         & ! static friction => slow decrease to v=0
                     &              ) * zmsk01x(ji,jj) + u_oce(ji,jj) * 0.01_wp * ( 1._wp - zmsk01x(ji,jj) )                  & ! v_ice = v_oce/100 if mass < zmmin & conc < zamin
                     &            )   * zmsk00x(ji,jj)
               ENDIF
!extra code for test case boundary conditions
               IF (mjg(jj)<25 .or. mjg(jj)>975 .or. mig(ji)<25 .or. mig(ji)>975) THEN
                  u_ice(ji,jj) = zinvw*(ztaux_ai(ji,jj) + zCorU(ji,jj) + zspgU(ji,jj) + ztaux_oi(ji,jj))
               END IF

            END_2D
            CALL lbc_lnk( 'icedyn_rhg_eap', u_ice, 'U', -1.0_wp )
            !
#if defined key_agrif
!!          CALL agrif_interp_ice( 'U', jter, nn_nevp )
            CALL agrif_interp_ice( 'U' )
#endif
            IF( ln_bdy )   CALL bdy_ice_dyn( 'U' )
            !
         ELSE ! odd iterations
            !
            DO_2D( 0, 0, 0, 0 )
               !                 !--- tau_io/(u_oce - u_ice)
               zTauO = zaU(ji,jj) * zrhoco * SQRT( ( u_ice (ji,jj) - u_oce (ji,jj) ) * ( u_ice (ji,jj) - u_oce (ji,jj) )  &
                  &                              + ( v_iceU(ji,jj) - v_oceU(ji,jj) ) * ( v_iceU(ji,jj) - v_oceU(ji,jj) ) )
               !                 !--- Ocean-to-Ice stress
               ztaux_oi(ji,jj) = zTauO * ( u_oce(ji,jj) - u_ice(ji,jj) )
               !
               !                 !--- tau_bottom/u_ice
               zvel  = 5.e-05_wp + SQRT( v_iceU(ji,jj) * v_iceU(ji,jj) + u_ice(ji,jj) * u_ice(ji,jj) )
               zTauB = ztaux_base(ji,jj) / zvel
               !                 !--- OceanBottom-to-Ice stress
               ztaux_bi(ji,jj) = zTauB * u_ice(ji,jj)
               !
               !                 !--- Coriolis at U-points (energy conserving formulation)
               zCorU(ji,jj)  =   0.25_wp * r1_e1u(ji,jj) *  &
                  &    ( zmf(ji  ,jj) * ( e1v(ji  ,jj) * v_ice(ji  ,jj) + e1v(ji  ,jj-1) * v_ice(ji  ,jj-1) )  &
                  &    + zmf(ji+1,jj) * ( e1v(ji+1,jj) * v_ice(ji+1,jj) + e1v(ji+1,jj-1) * v_ice(ji+1,jj-1) ) )
               !
               !                 !--- Sum of external forces (explicit solution) = F + tau_ia + Coriolis + spg + tau_io
               zRHS = zfU(ji,jj) + ztaux_ai(ji,jj) + zCorU(ji,jj) + zspgU(ji,jj) + ztaux_oi(ji,jj)
               !
               !                 !--- landfast switch => 0 = static  friction : TauB > RHS & sign(TauB) /= sign(RHS)
               !                                         1 = sliding friction : TauB < RHS
               rswitch = 1._wp - MIN( 1._wp, ABS( SIGN( 1._wp, zRHS + ztaux_base(ji,jj) ) - SIGN( 1._wp, zRHS ) ) )
               !
               IF( ln_aEVP ) THEN !--- ice velocity using aEVP (Kimmritz et al 2016 & 2017)
                  zbetau = MAX( zbeta(ji,jj), zbeta(ji+1,jj) )
                  u_ice(ji,jj) = ( (          rswitch   * ( zmU_t(ji,jj) * ( zbetau * u_ice(ji,jj) + u_ice_b(ji,jj) )         & ! previous velocity
                     &                                    + zRHS + zTauO * u_ice(ji,jj)                                       & ! F + tau_ia + Coriolis + spg + tau_io(only ocean part)
                     &                                    ) / MAX( zepsi, zmU_t(ji,jj) * ( zbetau + 1._wp ) + zTauO - zTauB ) & ! m/dt + tau_io(only ice part) + landfast
                     &            + ( 1._wp - rswitch ) * (  u_ice_b(ji,jj)                                                   &
                     &                                     + u_ice  (ji,jj) * MAX( 0._wp, zbetau - zdtevp * rn_lf_relax )     & ! static friction => slow decrease to v=0
                     &                                    ) / ( zbetau + 1._wp )                                              &
                     &             ) * zmsk01x(ji,jj) + u_oce(ji,jj) * 0.01_wp * ( 1._wp - zmsk01x(ji,jj) )                   & ! v_ice = v_oce/100 if mass < zmmin & conc < zamin
                     &           )   * zmsk00x(ji,jj)
               ELSE               !--- ice velocity using EVP implicit formulation (cf Madec doc & Bouillon 2009)
                  u_ice(ji,jj) = ( (         rswitch   * ( zmU_t(ji,jj)  * u_ice(ji,jj)                                       & ! previous velocity
                     &                                    + zRHS + zTauO * u_ice(ji,jj)                                       & ! F + tau_ia + Coriolis + spg + tau_io(only ocean part)
                     &                                    ) / MAX( zepsi, zmU_t(ji,jj) + zTauO - zTauB )                      & ! m/dt + tau_io(only ice part) + landfast
                     &            + ( 1._wp - rswitch ) *   u_ice(ji,jj) * MAX( 0._wp, 1._wp - zdtevp * rn_lf_relax )         & ! static friction => slow decrease to v=0
                     &              ) * zmsk01x(ji,jj) + u_oce(ji,jj) * 0.01_wp * ( 1._wp - zmsk01x(ji,jj) )                  & ! v_ice = v_oce/100 if mass < zmmin & conc < zamin
                     &            )   * zmsk00x(ji,jj)
               ENDIF
!extra code for test case boundary conditions
               IF (mjg(jj)<25 .or. mjg(jj)>975 .or. mig(ji)<25 .or. mig(ji)>975) THEN
                  u_ice(ji,jj) = zinvw*(ztaux_ai(ji,jj) + zCorU(ji,jj) + zspgU(ji,jj) + ztaux_oi(ji,jj))
               END IF
            END_2D
            CALL lbc_lnk( 'icedyn_rhg_eap', u_ice, 'U', -1.0_wp )
            !
#if defined key_agrif
!!          CALL agrif_interp_ice( 'U', jter, nn_nevp )
            CALL agrif_interp_ice( 'U' )
#endif
            IF( ln_bdy )   CALL bdy_ice_dyn( 'U' )
            !
            DO_2D( 0, 0, 0, 0 )
               !                 !--- tau_io/(v_oce - v_ice)
               zTauO = zaV(ji,jj) * zrhoco * SQRT( ( v_ice (ji,jj) - v_oce (ji,jj) ) * ( v_ice (ji,jj) - v_oce (ji,jj) )  &
                  &                              + ( u_iceV(ji,jj) - u_oceV(ji,jj) ) * ( u_iceV(ji,jj) - u_oceV(ji,jj) ) )
               !                 !--- Ocean-to-Ice stress
               ztauy_oi(ji,jj) = zTauO * ( v_oce(ji,jj) - v_ice(ji,jj) )
               !
               !                 !--- tau_bottom/v_ice
               zvel  = 5.e-05_wp + SQRT( v_ice(ji,jj) * v_ice(ji,jj) + u_iceV(ji,jj) * u_iceV(ji,jj) )
               zTauB = ztauy_base(ji,jj) / zvel
               !                 !--- OceanBottom-to-Ice stress
               ztauy_bi(ji,jj) = zTauB * v_ice(ji,jj)
               !
               !                 !--- Coriolis at v-points (energy conserving formulation)
               zCorV(ji,jj)  = - 0.25_wp * r1_e2v(ji,jj) *  &
                  &    ( zmf(ji,jj  ) * ( e2u(ji,jj  ) * u_ice(ji,jj  ) + e2u(ji-1,jj  ) * u_ice(ji-1,jj  ) )  &
                  &    + zmf(ji,jj+1) * ( e2u(ji,jj+1) * u_ice(ji,jj+1) + e2u(ji-1,jj+1) * u_ice(ji-1,jj+1) ) )
               !
               !                 !--- Sum of external forces (explicit solution) = F + tau_ia + Coriolis + spg + tau_io
               zRHS = zfV(ji,jj) + ztauy_ai(ji,jj) + zCorV(ji,jj) + zspgV(ji,jj) + ztauy_oi(ji,jj)
               !
               !                 !--- landfast switch => 0 = static  friction : TauB > RHS & sign(TauB) /= sign(RHS)
               !                                         1 = sliding friction : TauB < RHS
               rswitch = 1._wp - MIN( 1._wp, ABS( SIGN( 1._wp, zRHS + ztauy_base(ji,jj) ) - SIGN( 1._wp, zRHS ) ) )
               !
               IF( ln_aEVP ) THEN !--- ice velocity using aEVP (Kimmritz et al 2016 & 2017)
                  zbetav = MAX( zbeta(ji,jj), zbeta(ji,jj+1) )
                  v_ice(ji,jj) = ( (          rswitch   * ( zmV_t(ji,jj) * ( zbetav * v_ice(ji,jj) + v_ice_b(ji,jj) )         & ! previous velocity
                     &                                    + zRHS + zTauO * v_ice(ji,jj)                                       & ! F + tau_ia + Coriolis + spg + tau_io(only ocean part)
                     &                                    ) / MAX( zepsi, zmV_t(ji,jj) * ( zbetav + 1._wp ) + zTauO - zTauB ) & ! m/dt + tau_io(only ice part) + landfast
                     &            + ( 1._wp - rswitch ) * (  v_ice_b(ji,jj)                                                   &
                     &                                     + v_ice  (ji,jj) * MAX( 0._wp, zbetav - zdtevp * rn_lf_relax )     & ! static friction => slow decrease to v=0
                     &                                    ) / ( zbetav + 1._wp )                                              &
                     &             ) * zmsk01y(ji,jj) + v_oce(ji,jj) * 0.01_wp * ( 1._wp - zmsk01y(ji,jj) )                   & ! v_ice = v_oce/100 if mass < zmmin & conc < zamin
                     &           )   * zmsk00y(ji,jj)
               ELSE               !--- ice velocity using EVP implicit formulation (cf Madec doc & Bouillon 2009)
                  v_ice(ji,jj) = ( (         rswitch   * ( zmV_t(ji,jj)  * v_ice(ji,jj)                                       & ! previous velocity
                     &                                    + zRHS + zTauO * v_ice(ji,jj)                                       & ! F + tau_ia + Coriolis + spg + tau_io(only ocean part)
                     &                                    ) / MAX( zepsi, zmV_t(ji,jj) + zTauO - zTauB )                      & ! m/dt + tau_io(only ice part) + landfast
                     &            + ( 1._wp - rswitch ) *   v_ice(ji,jj) * MAX( 0._wp, 1._wp - zdtevp * rn_lf_relax )         & ! static friction => slow decrease to v=0
                     &              ) * zmsk01y(ji,jj) + v_oce(ji,jj) * 0.01_wp * ( 1._wp - zmsk01y(ji,jj) )                  & ! v_ice = v_oce/100 if mass < zmmin & conc < zamin
                     &            )   * zmsk00y(ji,jj)
               ENDIF
!extra code for test case boundary conditions
               IF (mjg(jj)<25 .or. mjg(jj)>975 .or. mig(ji)<25 .or. mig(ji)>975) THEN
                  v_ice(ji,jj) = zinvw*(ztauy_ai(ji,jj) + zCorV(ji,jj) + zspgV(ji,jj) + ztauy_oi(ji,jj))
               END IF
            END_2D
            CALL lbc_lnk( 'icedyn_rhg_eap', v_ice, 'V', -1.0_wp )
            !
#if defined key_agrif
!!          CALL agrif_interp_ice( 'V', jter, nn_nevp )
            CALL agrif_interp_ice( 'V' )
#endif
            IF( ln_bdy )   CALL bdy_ice_dyn( 'V' )
            !
         ENDIF

         ! convergence test
         IF( nn_rhg_chkcvg == 2 )   CALL rhg_cvg_eap( kt, jter, nn_nevp, u_ice, v_ice, zu_ice, zv_ice, zmsk15 )
         !
         !                                                ! ==================== !
      END DO                                              !  end loop over jter  !
      !                                                   ! ==================== !
      IF( ln_aEVP )   CALL iom_put( 'beta_evp' , zbeta )
      !
      CALL lbc_lnk( 'icedyn_rhg_eap', prdg_conv, 'T', 1.0_wp )  ! only need this in ridging module after rheology completed
      !
      !------------------------------------------------------------------------------!
      ! 4) Recompute delta, shear and div (inputs for mechanical redistribution)
      !------------------------------------------------------------------------------!
      DO_2D( 1, 0, 1, 0 )

         ! shear at F points
         zds(ji,jj) = ( ( u_ice(ji,jj+1) * r1_e1u(ji,jj+1) - u_ice(ji,jj) * r1_e1u(ji,jj) ) * e1f(ji,jj) * e1f(ji,jj)   &
            &         + ( v_ice(ji+1,jj) * r1_e2v(ji+1,jj) - v_ice(ji,jj) * r1_e2v(ji,jj) ) * e2f(ji,jj) * e2f(ji,jj)   &
            &         ) * r1_e1e2f(ji,jj) * fimask(ji,jj)

      END_2D

      DO_2D( 0, 0, 0, 0 )

         ! tension**2 at T points
         zdt  = ( ( u_ice(ji,jj) * r1_e2u(ji,jj) - u_ice(ji-1,jj) * r1_e2u(ji-1,jj) ) * e2t(ji,jj) * e2t(ji,jj)   &
            &   - ( v_ice(ji,jj) * r1_e1v(ji,jj) - v_ice(ji,jj-1) * r1_e1v(ji,jj-1) ) * e1t(ji,jj) * e1t(ji,jj)   &
            &   ) * r1_e1e2t(ji,jj)
         zdt2 = zdt * zdt

         zten_i(ji,jj) = zdt

         ! shear**2 at T points (doc eq. A16)
         zds2 = ( zds(ji,jj  ) * zds(ji,jj  ) * e1e2f(ji,jj  ) + zds(ji-1,jj  ) * zds(ji-1,jj  ) * e1e2f(ji-1,jj  )  &
            &   + zds(ji,jj-1) * zds(ji,jj-1) * e1e2f(ji,jj-1) + zds(ji-1,jj-1) * zds(ji-1,jj-1) * e1e2f(ji-1,jj-1)  &
            &   ) * 0.25_wp * r1_e1e2t(ji,jj)

         ! shear at T points
         pshear_i(ji,jj) = SQRT( zdt2 + zds2 )

         ! divergence at T points
         pdivu_i(ji,jj) = ( e2u(ji,jj) * u_ice(ji,jj) - e2u(ji-1,jj) * u_ice(ji-1,jj)   &
            &             + e1v(ji,jj) * v_ice(ji,jj) - e1v(ji,jj-1) * v_ice(ji,jj-1)   &
            &             ) * r1_e1e2t(ji,jj)

         ! delta at T points
         zfac            = SQRT( pdivu_i(ji,jj) * pdivu_i(ji,jj) + ( zdt2 + zds2 ) * z1_ecc2 ) ! delta
         rswitch         = 1._wp - MAX( 0._wp, SIGN( 1._wp, -zfac ) ) ! 0 if delta=0
         pdelta_i(ji,jj) = zfac + rn_creepl * rswitch ! delta+creepl

      END_2D
      CALL lbc_lnk( 'icedyn_rhg_eap', pshear_i, 'T', 1.0_wp, pdivu_i, 'T', 1.0_wp, pdelta_i, 'T', 1.0_wp, &
         &                              zten_i, 'T', 1.0_wp, zs1    , 'T', 1.0_wp, zs2     , 'T', 1.0_wp, &
         &                                zs12, 'F', 1.0_wp )

      ! --- Store the stress tensor for the next time step --- !
      pstress1_i (:,:) = zs1 (:,:)
      pstress2_i (:,:) = zs2 (:,:)
      pstress12_i(:,:) = zs12(:,:)
      !

      !------------------------------------------------------------------------------!
      ! 5) diagnostics
      !------------------------------------------------------------------------------!
      ! --- ice-ocean, ice-atm. & ice-oceanbottom(landfast) stresses --- !
      IF(  iom_use('utau_oi') .OR. iom_use('vtau_oi') .OR. iom_use('utau_ai') .OR. iom_use('vtau_ai') .OR. &
         & iom_use('utau_bi') .OR. iom_use('vtau_bi') ) THEN
         !
         CALL lbc_lnk( 'icedyn_rhg_eap', ztaux_oi, 'U', -1.0_wp, ztauy_oi, 'V', -1.0_wp, ztaux_ai, 'U', -1.0_wp, &
            &                            ztauy_ai, 'V', -1.0_wp, ztaux_bi, 'U', -1.0_wp, ztauy_bi, 'V', -1.0_wp )
         !
         CALL iom_put( 'utau_oi' , ztaux_oi * zmsk00 )
         CALL iom_put( 'vtau_oi' , ztauy_oi * zmsk00 )
         CALL iom_put( 'utau_ai' , ztaux_ai * zmsk00 )
         CALL iom_put( 'vtau_ai' , ztauy_ai * zmsk00 )
         CALL iom_put( 'utau_bi' , ztaux_bi * zmsk00 )
         CALL iom_put( 'vtau_bi' , ztauy_bi * zmsk00 )
      ENDIF

      ! --- divergence, shear and strength --- !
      IF( iom_use('icediv') )   CALL iom_put( 'icediv' , pdivu_i  * zmsk00 )   ! divergence
      IF( iom_use('iceshe') )   CALL iom_put( 'iceshe' , pshear_i * zmsk00 )   ! shear
      IF( iom_use('icedlt') )   CALL iom_put( 'icedlt' , pdelta_i * zmsk00 )   ! delta
      IF( iom_use('icestr') )   CALL iom_put( 'icestr' , strength * zmsk00 )   ! strength

      ! --- Stress tensor invariants (SIMIP diags) --- !
      IF( iom_use('normstr') .OR. iom_use('sheastr') ) THEN
         !
         ALLOCATE( zsig_I(jpi,jpj) , zsig_II(jpi,jpj) )
         !
         DO_2D( 1, 1, 1, 1 )

            ! Ice stresses
            ! sigma1, sigma2, sigma12 are some useful recombination of the stresses (Hunke and Dukowicz MWR 2002, Bouillon et al., OM2013)
            ! These are NOT stress tensor components, neither stress invariants, neither stress principal components
            ! I know, this can be confusing...
            zfac             =   strength(ji,jj) / ( pdelta_i(ji,jj) + rn_creepl )
            zsig1            =   zfac * ( pdivu_i(ji,jj) - pdelta_i(ji,jj) )
            zsig2            =   zfac * z1_ecc2 * zten_i(ji,jj)
            zsig12           =   zfac * z1_ecc2 * pshear_i(ji,jj)

            ! Stress invariants (sigma_I, sigma_II, Coon 1974, Feltham 2008)
            zsig_I (ji,jj)   =   zsig1 * 0.5_wp                                           ! 1st stress invariant, aka average normal stress, aka negative pressure
            zsig_II(ji,jj)   =   SQRT ( MAX( 0._wp, zsig2 * zsig2 * 0.25_wp + zsig12 ) )  ! 2nd  ''       '', aka maximum shear stress

         END_2D
         !
         ! Stress tensor invariants (normal and shear stress N/m) - SIMIP diags - definitions following Coon (1974) and Feltham (2008)
         IF( iom_use('normstr') )   CALL iom_put( 'normstr', zsig_I (:,:) * zmsk00(:,:) ) ! Normal stress
         IF( iom_use('sheastr') )   CALL iom_put( 'sheastr', zsig_II(:,:) * zmsk00(:,:) ) ! Maximum shear stress

         DEALLOCATE ( zsig_I, zsig_II )

      ENDIF

      ! --- Normalized stress tensor principal components --- !
      ! This are used to plot the normalized yield curve, see Lemieux & Dupont, 2020
      ! Recommendation 1 : we use ice strength, not replacement pressure
      ! Recommendation 2 : need to use deformations at PREVIOUS iterate for viscosities
      IF( iom_use('sig1_pnorm') .OR. iom_use('sig2_pnorm') ) THEN
         !
         ALLOCATE( zsig1_p(jpi,jpj) , zsig2_p(jpi,jpj) , zsig_I(jpi,jpj) , zsig_II(jpi,jpj) )
         !
         DO_2D( 1, 1, 1, 1 )

            ! Ice stresses computed with **viscosities** (delta, p/delta) at **previous** iterates
            !                        and **deformations** at current iterates
            !                        following Lemieux & Dupont (2020)
            zfac             =   zp_delt(ji,jj)
            zsig1            =   zfac * ( pdivu_i(ji,jj) - ( zdelta(ji,jj) + rn_creepl ) )
            zsig2            =   zfac * z1_ecc2 * zten_i(ji,jj)
            zsig12           =   zfac * z1_ecc2 * pshear_i(ji,jj)

            ! Stress invariants (sigma_I, sigma_II, Coon 1974, Feltham 2008), T-point
            zsig_I(ji,jj)    =   zsig1 * 0.5_wp                                           ! 1st stress invariant, aka average normal stress, aka negative pressure
            zsig_II(ji,jj)   =   SQRT ( MAX( 0._wp, zsig2 * zsig2 * 0.25_wp + zsig12 ) )  ! 2nd  ''       '', aka maximum shear stress

            ! Normalized  principal stresses (used to display the ellipse)
            z1_strength      =   1._wp / MAX( 1._wp, strength(ji,jj) )
            zsig1_p(ji,jj)   =   ( zsig_I(ji,jj) + zsig_II(ji,jj) ) * z1_strength
            zsig2_p(ji,jj)   =   ( zsig_I(ji,jj) - zsig_II(ji,jj) ) * z1_strength
         END_2D
         !
         CALL iom_put( 'sig1_pnorm' , zsig1_p )
         CALL iom_put( 'sig2_pnorm' , zsig2_p )

         DEALLOCATE( zsig1_p , zsig2_p , zsig_I, zsig_II )

      ENDIF

      ! --- yieldcurve --- !
      IF( iom_use('yield11') .OR. iom_use('yield12') .OR. iom_use('yield22')) THEN

         CALL lbc_lnk( 'icedyn_rhg_eap', zyield11, 'T', 1.0_wp, zyield22, 'T', 1.0_wp, zyield12, 'T', 1.0_wp )

         CALL iom_put( 'yield11', zyield11 * zmsk00 )
         CALL iom_put( 'yield22', zyield22 * zmsk00 )
         CALL iom_put( 'yield12', zyield12 * zmsk00 )
      ENDIF

      ! --- anisotropy tensor --- !
      IF( iom_use('aniso') ) THEN
         CALL lbc_lnk( 'icedyn_rhg_eap', paniso_11, 'T', 1.0_wp )
         CALL iom_put( 'aniso' , paniso_11 * zmsk00 )
      ENDIF

      ! --- SIMIP --- !
      IF(  iom_use('dssh_dx') .OR. iom_use('dssh_dy') .OR. &
         & iom_use('corstrx') .OR. iom_use('corstry') .OR. iom_use('intstrx') .OR. iom_use('intstry') ) THEN
         !
         CALL lbc_lnk( 'icedyn_rhg_eap', zspgU, 'U', -1.0_wp, zspgV, 'V', -1.0_wp, &
            &                            zCorU, 'U', -1.0_wp, zCorV, 'V', -1.0_wp, &
            &                              zfU, 'U', -1.0_wp,   zfV, 'V', -1.0_wp )

         CALL iom_put( 'dssh_dx' , zspgU * zmsk00 )   ! Sea-surface tilt term in force balance (x)
         CALL iom_put( 'dssh_dy' , zspgV * zmsk00 )   ! Sea-surface tilt term in force balance (y)
         CALL iom_put( 'corstrx' , zCorU * zmsk00 )   ! Coriolis force term in force balance (x)
         CALL iom_put( 'corstry' , zCorV * zmsk00 )   ! Coriolis force term in force balance (y)
         CALL iom_put( 'intstrx' , zfU   * zmsk00 )   ! Internal force term in force balance (x)
         CALL iom_put( 'intstry' , zfV   * zmsk00 )   ! Internal force term in force balance (y)
      ENDIF

      IF(  iom_use('xmtrpice') .OR. iom_use('ymtrpice') .OR. &
         & iom_use('xmtrpsnw') .OR. iom_use('ymtrpsnw') .OR. iom_use('xatrp') .OR. iom_use('yatrp') ) THEN
         !
         ALLOCATE( zdiag_xmtrp_ice(jpi,jpj) , zdiag_ymtrp_ice(jpi,jpj) , &
            &      zdiag_xmtrp_snw(jpi,jpj) , zdiag_ymtrp_snw(jpi,jpj) , zdiag_xatrp(jpi,jpj) , zdiag_yatrp(jpi,jpj) )
         !
         DO_2D( 0, 0, 0, 0 )
            ! 2D ice mass, snow mass, area transport arrays (X, Y)
            zfac_x = 0.5 * u_ice(ji,jj) * e2u(ji,jj) * zmsk00(ji,jj)
            zfac_y = 0.5 * v_ice(ji,jj) * e1v(ji,jj) * zmsk00(ji,jj)

            zdiag_xmtrp_ice(ji,jj) = rhoi * zfac_x * ( vt_i(ji+1,jj) + vt_i(ji,jj) ) ! ice mass transport, X-component
            zdiag_ymtrp_ice(ji,jj) = rhoi * zfac_y * ( vt_i(ji,jj+1) + vt_i(ji,jj) ) !        ''           Y-   ''

            zdiag_xmtrp_snw(ji,jj) = rhos * zfac_x * ( vt_s(ji+1,jj) + vt_s(ji,jj) ) ! snow mass transport, X-component
            zdiag_ymtrp_snw(ji,jj) = rhos * zfac_y * ( vt_s(ji,jj+1) + vt_s(ji,jj) ) !          ''          Y-   ''

            zdiag_xatrp(ji,jj)     = zfac_x * ( at_i(ji+1,jj) + at_i(ji,jj) )        ! area transport,      X-component
            zdiag_yatrp(ji,jj)     = zfac_y * ( at_i(ji,jj+1) + at_i(ji,jj) )        !        ''            Y-   ''

         END_2D

         CALL lbc_lnk( 'icedyn_rhg_eap', zdiag_xmtrp_ice, 'U', -1.0_wp, zdiag_ymtrp_ice, 'V', -1.0_wp, &
            &                            zdiag_xmtrp_snw, 'U', -1.0_wp, zdiag_ymtrp_snw, 'V', -1.0_wp, &
            &                            zdiag_xatrp    , 'U', -1.0_wp, zdiag_yatrp    , 'V', -1.0_wp )

         CALL iom_put( 'xmtrpice' , zdiag_xmtrp_ice )   ! X-component of sea-ice mass transport (kg/s)
         CALL iom_put( 'ymtrpice' , zdiag_ymtrp_ice )   ! Y-component of sea-ice mass transport
         CALL iom_put( 'xmtrpsnw' , zdiag_xmtrp_snw )   ! X-component of snow mass transport (kg/s)
         CALL iom_put( 'ymtrpsnw' , zdiag_ymtrp_snw )   ! Y-component of snow mass transport
         CALL iom_put( 'xatrp'    , zdiag_xatrp     )   ! X-component of ice area transport
         CALL iom_put( 'yatrp'    , zdiag_yatrp     )   ! Y-component of ice area transport

         DEALLOCATE( zdiag_xmtrp_ice , zdiag_ymtrp_ice , &
            &        zdiag_xmtrp_snw , zdiag_ymtrp_snw , zdiag_xatrp , zdiag_yatrp )

      ENDIF
      !
      ! --- convergence tests --- !
      IF( nn_rhg_chkcvg == 1 .OR. nn_rhg_chkcvg == 2 ) THEN
         IF( iom_use('uice_cvg') ) THEN
            IF( ln_aEVP ) THEN   ! output: beta * ( u(t=nn_nevp) - u(t=nn_nevp-1) )
               CALL iom_put( 'uice_cvg', MAX( ABS( u_ice(:,:) - zu_ice(:,:) ) * zbeta(:,:) * umask(:,:,1) , &
                  &                           ABS( v_ice(:,:) - zv_ice(:,:) ) * zbeta(:,:) * vmask(:,:,1) ) * zmsk15(:,:) )
            ELSE                 ! output: nn_nevp * ( u(t=nn_nevp) - u(t=nn_nevp-1) )
               CALL iom_put( 'uice_cvg', REAL( nn_nevp ) * MAX( ABS( u_ice(:,:) - zu_ice(:,:) ) * umask(:,:,1) , &
                  &                                             ABS( v_ice(:,:) - zv_ice(:,:) ) * vmask(:,:,1) ) * zmsk15(:,:) )
            ENDIF
         ENDIF
      ENDIF
      !
   END SUBROUTINE ice_dyn_rhg_eap


   SUBROUTINE rhg_cvg_eap( kt, kiter, kitermax, pu, pv, pub, pvb, pmsk15 )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE rhg_cvg_eap  ***
      !!
      !! ** Purpose :   check convergence of oce rheology
      !!
      !! ** Method  :   create a file ice_cvg.nc containing the convergence of ice velocity
      !!                during the sub timestepping of rheology so as:
      !!                  uice_cvg = MAX( u(t+1) - u(t) , v(t+1) - v(t) )
      !!                This routine is called every sub-iteration, so it is cpu expensive
      !!
      !! ** Note    :   for the first sub-iteration, uice_cvg is set to 0 (too large otherwise)
      !!----------------------------------------------------------------------
      INTEGER ,                 INTENT(in) ::   kt, kiter, kitermax       ! ocean time-step index
      REAL(wp), DIMENSION(:,:), INTENT(in) ::   pu, pv, pub, pvb          ! now and before velocities
      REAL(wp), DIMENSION(:,:), INTENT(in) ::   pmsk15
      !!
      INTEGER           ::   it, idtime, istatus
      INTEGER           ::   ji, jj          ! dummy loop indices
      REAL(wp)          ::   zresm           ! local real
      CHARACTER(len=20) ::   clname
      !!----------------------------------------------------------------------

      ! create file
      IF( kt == nit000 .AND. kiter == 1 ) THEN
         !
         IF( lwp ) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'rhg_cvg_eap : ice rheology convergence control'
            WRITE(numout,*) '~~~~~~~'
         ENDIF
         !
         IF( lwm ) THEN
            clname = 'ice_cvg.nc'
            IF( .NOT. Agrif_Root() )   clname = TRIM(Agrif_CFixed())//"_"//TRIM(clname)
            istatus = NF90_CREATE( TRIM(clname), NF90_CLOBBER, ncvgid )
            istatus = NF90_DEF_DIM( ncvgid, 'time'  , NF90_UNLIMITED, idtime )
            istatus = NF90_DEF_VAR( ncvgid, 'uice_cvg', NF90_DOUBLE   , (/ idtime /), nvarid )
            istatus = NF90_ENDDEF(ncvgid)
         ENDIF
         !
      ENDIF

      ! time
      it = ( kt - 1 ) * kitermax + kiter

      ! convergence
      IF( kiter == 1 ) THEN ! remove the first iteration for calculations of convergence (always very large)
         zresm = 0._wp
      ELSE
         zresm = 0._wp
         DO_2D( 0, 0, 0, 0 )
            ! cut of the boundary of the box (forced velocities)
            IF (mjg0(jj)>30 .AND. mjg0(jj)<=970 .AND. mig0(ji)>30 .AND. mig0(ji)<=970) THEN
               zresm = MAX( zresm, MAX( ABS( pu(ji,jj) - pub(ji,jj) ) * umask(ji,jj,1), &
                  &                     ABS( pv(ji,jj) - pvb(ji,jj) ) * vmask(ji,jj,1) ) * pmsk15(ji,jj) )
            ENDIF
         END_2D
         CALL mpp_max( 'icedyn_rhg_evp', zresm )   ! max over the global domain
      ENDIF

      IF( lwm ) THEN
         ! write variables
         istatus = NF90_PUT_VAR( ncvgid, nvarid, (/zresm/), (/it/), (/1/) )
         ! close file
         IF( kt == nitend )   istatus = NF90_CLOSE(ncvgid)
      ENDIF

   END SUBROUTINE rhg_cvg_eap


   SUBROUTINE update_stress_rdg( ksub, kndte, pdivu, ptension, pshear, pa11, pa12, &
      &                          pstressp,  pstressm, pstress12, pstrength, palphar, palphas )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE update_stress_rdg  ***
      !!
      !! ** Purpose :   Updates the stress depending on values of strain rate and structure
      !!                tensor and for the last subcycle step it computes closing and sliding rate
      !!---------------------------------------------------------------------
      INTEGER,  INTENT(in   ) ::   ksub, kndte
      REAL(wp), INTENT(in   ) ::   pstrength
      REAL(wp), INTENT(in   ) ::   pdivu, ptension, pshear
      REAL(wp), INTENT(in   ) ::   pa11, pa12
      REAL(wp), INTENT(  out) ::   pstressp, pstressm, pstress12
      REAL(wp), INTENT(  out) ::   palphar, palphas

      INTEGER  ::   kx ,ky, ka

      REAL(wp) ::   zstemp11r, zstemp12r, zstemp22r
      REAL(wp) ::   zstemp11s, zstemp12s, zstemp22s
      REAL(wp) ::   za22, zQd11Qd11, zQd11Qd12, zQd12Qd12
      REAL(wp) ::   zQ11Q11, zQ11Q12, zQ12Q12
      REAL(wp) ::   zdtemp11, zdtemp12, zdtemp22
      REAL(wp) ::   zrotstemp11r, zrotstemp12r, zrotstemp22r
      REAL(wp) ::   zrotstemp11s, zrotstemp12s, zrotstemp22s
      REAL(wp) ::   zsig11, zsig12, zsig22
      REAL(wp) ::   zsgprm11, zsgprm12, zsgprm22
      REAL(wp) ::   zAngle_denom_gamma,  zAngle_denom_alpha
      REAL(wp) ::   zTany_1, zTany_2
      REAL(wp) ::   zx, zy, zkxw, kyw, kaw
      REAL(wp) ::   zinvdx, zinvdy, zinvda
      REAL(wp) ::   zdtemp1, zdtemp2, zatempprime

      REAL(wp), PARAMETER ::   ppkfriction = 0.45_wp
      ! Factor to maintain the same stress as in EVP (see Section 3)
      ! Can be set to 1 otherwise
!     REAL(wp), PARAMETER ::   ppinvstressconviso = 1._wp/(1._wp+ppkfriction*ppkfriction)
      REAL(wp), PARAMETER ::   ppinvstressconviso = 1._wp

      ! next statement uses pphi set in main module (icedyn_rhg_eap)
      REAL(wp), PARAMETER ::   ppinvsin = 1._wp/sin(2._wp*pphi) * ppinvstressconviso

      ! compute eigenvalues, eigenvectors and angles for structure tensor, strain
      ! rates

      ! 1) structure tensor
      za22 = 1._wp-pa11
      zQ11Q11 = 1._wp
      zQ12Q12 = rsmall
      zQ11Q12 = rsmall

      ! gamma: angle between general coordiantes and principal axis of A
      ! here Tan2gamma = 2 a12 / (a11 - a22)
      IF((ABS(pa11 - za22) > rsmall).OR.(ABS(pa12) > rsmall)) THEN
         zAngle_denom_gamma = 1._wp/sqrt( ( pa11 - za22 )*( pa11 - za22) + &
                              4._wp*pa12*pa12 )

         zQ11Q11 = 0.5_wp + ( pa11 - za22 )*0.5_wp*zAngle_denom_gamma   !Cos^2
         zQ12Q12 = 0.5_wp - ( pa11 - za22 )*0.5_wp*zAngle_denom_gamma   !Sin^2
         zQ11Q12 = pa12*zAngle_denom_gamma                     !CosSin
      ENDIF

      ! rotation Q*atemp*Q^T
      zatempprime = zQ11Q11*pa11 + 2.0_wp*zQ11Q12*pa12 + zQ12Q12*za22

      ! make first principal value the largest
      zatempprime = max(zatempprime, 1.0_wp - zatempprime)

      ! 2) strain rate
      zdtemp11 = 0.5_wp*(pdivu + ptension)
      zdtemp12 = pshear*0.5_wp
      zdtemp22 = 0.5_wp*(pdivu - ptension)

      ! here Tan2alpha = 2 dtemp12 / (dtemp11 - dtemp22)

      zQd11Qd11 = 1.0_wp
      zQd12Qd12 = rsmall
      zQd11Qd12 = rsmall

      IF((ABS( zdtemp11 - zdtemp22) > rsmall).OR. (ABS(zdtemp12) > rsmall)) THEN

         zAngle_denom_alpha = 1.0_wp/sqrt( ( zdtemp11 - zdtemp22 )* &
                             ( zdtemp11 - zdtemp22 ) + 4.0_wp*zdtemp12*zdtemp12)

         zQd11Qd11 = 0.5_wp + ( zdtemp11 - zdtemp22 )*0.5_wp*zAngle_denom_alpha !Cos^2
         zQd12Qd12 = 0.5_wp - ( zdtemp11 - zdtemp22 )*0.5_wp*zAngle_denom_alpha !Sin^2
         zQd11Qd12 = zdtemp12*zAngle_denom_alpha !CosSin
      ENDIF

      zdtemp1 = zQd11Qd11*zdtemp11 + 2.0_wp*zQd11Qd12*zdtemp12 + zQd12Qd12*zdtemp22
      zdtemp2 = zQd12Qd12*zdtemp11 - 2.0_wp*zQd11Qd12*zdtemp12 + zQd11Qd11*zdtemp22
      ! In cos and sin values
      zx = 0._wp
      IF ((ABS(zdtemp1) > rsmall).OR.(ABS(zdtemp2) > rsmall)) THEN
         zx = atan2(zdtemp2,zdtemp1)
      ENDIF

      ! to ensure the angle lies between pi/4 and 9 pi/4
      IF (zx < rpi*0.25_wp) zx = zx + rpi*2.0_wp

      ! y: angle between major principal axis of strain rate and structure
      ! tensor
      ! y = gamma - alpha
      ! Expressesed componently with
      ! Tany = (Singamma*Cosgamma - Sinalpha*Cosgamma)/(Cos^2gamma - Sin^alpha)

      zTany_1 = zQ11Q12 - zQd11Qd12
      zTany_2 = zQ11Q11 - zQd12Qd12

      zy = 0._wp

      IF ((ABS(zTany_1) > rsmall).OR.(ABS(zTany_2) > rsmall)) THEN
         zy = atan2(zTany_1,zTany_2)
      ENDIF

      ! to make sure y is between 0 and pi
      IF (zy > rpi) zy = zy - rpi
      IF (zy < 0)  zy = zy + rpi

      ! 3) update anisotropic stress tensor
      zinvdx = real(nx_yield-1,kind=wp)/rpi
      zinvdy = real(ny_yield-1,kind=wp)/rpi
      zinvda = 2._wp*real(na_yield-1,kind=wp)

      ! % need 8 coords and 8 weights
      ! % range in kx
      kx  = int((zx-rpi*0.25_wp-rpi)*zinvdx) + 1
      !!clem kx  = MAX( 1, MIN( nx_yield-1, INT((zx-rpi*0.25_wp-rpi)*zinvdx) + 1  ) )
      zkxw = kx - (zx-rpi*0.25_wp-rpi)*zinvdx

      ky  = int(zy*zinvdy) + 1
      !!clem ky  = MAX( 1, MIN( ny_yield-1, INT(zy*zinvdy) + 1 ) )
      kyw = ky - zy*zinvdy

      ka  = int((zatempprime-0.5_wp)*zinvda) + 1
      !!clem ka  = MAX( 1, MIN( na_yield-1, INT((zatempprime-0.5_wp)*zinvda) + 1 ) )
      kaw = ka - (zatempprime-0.5_wp)*zinvda

      ! % Determine sigma_r(A1,Zeta,y) and sigma_s (see Section A1 of Tsamados 2013)
!!$      zstemp11r =  zkxw  * kyw         * kaw         * s11r(kx  ,ky  ,ka  ) &
!!$        & + (1._wp-zkxw) * kyw         * kaw         * s11r(kx+1,ky  ,ka  ) &
!!$        & + zkxw         * (1._wp-kyw) * kaw         * s11r(kx  ,ky+1,ka  ) &
!!$        & + zkxw         * kyw         * (1._wp-kaw) * s11r(kx  ,ky  ,ka+1) &
!!$        & + (1._wp-zkxw) * (1._wp-kyw) * kaw         * s11r(kx+1,ky+1,ka  ) &
!!$        & + (1._wp-zkxw) * kyw         * (1._wp-kaw) * s11r(kx+1,ky  ,ka+1) &
!!$        & + zkxw         * (1._wp-kyw) * (1._wp-kaw) * s11r(kx  ,ky+1,ka+1) &
!!$        & + (1._wp-zkxw) * (1._wp-kyw) * (1._wp-kaw) * s11r(kx+1,ky+1,ka+1)
!!$      zstemp12r =  zkxw  * kyw         * kaw         * s12r(kx  ,ky  ,ka  ) &
!!$        & + (1._wp-zkxw) * kyw         * kaw         * s12r(kx+1,ky  ,ka  ) &
!!$        & + zkxw         * (1._wp-kyw) * kaw         * s12r(kx  ,ky+1,ka  ) &
!!$        & + zkxw         * kyw         * (1._wp-kaw) * s12r(kx  ,ky  ,ka+1) &
!!$        & + (1._wp-zkxw) * (1._wp-kyw) * kaw         * s12r(kx+1,ky+1,ka  ) &
!!$        & + (1._wp-zkxw) * kyw         * (1._wp-kaw) * s12r(kx+1,ky  ,ka+1) &
!!$        & + zkxw         * (1._wp-kyw) * (1._wp-kaw) * s12r(kx  ,ky+1,ka+1) &
!!$        & + (1._wp-zkxw) * (1._wp-kyw) * (1._wp-kaw) * s12r(kx+1,ky+1,ka+1)
!!$      zstemp22r =  zkxw  * kyw         * kaw         * s22r(kx  ,ky  ,ka  ) &
!!$        & + (1._wp-zkxw) * kyw         * kaw         * s22r(kx+1,ky  ,ka  ) &
!!$        & + zkxw         * (1._wp-kyw) * kaw         * s22r(kx  ,ky+1,ka  ) &
!!$        & + zkxw         * kyw         * (1._wp-kaw) * s22r(kx  ,ky  ,ka+1) &
!!$        & + (1._wp-zkxw) * (1._wp-kyw) * kaw         * s22r(kx+1,ky+1,ka  ) &
!!$        & + (1._wp-zkxw) * kyw         * (1._wp-kaw) * s22r(kx+1,ky  ,ka+1) &
!!$        & + zkxw         * (1._wp-kyw) * (1._wp-kaw) * s22r(kx  ,ky+1,ka+1) &
!!$        & + (1._wp-zkxw) * (1._wp-kyw) * (1._wp-kaw) * s22r(kx+1,ky+1,ka+1)
!!$
!!$      zstemp11s =  zkxw  * kyw         * kaw         * s11s(kx  ,ky  ,ka  ) &
!!$        & + (1._wp-zkxw) * kyw         * kaw         * s11s(kx+1,ky  ,ka  ) &
!!$        & + zkxw         * (1._wp-kyw) * kaw         * s11s(kx  ,ky+1,ka  ) &
!!$        & + zkxw         * kyw         * (1._wp-kaw) * s11s(kx  ,ky  ,ka+1) &
!!$        & + (1._wp-zkxw) * (1._wp-kyw) * kaw         * s11s(kx+1,ky+1,ka  ) &
!!$        & + (1._wp-zkxw) * kyw         * (1._wp-kaw) * s11s(kx+1,ky  ,ka+1) &
!!$        & + zkxw         * (1._wp-kyw) * (1._wp-kaw) * s11s(kx  ,ky+1,ka+1) &
!!$        & + (1._wp-zkxw) * (1._wp-kyw) * (1._wp-kaw) * s11s(kx+1,ky+1,ka+1)
!!$      zstemp12s =  zkxw  * kyw         * kaw         * s12s(kx  ,ky  ,ka  ) &
!!$        & + (1._wp-zkxw) * kyw         * kaw         * s12s(kx+1,ky  ,ka  ) &
!!$        & + zkxw         * (1._wp-kyw) * kaw         * s12s(kx  ,ky+1,ka  ) &
!!$        & + zkxw         * kyw         * (1._wp-kaw) * s12s(kx  ,ky  ,ka+1) &
!!$        & + (1._wp-zkxw) * (1._wp-kyw) * kaw         * s12s(kx+1,ky+1,ka  ) &
!!$        & + (1._wp-zkxw) * kyw         * (1._wp-kaw) * s12s(kx+1,ky  ,ka+1) &
!!$        & + zkxw         * (1._wp-kyw) * (1._wp-kaw) * s12s(kx  ,ky+1,ka+1) &
!!$        & + (1._wp-zkxw) * (1._wp-kyw) * (1._wp-kaw) * s12s(kx+1,ky+1,ka+1)
!!$      zstemp22s =  zkxw  * kyw         * kaw         * s22s(kx  ,ky  ,ka  ) &
!!$        & + (1._wp-zkxw) * kyw         * kaw         * s22s(kx+1,ky  ,ka  ) &
!!$        & + zkxw         * (1._wp-kyw) * kaw         * s22s(kx  ,ky+1,ka  ) &
!!$        & + zkxw         * kyw         * (1._wp-kaw) * s22s(kx  ,ky  ,ka+1) &
!!$        & + (1._wp-zkxw) * (1._wp-kyw) * kaw         * s22s(kx+1,ky+1,ka  ) &
!!$        & + (1._wp-zkxw) * kyw         * (1._wp-kaw) * s22s(kx+1,ky  ,ka+1) &
!!$        & + zkxw         * (1._wp-kyw) * (1._wp-kaw) * s22s(kx  ,ky+1,ka+1) &
!!$        & + (1._wp-zkxw) * (1._wp-kyw) * (1._wp-kaw) * s22s(kx+1,ky+1,ka+1)

      zstemp11r = s11r(kx,ky,ka)
      zstemp12r = s12r(kx,ky,ka)
      zstemp22r = s22r(kx,ky,ka)
      zstemp11s = s11s(kx,ky,ka)
      zstemp12s = s12s(kx,ky,ka)
      zstemp22s = s22s(kx,ky,ka)


      ! Calculate mean ice stress over a collection of floes (Equation 3 in
      ! Tsamados 2013)

      zsig11  = pstrength*(zstemp11r + ppkfriction*zstemp11s) * ppinvsin
      zsig12  = pstrength*(zstemp12r + ppkfriction*zstemp12s) * ppinvsin
      zsig22  = pstrength*(zstemp22r + ppkfriction*zstemp22s) * ppinvsin

      ! Back - rotation of the stress from principal axes into general coordinates

      ! Update stress
      zsgprm11 = zQ11Q11*zsig11 + zQ12Q12*zsig22 -      2._wp*zQ11Q12 *zsig12
      zsgprm12 = zQ11Q12*zsig11 - zQ11Q12*zsig22 + (zQ11Q11 - zQ12Q12)*zsig12
      zsgprm22 = zQ12Q12*zsig11 + zQ11Q11*zsig22 +      2._wp*zQ11Q12 *zsig12

      pstressp  = zsgprm11 + zsgprm22
      pstress12 = zsgprm12
      pstressm  = zsgprm11 - zsgprm22

      ! Compute ridging and sliding functions in general coordinates
      ! (Equation 11 in Tsamados 2013)
      IF (ksub == kndte) THEN
         zrotstemp11r = zQ11Q11*zstemp11r - 2._wp*zQ11Q12* zstemp12r &
                      + zQ12Q12*zstemp22r
         zrotstemp12r = zQ11Q11*zstemp12r +       zQ11Q12*(zstemp11r-zstemp22r) &
                      - zQ12Q12*zstemp12r
         zrotstemp22r = zQ12Q12*zstemp11r + 2._wp*zQ11Q12* zstemp12r &
                      + zQ11Q11*zstemp22r

         zrotstemp11s = zQ11Q11*zstemp11s - 2._wp*zQ11Q12* zstemp12s &
                      + zQ12Q12*zstemp22s
         zrotstemp12s = zQ11Q11*zstemp12s +       zQ11Q12*(zstemp11s-zstemp22s) &
                      - zQ12Q12*zstemp12s
         zrotstemp22s = zQ12Q12*zstemp11s + 2._wp*zQ11Q12* zstemp12s &
                      + zQ11Q11*zstemp22s

         palphar =  zrotstemp11r*zdtemp11 + 2._wp*zrotstemp12r*zdtemp12 &
                  + zrotstemp22r*zdtemp22
         palphas =  zrotstemp11s*zdtemp11 + 2._wp*zrotstemp12s*zdtemp12 &
                  + zrotstemp22s*zdtemp22
      ENDIF
   END SUBROUTINE update_stress_rdg

!=======================================================================


   SUBROUTINE calc_ffrac( pstressp, pstressm, pstress12, pa11, pa12, &
      &                   pmresult11, pmresult12 )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE calc_ffrac  ***
      !!
      !! ** Purpose :   Computes term in evolution equation for structure tensor
      !!                which determines the ice floe re-orientation due to fracture
      !! ** Method :    Eq. 7: Ffrac = -kf(A-S) or = 0 depending on sigma_1 and sigma_2
      !!---------------------------------------------------------------------
      REAL (wp), INTENT(in)  ::   pstressp, pstressm, pstress12, pa11, pa12
      REAL (wp), INTENT(out) ::   pmresult11, pmresult12

      ! local variables

      REAL (wp) ::   zsigma11, zsigma12, zsigma22  ! stress tensor elements
      REAL (wp) ::   zAngle_denom        ! angle with principal component axis
      REAL (wp) ::   zsigma_1, zsigma_2   ! principal components of stress
      REAL (wp) ::   zQ11, zQ12, zQ11Q11, zQ11Q12, zQ12Q12

!!$   REAL (wp), PARAMETER ::   ppkfrac = 0.0001_wp   ! rate of fracture formation
      REAL (wp), PARAMETER ::   ppkfrac = 1.e-3_wp    ! rate of fracture formation
      REAL (wp), PARAMETER ::   ppthreshold = 0.3_wp  ! critical confinement ratio
      !!---------------------------------------------------------------
      !
      zsigma11 = 0.5_wp*(pstressp+pstressm)
      zsigma12 = pstress12
      zsigma22 = 0.5_wp*(pstressp-pstressm)

      ! Here's the change - no longer calculate gamma,
      ! use trig formulation, angle signs are kept correct, don't worry

      ! rotate tensor to get into sigma principal axis

      ! here Tan2gamma = 2 sig12 / (sig11 - sig12)
      ! add rsmall to the denominator to stop 1/0 errors, makes very little
      ! error to the calculated angles < rsmall

      zQ11Q11 = 0.1_wp
      zQ12Q12 = rsmall
      zQ11Q12 = rsmall

      IF((ABS( zsigma11 - zsigma22) > rsmall).OR.(ABS(zsigma12) > rsmall)) THEN

         zAngle_denom = 1.0_wp/sqrt( ( zsigma11 - zsigma22 )*( zsigma11 - &
                       zsigma22 ) + 4.0_wp*zsigma12*zsigma12)

         zQ11Q11 = 0.5_wp + ( zsigma11 - zsigma22 )*0.5_wp*zAngle_denom   !Cos^2
         zQ12Q12 = 0.5_wp - ( zsigma11 - zsigma22 )*0.5_wp*zAngle_denom   !Sin^2
         zQ11Q12 = zsigma12*zAngle_denom                      !CosSin
      ENDIF

      zsigma_1 = zQ11Q11*zsigma11 + 2.0_wp*zQ11Q12*zsigma12 + zQ12Q12*zsigma22 ! S(1,1)
      zsigma_2 = zQ12Q12*zsigma11 - 2.0_wp*zQ11Q12*zsigma12 + zQ11Q11*zsigma22 ! S(2,2)

      ! Pure divergence
      IF ((zsigma_1 >= 0.0_wp).AND.(zsigma_2 >= 0.0_wp))  THEN
         pmresult11 = 0.0_wp
         pmresult12 = 0.0_wp

      ! Unconfined compression: cracking of blocks not along the axial splitting
      ! direction
      ! which leads to the loss of their shape, so we again model it through diffusion
      ELSEIF ((zsigma_1 >= 0.0_wp).AND.(zsigma_2 < 0.0_wp))  THEN
         pmresult11 = - ppkfrac * (pa11 - zQ12Q12)
         pmresult12 = - ppkfrac * (pa12 + zQ11Q12)

      ! Shear faulting
      ELSEIF (zsigma_2 == 0.0_wp) THEN
         pmresult11 = 0.0_wp
         pmresult12 = 0.0_wp
      ELSEIF ((zsigma_1 <= 0.0_wp).AND.(zsigma_1/zsigma_2 <= ppthreshold)) THEN
         pmresult11 = - ppkfrac * (pa11 - zQ12Q12)
         pmresult12 = - ppkfrac * (pa12 + zQ11Q12)

      ! Horizontal spalling
      ELSE
         pmresult11 = 0.0_wp
         pmresult12 = 0.0_wp
      ENDIF

   END SUBROUTINE calc_ffrac


   SUBROUTINE rhg_eap_rst( cdrw, kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE rhg_eap_rst  ***
      !!
      !! ** Purpose :   Read or write RHG file in restart file
      !!
      !! ** Method  :   use of IOM library
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      INTEGER, OPTIONAL, INTENT(in) ::   kt     ! ice time-step
      !
      INTEGER  ::   iter                      ! local integer
      INTEGER  ::   id1, id2, id3, id4, id5   ! local integers
      INTEGER  ::   ix, iy, ip, iz, n, ia     ! local integers

      INTEGER, PARAMETER            ::    nz = 100

      REAL(wp) ::   ainit, xinit, yinit, pinit, zinit
      REAL(wp) ::   da, dx, dy, dp, dz, a1

      !!clem
      REAL(wp) ::   zw1, zw2, zfac, ztemp
      REAL(wp) ::   zidx, zidy, zidz
      REAL(wp) ::   zsaak(6)                  ! temporary array

      REAL(wp), PARAMETER           ::   eps6 = 1.0e-6_wp
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialize
         !                                   ! ---------------
         IF( ln_rstart ) THEN                   !* Read the restart file
            !
            id1 = iom_varid( numrir, 'stress1_i' , ldstop = .FALSE. )
            id2 = iom_varid( numrir, 'stress2_i' , ldstop = .FALSE. )
            id3 = iom_varid( numrir, 'stress12_i', ldstop = .FALSE. )
            id4 = iom_varid( numrir, 'aniso_11'  , ldstop = .FALSE. )
            id5 = iom_varid( numrir, 'aniso_12'  , ldstop = .FALSE. )
            !
            IF( MIN( id1, id2, id3, id4, id5 ) > 0 ) THEN      ! fields exist
               CALL iom_get( numrir, jpdom_auto, 'stress1_i' , stress1_i , cd_type = 'T' )
               CALL iom_get( numrir, jpdom_auto, 'stress2_i' , stress2_i , cd_type = 'T' )
               CALL iom_get( numrir, jpdom_auto, 'stress12_i', stress12_i, cd_type = 'F' )
               CALL iom_get( numrir, jpdom_auto, 'aniso_11'  , aniso_11  , cd_type = 'T' )
               CALL iom_get( numrir, jpdom_auto, 'aniso_12'  , aniso_12  , cd_type = 'T' )
            ELSE                                     ! start rheology from rest
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '   ==>>>   previous run without rheology, set stresses to 0'
               stress1_i (:,:) = 0._wp
               stress2_i (:,:) = 0._wp
               stress12_i(:,:) = 0._wp
               aniso_11  (:,:) = 0.5_wp
               aniso_12  (:,:) = 0._wp
            ENDIF
         ELSE                                   !* Start from rest
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '   ==>>>   start from rest: set stresses to 0'
            stress1_i (:,:) = 0._wp
            stress2_i (:,:) = 0._wp
            stress12_i(:,:) = 0._wp
            aniso_11  (:,:) = 0.5_wp
            aniso_12  (:,:) = 0._wp
         ENDIF
         !

         da = 0.5_wp/real(na_yield-1,kind=wp)
         ainit = 0.5_wp - da
         dx = rpi/real(nx_yield-1,kind=wp)
         xinit = rpi + 0.25_wp*rpi - dx
         dz = rpi/real(nz,kind=wp)
         zinit = -rpi*0.5_wp
         dy = rpi/real(ny_yield-1,kind=wp)
         yinit = -dy

         s11r(:,:,:) = 0._wp
         s12r(:,:,:) = 0._wp
         s22r(:,:,:) = 0._wp
         s11s(:,:,:) = 0._wp
         s12s(:,:,:) = 0._wp
         s22s(:,:,:) = 0._wp

!!$         DO ia=1,na_yield
!!$            DO ix=1,nx_yield
!!$               DO iy=1,ny_yield
!!$                  s11r(ix,iy,ia) = 0._wp
!!$                  s12r(ix,iy,ia) = 0._wp
!!$                  s22r(ix,iy,ia) = 0._wp
!!$                  s11s(ix,iy,ia) = 0._wp
!!$                  s12s(ix,iy,ia) = 0._wp
!!$                  s22s(ix,iy,ia) = 0._wp
!!$                  IF (ia <= na_yield-1) THEN
!!$                     DO iz=1,nz
!!$                        s11r(ix,iy,ia) = s11r(ix,iy,ia) + 1*w1(ainit+ia*da)* &
!!$                           exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
!!$                           s11kr(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz)*dz/sin(2._wp*pphi)
!!$                        s12r(ix,iy,ia) = s12r(ix,iy,ia) + 1*w1(ainit+ia*da)* &
!!$                           exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
!!$                           s12kr(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz)*dz/sin(2._wp*pphi)
!!$                        s22r(ix,iy,ia) = s22r(ix,iy,ia) + 1*w1(ainit+ia*da)* &
!!$                           exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
!!$                           s22kr(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz)*dz/sin(2._wp*pphi)
!!$                        s11s(ix,iy,ia) = s11s(ix,iy,ia) + 1*w1(ainit+ia*da)* &
!!$                           exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
!!$                           s11ks(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz)*dz/sin(2._wp*pphi)
!!$                        s12s(ix,iy,ia) = s12s(ix,iy,ia) + 1*w1(ainit+ia*da)* &
!!$                           exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
!!$                           s12ks(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz)*dz/sin(2._wp*pphi)
!!$                        s22s(ix,iy,ia) = s22s(ix,iy,ia) + 1*w1(ainit+ia*da)* &
!!$                           exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
!!$                           s22ks(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz)*dz/sin(2._wp*pphi)
!!$                     ENDDO
!!$                     IF (abs(s11r(ix,iy,ia)) < eps6) s11r(ix,iy,ia) = 0._wp
!!$                     IF (abs(s12r(ix,iy,ia)) < eps6) s12r(ix,iy,ia) = 0._wp
!!$                     IF (abs(s22r(ix,iy,ia)) < eps6) s22r(ix,iy,ia) = 0._wp
!!$                     IF (abs(s11s(ix,iy,ia)) < eps6) s11s(ix,iy,ia) = 0._wp
!!$                     IF (abs(s12s(ix,iy,ia)) < eps6) s12s(ix,iy,ia) = 0._wp
!!$                     IF (abs(s22s(ix,iy,ia)) < eps6) s22s(ix,iy,ia) = 0._wp
!!$                  ELSE
!!$                     s11r(ix,iy,ia) = 0.5_wp*s11kr(xinit+ix*dx,yinit+iy*dy,0._wp)/sin(2._wp*pphi)
!!$                     s12r(ix,iy,ia) = 0.5_wp*s12kr(xinit+ix*dx,yinit+iy*dy,0._wp)/sin(2._wp*pphi)
!!$                     s22r(ix,iy,ia) = 0.5_wp*s22kr(xinit+ix*dx,yinit+iy*dy,0._wp)/sin(2._wp*pphi)
!!$                     s11s(ix,iy,ia) = 0.5_wp*s11ks(xinit+ix*dx,yinit+iy*dy,0._wp)/sin(2._wp*pphi)
!!$                     s12s(ix,iy,ia) = 0.5_wp*s12ks(xinit+ix*dx,yinit+iy*dy,0._wp)/sin(2._wp*pphi)
!!$                     s22s(ix,iy,ia) = 0.5_wp*s22ks(xinit+ix*dx,yinit+iy*dy,0._wp)/sin(2._wp*pphi)
!!$                     IF (abs(s11r(ix,iy,ia)) < eps6) s11r(ix,iy,ia) = 0._wp
!!$                     IF (abs(s12r(ix,iy,ia)) < eps6) s12r(ix,iy,ia) = 0._wp
!!$                     IF (abs(s22r(ix,iy,ia)) < eps6) s22r(ix,iy,ia) = 0._wp
!!$                     IF (abs(s11s(ix,iy,ia)) < eps6) s11s(ix,iy,ia) = 0._wp
!!$                     IF (abs(s12s(ix,iy,ia)) < eps6) s12s(ix,iy,ia) = 0._wp
!!$                     IF (abs(s22s(ix,iy,ia)) < eps6) s22s(ix,iy,ia) = 0._wp
!!$                  ENDIF
!!$               ENDDO
!!$            ENDDO
!!$         ENDDO

         !! faster but still very slow => to be improved
         zfac = dz/sin(2._wp*pphi)
         DO ia = 1, na_yield-1
            zw1 = w1(ainit+ia*da)
            zw2 = w2(ainit+ia*da)
            DO iz = 1, nz
               zidz = zinit+iz*dz
               ztemp = zw1 * EXP(-zw2*(zinit+iz*dz)*(zinit+iz*dz))
               DO iy = 1, ny_yield
                  zidy = yinit+iy*dy
                  DO ix = 1, nx_yield
                     zidx = xinit+ix*dx
                     call all_skr_sks(zidx,zidy,zidz,zsaak)
                     zsaak = ztemp*zsaak*zfac
                     s11r(ix,iy,ia) = s11r(ix,iy,ia) + zsaak(1)
                     s12r(ix,iy,ia) = s12r(ix,iy,ia) + zsaak(2)
                     s22r(ix,iy,ia) = s22r(ix,iy,ia) + zsaak(3)
                     s11s(ix,iy,ia) = s11s(ix,iy,ia) + zsaak(4)
                     s12s(ix,iy,ia) = s12s(ix,iy,ia) + zsaak(5)
                     s22s(ix,iy,ia) = s22s(ix,iy,ia) + zsaak(6)
                  END DO
               END DO
            END DO
         END DO
         zfac = 1._wp/sin(2._wp*pphi)
         ia = na_yield
         DO iy = 1, ny_yield
            zidy = yinit+iy*dy
            DO ix = 1, nx_yield
               zidx = xinit+ix*dx
               call all_skr_sks(zidx,zidy,0._wp,zsaak)
               zsaak = 0.5_wp*zsaak*zfac
               s11r(ix,iy,ia) = zsaak(1)
               s12r(ix,iy,ia) = zsaak(2)
               s22r(ix,iy,ia) = zsaak(3)
               s11s(ix,iy,ia) = zsaak(4)
               s12s(ix,iy,ia) = zsaak(5)
               s22s(ix,iy,ia) = zsaak(6)
            ENDDO
         ENDDO
         WHERE (ABS(s11r(:,:,:)) < eps6) s11r(:,:,:) = 0._wp
         WHERE (ABS(s12r(:,:,:)) < eps6) s12r(:,:,:) = 0._wp
         WHERE (ABS(s22r(:,:,:)) < eps6) s22r(:,:,:) = 0._wp
         WHERE (ABS(s11s(:,:,:)) < eps6) s11s(:,:,:) = 0._wp
         WHERE (ABS(s12s(:,:,:)) < eps6) s12s(:,:,:) = 0._wp
         WHERE (ABS(s22s(:,:,:)) < eps6) s22s(:,:,:) = 0._wp


      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         IF(lwp) WRITE(numout,*) '---- rhg-rst ----'
         iter = kt + nn_fsbc - 1             ! ice restarts are written at kt == nitrst - nn_fsbc + 1
         !
         CALL iom_rstput( iter, nitrst, numriw, 'stress1_i' , stress1_i  )
         CALL iom_rstput( iter, nitrst, numriw, 'stress2_i' , stress2_i  )
         CALL iom_rstput( iter, nitrst, numriw, 'stress12_i', stress12_i )
         CALL iom_rstput( iter, nitrst, numriw, 'aniso_11'  , aniso_11   )
         CALL iom_rstput( iter, nitrst, numriw, 'aniso_12'  , aniso_12   )
         !
      ENDIF
      !
   END SUBROUTINE rhg_eap_rst


   FUNCTION w1(pa)
      !!-------------------------------------------------------------------
      !! Function : w1 (see Gaussian function psi in Tsamados et al 2013)
      !!-------------------------------------------------------------------
      REAL(wp), INTENT(in   ) ::   pa
      REAL(wp) ::   w1
      !!-------------------------------------------------------------------

      w1 = -   223.87569446_wp &
       &   +  2361.21986630_wp*pa &
       &   - 10606.56079975_wp*pa*pa &
       &   + 26315.50025642_wp*pa*pa*pa &
       &   - 38948.30444297_wp*pa*pa*pa*pa &
       &   + 34397.72407466_wp*pa*pa*pa*pa*pa &
       &   - 16789.98003081_wp*pa*pa*pa*pa*pa*pa &
       &   +  3495.82839237_wp*pa*pa*pa*pa*pa*pa*pa

   END FUNCTION w1


   FUNCTION w2(pa)
      !!-------------------------------------------------------------------
      !! Function : w2 (see Gaussian function psi in Tsamados et al 2013)
      !!-------------------------------------------------------------------
      REAL(wp), INTENT(in   ) ::   pa
      REAL(wp) ::   w2
      !!-------------------------------------------------------------------

      w2 = -    6670.68911883_wp &
       &   +   70222.33061536_wp*pa &
       &   -  314871.71525448_wp*pa*pa &
       &   +  779570.02793492_wp*pa*pa*pa &
       &   - 1151098.82436864_wp*pa*pa*pa*pa &
       &   + 1013896.59464498_wp*pa*pa*pa*pa*pa &
       &   -  493379.44906738_wp*pa*pa*pa*pa*pa*pa &
       &   +  102356.55151800_wp*pa*pa*pa*pa*pa*pa*pa

   END FUNCTION w2

   SUBROUTINE all_skr_sks( px, py, pz, allsk )
      REAL(wp), INTENT(in   ) ::   px,py,pz
      REAL(wp), INTENT(out  ) ::   allsk(6)

      REAL(wp) ::   zs12r0, zs21r0
      REAL(wp) ::   zs12s0, zs21s0

      REAL(wp) ::   zpih
      REAL(wp) ::   zn1t2i11, zn1t2i12, zn1t2i21, zn1t2i22
      REAL(wp) ::   zn2t1i11, zn2t1i12, zn2t1i21, zn2t1i22
      REAL(wp) ::   zt1t2i11, zt1t2i12, zt1t2i21, zt1t2i22
      REAL(wp) ::   zt2t1i11, zt2t1i12, zt2t1i21, zt2t1i22
      REAL(wp) ::   zd11, zd12, zd22
      REAL(wp) ::   zIIn1t2, zIIn2t1, zIIt1t2
      REAL(wp) ::   zHen1t2, zHen2t1
      !!-------------------------------------------------------------------

      zpih = 0.5_wp*rpi

      zn1t2i11 = cos(pz+zpih-pphi) * cos(pz+pphi)
      zn1t2i12 = cos(pz+zpih-pphi) * sin(pz+pphi)
      zn1t2i21 = sin(pz+zpih-pphi) * cos(pz+pphi)
      zn1t2i22 = sin(pz+zpih-pphi) * sin(pz+pphi)
      zn2t1i11 = cos(pz-zpih+pphi) * cos(pz-pphi)
      zn2t1i12 = cos(pz-zpih+pphi) * sin(pz-pphi)
      zn2t1i21 = sin(pz-zpih+pphi) * cos(pz-pphi)
      zn2t1i22 = sin(pz-zpih+pphi) * sin(pz-pphi)
      zt1t2i11 = cos(pz-pphi) * cos(pz+pphi)
      zt1t2i12 = cos(pz-pphi) * sin(pz+pphi)
      zt1t2i21 = sin(pz-pphi) * cos(pz+pphi)
      zt1t2i22 = sin(pz-pphi) * sin(pz+pphi)
      zt2t1i11 = cos(pz+pphi) * cos(pz-pphi)
      zt2t1i12 = cos(pz+pphi) * sin(pz-pphi)
      zt2t1i21 = sin(pz+pphi) * cos(pz-pphi)
      zt2t1i22 = sin(pz+pphi) * sin(pz-pphi)
   ! In expression of tensor d, with this formulatin d(x)=-d(x+pi)
   ! Solution, when diagonalizing always check sgn(a11-a22) if > then keep x else
   ! x=x-pi/2
      zd11 = cos(py)*cos(py)*(cos(px)+sin(px)*tan(py)*tan(py))
      zd12 = cos(py)*cos(py)*tan(py)*(-cos(px)+sin(px))
      zd22 = cos(py)*cos(py)*(sin(px)+cos(px)*tan(py)*tan(py))
      zIIn1t2 = zn1t2i11 * zd11 + (zn1t2i12 + zn1t2i21) * zd12 + zn1t2i22 * zd22
      zIIn2t1 = zn2t1i11 * zd11 + (zn2t1i12 + zn2t1i21) * zd12 + zn2t1i22 * zd22
      zIIt1t2 = zt1t2i11 * zd11 + (zt1t2i12 + zt1t2i21) * zd12 + zt1t2i22 * zd22

      IF (-zIIn1t2>=rsmall) THEN
      zHen1t2 = 1._wp
      ELSE
      zHen1t2 = 0._wp
      ENDIF

      IF (-zIIn2t1>=rsmall) THEN
      zHen2t1 = 1._wp
      ELSE
      zHen2t1 = 0._wp
      ENDIF

      !!-------------------------------------------------------------------
      !! Function : s11kr
      !!-------------------------------------------------------------------
      allsk(1) = (- zHen1t2 * zn1t2i11 - zHen2t1 * zn2t1i11)
      !!-------------------------------------------------------------------
      !! Function : s12kr
      !!-------------------------------------------------------------------
      zs12r0 = (- zHen1t2 * zn1t2i12 - zHen2t1 * zn2t1i12)
      zs21r0 = (- zHen1t2 * zn1t2i21 - zHen2t1 * zn2t1i21)
      allsk(2)=0.5_wp*(zs12r0+zs21r0)
      !!-------------------------------------------------------------------
      !! Function : s22kr
      !!-------------------------------------------------------------------
      allsk(3) = (- zHen1t2 * zn1t2i22 - zHen2t1 * zn2t1i22)
      !!-------------------------------------------------------------------
      !! Function : s11ks
      !!-------------------------------------------------------------------

      allsk(4) = sign(1._wp,zIIt1t2+rsmall)*(zHen1t2 * zt1t2i11 + zHen2t1 * zt2t1i11)
      !!-------------------------------------------------------------------
      !! Function : s12ks
      !!-------------------------------------------------------------------
      zs12s0 = sign(1._wp,zIIt1t2+rsmall)*(zHen1t2 * zt1t2i12 + zHen2t1 * zt2t1i12)
      zs21s0 = sign(1._wp,zIIt1t2+rsmall)*(zHen1t2 * zt1t2i21 + zHen2t1 * zt2t1i21)
      allsk(5)=0.5_wp*(zs12s0+zs21s0)
      !!-------------------------------------------------------------------
      !! Function : s22ks
      !!-------------------------------------------------------------------
      allsk(6) = sign(1._wp,zIIt1t2+rsmall)*(zHen1t2 * zt1t2i22 + zHen2t1 * zt2t1i22)
   END SUBROUTINE all_skr_sks

#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
   USE par_kind
   USE lib_mpp
   CONTAINS
   SUBROUTINE ice_dyn_rhg_eap( kt, Kmm, pstress1_i, pstress2_i, pstress12_i, pshear_i, pdivu_i, pdelta_i, paniso_11, paniso_12, prdg_conv )
      INTEGER                 , INTENT(in   ) ::   kt                                    ! time step
      INTEGER                 , INTENT(in   ) ::   Kmm                                   ! ocean time level index
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   pstress1_i, pstress2_i, pstress12_i   !
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   pshear_i  , pdivu_i   , pdelta_i      !
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   paniso_11 , paniso_12                 ! structure tensor components
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   prdg_conv                             ! for ridging
      CALL ctl_stop('EAP rheology is currently dsabled due to issues with AGRIF preprocessing')
   END SUBROUTINE ice_dyn_rhg_eap
   SUBROUTINE rhg_eap_rst( cdrw, kt )
      CHARACTER(len=*) , INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      INTEGER, OPTIONAL, INTENT(in) ::   kt     ! ice time-step
   END SUBROUTINE rhg_eap_rst
#endif

   !!==============================================================================
END MODULE icedyn_rhg_eap
