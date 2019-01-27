MODULE domvvl
   !!======================================================================
   !!                       ***  MODULE domvvl   ***
   !! Ocean : 
   !!======================================================================
   !! History :  2.0  !  2006-06  (B. Levier, L. Marie)  original code
   !!            3.1  !  2009-02  (G. Madec, M. Leclair, R. Benshila)  pure z* coordinate
   !!            3.3  !  2011-10  (M. Leclair) totally rewrote domvvl: vvl option includes z_star and z_tilde coordinates
   !!            3.6  !  2014-11  (P. Mathiot) add ice shelf capability
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_vvl_init     : define initial vertical scale factors, depths and column thickness
   !!   dom_vvl_sf_nxt   : Compute next vertical scale factors
   !!   dom_vvl_sf_swp   : Swap vertical scale factors and update the vertical grid
   !!   dom_vvl_interpol : Interpolate vertical scale factors from one grid point to another
   !!   dom_vvl_rst      : read/write restart file
   !!   dom_vvl_ctl      : Check the vvl options
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE phycst          ! physical constant
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! ocean surface boundary condition
   USE wet_dry         ! wetting and drying
   USE usrdef_istate   ! user defined initial state (wad only)
   USE restart         ! ocean restart
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O manager library
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC  dom_vvl_init       ! called by domain.F90
   PUBLIC  dom_vvl_sf_nxt     ! called by step.F90
   PUBLIC  dom_vvl_sf_swp     ! called by step.F90
   PUBLIC  dom_vvl_interpol   ! called by dynnxt.F90

   !                                                      !!* Namelist nam_vvl
   LOGICAL , PUBLIC :: ln_vvl_zstar           = .FALSE.    ! zstar  vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_ztilde          = .FALSE.    ! ztilde vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_layer           = .FALSE.    ! level  vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_ztilde_as_zstar = .FALSE.    ! ztilde vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_zstar_at_eqtor  = .FALSE.    ! ztilde vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_kepe            = .FALSE.    ! kinetic/potential energy transfer
   !                                                       ! conservation: not used yet
   REAL(wp)         :: rn_ahe3                             ! thickness diffusion coefficient
   REAL(wp)         :: rn_rst_e3t                          ! ztilde to zstar restoration timescale [days]
   REAL(wp)         :: rn_lf_cutoff                        ! cutoff frequency for low-pass filter  [days]
   REAL(wp)         :: rn_zdef_max                         ! maximum fractional e3t deformation
   LOGICAL , PUBLIC :: ln_vvl_dbg = .FALSE.                ! debug control prints

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: un_td, vn_td                ! thickness diffusion transport
   REAL(wp)        , ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: hdiv_lf                     ! low frequency part of hz divergence
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tilde_e3t_b, tilde_e3t_n    ! baroclinic scale factors
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tilde_e3t_a, dtilde_e3t_a   ! baroclinic scale factors
   REAL(wp)        , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: frq_rst_e3t                 ! retoring period for scale factors
   REAL(wp)        , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: frq_rst_hdv                 ! retoring period for low freq. divergence

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: domvvl.F90 10572 2019-01-24 15:37:13Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dom_vvl_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION dom_vvl_alloc  ***
      !!----------------------------------------------------------------------
      IF( ln_vvl_zstar )   dom_vvl_alloc = 0
      IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN
         ALLOCATE( tilde_e3t_b(jpi,jpj,jpk)  , tilde_e3t_n(jpi,jpj,jpk) , tilde_e3t_a(jpi,jpj,jpk) ,   &
            &      dtilde_e3t_a(jpi,jpj,jpk) , un_td  (jpi,jpj,jpk)     , vn_td  (jpi,jpj,jpk)     ,   &
            &      STAT = dom_vvl_alloc        )
         CALL mpp_sum ( 'domvvl', dom_vvl_alloc )
         IF( dom_vvl_alloc /= 0 )   CALL ctl_stop( 'STOP', 'dom_vvl_alloc: failed to allocate arrays' )
         un_td = 0._wp
         vn_td = 0._wp
      ENDIF
      IF( ln_vvl_ztilde ) THEN
         ALLOCATE( frq_rst_e3t(jpi,jpj) , frq_rst_hdv(jpi,jpj) , hdiv_lf(jpi,jpj,jpk) , STAT= dom_vvl_alloc )
         CALL mpp_sum ( 'domvvl', dom_vvl_alloc )
         IF( dom_vvl_alloc /= 0 )   CALL ctl_stop( 'STOP', 'dom_vvl_alloc: failed to allocate arrays' )
      ENDIF
      !
   END FUNCTION dom_vvl_alloc


   SUBROUTINE dom_vvl_init
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_vvl_init  ***
      !!                   
      !! ** Purpose :  Initialization of all scale factors, depths
      !!               and water column heights
      !!
      !! ** Method  :  - use restart file and/or initialize
      !!               - interpolate scale factors
      !!
      !! ** Action  : - e3t_(n/b) and tilde_e3t_(n/b)
      !!              - Regrid: e3(u/v)_n
      !!                        e3(u/v)_b       
      !!                        e3w_n           
      !!                        e3(u/v)w_b      
      !!                        e3(u/v)w_n      
      !!                        gdept_n, gdepw_n and gde3w_n
      !!              - h(t/u/v)_0
      !!              - frq_rst_e3t and frq_rst_hdv
      !!
      !! Reference  : Leclair, M., and G. Madec, 2011, Ocean Modelling.
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk
      INTEGER ::   ii0, ii1, ij0, ij1
      REAL(wp)::   zcoef
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dom_vvl_init : Variable volume activated'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      !
      CALL dom_vvl_ctl     ! choose vertical coordinate (z_star, z_tilde or layer)
      !
      !                    ! Allocate module arrays
      IF( dom_vvl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'dom_vvl_init : unable to allocate arrays' )
      !
      !                    ! Read or initialize e3t_(b/n), tilde_e3t_(b/n) and hdiv_lf
      CALL dom_vvl_rst( nit000, 'READ' )
      e3t_a(:,:,jpk) = e3t_0(:,:,jpk)  ! last level always inside the sea floor set one for all
      !
      !                    !== Set of all other vertical scale factors  ==!  (now and before)
      !                                ! Horizontal interpolation of e3t
      CALL dom_vvl_interpol( e3t_b(:,:,:), e3u_b(:,:,:), 'U' )    ! from T to U
      CALL dom_vvl_interpol( e3t_n(:,:,:), e3u_n(:,:,:), 'U' )
      CALL dom_vvl_interpol( e3t_b(:,:,:), e3v_b(:,:,:), 'V' )    ! from T to V 
      CALL dom_vvl_interpol( e3t_n(:,:,:), e3v_n(:,:,:), 'V' )
      CALL dom_vvl_interpol( e3u_n(:,:,:), e3f_n(:,:,:), 'F' )    ! from U to F
      !                                ! Vertical interpolation of e3t,u,v 
      CALL dom_vvl_interpol( e3t_n(:,:,:), e3w_n (:,:,:), 'W'  )  ! from T to W
      CALL dom_vvl_interpol( e3t_b(:,:,:), e3w_b (:,:,:), 'W'  )
      CALL dom_vvl_interpol( e3u_n(:,:,:), e3uw_n(:,:,:), 'UW' )  ! from U to UW
      CALL dom_vvl_interpol( e3u_b(:,:,:), e3uw_b(:,:,:), 'UW' )
      CALL dom_vvl_interpol( e3v_n(:,:,:), e3vw_n(:,:,:), 'VW' )  ! from V to UW
      CALL dom_vvl_interpol( e3v_b(:,:,:), e3vw_b(:,:,:), 'VW' )

      ! We need to define e3[tuv]_a for AGRIF initialisation (should not be a problem for the restartability...)
      e3t_a(:,:,:) = e3t_n(:,:,:)
      e3u_a(:,:,:) = e3u_n(:,:,:)
      e3v_a(:,:,:) = e3v_n(:,:,:)
      !
      !                    !==  depth of t and w-point  ==!   (set the isf depth as it is in the initial timestep)
      gdept_n(:,:,1) = 0.5_wp * e3w_n(:,:,1)       ! reference to the ocean surface (used for MLD and light penetration)
      gdepw_n(:,:,1) = 0.0_wp
      gde3w_n(:,:,1) = gdept_n(:,:,1) - sshn(:,:)  ! reference to a common level z=0 for hpg
      gdept_b(:,:,1) = 0.5_wp * e3w_b(:,:,1)
      gdepw_b(:,:,1) = 0.0_wp
      DO jk = 2, jpk                               ! vertical sum
         DO jj = 1,jpj
            DO ji = 1,jpi
               !    zcoef = tmask - wmask    ! 0 everywhere tmask = wmask, ie everywhere expect at jk = mikt
               !                             ! 1 everywhere from mbkt to mikt + 1 or 1 (if no isf)
               !                             ! 0.5 where jk = mikt     
!!gm ???????   BUG ?  gdept_n as well as gde3w_n  does not include the thickness of ISF ??
               zcoef = ( tmask(ji,jj,jk) - wmask(ji,jj,jk) )
               gdepw_n(ji,jj,jk) = gdepw_n(ji,jj,jk-1) + e3t_n(ji,jj,jk-1)
               gdept_n(ji,jj,jk) =      zcoef  * ( gdepw_n(ji,jj,jk  ) + 0.5 * e3w_n(ji,jj,jk))  &
                  &                + (1-zcoef) * ( gdept_n(ji,jj,jk-1) +       e3w_n(ji,jj,jk)) 
               gde3w_n(ji,jj,jk) = gdept_n(ji,jj,jk) - sshn(ji,jj)
               gdepw_b(ji,jj,jk) = gdepw_b(ji,jj,jk-1) + e3t_b(ji,jj,jk-1)
               gdept_b(ji,jj,jk) =      zcoef  * ( gdepw_b(ji,jj,jk  ) + 0.5 * e3w_b(ji,jj,jk))  &
                  &                + (1-zcoef) * ( gdept_b(ji,jj,jk-1) +       e3w_b(ji,jj,jk)) 
            END DO
         END DO
      END DO
      !
      !                    !==  thickness of the water column  !!   (ocean portion only)
      ht_n(:,:) = e3t_n(:,:,1) * tmask(:,:,1)   !!gm  BUG  :  this should be 1/2 * e3w(k=1) ....
      hu_b(:,:) = e3u_b(:,:,1) * umask(:,:,1)
      hu_n(:,:) = e3u_n(:,:,1) * umask(:,:,1)
      hv_b(:,:) = e3v_b(:,:,1) * vmask(:,:,1)
      hv_n(:,:) = e3v_n(:,:,1) * vmask(:,:,1)
      DO jk = 2, jpkm1
         ht_n(:,:) = ht_n(:,:) + e3t_n(:,:,jk) * tmask(:,:,jk)
         hu_b(:,:) = hu_b(:,:) + e3u_b(:,:,jk) * umask(:,:,jk)
         hu_n(:,:) = hu_n(:,:) + e3u_n(:,:,jk) * umask(:,:,jk)
         hv_b(:,:) = hv_b(:,:) + e3v_b(:,:,jk) * vmask(:,:,jk)
         hv_n(:,:) = hv_n(:,:) + e3v_n(:,:,jk) * vmask(:,:,jk)
      END DO
      !
      !                    !==  inverse of water column thickness   ==!   (u- and v- points)
      r1_hu_b(:,:) = ssumask(:,:) / ( hu_b(:,:) + 1._wp - ssumask(:,:) )    ! _i mask due to ISF
      r1_hu_n(:,:) = ssumask(:,:) / ( hu_n(:,:) + 1._wp - ssumask(:,:) )
      r1_hv_b(:,:) = ssvmask(:,:) / ( hv_b(:,:) + 1._wp - ssvmask(:,:) )
      r1_hv_n(:,:) = ssvmask(:,:) / ( hv_n(:,:) + 1._wp - ssvmask(:,:) )

      !                    !==   z_tilde coordinate case  ==!   (Restoring frequencies)
      IF( ln_vvl_ztilde ) THEN
!!gm : idea: add here a READ in a file of custumized restoring frequency
         !                                   ! Values in days provided via the namelist
         !                                   ! use rsmall to avoid possible division by zero errors with faulty settings
         frq_rst_e3t(:,:) = 2._wp * rpi / ( MAX( rn_rst_e3t  , rsmall ) * 86400.0_wp )
         frq_rst_hdv(:,:) = 2._wp * rpi / ( MAX( rn_lf_cutoff, rsmall ) * 86400.0_wp )
         !
         IF( ln_vvl_ztilde_as_zstar ) THEN   ! z-star emulation using z-tile
            frq_rst_e3t(:,:) = 0._wp               !Ignore namelist settings
            frq_rst_hdv(:,:) = 1._wp / rdt
         ENDIF
         IF ( ln_vvl_zstar_at_eqtor ) THEN   ! use z-star in vicinity of the Equator
            DO jj = 1, jpj
               DO ji = 1, jpi
!!gm  case |gphi| >= 6 degrees is useless   initialized just above by default
                  IF( ABS(gphit(ji,jj)) >= 6.) THEN
                     ! values outside the equatorial band and transition zone (ztilde)
                     frq_rst_e3t(ji,jj) =  2.0_wp * rpi / ( MAX( rn_rst_e3t  , rsmall ) * 86400.e0_wp )
                     frq_rst_hdv(ji,jj) =  2.0_wp * rpi / ( MAX( rn_lf_cutoff, rsmall ) * 86400.e0_wp )
                  ELSEIF( ABS(gphit(ji,jj)) <= 2.5) THEN    ! Equator strip ==> z-star
                     ! values inside the equatorial band (ztilde as zstar)
                     frq_rst_e3t(ji,jj) =  0.0_wp
                     frq_rst_hdv(ji,jj) =  1.0_wp / rdt
                  ELSE                                      ! transition band (2.5 to 6 degrees N/S)
                     !                                      ! (linearly transition from z-tilde to z-star)
                     frq_rst_e3t(ji,jj) = 0.0_wp + (frq_rst_e3t(ji,jj)-0.0_wp)*0.5_wp   &
                        &            * (  1.0_wp - COS( rad*(ABS(gphit(ji,jj))-2.5_wp)  &
                        &                                          * 180._wp / 3.5_wp ) )
                     frq_rst_hdv(ji,jj) = (1.0_wp / rdt)                                &
                        &            + (  frq_rst_hdv(ji,jj)-(1.e0_wp / rdt) )*0.5_wp   &
                        &            * (  1._wp  - COS( rad*(ABS(gphit(ji,jj))-2.5_wp)  &
                        &                                          * 180._wp / 3.5_wp ) )
                  ENDIF
               END DO
            END DO
            IF( cn_cfg == "orca" .OR. cn_cfg == "ORCA" ) THEN
               IF( nn_cfg == 3 ) THEN   ! ORCA2: Suppress ztilde in the Foxe Basin for ORCA2
                  ii0 = 103   ;   ii1 = 111       
                  ij0 = 128   ;   ij1 = 135   ;   
                  frq_rst_e3t( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) =  0.0_wp
                  frq_rst_hdv( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) =  1.e0_wp / rdt
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      !
      IF(lwxios) THEN
! define variables in restart file when writing with XIOS
         CALL iom_set_rstw_var_active('e3t_b')
         CALL iom_set_rstw_var_active('e3t_n')
         !                                           ! ----------------------- !
         IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN  ! z_tilde and layer cases !
            !                                        ! ----------------------- !
            CALL iom_set_rstw_var_active('tilde_e3t_b')
            CALL iom_set_rstw_var_active('tilde_e3t_n')
         END IF
         !                                           ! -------------!    
         IF( ln_vvl_ztilde ) THEN                    ! z_tilde case !
            !                                        ! ------------ !
            CALL iom_set_rstw_var_active('hdiv_lf')
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE dom_vvl_init


   SUBROUTINE dom_vvl_sf_nxt( kt, kcall ) 
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_vvl_sf_nxt  ***
      !!                   
      !! ** Purpose :  - compute the after scale factors used in tra_zdf, dynnxt,
      !!                 tranxt and dynspg routines
      !!
      !! ** Method  :  - z_star case:  Repartition of ssh INCREMENT proportionnaly to the level thickness.
      !!               - z_tilde_case: after scale factor increment = 
      !!                                    high frequency part of horizontal divergence
      !!                                  + retsoring towards the background grid
      !!                                  + thickness difusion
      !!                               Then repartition of ssh INCREMENT proportionnaly
      !!                               to the "baroclinic" level thickness.
      !!
      !! ** Action  :  - hdiv_lf    : restoring towards full baroclinic divergence in z_tilde case
      !!               - tilde_e3t_a: after increment of vertical scale factor 
      !!                              in z_tilde case
      !!               - e3(t/u/v)_a
      !!
      !! Reference  : Leclair, M., and Madec, G. 2011, Ocean Modelling.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )           ::   kt      ! time step
      INTEGER, INTENT( in ), OPTIONAL ::   kcall   ! optional argument indicating call sequence
      !
      INTEGER                ::   ji, jj, jk            ! dummy loop indices
      INTEGER , DIMENSION(3) ::   ijk_max, ijk_min      ! temporary integers
      REAL(wp)               ::   z2dt, z_tmin, z_tmax  ! local scalars
      LOGICAL                ::   ll_do_bclinic         ! local logical
      REAL(wp), DIMENSION(jpi,jpj)     ::   zht, z_scale, zwu, zwv, zhdiv
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ze3t
      !!----------------------------------------------------------------------
      !
      IF( ln_linssh )   RETURN      ! No calculation in linear free surface
      !
      IF( ln_timing )   CALL timing_start('dom_vvl_sf_nxt')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dom_vvl_sf_nxt : compute after scale factors'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
      ENDIF

      ll_do_bclinic = .TRUE.
      IF( PRESENT(kcall) ) THEN
         IF( kcall == 2 .AND. ln_vvl_ztilde )   ll_do_bclinic = .FALSE.
      ENDIF

      ! ******************************* !
      ! After acale factors at t-points !
      ! ******************************* !
      !                                                ! --------------------------------------------- !
      !                                                ! z_star coordinate and barotropic z-tilde part !
      !                                                ! --------------------------------------------- !
      !
      z_scale(:,:) = ( ssha(:,:) - sshb(:,:) ) * ssmask(:,:) / ( ht_0(:,:) + sshn(:,:) + 1. - ssmask(:,:) )
      DO jk = 1, jpkm1
         ! formally this is the same as e3t_a = e3t_0*(1+ssha/ht_0)
         e3t_a(:,:,jk) = e3t_b(:,:,jk) + e3t_n(:,:,jk) * z_scale(:,:) * tmask(:,:,jk)
      END DO
      !
      IF( ln_vvl_ztilde .OR. ln_vvl_layer .AND. ll_do_bclinic ) THEN   ! z_tilde or layer coordinate !
         !                                                            ! ------baroclinic part------ !
         ! I - initialization
         ! ==================

         ! 1 - barotropic divergence
         ! -------------------------
         zhdiv(:,:) = 0._wp
         zht(:,:)   = 0._wp
         DO jk = 1, jpkm1
            zhdiv(:,:) = zhdiv(:,:) + e3t_n(:,:,jk) * hdivn(:,:,jk)
            zht  (:,:) = zht  (:,:) + e3t_n(:,:,jk) * tmask(:,:,jk)
         END DO
         zhdiv(:,:) = zhdiv(:,:) / ( zht(:,:) + 1. - tmask_i(:,:) )

         ! 2 - Low frequency baroclinic horizontal divergence  (z-tilde case only)
         ! --------------------------------------------------
         IF( ln_vvl_ztilde ) THEN
            IF( kt > nit000 ) THEN
               DO jk = 1, jpkm1
                  hdiv_lf(:,:,jk) = hdiv_lf(:,:,jk) - rdt * frq_rst_hdv(:,:)   &
                     &          * ( hdiv_lf(:,:,jk) - e3t_n(:,:,jk) * ( hdivn(:,:,jk) - zhdiv(:,:) ) )
               END DO
            ENDIF
         ENDIF

         ! II - after z_tilde increments of vertical scale factors
         ! =======================================================
         tilde_e3t_a(:,:,:) = 0._wp  ! tilde_e3t_a used to store tendency terms

         ! 1 - High frequency divergence term
         ! ----------------------------------
         IF( ln_vvl_ztilde ) THEN     ! z_tilde case
            DO jk = 1, jpkm1
               tilde_e3t_a(:,:,jk) = tilde_e3t_a(:,:,jk) - ( e3t_n(:,:,jk) * ( hdivn(:,:,jk) - zhdiv(:,:) ) - hdiv_lf(:,:,jk) )
            END DO
         ELSE                         ! layer case
            DO jk = 1, jpkm1
               tilde_e3t_a(:,:,jk) = tilde_e3t_a(:,:,jk) -   e3t_n(:,:,jk) * ( hdivn(:,:,jk) - zhdiv(:,:) ) * tmask(:,:,jk)
            END DO
         ENDIF

         ! 2 - Restoring term (z-tilde case only)
         ! ------------------
         IF( ln_vvl_ztilde ) THEN
            DO jk = 1, jpk
               tilde_e3t_a(:,:,jk) = tilde_e3t_a(:,:,jk) - frq_rst_e3t(:,:) * tilde_e3t_b(:,:,jk)
            END DO
         ENDIF

         ! 3 - Thickness diffusion term
         ! ----------------------------
         zwu(:,:) = 0._wp
         zwv(:,:) = 0._wp
         DO jk = 1, jpkm1        ! a - first derivative: diffusive fluxes
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  un_td(ji,jj,jk) = rn_ahe3 * umask(ji,jj,jk) * e2_e1u(ji,jj)           &
                     &            * ( tilde_e3t_b(ji,jj,jk) - tilde_e3t_b(ji+1,jj  ,jk) )
                  vn_td(ji,jj,jk) = rn_ahe3 * vmask(ji,jj,jk) * e1_e2v(ji,jj)           & 
                     &            * ( tilde_e3t_b(ji,jj,jk) - tilde_e3t_b(ji  ,jj+1,jk) )
                  zwu(ji,jj) = zwu(ji,jj) + un_td(ji,jj,jk)
                  zwv(ji,jj) = zwv(ji,jj) + vn_td(ji,jj,jk)
               END DO
            END DO
         END DO
         DO jj = 1, jpj          ! b - correction for last oceanic u-v points
            DO ji = 1, jpi
               un_td(ji,jj,mbku(ji,jj)) = un_td(ji,jj,mbku(ji,jj)) - zwu(ji,jj)
               vn_td(ji,jj,mbkv(ji,jj)) = vn_td(ji,jj,mbkv(ji,jj)) - zwv(ji,jj)
            END DO
         END DO
         DO jk = 1, jpkm1        ! c - second derivative: divergence of diffusive fluxes
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  tilde_e3t_a(ji,jj,jk) = tilde_e3t_a(ji,jj,jk) + (   un_td(ji-1,jj  ,jk) - un_td(ji,jj,jk)    &
                     &                                          +     vn_td(ji  ,jj-1,jk) - vn_td(ji,jj,jk)    &
                     &                                            ) * r1_e1e2t(ji,jj)
               END DO
            END DO
         END DO
         !                       ! d - thickness diffusion transport: boundary conditions
         !                             (stored for tracer advction and continuity equation)
         CALL lbc_lnk_multi( 'domvvl', un_td , 'U' , -1._wp, vn_td , 'V' , -1._wp)

         ! 4 - Time stepping of baroclinic scale factors
         ! ---------------------------------------------
         ! Leapfrog time stepping
         ! ~~~~~~~~~~~~~~~~~~~~~~
         IF( neuler == 0 .AND. kt == nit000 ) THEN
            z2dt =  rdt
         ELSE
            z2dt = 2.0_wp * rdt
         ENDIF
         CALL lbc_lnk( 'domvvl', tilde_e3t_a(:,:,:), 'T', 1._wp )
         tilde_e3t_a(:,:,:) = tilde_e3t_b(:,:,:) + z2dt * tmask(:,:,:) * tilde_e3t_a(:,:,:)

         ! Maximum deformation control
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ze3t(:,:,jpk) = 0._wp
         DO jk = 1, jpkm1
            ze3t(:,:,jk) = tilde_e3t_a(:,:,jk) / e3t_0(:,:,jk) * tmask(:,:,jk) * tmask_i(:,:)
         END DO
         z_tmax = MAXVAL( ze3t(:,:,:) )
         CALL mpp_max( 'domvvl', z_tmax )                 ! max over the global domain
         z_tmin = MINVAL( ze3t(:,:,:) )
         CALL mpp_min( 'domvvl', z_tmin )                 ! min over the global domain
         ! - ML - test: for the moment, stop simulation for too large e3_t variations
         IF( ( z_tmax >  rn_zdef_max ) .OR. ( z_tmin < - rn_zdef_max ) ) THEN
            IF( lk_mpp ) THEN
               CALL mpp_maxloc( 'domvvl', ze3t, tmask, z_tmax, ijk_max )
               CALL mpp_minloc( 'domvvl', ze3t, tmask, z_tmin, ijk_min )
            ELSE
               ijk_max = MAXLOC( ze3t(:,:,:) )
               ijk_max(1) = ijk_max(1) + nimpp - 1
               ijk_max(2) = ijk_max(2) + njmpp - 1
               ijk_min = MINLOC( ze3t(:,:,:) )
               ijk_min(1) = ijk_min(1) + nimpp - 1
               ijk_min(2) = ijk_min(2) + njmpp - 1
            ENDIF
            IF (lwp) THEN
               WRITE(numout, *) 'MAX( tilde_e3t_a(:,:,:) / e3t_0(:,:,:) ) =', z_tmax
               WRITE(numout, *) 'at i, j, k=', ijk_max
               WRITE(numout, *) 'MIN( tilde_e3t_a(:,:,:) / e3t_0(:,:,:) ) =', z_tmin
               WRITE(numout, *) 'at i, j, k=', ijk_min            
               CALL ctl_stop( 'STOP', 'MAX( ABS( tilde_e3t_a(:,:,: ) ) / e3t_0(:,:,:) ) too high')
            ENDIF
         ENDIF
         ! - ML - end test
         ! - ML - Imposing these limits will cause a baroclinicity error which is corrected for below
         tilde_e3t_a(:,:,:) = MIN( tilde_e3t_a(:,:,:),   rn_zdef_max * e3t_0(:,:,:) )
         tilde_e3t_a(:,:,:) = MAX( tilde_e3t_a(:,:,:), - rn_zdef_max * e3t_0(:,:,:) )

         !
         ! "tilda" change in the after scale factor
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         DO jk = 1, jpkm1
            dtilde_e3t_a(:,:,jk) = tilde_e3t_a(:,:,jk) - tilde_e3t_b(:,:,jk)
         END DO
         ! III - Barotropic repartition of the sea surface height over the baroclinic profile
         ! ==================================================================================
         ! add ( ssh increment + "baroclinicity error" ) proportionly to e3t(n)
         ! - ML - baroclinicity error should be better treated in the future
         !        i.e. locally and not spread over the water column.
         !        (keep in mind that the idea is to reduce Eulerian velocity as much as possible)
         zht(:,:) = 0.
         DO jk = 1, jpkm1
            zht(:,:)  = zht(:,:) + tilde_e3t_a(:,:,jk) * tmask(:,:,jk)
         END DO
         z_scale(:,:) =  - zht(:,:) / ( ht_0(:,:) + sshn(:,:) + 1. - ssmask(:,:) )
         DO jk = 1, jpkm1
            dtilde_e3t_a(:,:,jk) = dtilde_e3t_a(:,:,jk) + e3t_n(:,:,jk) * z_scale(:,:) * tmask(:,:,jk)
         END DO

      ENDIF

      IF( ln_vvl_ztilde .OR. ln_vvl_layer )  THEN   ! z_tilde or layer coordinate !
      !                                           ! ---baroclinic part--------- !
         DO jk = 1, jpkm1
            e3t_a(:,:,jk) = e3t_a(:,:,jk) + dtilde_e3t_a(:,:,jk) * tmask(:,:,jk)
         END DO
      ENDIF

      IF( ln_vvl_dbg .AND. .NOT. ll_do_bclinic ) THEN   ! - ML - test: control prints for debuging
         !
         IF( lwp ) WRITE(numout, *) 'kt =', kt
         IF ( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN
            z_tmax = MAXVAL( tmask(:,:,1) * tmask_i(:,:) * ABS( zht(:,:) ) )
            CALL mpp_max( 'domvvl', z_tmax )                             ! max over the global domain
            IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(SUM(tilde_e3t_a))) =', z_tmax
         END IF
         !
         zht(:,:) = 0.0_wp
         DO jk = 1, jpkm1
            zht(:,:) = zht(:,:) + e3t_n(:,:,jk) * tmask(:,:,jk)
         END DO
         z_tmax = MAXVAL( tmask(:,:,1) * tmask_i(:,:) * ABS( ht_0(:,:) + sshn(:,:) - zht(:,:) ) )
         CALL mpp_max( 'domvvl', z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(ht_0+sshn-SUM(e3t_n))) =', z_tmax
         !
         zht(:,:) = 0.0_wp
         DO jk = 1, jpkm1
            zht(:,:) = zht(:,:) + e3t_a(:,:,jk) * tmask(:,:,jk)
         END DO
         z_tmax = MAXVAL( tmask(:,:,1) * tmask_i(:,:) * ABS( ht_0(:,:) + ssha(:,:) - zht(:,:) ) )
         CALL mpp_max( 'domvvl', z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(ht_0+ssha-SUM(e3t_a))) =', z_tmax
         !
         zht(:,:) = 0.0_wp
         DO jk = 1, jpkm1
            zht(:,:) = zht(:,:) + e3t_b(:,:,jk) * tmask(:,:,jk)
         END DO
         z_tmax = MAXVAL( tmask(:,:,1) * tmask_i(:,:) * ABS( ht_0(:,:) + sshb(:,:) - zht(:,:) ) )
         CALL mpp_max( 'domvvl', z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(ht_0+sshb-SUM(e3t_b))) =', z_tmax
         !
         z_tmax = MAXVAL( tmask(:,:,1) *  ABS( sshb(:,:) ) )
         CALL mpp_max( 'domvvl', z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(sshb))) =', z_tmax
         !
         z_tmax = MAXVAL( tmask(:,:,1) *  ABS( sshn(:,:) ) )
         CALL mpp_max( 'domvvl', z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(sshn))) =', z_tmax
         !
         z_tmax = MAXVAL( tmask(:,:,1) *  ABS( ssha(:,:) ) )
         CALL mpp_max( 'domvvl', z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(ssha))) =', z_tmax
      END IF

      ! *********************************** !
      ! After scale factors at u- v- points !
      ! *********************************** !

      CALL dom_vvl_interpol( e3t_a(:,:,:), e3u_a(:,:,:), 'U' )
      CALL dom_vvl_interpol( e3t_a(:,:,:), e3v_a(:,:,:), 'V' )

      ! *********************************** !
      ! After depths at u- v points         !
      ! *********************************** !

      hu_a(:,:) = e3u_a(:,:,1) * umask(:,:,1)
      hv_a(:,:) = e3v_a(:,:,1) * vmask(:,:,1)
      DO jk = 2, jpkm1
         hu_a(:,:) = hu_a(:,:) + e3u_a(:,:,jk) * umask(:,:,jk)
         hv_a(:,:) = hv_a(:,:) + e3v_a(:,:,jk) * vmask(:,:,jk)
      END DO
      !                                        ! Inverse of the local depth
!!gm BUG ?  don't understand the use of umask_i here .....
      r1_hu_a(:,:) = ssumask(:,:) / ( hu_a(:,:) + 1._wp - ssumask(:,:) )
      r1_hv_a(:,:) = ssvmask(:,:) / ( hv_a(:,:) + 1._wp - ssvmask(:,:) )
      !
      IF( ln_timing )   CALL timing_stop('dom_vvl_sf_nxt')
      !
   END SUBROUTINE dom_vvl_sf_nxt


   SUBROUTINE dom_vvl_sf_swp( kt )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_vvl_sf_swp  ***
      !!                   
      !! ** Purpose :  compute time filter and swap of scale factors 
      !!               compute all depths and related variables for next time step
      !!               write outputs and restart file
      !!
      !! ** Method  :  - swap of e3t with trick for volume/tracer conservation
      !!               - reconstruct scale factor at other grid points (interpolate)
      !!               - recompute depths and water height fields
      !!
      !! ** Action  :  - e3t_(b/n), tilde_e3t_(b/n) and e3(u/v)_n ready for next time step
      !!               - Recompute:
      !!                    e3(u/v)_b       
      !!                    e3w_n           
      !!                    e3(u/v)w_b      
      !!                    e3(u/v)w_n      
      !!                    gdept_n, gdepw_n  and gde3w_n
      !!                    h(u/v) and h(u/v)r
      !!
      !! Reference  : Leclair, M., and G. Madec, 2009, Ocean Modelling.
      !!              Leclair, M., and G. Madec, 2011, Ocean Modelling.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! time step
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zcoef        ! local scalar
      !!----------------------------------------------------------------------
      !
      IF( ln_linssh )   RETURN      ! No calculation in linear free surface
      !
      IF( ln_timing )   CALL timing_start('dom_vvl_sf_swp')
      !
      IF( kt == nit000 )   THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dom_vvl_sf_swp : - time filter and swap of scale factors'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   - interpolate scale factors and compute depths for next time step'
      ENDIF
      !
      ! Time filter and swap of scale factors
      ! =====================================
      ! - ML - e3(t/u/v)_b are allready computed in dynnxt.
      IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN
         IF( neuler == 0 .AND. kt == nit000 ) THEN
            tilde_e3t_b(:,:,:) = tilde_e3t_n(:,:,:)
         ELSE
            tilde_e3t_b(:,:,:) = tilde_e3t_n(:,:,:) & 
            &         + atfp * ( tilde_e3t_b(:,:,:) - 2.0_wp * tilde_e3t_n(:,:,:) + tilde_e3t_a(:,:,:) )
         ENDIF
         tilde_e3t_n(:,:,:) = tilde_e3t_a(:,:,:)
      ENDIF
      gdept_b(:,:,:) = gdept_n(:,:,:)
      gdepw_b(:,:,:) = gdepw_n(:,:,:)

      e3t_n(:,:,:) = e3t_a(:,:,:)
      e3u_n(:,:,:) = e3u_a(:,:,:)
      e3v_n(:,:,:) = e3v_a(:,:,:)

      ! Compute all missing vertical scale factor and depths
      ! ====================================================
      ! Horizontal scale factor interpolations
      ! --------------------------------------
      ! - ML - e3u_b and e3v_b are allready computed in dynnxt
      ! - JC - hu_b, hv_b, hur_b, hvr_b also
      
      CALL dom_vvl_interpol( e3u_n(:,:,:), e3f_n(:,:,:), 'F'  )
      
      ! Vertical scale factor interpolations
      CALL dom_vvl_interpol( e3t_n(:,:,:),  e3w_n(:,:,:), 'W'  )
      CALL dom_vvl_interpol( e3u_n(:,:,:), e3uw_n(:,:,:), 'UW' )
      CALL dom_vvl_interpol( e3v_n(:,:,:), e3vw_n(:,:,:), 'VW' )
      CALL dom_vvl_interpol( e3t_b(:,:,:),  e3w_b(:,:,:), 'W'  )
      CALL dom_vvl_interpol( e3u_b(:,:,:), e3uw_b(:,:,:), 'UW' )
      CALL dom_vvl_interpol( e3v_b(:,:,:), e3vw_b(:,:,:), 'VW' )

      ! t- and w- points depth (set the isf depth as it is in the initial step)
      gdept_n(:,:,1) = 0.5_wp * e3w_n(:,:,1)
      gdepw_n(:,:,1) = 0.0_wp
      gde3w_n(:,:,1) = gdept_n(:,:,1) - sshn(:,:)
      DO jk = 2, jpk
         DO jj = 1,jpj
            DO ji = 1,jpi
              !    zcoef = (tmask(ji,jj,jk) - wmask(ji,jj,jk))   ! 0 everywhere tmask = wmask, ie everywhere expect at jk = mikt
                                                                 ! 1 for jk = mikt
               zcoef = (tmask(ji,jj,jk) - wmask(ji,jj,jk))
               gdepw_n(ji,jj,jk) = gdepw_n(ji,jj,jk-1) + e3t_n(ji,jj,jk-1)
               gdept_n(ji,jj,jk) =    zcoef  * ( gdepw_n(ji,jj,jk  ) + 0.5 * e3w_n(ji,jj,jk) )  &
                   &             + (1-zcoef) * ( gdept_n(ji,jj,jk-1) +       e3w_n(ji,jj,jk) ) 
               gde3w_n(ji,jj,jk) = gdept_n(ji,jj,jk) - sshn(ji,jj)
            END DO
         END DO
      END DO

      ! Local depth and Inverse of the local depth of the water
      ! -------------------------------------------------------
      hu_n(:,:) = hu_a(:,:)   ;   r1_hu_n(:,:) = r1_hu_a(:,:)
      hv_n(:,:) = hv_a(:,:)   ;   r1_hv_n(:,:) = r1_hv_a(:,:)
      !
      ht_n(:,:) = e3t_n(:,:,1) * tmask(:,:,1)
      DO jk = 2, jpkm1
         ht_n(:,:) = ht_n(:,:) + e3t_n(:,:,jk) * tmask(:,:,jk)
      END DO

      ! write restart file
      ! ==================
      IF( lrst_oce  )   CALL dom_vvl_rst( kt, 'WRITE' )
      !
      IF( ln_timing )   CALL timing_stop('dom_vvl_sf_swp')
      !
   END SUBROUTINE dom_vvl_sf_swp


   SUBROUTINE dom_vvl_interpol( pe3_in, pe3_out, pout )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dom_vvl__interpol  ***
      !!
      !! ** Purpose :   interpolate scale factors from one grid point to another
      !!
      !! ** Method  :   e3_out = e3_0 + interpolation(e3_in - e3_0)
      !!                - horizontal interpolation: grid cell surface averaging
      !!                - vertical interpolation: simple averaging
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  pe3_in    ! input e3 to be interpolated
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  pe3_out   ! output interpolated e3
      CHARACTER(LEN=*)                , INTENT(in   ) ::  pout      ! grid point of out scale factors
      !                                                             !   =  'U', 'V', 'W, 'F', 'UW' or 'VW'
      !
      INTEGER ::   ji, jj, jk                                       ! dummy loop indices
      REAL(wp) ::  zlnwd                                            ! =1./0. when ln_wd_il = T/F
      !!----------------------------------------------------------------------
      !
      IF(ln_wd_il) THEN
        zlnwd = 1.0_wp
      ELSE
        zlnwd = 0.0_wp
      END IF
      !
      SELECT CASE ( pout )    !==  type of interpolation  ==!
         !
      CASE( 'U' )                   !* from T- to U-point : hor. surface weighted mean
         DO jk = 1, jpk
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  pe3_out(ji,jj,jk) = 0.5_wp * (  umask(ji,jj,jk) * (1.0_wp - zlnwd) + zlnwd ) * r1_e1e2u(ji,jj)   &
                     &                       * (   e1e2t(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) )     &
                     &                           + e1e2t(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'domvvl', pe3_out(:,:,:), 'U', 1._wp )
         pe3_out(:,:,:) = pe3_out(:,:,:) + e3u_0(:,:,:)
         !
      CASE( 'V' )                   !* from T- to V-point : hor. surface weighted mean
         DO jk = 1, jpk
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  pe3_out(ji,jj,jk) = 0.5_wp * ( vmask(ji,jj,jk)  * (1.0_wp - zlnwd) + zlnwd ) * r1_e1e2v(ji,jj)   &
                     &                       * (   e1e2t(ji,jj  ) * ( pe3_in(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) )     &
                     &                           + e1e2t(ji,jj+1) * ( pe3_in(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'domvvl', pe3_out(:,:,:), 'V', 1._wp )
         pe3_out(:,:,:) = pe3_out(:,:,:) + e3v_0(:,:,:)
         !
      CASE( 'F' )                   !* from U-point to F-point : hor. surface weighted mean
         DO jk = 1, jpk
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  pe3_out(ji,jj,jk) = 0.5_wp * (  umask(ji,jj,jk) * umask(ji,jj+1,jk) * (1.0_wp - zlnwd) + zlnwd ) &
                     &                       *    r1_e1e2f(ji,jj)                                                  &
                     &                       * (   e1e2u(ji,jj  ) * ( pe3_in(ji,jj  ,jk) - e3u_0(ji,jj  ,jk) )     &
                     &                           + e1e2u(ji,jj+1) * ( pe3_in(ji,jj+1,jk) - e3u_0(ji,jj+1,jk) ) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'domvvl', pe3_out(:,:,:), 'F', 1._wp )
         pe3_out(:,:,:) = pe3_out(:,:,:) + e3f_0(:,:,:)
         !
      CASE( 'W' )                   !* from T- to W-point : vertical simple mean
         !
         pe3_out(:,:,1) = e3w_0(:,:,1) + pe3_in(:,:,1) - e3t_0(:,:,1)
         ! - ML - The use of mask in this formulea enables the special treatment of the last w-point without indirect adressing
!!gm BUG? use here wmask in case of ISF ?  to be checked
         DO jk = 2, jpk
            pe3_out(:,:,jk) = e3w_0(:,:,jk) + ( 1.0_wp - 0.5_wp * ( tmask(:,:,jk) * (1.0_wp - zlnwd) + zlnwd ) )   &
               &                            * ( pe3_in(:,:,jk-1) - e3t_0(:,:,jk-1) )                               &
               &                            +            0.5_wp * ( tmask(:,:,jk) * (1.0_wp - zlnwd) + zlnwd )     &
               &                            * ( pe3_in(:,:,jk  ) - e3t_0(:,:,jk  ) )
         END DO
         !
      CASE( 'UW' )                  !* from U- to UW-point : vertical simple mean
         !
         pe3_out(:,:,1) = e3uw_0(:,:,1) + pe3_in(:,:,1) - e3u_0(:,:,1)
         ! - ML - The use of mask in this formaula enables the special treatment of the last w- point without indirect adressing
!!gm BUG? use here wumask in case of ISF ?  to be checked
         DO jk = 2, jpk
            pe3_out(:,:,jk) = e3uw_0(:,:,jk) + ( 1.0_wp - 0.5_wp * ( umask(:,:,jk) * (1.0_wp - zlnwd) + zlnwd ) )  &
               &                             * ( pe3_in(:,:,jk-1) - e3u_0(:,:,jk-1) )                              &
               &                             +            0.5_wp * ( umask(:,:,jk) * (1.0_wp - zlnwd) + zlnwd )    &
               &                             * ( pe3_in(:,:,jk  ) - e3u_0(:,:,jk  ) )
         END DO
         !
      CASE( 'VW' )                  !* from V- to VW-point : vertical simple mean
         !
         pe3_out(:,:,1) = e3vw_0(:,:,1) + pe3_in(:,:,1) - e3v_0(:,:,1)
         ! - ML - The use of mask in this formaula enables the special treatment of the last w- point without indirect adressing
!!gm BUG? use here wvmask in case of ISF ?  to be checked
         DO jk = 2, jpk
            pe3_out(:,:,jk) = e3vw_0(:,:,jk) + ( 1.0_wp - 0.5_wp * ( vmask(:,:,jk) * (1.0_wp - zlnwd) + zlnwd ) )  &
               &                             * ( pe3_in(:,:,jk-1) - e3v_0(:,:,jk-1) )                              &
               &                             +            0.5_wp * ( vmask(:,:,jk) * (1.0_wp - zlnwd) + zlnwd )    &
               &                             * ( pe3_in(:,:,jk  ) - e3v_0(:,:,jk  ) )
         END DO
      END SELECT
      !
   END SUBROUTINE dom_vvl_interpol


   SUBROUTINE dom_vvl_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dom_vvl_rst  ***
      !!                     
      !! ** Purpose :   Read or write VVL file in restart file
      !!
      !! ** Method  :   use of IOM library
      !!                if the restart does not contain vertical scale factors,
      !!                they are set to the _0 values
      !!                if the restart does not contain vertical scale factors increments (z_tilde),
      !!                they are set to 0.
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt     ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      !
      INTEGER ::   ji, jj, jk
      INTEGER ::   id1, id2, id3, id4, id5     ! local integers
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise 
         !                                   ! ===============
         IF( ln_rstart ) THEN                   !* Read the restart file
            CALL rst_read_open                  !  open the restart file if necessary
            CALL iom_get( numror, jpdom_autoglo, 'sshn'   , sshn, ldxios = lrxios    )
            !
            id1 = iom_varid( numror, 'e3t_b', ldstop = .FALSE. )
            id2 = iom_varid( numror, 'e3t_n', ldstop = .FALSE. )
            id3 = iom_varid( numror, 'tilde_e3t_b', ldstop = .FALSE. )
            id4 = iom_varid( numror, 'tilde_e3t_n', ldstop = .FALSE. )
            id5 = iom_varid( numror, 'hdiv_lf', ldstop = .FALSE. )
            !                             ! --------- !
            !                             ! all cases !
            !                             ! --------- !
            IF( MIN( id1, id2 ) > 0 ) THEN       ! all required arrays exist
               CALL iom_get( numror, jpdom_autoglo, 'e3t_b', e3t_b(:,:,:), ldxios = lrxios )
               CALL iom_get( numror, jpdom_autoglo, 'e3t_n', e3t_n(:,:,:), ldxios = lrxios )
               ! needed to restart if land processor not computed 
               IF(lwp) write(numout,*) 'dom_vvl_rst : e3t_b and e3t_n found in restart files'
               WHERE ( tmask(:,:,:) == 0.0_wp ) 
                  e3t_n(:,:,:) = e3t_0(:,:,:)
                  e3t_b(:,:,:) = e3t_0(:,:,:)
               END WHERE
               IF( neuler == 0 ) THEN
                  e3t_b(:,:,:) = e3t_n(:,:,:)
               ENDIF
            ELSE IF( id1 > 0 ) THEN
               IF(lwp) write(numout,*) 'dom_vvl_rst WARNING : e3t_n not found in restart files'
               IF(lwp) write(numout,*) 'e3t_n set equal to e3t_b.'
               IF(lwp) write(numout,*) 'neuler is forced to 0'
               CALL iom_get( numror, jpdom_autoglo, 'e3t_b', e3t_b(:,:,:), ldxios = lrxios )
               e3t_n(:,:,:) = e3t_b(:,:,:)
               neuler = 0
            ELSE IF( id2 > 0 ) THEN
               IF(lwp) write(numout,*) 'dom_vvl_rst WARNING : e3t_b not found in restart files'
               IF(lwp) write(numout,*) 'e3t_b set equal to e3t_n.'
               IF(lwp) write(numout,*) 'neuler is forced to 0'
               CALL iom_get( numror, jpdom_autoglo, 'e3t_n', e3t_n(:,:,:), ldxios = lrxios )
               e3t_b(:,:,:) = e3t_n(:,:,:)
               neuler = 0
            ELSE
               IF(lwp) write(numout,*) 'dom_vvl_rst WARNING : e3t_n not found in restart file'
               IF(lwp) write(numout,*) 'Compute scale factor from sshn'
               IF(lwp) write(numout,*) 'neuler is forced to 0'
               DO jk = 1, jpk
                  e3t_n(:,:,jk) =  e3t_0(:,:,jk) * ( ht_0(:,:) + sshn(:,:) ) &
                      &                          / ( ht_0(:,:) + 1._wp - ssmask(:,:) ) * tmask(:,:,jk)   &
                      &          + e3t_0(:,:,jk)                               * (1._wp -tmask(:,:,jk))
               END DO
               e3t_b(:,:,:) = e3t_n(:,:,:)
               neuler = 0
            ENDIF
            !                             ! ----------- !
            IF( ln_vvl_zstar ) THEN       ! z_star case !
               !                          ! ----------- !
               IF( MIN( id3, id4 ) > 0 ) THEN
                  CALL ctl_stop( 'dom_vvl_rst: z_star cannot restart from a z_tilde or layer run' )
               ENDIF
               !                          ! ----------------------- !
            ELSE                          ! z_tilde and layer cases !
               !                          ! ----------------------- !
               IF( MIN( id3, id4 ) > 0 ) THEN  ! all required arrays exist
                  CALL iom_get( numror, jpdom_autoglo, 'tilde_e3t_b', tilde_e3t_b(:,:,:), ldxios = lrxios )
                  CALL iom_get( numror, jpdom_autoglo, 'tilde_e3t_n', tilde_e3t_n(:,:,:), ldxios = lrxios )
               ELSE                            ! one at least array is missing
                  tilde_e3t_b(:,:,:) = 0.0_wp
                  tilde_e3t_n(:,:,:) = 0.0_wp
               ENDIF
               !                          ! ------------ !
               IF( ln_vvl_ztilde ) THEN   ! z_tilde case !
                  !                       ! ------------ !
                  IF( id5 > 0 ) THEN  ! required array exists
                     CALL iom_get( numror, jpdom_autoglo, 'hdiv_lf', hdiv_lf(:,:,:), ldxios = lrxios )
                  ELSE                ! array is missing
                     hdiv_lf(:,:,:) = 0.0_wp
                  ENDIF
               ENDIF
            ENDIF
            !
         ELSE                                   !* Initialize at "rest"
            !

            IF( ll_wd ) THEN   ! MJB ll_wd edits start here - these are essential 
               !
               IF( cn_cfg == 'wad' ) THEN
                  ! Wetting and drying test case
                  CALL usr_def_istate( gdept_b, tmask, tsb, ub, vb, sshb  )
                  tsn  (:,:,:,:) = tsb (:,:,:,:)       ! set now values from to before ones
                  sshn (:,:)     = sshb(:,:)
                  un   (:,:,:)   = ub  (:,:,:)
                  vn   (:,:,:)   = vb  (:,:,:)
               ELSE
                  ! if not test case
                  sshn(:,:) = -ssh_ref
                  sshb(:,:) = -ssh_ref

                  DO jj = 1, jpj
                     DO ji = 1, jpi
                        IF( ht_0(ji,jj)-ssh_ref <  rn_wdmin1 ) THEN ! if total depth is less than min depth

                           sshb(ji,jj) = rn_wdmin1 - (ht_0(ji,jj) )
                           sshn(ji,jj) = rn_wdmin1 - (ht_0(ji,jj) )
                           ssha(ji,jj) = rn_wdmin1 - (ht_0(ji,jj) )
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF !If test case else

               ! Adjust vertical metrics for all wad
               DO jk = 1, jpk
                  e3t_n(:,:,jk) =  e3t_0(:,:,jk) * ( ht_0(:,:) + sshn(:,:)  ) &
                    &                            / ( ht_0(:,:) + 1._wp - ssmask(:,:) ) * tmask(:,:,jk)   &
                    &            + e3t_0(:,:,jk) * ( 1._wp - tmask(:,:,jk) )
               END DO
               e3t_b(:,:,:) = e3t_n(:,:,:)

               DO ji = 1, jpi
                  DO jj = 1, jpj
                     IF ( ht_0(ji,jj) .LE. 0.0 .AND. NINT( ssmask(ji,jj) ) .EQ. 1) THEN
                       CALL ctl_stop( 'dom_vvl_rst: ht_0 must be positive at potentially wet points' )
                     ENDIF
                  END DO 
               END DO 
               !
            ELSE
               !
               ! usr_def_istate called here only to get sshb, that is needed to initialize e3t_b and e3t_n
               CALL usr_def_istate( gdept_0, tmask, tsb, ub, vb, sshb  )  
               ! usr_def_istate will be called again in istate_init to initialize ts(bn), ssh(bn), u(bn) and v(bn)
               !
               DO jk=1,jpk
                  e3t_b(:,:,jk) =  e3t_0(:,:,jk) * ( ht_0(:,:) + sshb(:,:) ) &
                    &                            / ( ht_0(:,:) + 1._wp - ssmask(:,:) ) * tmask(:,:,jk)   &
                    &            + e3t_0(:,:,jk) * ( 1._wp - tmask(:,:,jk) )   ! make sure e3t_b != 0 on land points
               END DO
               e3t_n(:,:,:) = e3t_b(:,:,:)
               sshn(:,:) = sshb(:,:)   ! needed later for gde3w
!!$                e3t_n(:,:,:)=e3t_0(:,:,:)
!!$                e3t_b(:,:,:)=e3t_0(:,:,:)
               !
            END IF           ! end of ll_wd edits

            IF( ln_vvl_ztilde .OR. ln_vvl_layer) THEN
               tilde_e3t_b(:,:,:) = 0._wp
               tilde_e3t_n(:,:,:) = 0._wp
               IF( ln_vvl_ztilde ) hdiv_lf(:,:,:) = 0._wp
            END IF
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! ===================
         IF(lwp) WRITE(numout,*) '---- dom_vvl_rst ----'
         IF( lwxios ) CALL iom_swap(      cwxios_context          )
         !                                           ! --------- !
         !                                           ! all cases !
         !                                           ! --------- !
         CALL iom_rstput( kt, nitrst, numrow, 'e3t_b', e3t_b(:,:,:), ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'e3t_n', e3t_n(:,:,:), ldxios = lwxios )
         !                                           ! ----------------------- !
         IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN  ! z_tilde and layer cases !
            !                                        ! ----------------------- !
            CALL iom_rstput( kt, nitrst, numrow, 'tilde_e3t_b', tilde_e3t_b(:,:,:), ldxios = lwxios)
            CALL iom_rstput( kt, nitrst, numrow, 'tilde_e3t_n', tilde_e3t_n(:,:,:), ldxios = lwxios)
         END IF
         !                                           ! -------------!    
         IF( ln_vvl_ztilde ) THEN                    ! z_tilde case !
            !                                        ! ------------ !
            CALL iom_rstput( kt, nitrst, numrow, 'hdiv_lf', hdiv_lf(:,:,:), ldxios = lwxios)
         ENDIF
         !
         IF( lwxios ) CALL iom_swap(      cxios_context          )
      ENDIF
      !
   END SUBROUTINE dom_vvl_rst


   SUBROUTINE dom_vvl_ctl
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dom_vvl_ctl  ***
      !!                
      !! ** Purpose :   Control the consistency between namelist options
      !!                for vertical coordinate
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ios
      !!
      NAMELIST/nam_vvl/ ln_vvl_zstar, ln_vvl_ztilde, ln_vvl_layer, ln_vvl_ztilde_as_zstar, &
         &              ln_vvl_zstar_at_eqtor      , rn_ahe3     , rn_rst_e3t            , &
         &              rn_lf_cutoff               , rn_zdef_max , ln_vvl_dbg                ! not yet implemented: ln_vvl_kepe
      !!---------------------------------------------------------------------- 
      !
      REWIND( numnam_ref )              ! Namelist nam_vvl in reference namelist : 
      READ  ( numnam_ref, nam_vvl, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nam_vvl in reference namelist', lwp )
      REWIND( numnam_cfg )              ! Namelist nam_vvl in configuration namelist : Parameters of the run
      READ  ( numnam_cfg, nam_vvl, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'nam_vvl in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nam_vvl )
      !
      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'dom_vvl_ctl : choice/control of the variable vertical coordinate'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist nam_vvl : chose a vertical coordinate'
         WRITE(numout,*) '      zstar                      ln_vvl_zstar   = ', ln_vvl_zstar
         WRITE(numout,*) '      ztilde                     ln_vvl_ztilde  = ', ln_vvl_ztilde
         WRITE(numout,*) '      layer                      ln_vvl_layer   = ', ln_vvl_layer
         WRITE(numout,*) '      ztilde as zstar   ln_vvl_ztilde_as_zstar  = ', ln_vvl_ztilde_as_zstar
         WRITE(numout,*) '      ztilde near the equator    ln_vvl_zstar_at_eqtor  = ', ln_vvl_zstar_at_eqtor
         WRITE(numout,*) '      !'
         WRITE(numout,*) '      thickness diffusion coefficient                      rn_ahe3      = ', rn_ahe3
         WRITE(numout,*) '      maximum e3t deformation fractional change            rn_zdef_max  = ', rn_zdef_max
         IF( ln_vvl_ztilde_as_zstar ) THEN
            WRITE(numout,*) '      ztilde running in zstar emulation mode (ln_vvl_ztilde_as_zstar=T) '
            WRITE(numout,*) '         ignoring namelist timescale parameters and using:'
            WRITE(numout,*) '            hard-wired : z-tilde to zstar restoration timescale (days)'
            WRITE(numout,*) '                         rn_rst_e3t     = 0.e0'
            WRITE(numout,*) '            hard-wired : z-tilde cutoff frequency of low-pass filter (days)'
            WRITE(numout,*) '                         rn_lf_cutoff   = 1.0/rdt'
         ELSE
            WRITE(numout,*) '      z-tilde to zstar restoration timescale (days)        rn_rst_e3t   = ', rn_rst_e3t
            WRITE(numout,*) '      z-tilde cutoff frequency of low-pass filter (days)   rn_lf_cutoff = ', rn_lf_cutoff
         ENDIF
         WRITE(numout,*) '         debug prints flag                                 ln_vvl_dbg   = ', ln_vvl_dbg
      ENDIF
      !
      ioptio = 0                      ! Parameter control
      IF( ln_vvl_ztilde_as_zstar )   ln_vvl_ztilde = .true.
      IF( ln_vvl_zstar           )   ioptio = ioptio + 1
      IF( ln_vvl_ztilde          )   ioptio = ioptio + 1
      IF( ln_vvl_layer           )   ioptio = ioptio + 1
      !
      IF( ioptio /= 1 )   CALL ctl_stop( 'Choose ONE vertical coordinate in namelist nam_vvl' )
      IF( .NOT. ln_vvl_zstar .AND. ln_isf ) CALL ctl_stop( 'Only vvl_zstar has been tested with ice shelf cavity' )
      !
      IF(lwp) THEN                   ! Print the choice
         WRITE(numout,*)
         IF( ln_vvl_zstar           ) WRITE(numout,*) '      ==>>>   zstar vertical coordinate is used'
         IF( ln_vvl_ztilde          ) WRITE(numout,*) '      ==>>>   ztilde vertical coordinate is used'
         IF( ln_vvl_layer           ) WRITE(numout,*) '      ==>>>   layer vertical coordinate is used'
         IF( ln_vvl_ztilde_as_zstar ) WRITE(numout,*) '      ==>>>   to emulate a zstar coordinate'
      ENDIF
      !
#if defined key_agrif
      IF( (.NOT.Agrif_Root()).AND.(.NOT.ln_vvl_zstar) )   CALL ctl_stop( 'AGRIF is implemented with zstar coordinate only' )
#endif
      !
   END SUBROUTINE dom_vvl_ctl

   !!======================================================================
END MODULE domvvl
