MODULE zdfiwm
   !!========================================================================
   !!                       ***  MODULE  zdfiwm  ***
   !! Ocean physics: Internal gravity wave-driven vertical mixing
   !!========================================================================
   !! History :  1.0  !  2004-04  (L. Bessieres, G. Madec)  Original code
   !!             -   !  2006-08  (A. Koch-Larrouy)  Indonesian strait
   !!            3.3  !  2010-10  (C. Ethe, G. Madec)  reorganisation of initialisation phase
   !!            3.6  !  2016-03  (C. de Lavergne)  New param: internal wave-driven mixing 
   !!            4.0  !  2017-04  (G. Madec)  renamed module, remove the old param. and the CPP keys
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zdf_iwm       : global     momentum & tracer Kz with wave induced Kz
   !!   zdf_iwm_init  : global     momentum & tracer Kz with wave induced Kz
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables
   USE zdf_oce        ! ocean vertical physics variables
   USE zdfddm         ! ocean vertical physics: double diffusive mixing
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE eosbn2         ! ocean equation of state
   USE phycst         ! physical constants
   !
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE iom            ! I/O Manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_iwm        ! called in step module 
   PUBLIC   zdf_iwm_init   ! called in nemogcm module 

   !                      !!* Namelist  namzdf_iwm : internal wave-driven mixing *
   INTEGER ::  nn_zpyc     ! pycnocline-intensified mixing energy proportional to N (=1) or N^2 (=2)
   LOGICAL ::  ln_mevar    ! variable (=T) or constant (=F) mixing efficiency
   LOGICAL ::  ln_tsdiff   ! account for differential T/S wave-driven mixing (=T) or not (=F)

   REAL(wp)::  r1_6 = 1._wp / 6._wp

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ebot_iwm   ! power available from high-mode wave breaking (W/m2)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   epyc_iwm   ! power available from low-mode, pycnocline-intensified wave breaking (W/m2)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ecri_iwm   ! power available from low-mode, critical slope wave breaking (W/m2)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hbot_iwm   ! WKB decay scale for high-mode energy dissipation (m)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hcri_iwm   ! decay scale for low-mode critical slope dissipation (m)

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: zdfiwm.F90 8093 2017-05-30 08:13:14Z gm $
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_iwm_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION zdf_iwm_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( ebot_iwm(jpi,jpj),  epyc_iwm(jpi,jpj),  ecri_iwm(jpi,jpj) ,     &
      &         hbot_iwm(jpi,jpj),  hcri_iwm(jpi,jpj)                     , STAT=zdf_iwm_alloc )
      !
      CALL mpp_sum ( 'zdfiwm', zdf_iwm_alloc )
      IF( zdf_iwm_alloc /= 0 )   CALL ctl_stop( 'STOP', 'zdf_iwm_alloc: failed to allocate arrays' )
   END FUNCTION zdf_iwm_alloc


   SUBROUTINE zdf_iwm( kt, p_avm, p_avt, p_avs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_iwm  ***
      !!                   
      !! ** Purpose :   add to the vertical mixing coefficients the effect of
      !!              breaking internal waves.
      !!
      !! ** Method  : - internal wave-driven vertical mixing is given by:
      !!                  Kz_wave = min(  100 cm2/s, f(  Reb = zemx_iwm /( Nu * N^2 )  )
      !!              where zemx_iwm is the 3D space distribution of the wave-breaking 
      !!              energy and Nu the molecular kinematic viscosity.
      !!              The function f(Reb) is linear (constant mixing efficiency)
      !!              if the namelist parameter ln_mevar = F and nonlinear if ln_mevar = T.
      !!
      !!              - Compute zemx_iwm, the 3D power density that allows to compute
      !!              Reb and therefrom the wave-induced vertical diffusivity.
      !!              This is divided into three components:
      !!                 1. Bottom-intensified low-mode dissipation at critical slopes
      !!                     zemx_iwm(z) = ( ecri_iwm / rau0 ) * EXP( -(H-z)/hcri_iwm )
      !!                                   / ( 1. - EXP( - H/hcri_iwm ) ) * hcri_iwm
      !!              where hcri_iwm is the characteristic length scale of the bottom 
      !!              intensification, ecri_iwm a map of available power, and H the ocean depth.
      !!                 2. Pycnocline-intensified low-mode dissipation
      !!                     zemx_iwm(z) = ( epyc_iwm / rau0 ) * ( sqrt(rn2(z))^nn_zpyc )
      !!                                   / SUM( sqrt(rn2(z))^nn_zpyc * e3w(z) )
      !!              where epyc_iwm is a map of available power, and nn_zpyc
      !!              is the chosen stratification-dependence of the internal wave
      !!              energy dissipation.
      !!                 3. WKB-height dependent high mode dissipation
      !!                     zemx_iwm(z) = ( ebot_iwm / rau0 ) * rn2(z) * EXP(-z_wkb(z)/hbot_iwm)
      !!                                   / SUM( rn2(z) * EXP(-z_wkb(z)/hbot_iwm) * e3w(z) )
      !!              where hbot_iwm is the characteristic length scale of the WKB bottom 
      !!              intensification, ebot_iwm is a map of available power, and z_wkb is the
      !!              WKB-stretched height above bottom defined as
      !!                    z_wkb(z) = H * SUM( sqrt(rn2(z'>=z)) * e3w(z'>=z) )
      !!                                 / SUM( sqrt(rn2(z'))    * e3w(z')    )
      !!
      !!              - update the model vertical eddy viscosity and diffusivity: 
      !!                     avt  = avt  +    av_wave
      !!                     avm  = avm  +    av_wave
      !!
      !!              - if namelist parameter ln_tsdiff = T, account for differential mixing:
      !!                     avs  = avt  +    av_wave * diffusivity_ratio(Reb)
      !!
      !! ** Action  : - avt, avs, avm, increased by tide internal wave-driven mixing    
      !!
      !! References :  de Lavergne et al. 2015, JPO; 2016, in prep.
      !!----------------------------------------------------------------------
      INTEGER                    , INTENT(in   ) ::   kt             ! ocean time step
      REAL(wp), DIMENSION(:,:,:) , INTENT(inout) ::   p_avm          ! momentum Kz (w-points)
      REAL(wp), DIMENSION(:,:,:) , INTENT(inout) ::   p_avt, p_avs   ! tracer   Kz (w-points)
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zztmp        ! scalar workspace
      REAL(wp), DIMENSION(jpi,jpj)     ::   zfact       ! Used for vertical structure
      REAL(wp), DIMENSION(jpi,jpj)     ::   zhdep       ! Ocean depth
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zwkb        ! WKB-stretched height above bottom
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zweight     ! Weight for high mode vertical distribution
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   znu_t       ! Molecular kinematic viscosity (T grid)
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   znu_w       ! Molecular kinematic viscosity (W grid)
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zReb        ! Turbulence intensity parameter
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zemx_iwm    ! local energy density available for mixing (W/kg)
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zav_ratio   ! S/T diffusivity ratio (only for ln_tsdiff=T)
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zav_wave    ! Internal wave-induced diffusivity
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   z3d  ! 3D workspace used for iom_put 
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   z2d  ! 2D     -      -    -     -
      !!----------------------------------------------------------------------
      !
      !                       !* Set to zero the 1st and last vertical levels of appropriate variables
      zemx_iwm (:,:,1) = 0._wp   ;   zemx_iwm (:,:,jpk) = 0._wp
      zav_ratio(:,:,1) = 0._wp   ;   zav_ratio(:,:,jpk) = 0._wp
      zav_wave (:,:,1) = 0._wp   ;   zav_wave (:,:,jpk) = 0._wp
      !
      !                       ! ----------------------------- !
      !                       !  Internal wave-driven mixing  !  (compute zav_wave)
      !                       ! ----------------------------- !
      !                             
      !                       !* Critical slope mixing: distribute energy over the time-varying ocean depth,
      !                                                 using an exponential decay from the seafloor.
      DO jj = 1, jpj                ! part independent of the level
         DO ji = 1, jpi
            zhdep(ji,jj) = gdepw_0(ji,jj,mbkt(ji,jj)+1)       ! depth of the ocean
            zfact(ji,jj) = rau0 * (  1._wp - EXP( -zhdep(ji,jj) / hcri_iwm(ji,jj) )  )
            IF( zfact(ji,jj) /= 0._wp )   zfact(ji,jj) = ecri_iwm(ji,jj) / zfact(ji,jj)
         END DO
      END DO
!!gm gde3w ==>>>  check for ssh taken into account.... seem OK gde3w_n=gdept_n - sshn
      DO jk = 2, jpkm1              ! complete with the level-dependent part
         zemx_iwm(:,:,jk) = zfact(:,:) * (  EXP( ( gde3w_n(:,:,jk  ) - zhdep(:,:) ) / hcri_iwm(:,:) )                      &
            &                             - EXP( ( gde3w_n(:,:,jk-1) - zhdep(:,:) ) / hcri_iwm(:,:) )  ) * wmask(:,:,jk)   &
            &                          / ( gde3w_n(:,:,jk) - gde3w_n(:,:,jk-1) )

!!gm delta(gde3w_n) = e3t_n  !!  Please verify the grid-point position w versus t-point
!!gm it seems to me that only 1/hcri_iwm  is used ==>  compute it one for all

      END DO

      !                        !* Pycnocline-intensified mixing: distribute energy over the time-varying 
      !                        !* ocean depth as proportional to sqrt(rn2)^nn_zpyc
      !                                          ! (NB: N2 is masked, so no use of wmask here)
      SELECT CASE ( nn_zpyc )
      !
      CASE ( 1 )               ! Dissipation scales as N (recommended)
         !
         zfact(:,:) = 0._wp
         DO jk = 2, jpkm1              ! part independent of the level
            zfact(:,:) = zfact(:,:) + e3w_n(:,:,jk) * SQRT(  MAX( 0._wp, rn2(:,:,jk) )  ) * wmask(:,:,jk)
         END DO
         !
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( zfact(ji,jj) /= 0 )   zfact(ji,jj) = epyc_iwm(ji,jj) / ( rau0 * zfact(ji,jj) )
            END DO
         END DO
         !
         DO jk = 2, jpkm1              ! complete with the level-dependent part
            zemx_iwm(:,:,jk) = zemx_iwm(:,:,jk) + zfact(:,:) * SQRT(  MAX( 0._wp, rn2(:,:,jk) )  ) * wmask(:,:,jk)
         END DO
         !
      CASE ( 2 )               ! Dissipation scales as N^2
         !
         zfact(:,:) = 0._wp
         DO jk = 2, jpkm1              ! part independent of the level
            zfact(:,:) = zfact(:,:) + e3w_n(:,:,jk) * MAX( 0._wp, rn2(:,:,jk) ) * wmask(:,:,jk)
         END DO
         !
         DO jj= 1, jpj
            DO ji = 1, jpi
               IF( zfact(ji,jj) /= 0 )   zfact(ji,jj) = epyc_iwm(ji,jj) / ( rau0 * zfact(ji,jj) )
            END DO
         END DO
         !
         DO jk = 2, jpkm1              ! complete with the level-dependent part
            zemx_iwm(:,:,jk) = zemx_iwm(:,:,jk) + zfact(:,:) * MAX( 0._wp, rn2(:,:,jk) ) * wmask(:,:,jk)
         END DO
         !
      END SELECT

      !                        !* WKB-height dependent mixing: distribute energy over the time-varying 
      !                        !* ocean depth as proportional to rn2 * exp(-z_wkb/rn_hbot)
      !
      zwkb (:,:,:) = 0._wp
      zfact(:,:)   = 0._wp
      DO jk = 2, jpkm1
         zfact(:,:) = zfact(:,:) + e3w_n(:,:,jk) * SQRT(  MAX( 0._wp, rn2(:,:,jk) )  ) * wmask(:,:,jk)
         zwkb(:,:,jk) = zfact(:,:)
      END DO
!!gm even better:
!      DO jk = 2, jpkm1
!         zwkb(:,:) = zwkb(:,:) + e3w_n(:,:,jk) * SQRT(  MAX( 0._wp, rn2(:,:,jk) )  )
!      END DO
!      zfact(:,:) = zwkb(:,:,jpkm1)
!!gm or just use zwkb(k=jpk-1) instead of zfact...
!!gm
      !
      DO jk = 2, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( zfact(ji,jj) /= 0 )   zwkb(ji,jj,jk) = zhdep(ji,jj) * ( zfact(ji,jj) - zwkb(ji,jj,jk) )   &
                  &                                     * wmask(ji,jj,jk) / zfact(ji,jj)
            END DO
         END DO
      END DO
      zwkb(:,:,1) = zhdep(:,:) * wmask(:,:,1)
      !
      zweight(:,:,:) = 0._wp
      DO jk = 2, jpkm1
         zweight(:,:,jk) = MAX( 0._wp, rn2(:,:,jk) ) * hbot_iwm(:,:) * wmask(:,:,jk)                    &
            &   * (  EXP( -zwkb(:,:,jk) / hbot_iwm(:,:) ) - EXP( -zwkb(:,:,jk-1) / hbot_iwm(:,:) )  )
      END DO
      !
      zfact(:,:) = 0._wp
      DO jk = 2, jpkm1              ! part independent of the level
         zfact(:,:) = zfact(:,:) + zweight(:,:,jk)
      END DO
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( zfact(ji,jj) /= 0 )   zfact(ji,jj) = ebot_iwm(ji,jj) / ( rau0 * zfact(ji,jj) )
         END DO
      END DO
      !
      DO jk = 2, jpkm1              ! complete with the level-dependent part
         zemx_iwm(:,:,jk) = zemx_iwm(:,:,jk) + zweight(:,:,jk) * zfact(:,:) * wmask(:,:,jk)   &
            &                                / ( gde3w_n(:,:,jk) - gde3w_n(:,:,jk-1) )
!!gm  use of e3t_n just above?
      END DO
      !
!!gm  this is to be replaced by just a constant value znu=1.e-6 m2/s
      ! Calculate molecular kinematic viscosity
      znu_t(:,:,:) = 1.e-4_wp * (  17.91_wp - 0.53810_wp * tsn(:,:,:,jp_tem) + 0.00694_wp * tsn(:,:,:,jp_tem) * tsn(:,:,:,jp_tem)  &
         &                                  + 0.02305_wp * tsn(:,:,:,jp_sal)  ) * tmask(:,:,:) * r1_rau0
      DO jk = 2, jpkm1
         znu_w(:,:,jk) = 0.5_wp * ( znu_t(:,:,jk-1) + znu_t(:,:,jk) ) * wmask(:,:,jk)
      END DO
!!gm end
      !
      ! Calculate turbulence intensity parameter Reb
      DO jk = 2, jpkm1
         zReb(:,:,jk) = zemx_iwm(:,:,jk) / MAX( 1.e-20_wp, znu_w(:,:,jk) * rn2(:,:,jk) )
      END DO
      !
      ! Define internal wave-induced diffusivity
      DO jk = 2, jpkm1
         zav_wave(:,:,jk) = znu_w(:,:,jk) * zReb(:,:,jk) * r1_6   ! This corresponds to a constant mixing efficiency of 1/6
      END DO
      !
      IF( ln_mevar ) THEN              ! Variable mixing efficiency case : modify zav_wave in the
         DO jk = 2, jpkm1              ! energetic (Reb > 480) and buoyancy-controlled (Reb <10.224 ) regimes
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( zReb(ji,jj,jk) > 480.00_wp ) THEN
                     zav_wave(ji,jj,jk) = 3.6515_wp * znu_w(ji,jj,jk) * SQRT( zReb(ji,jj,jk) )
                  ELSEIF( zReb(ji,jj,jk) < 10.224_wp ) THEN
                     zav_wave(ji,jj,jk) = 0.052125_wp * znu_w(ji,jj,jk) * zReb(ji,jj,jk) * SQRT( zReb(ji,jj,jk) )
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
      !
      DO jk = 2, jpkm1                 ! Bound diffusivity by molecular value and 100 cm2/s
         zav_wave(:,:,jk) = MIN(  MAX( 1.4e-7_wp, zav_wave(:,:,jk) ), 1.e-2_wp  ) * wmask(:,:,jk)
      END DO
      !
      IF( kt == nit000 ) THEN        !* Control print at first time-step: diagnose the energy consumed by zav_wave
         zztmp = 0._wp
!!gm used of glosum 3D....
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zztmp = zztmp + e3w_n(ji,jj,jk) * e1e2t(ji,jj)   &
                     &          * MAX( 0._wp, rn2(ji,jj,jk) ) * zav_wave(ji,jj,jk) * wmask(ji,jj,jk) * tmask_i(ji,jj)
               END DO
            END DO
         END DO
         CALL mpp_sum( 'zdfiwm', zztmp )
         zztmp = rau0 * zztmp ! Global integral of rauo * Kz * N^2 = power contributing to mixing 
         !
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'zdf_iwm : Internal wave-driven mixing (iwm)'
            WRITE(numout,*) '~~~~~~~ '
            WRITE(numout,*)
            WRITE(numout,*) '      Total power consumption by av_wave =  ', zztmp * 1.e-12_wp, 'TW'
         ENDIF
      ENDIF

      !                          ! ----------------------- !
      !                          !   Update  mixing coefs  !                          
      !                          ! ----------------------- !
      !      
      IF( ln_tsdiff ) THEN          !* Option for differential mixing of salinity and temperature
         DO jk = 2, jpkm1              ! Calculate S/T diffusivity ratio as a function of Reb
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zav_ratio(ji,jj,jk) = ( 0.505_wp + 0.495_wp *                                                                  &
                      &   TANH(    0.92_wp * (   LOG10(  MAX( 1.e-20_wp, zReb(ji,jj,jk) * 5._wp * r1_6 )  ) - 0.60_wp   )    )   &
                      &                 ) * wmask(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL iom_put( "av_ratio", zav_ratio )
         DO jk = 2, jpkm1           !* update momentum & tracer diffusivity with wave-driven mixing
            p_avs(:,:,jk) = p_avs(:,:,jk) + zav_wave(:,:,jk) * zav_ratio(:,:,jk)
            p_avt(:,:,jk) = p_avt(:,:,jk) + zav_wave(:,:,jk)
            p_avm(:,:,jk) = p_avm(:,:,jk) + zav_wave(:,:,jk)
         END DO
         !
      ELSE                          !* update momentum & tracer diffusivity with wave-driven mixing
         DO jk = 2, jpkm1
            p_avs(:,:,jk) = p_avs(:,:,jk) + zav_wave(:,:,jk)
            p_avt(:,:,jk) = p_avt(:,:,jk) + zav_wave(:,:,jk)
            p_avm(:,:,jk) = p_avm(:,:,jk) + zav_wave(:,:,jk)
         END DO
      ENDIF

      !                             !* output internal wave-driven mixing coefficient
      CALL iom_put( "av_wave", zav_wave )
                                    !* output useful diagnostics: Kz*N^2 , 
!!gm Kz*N2 should take into account the ratio avs/avt if it is used.... (see diaar5)
                                    !  vertical integral of rau0 * Kz * N^2 , energy density (zemx_iwm)
      IF( iom_use("bflx_iwm") .OR. iom_use("pcmap_iwm") ) THEN
         ALLOCATE( z2d(jpi,jpj) , z3d(jpi,jpj,jpk) )
         z3d(:,:,:) = MAX( 0._wp, rn2(:,:,:) ) * zav_wave(:,:,:)
         z2d(:,:) = 0._wp
         DO jk = 2, jpkm1
            z2d(:,:) = z2d(:,:) + e3w_n(:,:,jk) * z3d(:,:,jk) * wmask(:,:,jk)
         END DO
         z2d(:,:) = rau0 * z2d(:,:)
         CALL iom_put( "bflx_iwm", z3d )
         CALL iom_put( "pcmap_iwm", z2d )
         DEALLOCATE( z2d , z3d )
      ENDIF
      CALL iom_put( "emix_iwm", zemx_iwm )
      
      IF(ln_ctl)   CALL prt_ctl(tab3d_1=zav_wave , clinfo1=' iwm - av_wave: ', tab3d_2=avt, clinfo2=' avt: ', kdim=jpk)
      !
   END SUBROUTINE zdf_iwm


   SUBROUTINE zdf_iwm_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_iwm_init  ***
      !!                     
      !! ** Purpose :   Initialization of the wave-driven vertical mixing, reading
      !!              of input power maps and decay length scales in netcdf files.
      !!
      !! ** Method  : - Read the namzdf_iwm namelist and check the parameters
      !!
      !!              - Read the input data in NetCDF files :
      !!              power available from high-mode wave breaking (mixing_power_bot.nc)
      !!              power available from pycnocline-intensified wave-breaking (mixing_power_pyc.nc)
      !!              power available from critical slope wave-breaking (mixing_power_cri.nc)
      !!              WKB decay scale for high-mode wave-breaking (decay_scale_bot.nc)
      !!              decay scale for critical slope wave-breaking (decay_scale_cri.nc)
      !!
      !! ** input   : - Namlist namzdf_iwm
      !!              - NetCDF files : mixing_power_bot.nc, mixing_power_pyc.nc, mixing_power_cri.nc,
      !!              decay_scale_bot.nc decay_scale_cri.nc
      !!
      !! ** Action  : - Increase by 1 the nstop flag is setting problem encounter
      !!              - Define ebot_iwm, epyc_iwm, ecri_iwm, hbot_iwm, hcri_iwm
      !!
      !! References : de Lavergne et al. JPO, 2015 ; de Lavergne PhD 2016
      !!              de Lavergne et al. in prep., 2017
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   inum         ! local integer
      INTEGER  ::   ios
      REAL(wp) ::   zbot, zpyc, zcri   ! local scalars
      !!
      NAMELIST/namzdf_iwm/ nn_zpyc, ln_mevar, ln_tsdiff
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )              ! Namelist namzdf_iwm in reference namelist : Wave-driven mixing
      READ  ( numnam_ref, namzdf_iwm, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namzdf_iwm in reference namelist', lwp )
      !
      REWIND( numnam_cfg )              ! Namelist namzdf_iwm in configuration namelist : Wave-driven mixing
      READ  ( numnam_cfg, namzdf_iwm, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namzdf_iwm in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namzdf_iwm )
      !
      IF(lwp) THEN                  ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_iwm_init : internal wave-driven mixing'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namzdf_iwm : set wave-driven mixing parameters'
         WRITE(numout,*) '      Pycnocline-intensified diss. scales as N (=1) or N^2 (=2) = ', nn_zpyc
         WRITE(numout,*) '      Variable (T) or constant (F) mixing efficiency            = ', ln_mevar
         WRITE(numout,*) '      Differential internal wave-driven mixing (T) or not (F)   = ', ln_tsdiff
      ENDIF
      
      ! The new wave-driven mixing parameterization elevates avt and avm in the interior, and
      ! ensures that avt remains larger than its molecular value (=1.4e-7). Therefore, avtb should 
      ! be set here to a very small value, and avmb to its (uniform) molecular value (=1.4e-6).
      avmb(:) = 1.4e-6_wp        ! viscous molecular value
      avtb(:) = 1.e-10_wp        ! very small diffusive minimum (background avt is specified in zdf_iwm)    
      avtb_2d(:,:) = 1.e0_wp     ! uniform 
      IF(lwp) THEN                  ! Control print
         WRITE(numout,*)
         WRITE(numout,*) '   Force the background value applied to avm & avt in TKE to be everywhere ',   &
            &               'the viscous molecular value & a very small diffusive value, resp.'
      ENDIF
            
      !                             ! allocate iwm arrays
      IF( zdf_iwm_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'zdf_iwm_init : unable to allocate iwm arrays' )
      !
      !                             ! read necessary fields
!!$      CALL iom_open('mixing_power_bot',inum)       ! energy flux for high-mode wave breaking [W/m2]
!!$      CALL iom_get  (inum, jpdom_data, 'field', ebot_iwm, 1 ) 
!!$      CALL iom_close(inum)
      ebot_iwm(:,:) = 1.e-6
      !
!!$      CALL iom_open('mixing_power_pyc',inum)       ! energy flux for pynocline-intensified wave breaking [W/m2]
!!$      CALL iom_get  (inum, jpdom_data, 'field', epyc_iwm, 1 )
!!$      CALL iom_close(inum)
      epyc_iwm(:,:) = 1.e-6
      !
!!$      CALL iom_open('mixing_power_cri',inum)       ! energy flux for critical slope wave breaking [W/m2]
!!$      CALL iom_get  (inum, jpdom_data, 'field', ecri_iwm, 1 )
!!$      CALL iom_close(inum)
      ecri_iwm(:,:) = 1.e-10
      !
!!$      CALL iom_open('decay_scale_bot',inum)        ! spatially variable decay scale for high-mode wave breaking [m]
!!$      CALL iom_get  (inum, jpdom_data, 'field', hbot_iwm, 1 )
!!$      CALL iom_close(inum)
      hbot_iwm(:,:) = 100.
      !
!!$      CALL iom_open('decay_scale_cri',inum)        ! spatially variable decay scale for critical slope wave breaking [m]
!!$      CALL iom_get  (inum, jpdom_data, 'field', hcri_iwm, 1 )
!!$      CALL iom_close(inum)
      hcri_iwm(:,:) = 100.

      ebot_iwm(:,:) = ebot_iwm(:,:) * ssmask(:,:)
      epyc_iwm(:,:) = epyc_iwm(:,:) * ssmask(:,:)
      ecri_iwm(:,:) = ecri_iwm(:,:) * ssmask(:,:)

      zbot = glob_sum( 'zdfiwm', e1e2t(:,:) * ebot_iwm(:,:) )
      zpyc = glob_sum( 'zdfiwm', e1e2t(:,:) * epyc_iwm(:,:) )
      zcri = glob_sum( 'zdfiwm', e1e2t(:,:) * ecri_iwm(:,:) )
      IF(lwp) THEN
         WRITE(numout,*) '      High-mode wave-breaking energy:             ', zbot * 1.e-12_wp, 'TW'
         WRITE(numout,*) '      Pycnocline-intensifed wave-breaking energy: ', zpyc * 1.e-12_wp, 'TW'
         WRITE(numout,*) '      Critical slope wave-breaking energy:        ', zcri * 1.e-12_wp, 'TW'
      ENDIF
      !
   END SUBROUTINE zdf_iwm_init

   !!======================================================================
END MODULE zdfiwm
