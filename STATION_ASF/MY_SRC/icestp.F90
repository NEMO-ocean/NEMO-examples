MODULE icestp
   !!======================================================================
   !!                       ***  MODULE  icestp  ***
   !! sea ice : Master routine for all the sea ice model
   !!=====================================================================
   !!
   !! The sea ice model SI3 (Sea Ice modelling Integrated Initiative),
   !!                        aka Sea Ice cube for its nickname
   !!
   !!    is originally based on LIM3, developed in Louvain-la-Neuve by:
   !!       * Martin Vancoppenolle (UCL-ASTR, Belgium)
   !!       * Sylvain Bouillon (UCL-ASTR, Belgium)
   !!       * Miguel Angel Morales Maqueda (NOC-L, UK)
   !!      thanks to valuable earlier work by
   !!       * Thierry Fichefet
   !!       * Hugues Goosse
   !!      thanks also to the following persons who contributed
   !!       * Gurvan Madec, Claude Talandier, Christian Ethe (LOCEAN, France)
   !!       * Xavier Fettweis (UCL-ASTR), Ralph Timmermann (AWI, Germany)
   !!       * Bill Lipscomb (LANL), Cecilia Bitz (UWa) and Elisabeth Hunke (LANL), USA.
   !!
   !! SI3 has been made possible by a handful of persons who met as working group
   !!      (from France, Belgium, UK and Italy)
   !!    * Clement Rousset, Martin Vancoppenolle & Gurvan Madec (LOCEAN, France)
   !!    * Matthieu Chevalier & David Salas (Meteo France, France)
   !!    * Gilles Garric (Mercator Ocean, France)
   !!    * Thierry Fichefet & Francois Massonnet (UCL, Belgium)
   !!    * Ed Blockley & Jeff Ridley (Met Office, UK)
   !!    * Danny Feltham & David Schroeder (CPOM, UK)
   !!    * Yevgeny Aksenov (NOC, UK)
   !!    * Paul Holland (BAS, UK)
   !!    * Dorotea Iovino (CMCC, Italy)
   !!======================================================================
   !! History :  4.0  !  2018     (C. Rousset)      Original code SI3
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_stp       : sea-ice model time-stepping and update ocean SBC over ice-covered area
   !!   ice_init      : initialize sea-ice
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE c1d            ! 1D vertical configuration
   USE ice            ! sea-ice: variables
   USE ice1D          ! sea-ice: thermodynamical 1D variables
   !
   USE phycst         ! Define parameters for the routines
   USE eosbn2         ! equation of state
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE sbc_ice        ! Surface boundary condition: ice   fields
   !
   USE icesbc         ! sea-ice: Surface boundary conditions
   USE iceupdate      ! sea-ice: sea surface boundary condition update
   USE icewri         ! sea-ice: outputs
   !
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE timing         ! Timing
   USE prtctl         ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_stp    ! called by sbcmod.F90
   PUBLIC   ice_init   ! called by sbcmod.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icestp.F90 13655 2020-10-21 14:15:13Z laurent $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_stp( kt, Kbb, Kmm, ksbc )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE ice_stp  ***
      !!
      !! ** Purpose :   sea-ice model time-stepping and update ocean surface
      !!              boundary condition over ice-covered area
      !!
      !! ** Method  :   ice model time stepping
      !!              - call the ice dynamics routine
      !!              - call the ice advection/diffusion routine
      !!              - call the ice thermodynamics routine
      !!              - call the routine that computes mass and
      !!                heat fluxes at the ice/ocean interface
      !!              - save the outputs
      !!              - save the outputs for restart when necessary
      !!
      !! ** Action  : - time evolution of the LIM sea-ice model
      !!              - update all sbc variables below sea-ice:
      !!                utau, vtau, taum, wndm, qns , qsr, emp , sfx
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm ! ocean time level indices
      INTEGER, INTENT(in) ::   ksbc     ! flux formulation (user defined, bulk, or Pure Coupled)
      !
      INTEGER ::   jl   ! dummy loop index
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('ice_stp')
      !
      !                                      !-----------------------!
      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN   ! --- Ice time step --- !
         !                                   !-----------------------!
         !
         kt_ice = kt                              ! -- Ice model time step
         !
         u_oce(:,:) = ssu_m(:,:)                  ! -- mean surface ocean current
         v_oce(:,:) = ssv_m(:,:)
         !
         CALL eos_fzp( sss_m(:,:) , t_bo(:,:) )   ! -- freezing temperature [Kelvin] (set to rt0 over land)
         t_bo(:,:) = ( t_bo(:,:) + rt0 ) * tmask(:,:,1) + rt0 * ( 1._wp - tmask(:,:,1) )
         !
         !
         !------------------------------------------------!
         ! --- Dynamical coupling with the atmosphere --- !
         !------------------------------------------------!
         ! It provides the following fields used in sea ice model:
         !    utau_ice, vtau_ice = surface ice stress [N/m2]
         !------------------------------------------------!
                                        CALL ice_sbc_tau( kt, ksbc, utau_ice, vtau_ice )
         !
         !------------------------------------------------------!
         ! --- Thermodynamical coupling with the atmosphere --- !
         !------------------------------------------------------!
         ! It provides the following fields used in sea ice model:
         !    emp_oce , emp_ice    = E-P over ocean and sea ice                    [Kg/m2/s]
         !    sprecip              = solid precipitation                           [Kg/m2/s]
         !    evap_ice             = sublimation                                   [Kg/m2/s]
         !    qsr_tot , qns_tot    = solar & non solar heat flux (total)           [W/m2]
         !    qsr_ice , qns_ice    = solar & non solar heat flux over ice          [W/m2]
         !    dqns_ice             = non solar  heat sensistivity                  [W/m2]
         !    qemp_oce, qemp_ice,  = sensible heat (associated with evap & precip) [W/m2]
         !    qprec_ice, qevap_ice
         !------------------------------------------------------!
                                        CALL ice_sbc_flx( kt, ksbc )
         !----------------------------!
         ! --- ice thermodynamics --- !
         !----------------------------!

         fr_i  (:,:)   = at_i(:,:) !#LB...

         !!#LB: lines stolen from "ice_update_flx()@iceupdate.F90"
         !!=============================================================

         ! --- case we bypass ice thermodynamics --- !
         !IF( .NOT. ln_icethd ) THEN   ! we suppose ice is impermeable => ocean is isolated from atmosphere
         qt_atm_oi  (:,:)   = ( 1._wp - at_i_b(:,:) ) * ( qns_oce(:,:) + qsr_oce(:,:) ) + qemp_oce(:,:)
         qt_oce_ai  (:,:)   = ( 1._wp - at_i_b(:,:) ) *   qns_oce(:,:)                  + qemp_oce(:,:)
         emp_ice    (:,:)   = 0._wp
         qemp_ice   (:,:)   = 0._wp
         qevap_ice  (:,:,:) = 0._wp
         !ENDIF

         ! output all fluxes
         !------------------

         ! --- mass fluxes [kg/m2/s] --- !
         IF( iom_use('emp_oce'    ) )   CALL iom_put( 'emp_oce', emp_oce )   ! emp over ocean (taking into account the snow blown away from the ice)
         IF( iom_use('emp_ice'    ) )   CALL iom_put( 'emp_ice', emp_ice )   ! emp over ice   (taking into account the snow blown away from the ice)
         ! --- heat fluxes [W/m2] --- !
         IF( iom_use('qsr_oce'    ) )   CALL iom_put( 'qsr_oce'    , qsr_oce * ( 1._wp - at_i_b )                               )   !     solar flux at ocean surface
         IF( iom_use('qns_oce'    ) )   CALL iom_put( 'qns_oce'    , qns_oce * ( 1._wp - at_i_b ) + qemp_oce                    )   ! non-solar flux at ocean surface
         IF( iom_use('qns_ice'    ) )   CALL iom_put( 'qns_ice'    , SUM( qns_ice * a_i_b, dim=3 ) + qemp_ice                   )   ! non-solar flux at ice surface
         IF( iom_use('qsr_ice'    ) )   CALL iom_put( 'qsr_ice'    , SUM( qsr_ice * a_i_b, dim=3 )                              )
         IF( iom_use('qt_oce'     ) )   CALL iom_put( 'qt_oce'     ,      ( qsr_oce + qns_oce ) * ( 1._wp - at_i_b ) + qemp_oce )
         IF( iom_use('qt_ice'     ) )   CALL iom_put( 'qt_ice'     , SUM( ( qns_ice + qsr_ice ) * a_i_b, dim=3 )     + qemp_ice )
         IF( iom_use('qemp_oce'   ) )   CALL iom_put( 'qemp_oce'   , qemp_oce                                                   )   ! Downward Heat Flux from E-P over ocean
         IF( iom_use('qemp_ice'   ) )   CALL iom_put( 'qemp_ice'   , qemp_ice                                                   )   ! Downward Heat Flux from E-P over ice


         !!=============================================================

         !
                                        CALL ice_wri( kt )            ! -- Ice outputs
         !
         !IF( lrst_ice )                 CALL ice_rst_write( kt )      ! -- Ice restart file
         !
         !IF( ln_icectl )                CALL ice_ctl( kt )            ! -- Control checks
         !
      ENDIF   ! End sea-ice time step only

      !-------------------------!
      ! --- Ocean time step --- !
      !-------------------------!
      !
      IF( ln_timing )   CALL timing_stop('ice_stp')
      !
   END SUBROUTINE ice_stp


   SUBROUTINE ice_init( Kbb, Kmm, Kaa )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ice_init  ***
      !!
      !! ** purpose :   Initialize sea-ice parameters
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: Kbb, Kmm, Kaa
      !
      INTEGER ::   ierr
      !!----------------------------------------------------------------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'Sea Ice Model: SI3 (Sea Ice modelling Integrated Initiative)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'ice_init: Arrays allocation & Initialization of all routines & init state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~'
      !
      !                                ! Load the reference and configuration namelist files and open namelist output file
      CALL load_nml( numnam_ice_ref, 'namelist_ice_ref',    numout, lwm )
      CALL load_nml( numnam_ice_cfg, 'namelist_ice_cfg',    numout, lwm )
      IF(lwm) CALL ctl_opn( numoni , 'output.namelist.ice', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, 1 )
      !
      CALL par_init                ! set some ice run parameters
      !
      !                                ! Allocate the ice arrays (sbc_ice already allocated in sbc_init)
      ierr =        ice_alloc        ()      ! ice variables
      ierr = ierr + sbc_ice_alloc    ()      ! surface boundary conditions
      ierr = ierr + ice1D_alloc      ()      ! thermodynamics
      !
      CALL mpp_sum( 'icestp', ierr )
      IF( ierr /= 0 )   CALL ctl_stop('STOP', 'ice_init : unable to allocate ice arrays')
      !
      !
      CALL ice_sbc_init                ! set ice-ocean and ice-atm. coupling parameters
      !
      fr_i  (:,:)   = at_i(:,:)        ! initialisation of sea-ice fraction
      !
      IF( ln_rstart )   CALL iom_close( numrir )  ! close input ice restart file
      !
   END SUBROUTINE ice_init


   SUBROUTINE par_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE par_init ***
      !!
      !! ** Purpose :   Definition generic parameters for ice model
      !!
      !! ** Method  :   Read namelist and check the parameter
      !!                values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist nampar
      !!-------------------------------------------------------------------
      INTEGER  ::   ios                 ! Local integer
      !!
      NAMELIST/nampar/ jpl, nlay_i, nlay_s, ln_virtual_itd, ln_icedyn, ln_icethd, rn_amax_n, rn_amax_s,  &
         &             cn_icerst_in, cn_icerst_indir, cn_icerst_out, cn_icerst_outdir
      !!-------------------------------------------------------------------
      !
      READ  ( numnam_ice_ref, nampar, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampar in reference namelist' )
      READ  ( numnam_ice_cfg, nampar, IOSTAT = ios, ERR = 902 )
902   IF( ios > 0 )   CALL ctl_nam ( ios , 'nampar in configuration namelist' )
      IF(lwm) WRITE( numoni, nampar )
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) '   par_init: ice parameters shared among all the routines'
         WRITE(numout,*) '   ~~~~~~~~'
         WRITE(numout,*) '      Namelist nampar: '
         WRITE(numout,*) '         number of ice  categories                           jpl       = ', jpl
         WRITE(numout,*) '         number of ice  layers                               nlay_i    = ', nlay_i
         WRITE(numout,*) '         number of snow layers                               nlay_s    = ', nlay_s
         WRITE(numout,*) '         virtual ITD param for jpl=1 (T) or not (F)     ln_virtual_itd = ', ln_virtual_itd
         WRITE(numout,*) '         Ice dynamics       (T) or not (F)                   ln_icedyn = ', ln_icedyn
         WRITE(numout,*) '         Ice thermodynamics (T) or not (F)                   ln_icethd = ', ln_icethd
         WRITE(numout,*) '         maximum ice concentration for NH                              = ', rn_amax_n
         WRITE(numout,*) '         maximum ice concentration for SH                              = ', rn_amax_s
      ENDIF
      !                                        !--- change max ice concentration for roundoff errors
      rn_amax_n = MIN( rn_amax_n, 1._wp - epsi10 )
      rn_amax_s = MIN( rn_amax_s, 1._wp - epsi10 )
      !                                        !--- check consistency
      IF ( jpl > 1 .AND. ln_virtual_itd ) THEN
         ln_virtual_itd = .FALSE.
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ln_virtual_itd forced to false as jpl>1, no need with multiple categories to emulate them'
      ENDIF
      !
      IF( ln_cpl .AND. nn_cats_cpl /= 1 .AND. nn_cats_cpl /= jpl ) THEN
         CALL ctl_stop( 'STOP', 'par_init: in coupled mode, nn_cats_cpl should be either 1 or jpl' )
      ENDIF
      !
      rDt_ice   = REAL(nn_fsbc) * rn_Dt          !--- sea-ice timestep and its inverse
      r1_Dt_ice = 1._wp / rDt_ice
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '      ice timestep rDt_ice = nn_fsbc*rn_Dt = ', rDt_ice
      !
      r1_nlay_i = 1._wp / REAL( nlay_i, wp )   !--- inverse of nlay_i and nlay_s
      r1_nlay_s = 1._wp / REAL( nlay_s, wp )
      !
   END SUBROUTINE par_init

#else
   !!----------------------------------------------------------------------
   !!   Default option           Dummy module         NO SI3 sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE ice_stp ( kt, ksbc )     ! Dummy routine
      INTEGER, INTENT(in) ::   kt, ksbc
      WRITE(*,*) 'ice_stp: You should not have seen this print! error?', kt
   END SUBROUTINE ice_stp
   SUBROUTINE ice_init                 ! Dummy routine
   END SUBROUTINE ice_init
#endif

   !!======================================================================
END MODULE icestp
