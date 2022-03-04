MODULE usrdef_sbc
   !!======================================================================
   !!                       ***  MODULE  usrdef_sbc  ***
   !! 
   !!                      ===  BENCH configuration  ===
   !!
   !! User defined :   surface forcing of a user configuration
   !!======================================================================
   !! History :  4.0   ! 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_sbc    : user defined surface bounday conditions in BENCH case
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE dom_oce        
   USE oce             ! ocean dynamics and tracers
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbc_ice         ! Surface boundary condition: ocean fields
   USE in_out_manager  ! I/O manager
   USE phycst          ! physical constants
   USE lib_mpp         ! MPP library
   USE lbclnk          ! lateral boundary conditions - mpp exchanges

#if defined key_si3
   USE ice, ONLY       : at_i_b, a_i_b
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usrdef_sbc_oce      ! routine called by sbcmod.F90 for sbc ocean
   PUBLIC   usrdef_sbc_ice_tau  ! routine called by icestp.F90 for ice dynamics
   PUBLIC   usrdef_sbc_ice_flx  ! routine called by icestp.F90 for ice thermo

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2016)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usrdef_sbc_oce( kt, Kbb )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE usr_def_sbc  ***
      !!              
      !! ** Purpose :   provide at each time-step the surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   all 0 fields, for BENCH case
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - set to ZERO all the ocean surface boundary condition, i.e.   
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, sfx
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      INTEGER, INTENT(in) ::   Kbb  ! ocean time index
      !!---------------------------------------------------------------------
      !     
      IF( kt == nit000 ) THEN
         !
         IF(lwp) WRITE(numout,*)' usr_sbc : BENCH case: surface forcing'
         IF(lwp) WRITE(numout,*)' ~~~~~~~~~~~   vtau = taum = wndm = qns = qsr = emp = sfx = 0'
         !
         utau(:,:) = 0._wp
         vtau(:,:) = 0._wp
         taum(:,:) = 0._wp
         wndm(:,:) = 0._wp
         !
         emp (:,:) = 0._wp
         sfx (:,:) = 0._wp
         qns (:,:) = 0._wp
         qsr (:,:) = 0._wp
         !
         utau_b(:,:) = 0._wp 
         vtau_b(:,:) = 0._wp
         emp_b (:,:) = 0._wp
         sfx_b (:,:) = 0._wp
         qns_b (:,:) = 0._wp
         !
      ENDIF
      !
   END SUBROUTINE usrdef_sbc_oce

   
   SUBROUTINE usrdef_sbc_ice_tau( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE usrdef_sbc_ice_tau  ***
      !!
      !! ** Purpose :   provide the surface boundary (momentum) condition over
      !sea-ice
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !
      REAL(wp) ::   zztmp
      INTEGER  ::   ji, jj
      !!---------------------------------------------------------------------
#if defined key_si3
      IF( kt==nit000 .AND. lwp)   WRITE(numout,*)' usrdef_sbc_ice : BENCH case: constant stress forcing'
      !
      ! define unique value on each point. z2d ranging from 0.05 to -0.05
      !
      DO_2D( 0, 0, 0, 0 )
         zztmp = 0.1 * ( 0.5 - REAL( mig0(ji) + (mjg0(jj)-1) * Ni0glo, wp ) / REAL( Ni0glo * Nj0glo, wp ) )
         utau_ice(ji,jj) = 0.1_wp + zztmp
         vtau_ice(ji,jj) = 0.1_wp + zztmp
      END_2D

      CALL lbc_lnk( 'usrdef_sbc', utau_ice, 'U', -1., vtau_ice, 'V', -1. )
#endif
      !
   END SUBROUTINE usrdef_sbc_ice_tau

   
   SUBROUTINE usrdef_sbc_ice_flx( kt, phs, phi )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE usrdef_sbc_ice_flx  ***
      !!
      !! ** Purpose :   provide the surface boundary (flux) condition over sea-ice
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phs    ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phi    ! ice thickness
      !!
      REAL(wp), DIMENSION(jpi,jpj) ::   zsnw   ! snw distribution after wind blowing
      !!---------------------------------------------------------------------
#if defined key_si3
      !
      IF( kt==nit000 .AND. lwp)   WRITE(numout,*)' usrdef_sbc_ice : BENCH case: NO flux forcing'
      !
      ! ocean variables (renaming)
      emp_oce (:,:)   = 0._wp   ! uniform value for freshwater budget (E-P)
      qsr_oce (:,:)   = 0._wp   ! uniform value for     solar radiation
      qns_oce (:,:)   = 0._wp   ! uniform value for non-solar heat flux

      ! ice variables
      alb_ice (:,:,:) = 0.7_wp  ! useless
      qsr_ice (:,:,:) = 0._wp   ! uniform value for     solar radiation
      qns_ice (:,:,:) = 0._wp   ! uniform value for non-solar heat flux
      dqns_ice(:,:,:) = 0._wp   ! uniform value for non solar heat flux sensitivity for ice
      sprecip (:,:)   = 0._wp   ! uniform value for snow precip
      evap_ice(:,:,:) = 0._wp   ! uniform value for sublimation

      ! ice fields deduced from above
      zsnw(:,:) = 1._wp
      !!CALL lim_thd_snwblow( at_i_b, zsnw )  ! snow distribution over ice after wind blowing 
      emp_ice  (:,:)   = SUM( a_i_b(:,:,:) * evap_ice(:,:,:), dim=3 ) - sprecip(:,:) * zsnw(:,:)
      emp_oce  (:,:)   = emp_oce(:,:) - sprecip(:,:) * (1._wp - zsnw(:,:) )
      qevap_ice(:,:,:) =   0._wp
      qprec_ice(:,:)   =   rhos * ( sst_m(:,:) * rcpi - rLfus ) * tmask(:,:,1) !  in J/m3
      qemp_oce (:,:)   = - emp_oce(:,:) * sst_m(:,:) * rcp
      qemp_ice (:,:)   =   sprecip(:,:) * zsnw * ( sst_m(:,:) * rcpi - rLfus ) * tmask(:,:,1) ! solid precip (only)

      ! total fluxes
      emp_tot (:,:) = emp_ice  + emp_oce
      qns_tot (:,:) = at_i_b(:,:) * qns_oce(:,:) + SUM( a_i_b(:,:,:) * qns_ice(:,:,:), dim=3 ) + qemp_ice(:,:) + qemp_oce(:,:)
      qsr_tot (:,:) = at_i_b(:,:) * qsr_oce(:,:) + SUM( a_i_b(:,:,:) * qsr_ice(:,:,:), dim=3 )

      ! --- shortwave radiation transmitted thru the surface scattering layer (W/m2) --- !
      qtr_ice_top(:,:,:) = 0._wp
#endif

   END SUBROUTINE usrdef_sbc_ice_flx

   !!======================================================================
END MODULE usrdef_sbc
