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

   PUBLIC   usrdef_sbc_oce      ! routine called in sbcmod module
   PUBLIC   usrdef_sbc_ice_tau  ! routine called by sbcice_lim.F90 for ice dynamics
   PUBLIC   usrdef_sbc_ice_flx  ! routine called by sbcice_lim.F90 for ice thermo

   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2016)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usrdef_sbc_oce( kt )
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
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! 2D workspace
      INTEGER  ::   ji, jj
      !!---------------------------------------------------------------------
#if defined key_si3
      IF( kt==nit000 .AND. lwp)   WRITE(numout,*)' usrdef_sbc_ice : BENCH case: constant stress forcing'
      !
      ! define unique value on each point. z2d ranging from 0.05 to -0.05
      DO jj = 1, jpj
         DO ji = 1, jpi
            z2d(ji,jj) = 0.1 * ( 0.5 - REAL( nimpp + ji - 1 + ( njmpp + jj - 2 ) * jpiglo, wp ) / REAL( jpiglo * jpjglo, wp ) )
         ENDDO
      ENDDO
      utau_ice(:,:) = 0.1_wp +  z2d(:,:)
      vtau_ice(:,:) = 0.1_wp +  z2d(:,:)

      CALL lbc_lnk_multi( 'usrdef_sbc', utau_ice, 'U', -1., vtau_ice, 'V', -1. )
#endif
      !
   END SUBROUTINE usrdef_sbc_ice_tau

   
   SUBROUTINE usrdef_sbc_ice_flx( kt, phs, phi )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE usrdef_sbc_ice_flx  ***
      !!
      !! ** Purpose :   provide the surface boundary (flux) condition over
      !sea-ice
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phs    ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phi    ! ice thickness
      !!
      REAL(wp) ::   zfr1, zfr2                 ! local variables
      REAL(wp), DIMENSION(jpi,jpj) ::   zsnw   ! snw distribution after wind blowing
      !!---------------------------------------------------------------------
      !
#if defined key_si3
      IF( kt==nit000 .AND. lwp)   WRITE(numout,*)' usrdef_sbc_ice : BENCH case: NO flux forcing'
      !
      ! ocean variables (renaming)
      emp_oce (:,:)   = 0._wp   ! uniform value for freshwater budget (E-P)
      qsr_oce (:,:)   = 0._wp   ! uniform value for     solar radiation
      qns_oce (:,:)   = 0._wp   ! uniform value for non-solar radiation

      ! ice variables
      alb_ice (:,:,:) = 0.7_wp  ! useless
      qsr_ice (:,:,:) = 0._wp   ! uniform value for     solar radiation
      qns_ice (:,:,:) = 0._wp   ! uniform value for non-solar radiation
      sprecip (:,:)   = 0._wp   ! uniform value for snow precip
      evap_ice(:,:,:) = 0._wp   ! uniform value for sublimation

      ! ice fields deduced from above
      zsnw(:,:) = 1._wp
      !!CALL lim_thd_snwblow( at_i_b, zsnw )  ! snow distribution over ice after
      !wind blowing 
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

      ! --- shortwave radiation transmitted below the surface (W/m2, see Grenfell Maykut 77) --- !
      zfr1 = ( 0.18 * ( 1.0 - cldf_ice ) + 0.35 * cldf_ice )            ! transmission when hi>10cm
      zfr2 = ( 0.82 * ( 1.0 - cldf_ice ) + 0.65 * cldf_ice )            ! zfr2 such that zfr1 + zfr2 to equal 1
      !
      WHERE    ( phs(:,:,:) <= 0._wp .AND. phi(:,:,:) <  0.1_wp )       ! linear decrease from hi=0 to 10cm  
         qtr_ice_top(:,:,:) = qsr_ice(:,:,:) * ( zfr1 + zfr2 * ( 1._wp - phi(:,:,:) * 10._wp ) )
      ELSEWHERE( phs(:,:,:) <= 0._wp .AND. phi(:,:,:) >= 0.1_wp )       ! constant (zfr1) when hi>10cm
         qtr_ice_top(:,:,:) = qsr_ice(:,:,:) * zfr1
      ELSEWHERE                                                         ! zero when hs>0
         qtr_ice_top(:,:,:) = 0._wp 
      END WHERE
#endif

   END SUBROUTINE usrdef_sbc_ice_flx

   !!======================================================================
END MODULE usrdef_sbc
