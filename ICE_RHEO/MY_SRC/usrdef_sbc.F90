MODULE usrdef_sbc
   !!======================================================================
   !!                       ***  MODULE  usrdef_sbc  ***
   !! 
   !!                      ===  ICE_RHEO configuration  ===
   !!
   !! User defined :   surface forcing of a user configuration
   !!======================================================================
   !! History :  4.0   ! 2016-03  (S. Flavoni, G. Madec)  user defined interface
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_sbc    : user defined surface bounday conditions in ICE_RHEO case
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE phycst          ! physical constants
   USE ice, ONLY       : at_i_b, a_i_b, at_i, u_ice, v_ice
   USE icethd_dh       ! for CALL ice_thd_snwblow
   USE sbc_phy, ONLY : pp_cldf
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined) 

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usrdef_sbc_oce      ! routine called by sbcmod.F90 for sbc ocean
   PUBLIC   usrdef_sbc_ice_tau  ! routine called by icestp.F90 for ice dynamics
   PUBLIC   usrdef_sbc_ice_flx  ! routine called by icestp.F90 for ice thermo

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_sbc.F90 10074 2018-08-28 16:15:49Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usrdef_sbc_oce( kt, Kbb )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE usr_def_sbc  ***
      !!              
      !! ** Purpose :   provide at each time-step the surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   all 0 fields, for ICE_RHEO case
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - set to ZERO all the ocean surface boundary condition, i.e.   
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, sfx
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      INTEGER, INTENT(in) ::   Kbb  ! ocean time index
      INTEGER ::    ij0, ij1, ii0, ii1, jj, ji   ! loop indices
      REAL(wp) ::   zrhoco          ! ocean density and drag coefficient product
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         !
         !IF(lwp)   WRITE(numout,*)' usrdef_sbc_oce : ICE_RHEO case: ocean boudary conditions'

         utau(:,:) = 0._wp
         utau(:,:) = 0._wp

         !ij0 = 1   ;   ij1 = 25                       ! set boundary condition
         !ii0 = 975   ;   ii1 = 1000
         !DO jj = mj0(ij0), mj1(ij1)
         !   DO ji = mi0(ii0), mi1(ii1)
         !      utau(ji,jj) = -utau_ice(ji,jj) 
         !      vtau(ji,jj) = -vtau_ice(ji,jj)
         !   END DO
         !END DO

         taum(:,:) = 0._wp   ! assume these are not used
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
      INTEGER ::   jj, ji   ! loop indices

      REAL(wp) ::   zwndi_f , zwndj_f, zwnorm_f       ! relative wind module and components at F-point
      REAL(wp) ::   zwndi_t , zwndj_t                 ! relative wind components at T-point
      REAL(wp), DIMENSION(jpi,jpj) ::   windu, windv  ! wind components (idealised forcing)
      REAL(wp), PARAMETER ::   r_vfac = 1._wp         ! relative velocity (make 0 for absolute velocity) 
      REAL(wp), PARAMETER ::   Rwind = -0.8_wp        ! ratio of wind components
      REAL(wp), PARAMETER ::   Umax = 15._wp          ! maximum wind speed (m/s)
      REAL(wp), PARAMETER ::   d = 2000._wp           ! size of the domain (km)
      REAL(wp), PARAMETER ::   res = 2._wp            ! gridcell size
      REAL(wp), PARAMETER ::   zrhoa  = 1.22          ! Air density kg/m3
      REAL(wp), PARAMETER ::   Cd_atm =  1.4e-3       ! transfer coefficient over ice
      !!---------------------------------------------------------------------
      ! extra code for test case
      IF( kt==nit000 .AND. lwp)   WRITE(numout,*)' usrdef_sbc_ice : ICE_RHEO case: analytical stress forcing'

      DO_2D( 0, 0, 0, 0 )
         ! wind spins up over 6 hours, factor 1000 to balance the units
         windu(ji,jj) = Umax/sqrt(d*1000)*(d-2*mig(ji)*res)/((d-2*mig(ji)*res)**2+(d-2*mjg(jj)*res)**2*Rwind**2)**(1/4)*min(kt*30./21600,1.)
         windv(ji,jj) = Umax/sqrt(d*1000)*(d-2*mjg(jj)*res)/((d-2*mig(ji)*res)**2+(d-2*mjg(jj)*res)**2*Rwind**2)**(1/4)*Rwind*min(kt*30./21600,1.)
      END_2D
      CALL lbc_lnk( 'usrdef_sbc', windu, 'U', -1., windv, 'V', -1. )

      wndm_ice(:,:) = 0._wp      !!gm brutal....

      ! ------------------------------------------------------------ !
      !    Wind module relative to the moving ice ( U10m - U_ice )   !
      ! ------------------------------------------------------------ !
      ! C-grid ice dynamics :   U & V-points (same as ocean)
      DO_2D( 0, 0, 0, 0 )
         zwndi_t = (  windu(ji,jj) - r_vfac * 0.5 * ( u_ice(ji-1,jj  ) + u_ice(ji,jj) )  )
         zwndj_t = (  windv(ji,jj) - r_vfac * 0.5 * ( v_ice(ji,jj-1) + v_ice(ji,jj) )  )
         wndm_ice(ji,jj) = SQRT( zwndi_t * zwndi_t + zwndj_t * zwndj_t ) * tmask(ji,jj,1)
      END_2D
      CALL lbc_lnk( 'usrdef_sbc', wndm_ice, 'T',  1. )

      !!gm brutal....
      utau_ice  (:,:) = 0._wp
      vtau_ice  (:,:) = 0._wp
      !!gm end

      ! ------------------------------------------------------------ !
      !    Wind stress relative to the moving ice ( U10m - U_ice )   !
      ! ------------------------------------------------------------ !
      ! C-grid ice dynamics :   U & V-points (same as ocean)
      DO_2D( 0, 0, 0, 0 )
         utau_ice(ji,jj) = 0.5 * zrhoa * Cd_atm * ( wndm_ice(ji+1,jj  ) + wndm_ice(ji,jj) )            &
            &          * ( 0.5 * (windu(ji+1,jj) + windu(ji,jj) ) - r_vfac * u_ice(ji,jj) )
         vtau_ice(ji,jj) = 0.5 * zrhoa * Cd_atm * ( wndm_ice(ji,jj+1  ) + wndm_ice(ji,jj) )            &
            &          * ( 0.5 * (windv(ji,jj+1) + windv(ji,jj) ) - r_vfac * v_ice(ji,jj) )
      END_2D
      CALL lbc_lnk( 'usrdef_sbc', utau_ice, 'U', -1., vtau_ice, 'V', -1. )
      !
   END SUBROUTINE usrdef_sbc_ice_tau

   SUBROUTINE usrdef_sbc_ice_flx( kt, phs, phi )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE usrdef_sbc_ice_flx  ***
      !!
      !! ** Purpose :   provide the surface boundary (flux) condition over sea-ice
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   phs    ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   phi    ! ice thickness
      !!
      REAL(wp) ::   zfr1, zfr2                 ! local variables
      REAL(wp), DIMENSION(jpi,jpj) ::   zsnw   ! snw distribution after wind blowing
      !!---------------------------------------------------------------------
      !
      IF( kt==nit000 .AND. lwp)   WRITE(numout,*)' usrdef_sbc_ice : ICE_RHEO case: NO flux forcing'
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

      ! --- shortwave radiation transmitted below the surface (W/m2, see Grenfell Maykut 77) --- !
      zfr1 = ( 0.18 * ( 1.0 - pp_cldf ) + 0.35 * pp_cldf )            ! transmission when hi>10cm
      zfr2 = ( 0.82 * ( 1.0 - pp_cldf ) + 0.65 * pp_cldf )            ! zfr2 such that zfr1 + zfr2 to equal 1
      !
      WHERE    ( phs(:,:,:) <= 0._wp .AND. phi(:,:,:) <  0.1_wp )       ! linear decrease from hi=0 to 10cm  
         qtr_ice_top(:,:,:) = qsr_ice(:,:,:) * ( zfr1 + zfr2 * ( 1._wp - phi(:,:,:) * 10._wp ) )
      ELSEWHERE( phs(:,:,:) <= 0._wp .AND. phi(:,:,:) >= 0.1_wp )       ! constant (zfr1) when hi>10cm
         qtr_ice_top(:,:,:) = qsr_ice(:,:,:) * zfr1
      ELSEWHERE                                                         ! zero when hs>0
         qtr_ice_top(:,:,:) = 0._wp 
      END WHERE
          
   END SUBROUTINE usrdef_sbc_ice_flx


   !!======================================================================
END MODULE usrdef_sbc
