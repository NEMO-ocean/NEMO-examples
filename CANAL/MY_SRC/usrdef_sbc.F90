MODULE usrdef_sbc
   !!======================================================================
   !!                       ***  MODULE  usrdef_sbc  ***
   !! 
   !!                      ===  CANAL configuration  ===
   !!
   !! User defined :   surface forcing of a user configuration
   !!======================================================================
   !! History :  4.0   ! 2017-11  (J.Chanut)  user defined interface
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_sbc    : user defined surface bounday conditions in OVERFLOW case
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE phycst          ! physical constants
   USE usrdef_nam, ONLY : rn_u10, rn_uofac, rn_windszy, rn_windszx 
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined) 

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usrdef_sbc_oce      ! routine called in sbcmod module
   PUBLIC   usrdef_sbc_ice_tau  ! routine called by icestp.F90 for ice dynamics
   PUBLIC   usrdef_sbc_ice_flx  ! routine called by icestp.F90 for ice thermo

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_sbc.F90 13472 2020-09-16 13:05:19Z smasson $
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
      !! ** Method  :   all 0 fields, for CANAL case
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - set to ZERO all the ocean surface boundary condition, i.e.   
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, sfx
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt        ! ocean time step
      INTEGER, INTENT(in) ::   Kbb       ! ocean time index
      INTEGER  ::   ji, jj               ! dummy loop indices
      REAL(wp) :: zrhoair = 1.22     ! approximate air density [Kg/m3]
      REAL(wp) :: zcd = 1.13e-3      ! approximate drag coefficient
      REAL(wp) :: zrhocd             ! Rho * Cd
      REAL(wp), DIMENSION(jpi,jpj) :: zwndrel   ! relative wind
      !!---------------------------------------------------------------------
      !
      zrhocd = zrhoair * zcd
      
      IF( kt == nit000 ) THEN
         !
         IF(lwp) WRITE(numout,*)' usr_sbc : EW_CANAL case: surface forcing'
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
      ENDIF

      IF( rn_u10 /= 0. .AND. rn_windszy > 0. ) THEN
         IF( nyear == 1 .AND. nmonth == 1 .AND. nday <= 10 ) THEN
            WHERE( ABS(gphit) <= rn_windszy/2. .AND. ABS(glamt) <= rn_windszx/2. ) utau(:,:) = zrhocd * rn_u10 * rn_u10
         ELSE
            utau(:,:) = 0.
         ENDIF
      ENDIF

      IF( rn_uofac /= 0. ) THEN
         
         WHERE( ABS(gphit) <= rn_windszy/2. )
            zwndrel(:,:) = rn_u10 - rn_uofac * uu(:,:,1,Kbb)
         ELSEWHERE
            zwndrel(:,:) =        - rn_uofac * uu(:,:,1,Kbb)
         END WHERE
         utau(:,:) = zrhocd * zwndrel(:,:) * zwndrel(:,:)

         zwndrel(:,:) = - rn_uofac * vv(:,:,1,Kbb)
         vtau(:,:) = zrhocd * zwndrel(:,:) * zwndrel(:,:)

      ENDIF
      !
   END SUBROUTINE usrdef_sbc_oce

   SUBROUTINE usrdef_sbc_ice_tau( kt )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
   END SUBROUTINE usrdef_sbc_ice_tau


   SUBROUTINE usrdef_sbc_ice_flx( kt, phs, phi )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phs    ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phi    ! ice thickness
   END SUBROUTINE usrdef_sbc_ice_flx

   !!======================================================================
END MODULE usrdef_sbc
