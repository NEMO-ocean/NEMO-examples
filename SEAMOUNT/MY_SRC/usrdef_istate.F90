MODULE usrdef_istate
   !!======================================================================
   !!                   ***  MODULE  usrdef_istate   ***
   !!
   !!                     ===  GYRE configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-03  (S. Flavoni) Original code
   !!                 ! 2020-11  (S. Techene, G. Madec) separate tsuv from ssh
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE dom_oce , ONLY : mi0, mig, mjg, glamt, gphit, ht_0
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE eosbn2, ONLY: rn_a0, rho0 
   USE usrdef_nam
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate       ! called in istate.F90
   PUBLIC   usr_def_istate_ssh   ! called by domqco.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_istate.F90 14053 2020-12-03 13:48:38Z techene $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here GYRE configuration example : (double gyre with rotated domain)
      !!
      !! ** Method  : - set temprature field
      !!              - set salinity   field
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
      !
      INTEGER :: ji, jj, jk  ! dummy loop indices
      !!----------------------------------------------------------------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : SEAMOUNT_TEST_CASE configuration, analytical definition of initial state.'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   Ocean at rest, with an initial density profile                           '
      IF(lwp) WRITE(numout,*) '                                rho_ini = rn_drho * exp(z / rn_delta)'
      IF(lwp) WRITE(numout,*) '                 This test case uses a linear EOS only function of temperature.'
      !
      IF(lwp) WRITE(numout,*) ''
      IF(lwp) WRITE(numout,*) '                 Initial condition settings:'
      IF(lwp) WRITE(numout,*) '                     *) Surface-bottom density difference    rn_drho  = ', rn_drho
      IF(lwp) WRITE(numout,*) '                     *) Steepness of the stratification      rn_delta = ', rn_delta
      IF(lwp) WRITE(numout,*) '                     *) Estimated Burger number              rn_Snum  = ', rn_Snum
      !
      pu  (:,:,:) = 0._wp        ! ocean at rest
      pv  (:,:,:) = 0._wp
      !
      !                          ! T & S profiles
      T0 = 10._wp
      S0 = 35._wp
      DO jk = 1, jpk
         rho(:,:,k) = rho0 - rn_drho * EXP(pdept(:,:,jk) / rn_delta)
         pts(:,:,jk,jp_tem) = ( T0 + (rho0 - rho) / rn_a0 ) * ptmask(:,:,jk)
      END DO
      !
      pts(:,:,:,jp_sal) = 35._wp * ptmask(:,:,:)
      !!----------------------------------------------------------------------
      !
   END SUBROUTINE usr_def_istate

   
   SUBROUTINE usr_def_istate_ssh( ptmask, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate_ssh  ***
      !! 
      !! ** Purpose :   Initialization of ssh
      !!
      !! ** Method  :   Set ssh as null, ptmask is required for test cases
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask   [m]
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height   [m]
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate_ssh : GYRE configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~   Ocean at rest, ssh is zero'
      !
      ! Sea level:
      pssh(:,:) = 0._wp
      !
   END SUBROUTINE usr_def_istate_ssh

   !!======================================================================
END MODULE usrdef_istate
