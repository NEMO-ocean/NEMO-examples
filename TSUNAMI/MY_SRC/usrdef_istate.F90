MODULE usrdef_istate
   !!======================================================================
   !!                     ***  MODULE usrdef_istate   ***
   !!
   !!                      ===  TSUNAMI configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  NEMO !
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE dom_oce        
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lbclnk          ! lateral boundary conditions - mpp exchanges
   !   
   USE usrdef_nam
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate       ! called by istate.F90
   PUBLIC   usr_def_istate_ssh   ! called by domqco.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2016)
   !! $Id$ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here TSUNAMI configuration 
      !!
      !! ** Method  :   Set a gaussian anomaly of pressure and associated
      !!                geostrophic velocities
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! 2D workspace
      REAL(wp) ::   zfact
      INTEGER  ::   ji, jj, jk
      INTEGER  ::   igloi, igloj   ! to be removed in the future, see comment bellow
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : TSUNAMI configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   '
      !
      pts(:,:,:,jp_tem) = 20._wp
      pts(:,:,:,jp_sal) = 30._wp
      pu( :,:,:       ) =  0._wp
      pv( :,:,:       ) =  0._wp
      
   END SUBROUTINE usr_def_istate


   SUBROUTINE usr_def_istate_ssh( ptmask, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate_ssh  ***
      !! 
      !! ** Purpose :   Initialization of ssh
      !!                Here TSUNAMI configuration 
      !!
      !! ** Method  :   Set ssh
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask   [m]
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height   [m]
      !
      INTEGER  ::   ji, jj
      REAL(wp), DIMENSION(jpi,jpj) ::  zdist
      REAL(wp)                     ::  zmax
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate_ssh : TSUNAMI configuration, analytical definition of initial ssh'
      !
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         zdist(ji,jj) = SQRT( glamt(ji,jj)**2 + gphit(ji,jj)**2 )
      END_2D
      zmax = MAXVAL( zdist ) / 20._wp
      CALL mpp_max( 'usrdef_istate', zmax )

      pssh(:,:) = 0 
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         IF( zdist(ji,jj) <= zmax )   pssh(ji,jj) = 0.1 * COS( zdist(ji,jj) / zmax * rpi * 0.5_wp )
      END_2D

      !
      CALL lbc_lnk('usrdef_istate', pssh, 'T',  1. )            ! apply boundary conditions
      !
   END SUBROUTINE usr_def_istate_ssh
   
   !!======================================================================
END MODULE usrdef_istate
