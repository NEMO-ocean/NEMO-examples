MODULE usrdef_istate
   !!======================================================================
   !!                     ***  MODULE usrdef_istate   ***
   !!
   !!                      ===  BENCH configuration  ===
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
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv ) !!st, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here BENCH configuration 
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
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : BENCH configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   '
      !
      ! define unique value on each point of the inner global domain. z2d ranging from 0.05 to -0.05
      !
      DO_2D( 0, 0, 0, 0 )   !  +/- 0.05
         z2d(ji,jj) = 0.1 * ( 0.5 - REAL( mig0(ji) + (mjg0(jj)-1) * Ni0glo, wp ) / REAL( Ni0glo * Nj0glo, wp ) )
      END_2D
      !
      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         zfact = REAL(jk-1,wp) / REAL(jpk-1,wp)   ! 0 to 1 to add a basic stratification
         ! temperature choosen to lead to ~50% ice at the beginning if rn_thres_sst = 0.5
         pts(ji,jj,jk,jp_tem) = 20._wp*z2d(ji,jj) - 1._wp - 0.5_wp * zfact    ! -1 to -1.5 +/- 1.0 degG
         ! salinity:  
         pts(ji,jj,jk,jp_sal) = 30._wp + 1._wp * zfact + z2d(ji,jj)           ! 30 to 31   +/- 0.05 psu
         ! velocities:
         pu(ji,jj,jk) = z2d(ji,jj) *  0.1_wp * umask(ji,jj,jk)                ! +/- 0.005  m/s
         pv(ji,jj,jk) = z2d(ji,jj) * 0.01_wp * vmask(ji,jj,jk)                ! +/- 0.0005 m/s
      END_3D
      pts(:,:,jpk,:) = 0._wp
      pu( :,:,jpk  ) = 0._wp
      pv( :,:,jpk  ) = 0._wp
      !
      CALL lbc_lnk('usrdef_istate',  pts, 'T',  1. )            ! apply boundary conditions
      CALL lbc_lnk('usrdef_istate',   pu, 'U', -1. )            ! apply boundary conditions
      CALL lbc_lnk('usrdef_istate',   pv, 'V', -1. )            ! apply boundary conditions
      
   END SUBROUTINE usr_def_istate


   SUBROUTINE usr_def_istate_ssh( ptmask, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate_ssh  ***
      !! 
      !! ** Purpose :   Initialization of ssh
      !!                Here BENCH configuration 
      !!
      !! ** Method  :   Set ssh
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask   [m]
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height   [m]
      !
      INTEGER  ::   ji, jj
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate_ssh : BENCH configuration, analytical definition of initial ssh'
      !
      DO_2D( 0, 0, 0, 0 )                              ! sea level:  +/- 0.05 m
         pssh(ji,jj) = 0.1 * ( 0.5 - REAL( mig0(ji) + (mjg0(jj)-1) * Ni0glo, wp ) / REAL( Ni0glo * Nj0glo, wp ) )
      END_2D
      !
      CALL lbc_lnk('usrdef_istate', pssh, 'T',  1. )   ! apply boundary conditions
      !
   END SUBROUTINE usr_def_istate_ssh
   
   !!======================================================================
END MODULE usrdef_istate
