MODULE usrdef_istate
   !!======================================================================
   !!                     ***  MODULE usrdef_istate   ***
   !!
   !!                      ===  VORTEX configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2017-11  (J. Chanut) Original code
   !!                 ! 2020-11  (S. Techene, G. Madec) separate tsuv from ssh
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE dom_oce , ONLY : glamt, gphit, glamu, gphiu, glamv, gphiv  
   USE phycst         ! physical constants
   USE eosbn2  , ONLY : rn_a0
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   !   
   USE usrdef_nam, ONLY : rn_ppgphi0  ! Reference latitude   
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate       ! called by istate.F90
   PUBLIC   usr_def_istate_ssh   ! called by domqco.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_istate.F90 14857 2021-05-12 16:47:25Z hadcv $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here VORTEX configuration 
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
      INTEGER  :: ji, jj, jk  ! dummy loop indices
      REAL(wp) :: zx, zy, zP0, zumax, zlambda, zn2, zf0, zH, zrho1, za, zf
      REAL(wp) :: zdt, zdu, zdv
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : VORTEX configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   '
      !
      !
      !
      zf0   = 2._wp * omega * SIN( rad * rn_ppgphi0 )
      zumax = 1._wp * SIGN(1._wp, zf0) ! Here Anticyclonic: set zumax=-1 for cyclonic
      zlambda = SQRT(2._wp)*60.e3      ! Horizontal scale in meters
      zn2 = 3.e-3**2
      zH = 0.5_wp * 5000._wp
      !
      zP0 = rho0 * zf0 * zumax * zlambda * SQRT(EXP(1._wp)/2._wp)
      !
      ! temperature:         
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         zx = glamt(ji,jj) * 1.e3
         zy = gphit(ji,jj) * 1.e3
         DO jk=1,jpk
            zdt =  pdept(ji,jj,jk) 
            zrho1 = rho0 * (1._wp + zn2*zdt/grav)
            IF (zdt < zH) THEN
               zrho1 = zrho1 - zP0 * (1._wp-EXP(zdt-zH)) &
                  & * EXP(-(zx**2+zy**2)/zlambda**2) / (grav*(zH -1._wp + EXP(-zH)));
            ENDIF
            pts(ji,jj,jk,jp_tem) = (20._wp + (rho0-zrho1) / rn_a0 ) * ptmask(ji,jj,jk)
         END DO
      END_2D
      !
      ! salinity:  
      pts(:,:,:,jp_sal) = 35._wp * ptmask(:,:,:) 
      !
      ! velocities:
      za = 2._wp * zP0 / (zf0 * rho0 * zlambda**2)
      DO_2D( 0, 0, 0, 0 )
         zx = glamu(ji,jj) * 1.e3
         zy = gphiu(ji,jj) * 1.e3
         DO jk=1, jpk
            zdu = 0.5_wp * (pdept(ji  ,jj,jk) + pdept(ji+1,jj,jk))
            IF (zdu < zH) THEN
               zf = (zH-1._wp-zdu+EXP(zdu-zH)) / (zH-1._wp+EXP(-zH))
               pu(ji,jj,jk) = (za * zf * zy * EXP(-(zx**2+zy**2)/zlambda**2)) * ptmask(ji,jj,jk) * ptmask(ji+1,jj,jk)
            ELSE
               pu(ji,jj,jk) = 0._wp
            ENDIF
         END DO
      END_2D
      !
      DO_2D( 0, 0, 0, 0 )
         zx = glamv(ji,jj) * 1.e3
         zy = gphiv(ji,jj) * 1.e3
         DO jk=1, jpk
            zdv = 0.5_wp * (pdept(ji  ,jj,jk) + pdept(ji,jj+1,jk))
            IF (zdv < zH) THEN
               zf = (zH-1._wp-zdv+EXP(zdv-zH)) / (zH-1._wp+EXP(-zH))
               pv(ji,jj,jk) = -(za * zf * zx * EXP(-(zx**2+zy**2)/zlambda**2)) * ptmask(ji,jj,jk) * ptmask(ji,jj+1,jk)
            ELSE
               pv(ji,jj,jk) = 0._wp
            ENDIF
         END DO
      END_2D
      !
      CALL lbc_lnk( 'usrdef_istate', pu, 'U', -1., pv, 'V', -1. )
      !   
   END SUBROUTINE usr_def_istate


   SUBROUTINE usr_def_istate_ssh( ptmask, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of ssh
      !!                Here VORTEX configuration 
      !!
      !! ** Method  :   Set ssh according to a gaussian anomaly of pressure and associated
      !!                geostrophic velocities
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask   [m]
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height   [m]
      !
      INTEGER  :: ji, jj ! dummy loop indices
      REAL(wp) :: zx, zy, zP0, zumax, zlambda, zf0, zH, zrho1, za
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate_ssh : VORTEX configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   '
      !
      !
      !
      zf0   = 2._wp * omega * SIN( rad * rn_ppgphi0 )
      zumax = 1._wp * SIGN(1._wp, zf0) ! Here Anticyclonic: set zumax=-1 for cyclonic
      zlambda = SQRT(2._wp)*60.e3      ! Horizontal scale in meters 
      zH = 0.5_wp * 5000._wp
      !
      zP0 = rho0 * zf0 * zumax * zlambda * SQRT(EXP(1._wp)/2._wp)
      !
      ! Sea level:
      za = -zP0 * (1._wp-EXP(-zH)) / (grav*(zH-1._wp + EXP(-zH)))
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         zx = glamt(ji,jj) * 1.e3
         zy = gphit(ji,jj) * 1.e3
         zrho1 = rho0 + za * EXP(-(zx**2+zy**2)/zlambda**2)
         pssh(ji,jj) = zP0 * EXP(-(zx**2+zy**2)/zlambda**2)/(zrho1*grav) * ptmask(ji,jj,1)
      END_2D
      
   END SUBROUTINE usr_def_istate_ssh

   !!======================================================================
END MODULE usrdef_istate
