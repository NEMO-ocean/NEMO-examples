MODULE usrdef_istate
   !!======================================================================
   !!                     ***  MODULE usrdef_istate   ***
   !!
   !!                      ===  VORTEX configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2017-11  (J. Chanut) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE dom_oce , ONLY : glamt, gphit, glamu, gphiu, glamv, gphiv  
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   !   
   USE usrdef_nam, ONLY : rn_ppgphi0  ! Reference latitude   
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate   ! called by istate.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_istate.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv, pssh )
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
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height
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
      zP0 = rau0 * zf0 * zumax * zlambda * SQRT(EXP(1._wp)/2._wp)
      !
      ! Sea level:
      za = -zP0 * (1._wp-EXP(-zH)) / (grav*(zH-1._wp + EXP(-zH)))
      DO ji=1, jpi
         DO jj=1, jpj
            zx = glamt(ji,jj) * 1.e3
            zy = gphit(ji,jj) * 1.e3
            zrho1 = rau0 + za * EXP(-(zx**2+zy**2)/zlambda**2)
            pssh(ji,jj) = zP0 * EXP(-(zx**2+zy**2)/zlambda**2)/(zrho1*grav) * ptmask(ji,jj,1)
         END DO
      END DO
      !
      ! temperature:         
      DO ji=1, jpi
         DO jj=1, jpj
            zx = glamt(ji,jj) * 1.e3
            zy = gphit(ji,jj) * 1.e3
            DO jk=1,jpk
               zdt =  pdept(ji,jj,jk) 
               zrho1 = rau0 * (1._wp + zn2*zdt/grav)
               IF (zdt < zH) THEN
                  zrho1 = zrho1 - zP0 * (1._wp-EXP(zdt-zH)) &
                          & * EXP(-(zx**2+zy**2)/zlambda**2) / (grav*(zH -1._wp + exp(-zH)));
               ENDIF
               pts(ji,jj,jk,jp_tem) = (20._wp + (rau0-zrho1) / 0.28_wp) * ptmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      ! salinity:  
      pts(:,:,:,jp_sal) = 35._wp 
      !
      ! velocities:
      za = 2._wp * zP0 / (zf0 * rau0 * zlambda**2)
      DO ji=1, jpim1
         DO jj=1, jpj
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
         END DO
      END DO
      !
      DO ji=1, jpi
         DO jj=1, jpjm1
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
         END DO
      END DO

      CALL lbc_lnk( 'usrdef_istate', pu, 'U', -1. )
      CALL lbc_lnk( 'usrdef_istate', pv, 'V', -1. )
      !   
   END SUBROUTINE usr_def_istate

   !!======================================================================
END MODULE usrdef_istate
