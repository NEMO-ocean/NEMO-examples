MODULE usrdef_istate
   !!======================================================================
   !!                     ***  MODULE usrdef_istate   ***
   !!
   !!                  ===  WAD_TEST_CASES configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  4.0 ! 2016-03  (S. Flavoni) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE dom_oce , ONLY : mi0, mig, mjg, glamt, gphit, ht_0
   USE phycst         ! physical constants
   USE wet_dry        ! Wetting and drying
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate   ! called in istate.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_istate.F90 10074 2018-08-28 16:15:49Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here WAD_TEST_CASES configuration 
      !!
      !! ** Method  : - set temprature field
      !!              - set salinity   field
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height
      INTEGER  ::   ji, jj            ! dummy loop indices
      REAL(wp) ::   zi, zj
      !
      INTEGER  ::   jk     ! dummy loop indices
      REAL(wp) ::   zdam   ! location of dam [Km]
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : WAD_TEST_CASES configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   Ocean at rest, with a constant temperature                    '
      IF(lwp) WRITE(numout,*) '                 and  constant salinity (not used as rho=F(T) '
      !
      !
      pu  (:,:,:) = 0._wp        ! ocean at rest
      pv  (:,:,:) = 0._wp
      pssh(:,:)   = 0._wp
      !
      !                          ! T & S profiles
      pts(:,:,:,jp_tem) = 10._wp * ptmask(:,:,:)
      !
      pts(:,:,:,jp_sal) = 35._wp * ptmask(:,:,:)
      !!----------------------------------------------------------------------
      !
      !!----------------------------------------------------------------------
      !
      ! Uniform T & S in most test cases
      pts(:,:,:,jp_tem) = 10._wp
      pts(:,:,:,jp_sal) = 35._wp
      SELECT CASE ( nn_cfg ) 
         !                                        ! ====================
         CASE ( 1 )                               ! WAD 1 configuration
            !                                     ! ====================
            !
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'usr_def_istate : WAD Closed box with EW linear bottom slope'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
            !
            do ji = 1,jpi
             pssh(ji,:) = ( -5.5_wp + 7.4_wp*glamt(ji,1)/50._wp)*ptmask(ji,:,1)
            end do
            !                                     ! ====================
         CASE ( 2, 8 )                            ! WAD 2 configuration
            !                                     ! ====================
            !
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'usr_def_istate : WAD Parobolic EW channel, mid-range initial ssh slope'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
            !
            do ji = 1,jpi
             pssh(ji,:) = ( -1.5_wp + 5.0_wp*glamt(ji,1)/50._wp)*ptmask(ji,:,1)
            end do
            !                                     ! ====================
         CASE ( 3 )                               ! WAD 3 configuration
            !                                     ! ====================
            !
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'usr_def_istate : WAD Parobolic EW channel, extreme initial ssh slope' 
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
            !
            do ji = 1,jpi
             pssh(ji,:) = ( -4.5_wp + 6.8_wp*glamt(ji,1)/50._wp)*ptmask(ji,:,1)
            end do

            !
            !                                     ! ====================
         CASE ( 4 )                               ! WAD 4 configuration
            !                                     ! ====================
            !
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'usr_def_istate : WAD Parobolic bowl, mid-range initial ssh slope' 
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
            !
            DO ji = 1, jpi
               zi = MAX(1.0-((glamt(ji,1)-25._wp)**2)/400.0, 0.0 )
               DO jj = 1, jpj
                  zj = MAX(1.0-((gphit(1,jj)-17._wp)**2)/144.0, 0.0 )
                  pssh(ji,jj) = -2.5_wp + 5.4_wp*zi*zj
               END DO
            END DO

            !
            !                                    ! ===========================
         CASE ( 5, 7 )                           ! WAD 5 and 7 configurations
            !                                    ! ===========================
            !
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'usr_def_istate : WAD Double slope with shelf'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
            !
            do ji = 1,jpi
             pssh(ji,:) = ( -2.5_wp + 5.5_wp*glamt(ji,1)/50._wp)*ptmask(ji,:,1)
            end do

            !
            !                                     ! ====================
         CASE ( 6 )                               ! WAD 6 configuration
            !                                     ! ====================
            !
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'usr_def_istate : WAD Parobolic EW channel with gaussian ridge' 
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
            !
            do ji = 1,jpi
             pssh(ji,:) = ( -2.5_wp + 5.5_wp*(50._wp-glamt(ji,1))/50._wp)*ptmask(ji,:,1)
            end do
            !
            do ji = mi0(jpiglo/2), mi0(jpiglo)
             pts(ji,:,:,jp_sal) = 30._wp
             pssh(ji,:) = -0.1*ptmask(ji,:,1)
            end do
            !
            !
            !                                    ! ===========================
         CASE DEFAULT                            ! NONE existing configuration
            !                                    ! ===========================
            WRITE(ctmp1,*) 'WAD test with a ', nn_cfg,' option is not coded'
            !
            CALL ctl_stop( ctmp1 )
            !
      END SELECT


      !
      ! Apply minimum wetdepth criterion
      !
      do jj = 1,jpj
         do ji = 1,jpi
            IF( ht_0(ji,jj) + pssh(ji,jj) < rn_wdmin1 ) THEN
               pssh(ji,jj) = ptmask(ji,jj,1)*( rn_wdmin1 - ht_0(ji,jj) )
            ENDIF
         end do
      end do
      !
   END SUBROUTINE usr_def_istate

   !!======================================================================
END MODULE usrdef_istate
