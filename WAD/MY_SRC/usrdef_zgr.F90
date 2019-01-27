MODULE usrdef_zgr
   !!======================================================================
   !!                   ***  MODULE  usrdef_zgr  ***
   !!
   !!                   ===  WAD_TEST_CASES case  ===
   !!
   !! Ocean domain : user defined vertical coordinate system 
   !!======================================================================
   !! History :  4.0  ! 2016-06  (G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_zgr   : user defined vertical coordinate system (required)
   !!       zgr_z     : reference 1D z-coordinate 
   !!---------------------------------------------------------------------
   USE oce            ! ocean variables
   USE dom_oce ,  ONLY: ht_0, mi0, mi1, nimpp, njmpp,  &
                      & mj0, mj1, glamt, gphit         ! ocean space and time domain
   USE usrdef_nam     ! User defined : namelist variables
   USE wet_dry ,  ONLY: rn_wdmin1, rn_wdmin2, rn_wdld  ! Wetting and drying
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_zgr        ! called by domzgr.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_zgr.F90 10425 2018-12-19 21:54:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS             

   SUBROUTINE usr_def_zgr( ld_zco  , ld_zps  , ld_sco  , ld_isfcav,    &   ! type of vertical coordinate
      &                    pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d  ,    &   ! 1D reference vertical coordinate
      &                    pdept , pdepw ,                             &   ! 3D t & w-points depth
      &                    pe3t  , pe3u  , pe3v   , pe3f ,             &   ! vertical scale factors
      &                    pe3w  , pe3uw , pe3vw         ,             &   !     -      -      -
      &                    k_top  , k_bot    )                             ! top & bottom ocean level
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE usr_def_zgr  ***
      !!
      !! ** Purpose :   User defined the vertical coordinates
      !!
      !!----------------------------------------------------------------------
      LOGICAL                   , INTENT(out) ::   ld_zco, ld_zps, ld_sco      ! vertical coordinate flags
      LOGICAL                   , INTENT(out) ::   ld_isfcav                   ! under iceshelf cavity flag
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pe3t_1d , pe3w_1d           ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pdept, pdepw                ! grid-point depth        [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors  [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pe3w , pe3uw, pe3vw         ! i-scale factors 
      INTEGER , DIMENSION(:,:)  , INTENT(out) ::   k_top, k_bot                ! first & last ocean level
      !
      INTEGER  ::   ji, jj, jk        ! dummy indices
      INTEGER  ::   ik                ! local integers
      REAL(wp) ::   zfact, z1_jpkm1   ! local scalar
      REAL(wp) ::   ze3min            ! local scalar
      REAL(wp) ::   zi, zj, zbathy    ! local scalar
      REAL(wp) ::   ztmpu, ztmpv, ztmpf, ztmpu1, ztmpv1, ztmpf1, zwet
      REAL(wp), DIMENSION(jpi,jpj) ::   zht, zhu, zhv, z2d   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr : WAD_TEST_CASES configuration (s-coordinate closed box ocean without cavities)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      !
      ! type of vertical coordinate ==>>>   here WAD_TEST_CASES : s-coordinate always
      ! ---------------------------
      ld_zco    = .FALSE.      ! z-partial-step coordinate
      ld_zps    = .FALSE.      ! z-partial-step coordinate
      ld_sco    = .TRUE.       ! s-coordinate
      ld_isfcav = .FALSE.      ! ISF Ice Shelves Flag
      !
      !
      ! Build the vertical coordinate system
      ! ------------------------------------
      !
      !                       !==  UNmasked meter bathymetry  ==!
      !
!
      zbathy=10.0
      IF( cn_cfg == 'wad' ) THEN
         SELECT CASE ( nn_cfg )
            !                                        ! ====================
            CASE ( 1 )                               ! WAD 1 configuration
               !                                     ! ====================
               !
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) 'usr_def_zgr (WAD) : Closed box with EW linear bottom slope'
               IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
               !
               zht = 1.5_wp
               DO ji = 1, jpi
                 zi = MIN((glamt(ji,1) - 10.0)/40.0, 1.0 )
                 zht(ji,:) = MAX(zbathy*zi, -2.0) 
               END DO
               zht(mi0(1):mi1(1),:) = -4._wp
               zht(mi0(jpiglo):mi1(jpiglo),:) = -4._wp
               zht(:,mj0(1):mj1(1)) = -4._wp
               zht(:,mj0(jpjglo):mj1(jpjglo)) = -4._wp
               !                                     ! ====================
            CASE ( 2, 3, 8 )                         ! WAD 2 or 3  configuration
               !                                     ! ====================
               !
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) 'usr_def_zgr (WAD) : Parobolic EW channel'
               IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
               !
               DO ji = 1, jpi
                 zi = MAX(1.0-((glamt(ji,1)-25.0)**2)/484.0, -0.3 )
                 zht(ji,:) = MAX(zbathy*zi, -2.0)
               END DO
               zht(mi0(1):mi1(1),:) = -4._wp
               zht(mi0(jpiglo):mi1(jpiglo),:) = -4._wp
               IF( nn_cfg /= 8 ) THEN
                  zht(:,mj0(1):mj1(1)) = -4._wp
                  zht(:,mj0(jpjglo):mj1(jpjglo)) = -4._wp
               ENDIF
               !                                     ! ====================
            CASE ( 4 )                               ! WAD 4 configuration
               !                                     ! ====================
               !
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) 'usr_def_zgr (WAD) : Parobolic bowl'
               IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
               !
               DO ji = 1, jpi
                 zi = MAX(1.0-((glamt(ji,1)-25.0)**2)/484.0, -2.0 )
               DO jj = 1, jpj
                 zj = MAX(1.0-((gphit(1,jj)-17.0)**2)/196.0, -2.0 )
                 zht(ji,jj) = MAX(zbathy*zi*zj, -2.0)
               END DO
               END DO
               zht(mi0(1):mi1(1),:) = -4._wp
               zht(mi0(jpiglo):mi1(jpiglo),:) = -4._wp
               zht(:,mj0(1):mj1(1)) = -4._wp
               zht(:,mj0(jpjglo):mj1(jpjglo)) = -4._wp
               !                                    ! ===========================
            CASE ( 5 )                              ! WAD 5 configuration
               !                                    ! ====================
               !
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) 'usr_def_zgr (WAD) : Double slope with shelf'
               IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
               !
               DO ji = 1, jpi
                 zi = MIN(glamt(ji,1)/45.0, 1.0 )
                 zht(ji,:) = MAX(zbathy*zi, -2.0)
                 IF( glamt(ji,1) >= 46.0 ) THEN
                   zht(ji,:) = 10.0
                 ELSE IF( glamt(ji,1) >= 20.0 .AND. glamt(ji,1) < 46.0 ) THEN
                   zi = 7.5/25.
                   zht(ji,:) = MAX(10. - zi*(47.-glamt(ji,1)),2.5)
                 ELSE IF( glamt(ji,1) >= 15.0 .AND. glamt(ji,1) < 20.0 ) THEN
                   zht(ji,:) = 2.5
                 ELSE IF( glamt(ji,1) >= 4.0 .AND. glamt(ji,1) < 15.0 ) THEN
                   zi = 4.5/11.0
                   zht(ji,:) = MAX(2.5 - zi*(16.0-glamt(ji,1)), -2.0)
                 ELSE IF( glamt(ji,1) >= 0.0 .AND. glamt(ji,1) < 4.0 ) THEN
                   zht(ji,:) = -2.0
                 ENDIF
               END DO
               !                                    ! ===========================
               zht(mi0(1):mi1(1),:) = -4._wp
               zht(mi0(jpiglo):mi1(jpiglo),:) = -4._wp
               zht(:,mj0(1):mj1(1)) = -4._wp
               zht(:,mj0(jpjglo):mj1(jpjglo)) = -4._wp
               !                                    ! ===========================
            CASE ( 6 )                              ! WAD 6 configuration
               !                                    ! ====================
               !
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) 'usr_def_zgr (WAD) : Parabolic channel with gaussian ridge'
               IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
               !
               DO ji = 1, jpi
                 zi = MAX(1.0-((glamt(ji,1)-25.0)**2)/484.0, -2.0 )
                 zj = 1.075*MAX(EXP(-1.0*((glamt(ji,1)-25.0)**2)/32.0) , 0.0 )
                 zht(ji,:) = MAX(zbathy*(zi-zj), -2.0)
               END DO
               zht(mi0(1):mi1(1),:) = -4._wp
               zht(mi0(jpiglo):mi1(jpiglo),:) = -4._wp
               zht(:,mj0(1):mj1(1)) = -4._wp
               zht(:,mj0(jpjglo):mj1(jpjglo)) = -4._wp
               !                                    ! ===========================
            CASE ( 7 )                              ! WAD 7 configuration
               !                                    ! ====================
               !
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) 'usr_def_zgr (WAD) : Double slope with open boundary'
               IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
               !
               DO ji = 1, jpi
                 zi = MIN(glamt(ji,1)/45.0, 1.0 )
                 zht(ji,:) = MAX(zbathy*zi, -2.0)
                 IF( glamt(ji,1) >= 46.0 ) THEN
                   zht(ji,:) = 10.0
                 ELSE IF( glamt(ji,1) >= 20.0 .AND. glamt(ji,1) < 46.0 ) THEN
                   zi = 7.5/25.
                   zht(ji,:) = MAX(10. - zi*(47.-glamt(ji,1)),2.5)
                 ELSE IF( glamt(ji,1) >= 15.0 .AND. glamt(ji,1) < 20.0 ) THEN
                   zht(ji,:) = 2.5
                 ELSE IF( glamt(ji,1) >= 4.0 .AND. glamt(ji,1) < 15.0 ) THEN
                   zi = 4.5/11.0
                   zht(ji,:) = MAX(2.5 - zi*(16.0-glamt(ji,1)), -2.0)
                 ELSE IF( glamt(ji,1) >= 0.0 .AND. glamt(ji,1) < 4.0 ) THEN
                   zht(ji,:) = -2.0
                 ENDIF
               END DO
               !                                    ! ===========================
               zht(mi0(1):mi1(1),:) = -4._wp
               zht(:,mj0(1):mj1(1)) = -4._wp
               zht(:,mj0(jpjglo):mj1(jpjglo)) = -4._wp
            CASE DEFAULT
               !                                    ! ===========================
               WRITE(ctmp1,*) 'WAD test with a ', nn_cfg,' option is not coded'
               !
               CALL ctl_stop( ctmp1 )
               !
         END SELECT
      END IF
              

      ! at u-point: averaging zht
      DO ji = 1, jpim1
         zhu(ji,:) = 0.5_wp * ( zht(ji,:) + zht(ji+1,:) )
      END DO
      CALL lbc_lnk( 'usrdef_zgr', zhu, 'U', 1. )     ! boundary condition: this mask the surrounding grid-points
      !                                ! ==>>>  set by hand non-zero value on first/last columns & rows 
      DO ji = mi0(1), mi1(1)              ! first row of global domain only
         zhu(ji,:) = zht(1,:)
      END DO
       DO ji = mi0(jpiglo), mi1(jpiglo)   ! last  row of global domain only
         zhu(ji,:) = zht(jpi,:)
      END DO
      ! at v-point: averaging zht
      zhv = 0._wp
      DO jj = 1, jpjm1
         zhv(:,jj) = 0.5_wp * ( zht(:,jj) + zht(:,jj+1) )
      END DO
      CALL lbc_lnk( 'usrdef_zgr', zhv, 'V', 1. )     ! boundary condition: this mask the surrounding grid-points
      DO jj = mj0(1), mj1(1)   ! first  row of global domain only
         zhv(:,jj) = zht(:,jj)
      END DO
      DO jj = mj0(jpjglo), mj1(jpjglo)   ! last  row of global domain only
         zhv(:,jj) = zht(:,jj)
      END DO
      !     
      CALL zgr_z( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! Reference z-coordinate system
      !
      !
      !                       !==  top masked level bathymetry  ==!  (all coordinates)
      !
      ! no ocean cavities : top ocean level is ONE, except over land
      ! the ocean basin surrounnded by land (1 grid-point) set through lbc_lnk call as jperio=0 
      z2d(:,:) = 1._wp                    ! surface ocean is the 1st level
      z2d(mi0(1):mi1(1),:) = 0._wp
      z2d(mi0(jpiglo):mi1(jpiglo),:) = 0._wp
      z2d(:,mj0(1):mj1(1)) = 0._wp
      z2d(:,mj0(jpjglo):mj1(jpjglo)) = 0._wp





      CALL lbc_lnk( 'usrdef_zgr', z2d, 'T', 1. )        ! closed basin since jperio = 0 (see userdef_nam.F90)
      k_top(:,:) = NINT( z2d(:,:) )
      !
      !                              
      !
      IF ( ld_sco ) THEN      !==  s-coordinate  ==!   (terrain-following coordinate)
         !
         ht_0 = zht
         k_bot(:,:) = jpkm1 * k_top(:,:)  !* bottom ocean = jpk-1 (here use k_top as a land mask)
         DO jj = 1, jpj
            DO ji = 1, jpi
              IF( zht(ji,jj) <= -(rn_wdld - rn_wdmin2)) THEN
                k_bot(ji,jj) = 0
                k_top(ji,jj) = 0
              ENDIF
           END DO
         END DO
         !
         !                                !* terrain-following coordinate with e3.(k)=cst)
         !                                !  OVERFLOW case : identical with j-index (T=V, U=F)
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
              z1_jpkm1 = 1._wp / REAL( k_bot(ji,jj) - k_top(ji,jj) + 1 , wp)
              DO jk = 1, jpk
                  zwet = MAX( zht(ji,jj), rn_wdmin1 )
                  pdept(ji,jj,jk) = zwet * z1_jpkm1 * ( REAL( jk   , wp ) - 0.5_wp )
                  pdepw(ji,jj,jk) = zwet * z1_jpkm1 * ( REAL( jk-1 , wp )          )
                  pe3t (ji,jj,jk) = zwet * z1_jpkm1
                  pe3w (ji,jj,jk) = zwet * z1_jpkm1
                  zwet = MAX( zhu(ji,jj), rn_wdmin1 )
                  pe3u (ji,jj,jk) = zwet * z1_jpkm1
                  pe3uw(ji,jj,jk) = zwet * z1_jpkm1
                  pe3f (ji,jj,jk) = zwet * z1_jpkm1
                  zwet = MAX( zhv(ji,jj), rn_wdmin1 )
                  pe3v (ji,jj,jk) = zwet * z1_jpkm1
                  pe3vw(ji,jj,jk) = zwet * z1_jpkm1
              END DO      
           END DO      
         END DO      
         CALL lbc_lnk( 'usrdef_zgr', pdept, 'T', 1. )
         CALL lbc_lnk( 'usrdef_zgr', pdepw, 'T', 1. )
         CALL lbc_lnk( 'usrdef_zgr', pe3t , 'T', 1. )
         CALL lbc_lnk( 'usrdef_zgr', pe3w , 'T', 1. )
         CALL lbc_lnk( 'usrdef_zgr', pe3u , 'U', 1. )
         CALL lbc_lnk( 'usrdef_zgr', pe3uw, 'U', 1. )
         CALL lbc_lnk( 'usrdef_zgr', pe3f , 'F', 1. )
         CALL lbc_lnk( 'usrdef_zgr', pe3v , 'V', 1. )
         CALL lbc_lnk( 'usrdef_zgr', pe3vw, 'V', 1. )
         WHERE( pe3t (:,:,:) == 0._wp )   pe3t (:,:,:) = 1._wp
         WHERE( pe3u (:,:,:) == 0._wp )   pe3u (:,:,:) = 1._wp
         WHERE( pe3v (:,:,:) == 0._wp )   pe3v (:,:,:) = 1._wp
         WHERE( pe3f (:,:,:) == 0._wp )   pe3f (:,:,:) = 1._wp
         WHERE( pe3w (:,:,:) == 0._wp )   pe3w (:,:,:) = 1._wp
         WHERE( pe3uw(:,:,:) == 0._wp )   pe3uw(:,:,:) = 1._wp
         WHERE( pe3vw(:,:,:) == 0._wp )   pe3vw(:,:,:) = 1._wp
      ENDIF
      !
      !
      !
   END SUBROUTINE usr_def_zgr


   SUBROUTINE zgr_z( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! 1D reference vertical coordinate
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_z  ***
      !!
      !! ** Purpose :   set the depth of model levels and the resulting 
      !!      vertical scale factors.
      !!
      !! ** Method  :   1D z-coordinate system (use in all type of coordinate)
      !!       The depth of model levels is set from dep(k), an analytical function:
      !!                   w-level: depw_1d  = dep(k)
      !!                   t-level: dept_1d  = dep(k+0.5)
      !!       The scale factors are the discrete derivative of the depth:
      !!                   e3w_1d(jk) = dk[ dept_1d ] 
      !!                   e3t_1d(jk) = dk[ depw_1d ]
      !!
      !!            ===    Here constant vertical resolution   ===
      !!
      !! ** Action  : - pdept_1d, pdepw_1d : depth of T- and W-point (m)
      !!              - pe3t_1d , pe3w_1d  : scale factors at T- and W-levels (m)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:), INTENT(out) ::   pdept_1d, pdepw_1d   ! 1D grid-point depth        [m]
      REAL(wp), DIMENSION(:), INTENT(out) ::   pe3t_1d , pe3w_1d    ! 1D vertical scale factors  [m]
      !
      INTEGER  ::   jk       ! dummy loop indices
      REAL(wp) ::   zt, zw   ! local scalar
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                         ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) '    zgr_z : Reference vertical z-coordinates: uniform dz = ', rn_dz
         WRITE(numout,*) '    ~~~~~~~'
      ENDIF
      !
      ! Reference z-coordinate (depth - scale factor at T- and W-points)   ! Madec & Imbard 1996 function
      ! ----------------------
      DO jk = 1, jpk
         zw = REAL( jk , wp )
         zt = REAL( jk , wp ) + 0.5_wp
         pdepw_1d(jk) =    rn_dz *   REAL( jk-1 , wp )
         pdept_1d(jk) =    rn_dz * ( REAL( jk-1 , wp ) + 0.5_wp )
         pe3w_1d (jk) =    rn_dz
         pe3t_1d (jk) =    rn_dz
      END DO
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) '              Reference 1D z-coordinate depth and scale factors:'
         WRITE(numout, "(9x,' level  gdept_1d  gdepw_1d  e3t_1d   e3w_1d  ')" )
         WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, pdept_1d(jk), pdepw_1d(jk), pe3t_1d(jk), pe3w_1d(jk), jk = 1, jpk )
      ENDIF
      !
   END SUBROUTINE zgr_z
   
   !!======================================================================
END MODULE usrdef_zgr
