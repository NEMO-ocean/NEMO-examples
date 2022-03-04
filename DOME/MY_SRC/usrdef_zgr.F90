MODULE usrdef_zgr
   !!======================================================================
   !!                     ***  MODULE usrdef_zgr  ***
   !!
   !!                       ===  DOME case  ===
   !!
   !! user defined :  vertical coordinate system of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2020-12  (J. Chanut)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_zgr   : user defined vertical coordinate system (required)
   !!       zgr_z1d   : reference 1D z-coordinate 
   !!---------------------------------------------------------------------
   USE oce            ! ocean variables
   USE dom_oce ,  ONLY: mi0, mi1          ! ocean space and time domain
   USE dom_oce ,  ONLY: glamt, gphit      ! ocean space and time domain
   USE usrdef_nam     ! User defined : namelist variables
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_zgr   ! called by domzgr.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_zgr.F90 14053 2020-12-03 13:48:38Z techene $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS             

   SUBROUTINE usr_def_zgr( ld_zco  , ld_zps  , ld_sco  , ld_isfcav,    &   ! type of vertical coordinate
      &                    pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d  ,    &   ! 1D reference vertical coordinate
      &                    pdept , pdepw ,                             &   ! 3D t & w-points depth
      &                    pe3t  , pe3u  , pe3v , pe3f ,               &   ! vertical scale factors
      &                    pe3w  , pe3uw , pe3vw,                      &   !     -      -      -
      &                    k_top  , k_bot    )                             ! top & bottom ocean level
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE usr_def_zgr  ***
      !!
      !! ** Purpose :   User defined the vertical coordinates
      !!
      !!----------------------------------------------------------------------
      LOGICAL                   , INTENT(in   ) ::   ld_zco, ld_zps, ld_sco      ! vertical coordinate flags ( read in namusr_def )
      LOGICAL                   , INTENT(  out) ::   ld_isfcav                   ! under iceshelf cavity flag
      REAL(wp), DIMENSION(:)    , INTENT(  out) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:)    , INTENT(  out) ::   pe3t_1d , pe3w_1d           ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! grid-point depth        [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors  [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         ! i-scale factors 
      INTEGER , DIMENSION(:,:)  , INTENT(  out) ::   k_top, k_bot                ! first & last ocean level
      !
      INTEGER  ::   ji, jj, jk        ! dummy indices
      INTEGER  ::   ik                ! local integers
      REAL(wp) ::   zfact, z1_jpkm1   ! local scalar
      REAL(wp) ::   ze3min            ! local scalar
      REAL(wp), DIMENSION(jpi,jpj) ::   zht, zhu, zhv, zhf, z2d   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr : DOME configuration (z(ps)- or s-coordinate closed box ocean without cavities)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      !
      ! type of vertical coordinate
      ! ---------------------------
      ! already set in usrdef_nam.F90 by reading the namusr_def namelist except for ISF
      ld_isfcav = .FALSE.
      !
      !
      ! Build the vertical coordinate system
      ! ------------------------------------
      !
      !                       !==  UNmasked meter bathymetry  ==!
      !
      zht(:,:) = MIN(3600._wp, MAX( 600._wp,  600._wp - gphit(:,:)*1.e3*0.01 ))
      !
      ! at u/v/f-point: averaging zht
      zhu(:,:) = 600_wp ; zhv(:,:) = 600_wp ; zhf(:,:) = 600_wp
      DO_2D( 0, 0, 0, 0 )
         zhu(ji,jj) =  0.5_wp * ( zht(ji,jj  ) + zht(ji+1,jj  ) )
         zhv(jj,jj) =  0.5_wp * ( zht(ji,jj  ) + zht(ji  ,jj+1) )
         zhf(ji,jj) = 0.25_wp * ( zht(ji,jj  ) + zht(ji+1,jj  ) &
                       &        + zht(ji,jj+1) + zht(ji+1,jj+1) ) 
      END_2D
      CALL lbc_lnk( 'usrdef_zgr', zhu, 'U', 1.0_wp, zhv, 'V', 1.0_wp, zhf, 'F', 1.0_wp)      
      !     
      CALL zgr_z1d( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! Reference z-coordinate system
      !
      !
      !                       !==  top masked level bathymetry  ==!  (all coordinates)
      !
      ! no ocean cavities : top ocean level is ONE, except over land
      ! the ocean basin surrounded by land (1+nn_hls grid-point) set through lbc_lnk call
      z2d(:,:) = 1._wp                    ! surface ocean is the 1st level
      WHERE (gphit(:,:)>0._wp) z2d(:,:) = 0._wp
      ! Dig inlet:
      WHERE ((gphit(:,:)>0._wp).AND.(glamt(:,:)>-50._wp).AND.(glamt(:,:)<50._wp)) z2d(:,:) = 1._wp
      CALL lbc_lnk( 'usrdef_zgr', z2d, 'T', 1. )        ! closed basin, see userdef_nam.F90
      k_top(:,:) = NINT( z2d(:,:) )
      !
      !                              
      !
      IF ( ld_sco ) THEN      !==  s-coordinate  ==!   (terrain-following coordinate)
         !
         k_bot(:,:) = jpkm1 * k_top(:,:)  !* bottom ocean = jpk-1 (here use k_top as a land mask)
         !
         !                                !* terrain-following coordinate with e3.(k)=cst)
         z1_jpkm1 = 1._wp / REAL( jpkm1 , wp)
         DO jk = 1, jpk
            pdept(:,:,jk) = zht(:,:) * z1_jpkm1 * ( REAL( jk   , wp ) - 0.5_wp )
            pdepw(:,:,jk) = zht(:,:) * z1_jpkm1 * ( REAL( jk-1 , wp )          )
            pe3t (:,:,jk) = zht(:,:) * z1_jpkm1
            pe3u (:,:,jk) = zhu(:,:) * z1_jpkm1
            pe3v (:,:,jk) = zhv(:,:) * z1_jpkm1
            pe3f (:,:,jk) = zhf(:,:) * z1_jpkm1
            pe3w (:,:,jk) = zht(:,:) * z1_jpkm1
            pe3uw(:,:,jk) = zhu(:,:) * z1_jpkm1
            pe3vw(:,:,jk) = zhv(:,:) * z1_jpkm1
         END DO      
      ENDIF
      !
      !
      IF ( ld_zco ) THEN      !==  z-coordinate  ==!   (step-like topography)
         !
         !                                !* bottom ocean compute from the depth of grid-points
         k_bot(:,:) = jpkm1 * k_top(:,:)     ! here use k_top as a land mask
         DO jk = 1, jpkm1
            WHERE( pdept_1d(jk) < zht(:,:) .AND. zht(:,:) <= pdept_1d(jk+1) )   k_bot(:,:) = jk * k_top(:,:)
         END DO
         !                                !* horizontally uniform coordinate (reference z-co everywhere)
         DO jk = 1, jpk
            pdept(:,:,jk) = pdept_1d(jk)
            pdepw(:,:,jk) = pdepw_1d(jk)
            pe3t (:,:,jk) = pe3t_1d (jk)
            pe3u (:,:,jk) = pe3t_1d (jk)
            pe3v (:,:,jk) = pe3t_1d (jk)
            pe3f (:,:,jk) = pe3t_1d (jk)
            pe3w (:,:,jk) = pe3w_1d (jk)
            pe3uw(:,:,jk) = pe3w_1d (jk)
            pe3vw(:,:,jk) = pe3w_1d (jk)
         END DO
      ENDIF
      !
      !
      IF ( ld_zps ) THEN      !==  zps-coordinate  ==!   (partial bottom-steps)
         !
         ze3min = 0.1_wp * rn_dz
         IF(lwp) WRITE(numout,*) '   minimum thickness of the partial cells = 10 % of e3 = ', ze3min
         !
         !
         !                                !* bottom ocean compute from the depth of grid-points
         k_bot(:,:) = jpkm1
         DO jk = jpkm1, 1, -1
            WHERE( zht(:,:) < pdepw_1d(jk) + ze3min )   k_bot(:,:) = jk-1
         END DO
         !
         !                                !* vertical coordinate system
         DO jk = 1, jpk                      ! initialization to the reference z-coordinate
            pdept(:,:,jk) = pdept_1d(jk)
            pdepw(:,:,jk) = pdepw_1d(jk)
            pe3t (:,:,jk) = pe3t_1d (jk)
            pe3u (:,:,jk) = pe3t_1d (jk)
            pe3v (:,:,jk) = pe3t_1d (jk)
            pe3f (:,:,jk) = pe3t_1d (jk)
            pe3w (:,:,jk) = pe3w_1d (jk)
            pe3uw(:,:,jk) = pe3w_1d (jk)
            pe3vw(:,:,jk) = pe3w_1d (jk)
         END DO
         DO_2D( 1, 1, 1, 1 )
            ik = k_bot(ji,jj)
            pdepw(ji,jj,ik+1) = MIN( zht(ji,jj) , pdepw_1d(ik+1) )
            pe3t (ji,jj,ik  ) = pdepw(ji,jj,ik+1) - pdepw(ji,jj,ik)
            pe3t (ji,jj,ik+1) = pe3t (ji,jj,ik  ) 
            !
            pdept(ji,jj,ik  ) = pdepw(ji,jj,ik  ) + pe3t (ji,jj,ik  ) * 0.5_wp
            pdept(ji,jj,ik+1) = pdepw(ji,jj,ik+1) + pe3t (ji,jj,ik+1) * 0.5_wp
            pe3w (ji,jj,ik+1) = pdept(ji,jj,ik+1) - pdept(ji,jj,ik)              ! = pe3t (ji,jj,ik  )
            pe3w (ji,jj,ik  ) = pdept(ji,jj,ik  ) - pdept(ji,jj,ik-1)            ! st caution ik > 1
         END_2D         
         !
         DO_3D( 0, 0, 0, 0, 1, jpk ) 
               pe3u (ji,jj,jk) = MIN( pe3t(ji,jj,jk), pe3t(ji+1,jj,jk) )
               pe3v (ji,jj,jk) = MIN( pe3t(ji,jj,jk), pe3t(ji,jj+1,jk) )
               pe3uw(ji,jj,jk) = MIN( pe3w(ji,jj,jk), pe3w(ji+1,jj,jk) )
               pe3vw(ji,jj,jk) = MIN( pe3w(ji,jj,jk), pe3w(ji,jj+1,jk) )
         END_3D 
         !
         CALL lbc_lnk('usrdef_zgr', pe3u , 'U', 1._wp, pe3uw, 'U', 1._wp )   
         CALL lbc_lnk('usrdef_zgr', pe3v , 'V', 1._wp, pe3vw, 'V', 1._wp ) 
         !
         DO jk = 1, jpk                 
            WHERE( pe3u (:,:,jk) == 0._wp )   pe3u (:,:,jk) = pe3t_1d(jk)
            WHERE( pe3v (:,:,jk) == 0._wp )   pe3v (:,:,jk) = pe3t_1d(jk)
            WHERE( pe3uw(:,:,jk) == 0._wp )   pe3uw(:,:,jk) = pe3w_1d(jk)
            WHERE( pe3vw(:,:,jk) == 0._wp )   pe3vw(:,:,jk) = pe3w_1d(jk)
         END DO

         DO_3D( 0, 0, 0, 0, 1, jpk )
               pe3f(ji,jj,jk) = MIN( pe3v(ji,jj,jk), pe3v(ji+1,jj,jk) )
         END_3D
         CALL lbc_lnk('usrdef_zgr', pe3f, 'F', 1._wp )      
         !      
      ENDIF
      !
   END SUBROUTINE usr_def_zgr


   SUBROUTINE zgr_z1d( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! 1D reference vertical coordinate
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_z1d  ***
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
         WRITE(numout,*) '    zgr_z1d : Reference vertical z-coordinates: uniform dz = ', rn_dz
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
   END SUBROUTINE zgr_z1d
   
   !!======================================================================
END MODULE usrdef_zgr
