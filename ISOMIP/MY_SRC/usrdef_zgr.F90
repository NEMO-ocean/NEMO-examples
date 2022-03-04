MODULE usrdef_zgr
   !!======================================================================
   !!                     ***  MODULE usrdef_zgr  ***
   !!
   !!                       ===  ISOMIP case  ===
   !!
   !! user defined :  vertical coordinate system of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-08  (G. Madec,   S. Flavoni)  Original code
   !!                 ! 2017-02  (P. Mathiot, S. Flavoni)  Adapt code to ISOMIP case
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_zgr   : user defined vertical coordinate system (required)
   !!       zgr_z1d   : reference 1D z-coordinate 
   !!---------------------------------------------------------------------
   USE oce            ! ocean variables
   USE dom_oce ,  ONLY: mj0   , mj1    ! ocean space and time domain
   USE dom_oce ,  ONLY: glamt , gphit  ! ocean space and time domain
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
   !! $Id: usrdef_zgr.F90 13295 2020-07-10 18:24:21Z acc $
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
      INTEGER  ::   ji , jj, jk       ! dummy indices
      INTEGER  ::   ij0, ij1          ! dummy indices  
      INTEGER  ::   ik                ! local integers
      REAL(wp) ::   zfact, z1_jpkm1   ! local scalar
      REAL(wp) ::   ze3min, zdepth    ! local scalar
      REAL(wp), DIMENSION(jpi,jpj) ::   zht  , zhu         ! bottom depth
      REAL(wp), DIMENSION(jpi,jpj) ::   zhisf, zhisfu      ! top depth
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr : ISOMIP configuration (z(ps)- or s-coordinate closed box ocean without cavities)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      !
      ! type of vertical coordinate
      ! ---------------------------
      ! set in usrdef_nam.F90 by reading the namusr_def namelist except for ISF
      ld_isfcav = .TRUE.       ! ISF Ice Shelves Flag
      !
      !
      ! Build the vertical coordinate system
      ! ------------------------------------
      !
      !                       !==  isfdraft  ==!
      !
      zht  (:,:) = rbathy 
      zhisf(:,:) = 200._wp
      ij0 = 1   ;   ij1 = 40+nn_hls
      DO jj = mj0(ij0), mj1(ij1)
         zhisf(:,jj)=700.0_wp-(gphit(:,jj)+80.0_wp)*125.0_wp
      END DO
      !
      CALL zgr_z1d( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! Reference z-coordinate system
      !
      !                       !==  top masked level bathymetry  ==!  (all coordinates)
      !
      IF ( ld_zps ) THEN      !==  zps-coordinate  ==!   (partial bottom-steps)
         !
         ze3min = 0.1_wp * rn_e3
         IF(lwp) WRITE(numout,*) '   minimum thickness of the partial cells = 10 % of e3 = ', ze3min
         !
         !                                !* bottom ocean compute from the depth of grid-points
         k_bot(:,:) = jpkm1
         DO jk = jpkm1, 1, -1
            WHERE( zht(:,:) < pdepw_1d(jk) + ze3min )   k_bot(:,:) = jk-1
         END DO
         !                                !* top ocean compute from the depth of grid-points
         k_top(:,:) = 1                   ! 
         DO jk = 2, jpkm1
            zdepth = pdepw_1d(jk+1) - ze3min
            WHERE( zhisf(:,:) > 0.0 .AND. zhisf(:,:) >= zdepth )   k_top(:,:) = (jk + 1) 
         END DO
         !
         !                                   !* vertical coordinate system
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
         ! top scale factors and depth at T- and W-points
         DO_2D( 1, 1, 1, 1 )
            ik = k_top(ji,jj)
            IF ( ik > 2 ) THEN
               ! pdeptw at the interface
               pdepw(ji,jj,ik  ) = MAX( zhisf(ji,jj) , pdepw(ji,jj,ik) )
               ! e3t in both side of the interface
               pe3t (ji,jj,ik  ) = pdepw(ji,jj,ik+1) - pdepw(ji,jj,ik)
               ! pdept in both side of the interface (from previous e3t)
               pdept(ji,jj,ik  ) = pdepw(ji,jj,ik  ) + pe3t (ji,jj,ik  ) * 0.5_wp
               pdept(ji,jj,ik-1) = pdepw(ji,jj,ik  ) - pe3t (ji,jj,ik  ) * 0.5_wp
               ! pe3w on both side of the interface
               pe3w (ji,jj,ik+1) = pdept(ji,jj,ik+1) - pdept(ji,jj,ik  )
               pe3w (ji,jj,ik  ) = pdept(ji,jj,ik  ) - pdept(ji,jj,ik-1)
               ! e3t into the ice shelf
               pe3t (ji,jj,ik-1) = pdepw(ji,jj,ik  ) - pdepw(ji,jj,ik-1)
               pe3w (ji,jj,ik-1) = pdept(ji,jj,ik-1) - pdept(ji,jj,ik-2)
            END IF
         END_2D
         ! bottom scale factors and depth at T- and W-points
         DO_2D( 1, 1, 1, 1 )
            ik = k_bot(ji,jj)
            pdepw(ji,jj,ik+1) = MIN( zht(ji,jj) , pdepw_1d(ik+1) )
            pe3t (ji,jj,ik  ) = pdepw(ji,jj,ik+1) - pdepw(ji,jj,ik)
            pe3t (ji,jj,ik+1) = pe3t (ji,jj,ik  ) 
            !
            pdept(ji,jj,ik  ) = pdepw(ji,jj,ik  ) + pe3t (ji,jj,ik  ) * 0.5_wp
            pdept(ji,jj,ik+1) = pdepw(ji,jj,ik+1) + pe3t (ji,jj,ik+1) * 0.5_wp
            pe3w (ji,jj,ik+1) = pdept(ji,jj,ik+1) - pdept(ji,jj,ik)
         END_2D       
         !                                   ! bottom scale factors and depth at  U-, V-, UW and VW-points
         pe3u (:,:,:) = pe3t(:,:,:)
         pe3uw(:,:,:) = pe3w(:,:,:)
         DO_3D( 0, 0, 0, 0, 1, jpk )
         !                                   ! Computed as the minimum of neighbooring scale factors
            pe3v (ji,jj,jk) = MIN( pe3t(ji,jj,jk), pe3t(ji,jj+1,jk) )
            pe3vw(ji,jj,jk) = MIN( pe3w(ji,jj,jk), pe3w(ji,jj+1,jk) )
            pe3f (ji,jj,jk) = pe3v(ji,jj,jk)
         END_3D
         CALL lbc_lnk( 'usrdef_zgr', pe3v , 'V', 1._wp )   ;   CALL lbc_lnk( 'usrdef_zgr', pe3vw, 'V', 1._wp )
         CALL lbc_lnk( 'usrdef_zgr', pe3f , 'F', 1._wp )
         DO jk = 1, jpk
            ! set to z-scale factor if zero (i.e. along closed boundaries) because of lbclnk
            WHERE( pe3u (:,:,jk) == 0._wp )   pe3u (:,:,jk) = pe3t_1d(jk)
            WHERE( pe3v (:,:,jk) == 0._wp )   pe3v (:,:,jk) = pe3t_1d(jk)
            WHERE( pe3f (:,:,jk) == 0._wp )   pe3f (:,:,jk) = pe3t_1d(jk)
            WHERE( pe3uw(:,:,jk) == 0._wp )   pe3uw(:,:,jk) = pe3w_1d(jk)
            WHERE( pe3vw(:,:,jk) == 0._wp )   pe3vw(:,:,jk) = pe3w_1d(jk)
         END DO
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
         WRITE(numout,*) '    zgr_z1d : Reference vertical z-coordinates: uniform dz = ', rn_e3
         WRITE(numout,*) '    ~~~~~~~'
      ENDIF
      !
      ! Reference z-coordinate (depth - scale factor at T- and W-points)   ! Madec & Imbard 1996 function
      ! ----------------------
      DO jk = 1, jpk
         zw = REAL( jk , wp )
         zt = REAL( jk , wp ) + 0.5_wp
         pdepw_1d(jk) =    rn_e3 *   REAL( jk-1 , wp )
         pdept_1d(jk) =    rn_e3 * ( REAL( jk-1 , wp ) + 0.5_wp )
         pe3w_1d (jk) =    rn_e3
         pe3t_1d (jk) =    rn_e3
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
