MODULE usrdef_zgr
   !!======================================================================
   !!                       ***  MODULE  usrdef_zgr  ***
   !!
   !!                       ===  SEAMOUNT configuration  ===
   !!
   !! User defined : vertical coordinate system of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-06  (G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_zgr   : user defined vertical coordinate system
   !!      zgr_z      : reference 1D z-coordinate 
   !!      zgr_top_bot: ocean top and bottom level indices
   !!      zgr_zco    : 3D verticl coordinate in pure z-coordinate case
   !!---------------------------------------------------------------------
   USE oce            ! ocean variables
   USE dom_oce        ! ocean domain
   USE depth_e3       ! depth <=> e3
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library
   USE usrdef_nam

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_zgr        ! called by domzgr.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_zgr.F90 13286 2020-07-09 15:48:29Z smasson $
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
      LOGICAL                   , INTENT(in ) ::   ld_zco, ld_zps, ld_sco      ! vertical coordinate flags
      LOGICAL                   , INTENT(in ) ::   ld_isfcav                   ! under iceshelf cavity flag
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
      REAL(wp) ::   zlam_mid, zphi_mid! local scalar
      REAL(wp), DIMENSION(jpi,jpj) ::   zht, zhu, zhv, z2d   ! 2D workspace
      !
      !!----------------------------------------------------------------------
      !
      !
      ! ------------------------------------
      ! Build the vertical coordinate system
      ! ------------------------------------
      !
      !                       !==  Unmasked meter bathymetry  ==!
      !
      !
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr (SEAMOUNT) : Isolated Gaussian bump in E-W periodic channel'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      !
      zlam_mid = 0.5_wp * 1000._wp * rn_length
      zphi_mid = 0.5_wp * 1000._wp * rn_width
      zht = 0._wp
      DO jj = 1, jpj
         DO ji = 1, jpi
            zht(ji,jj) = rn_bathy - rn_seamountheight * EXP( &
                    &  - ( ( 1000._wp * glamt(ji,jj) - zlam_mid) ** 2 + ( 1000._wp * gphit(ji,jj) - zphi_mid) ** 2) & 
                    &  / rn_l ** 2)
         END DO
      END DO
      ! 
      ! ------------------------------------
      CALL zgr_z( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! Reference z-coordinate system
      !
      CALL zgr_msk_top_bot( k_top , k_bot )                 ! masked top and bottom ocean t-level indices
      !
      !                                                     ! z-coordinate (3D arrays) from the 1D z-coord.
      IF ( ln_zco ) CALL zgr_zco( zht, pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d,   &   ! in  : 1D reference vertical coordinate
            &          pdept   , pdepw   ,                     &   ! out : 3D t & w-points depth
            &          pe3t    , pe3u    , pe3v   , pe3f   ,   &   !       vertical scale factors
            &          pe3w    , pe3uw   , pe3vw, k_bot, k_top             )     !           -      -      -
      !                                                     ! z-coordinate (3D arrays) from the 1D z-coord.
      IF ( ln_zps ) CALL zgr_zps( zht, pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d,   &   ! in  : 1D reference vertical coordinate
            &          pdept   , pdepw   ,                     &   ! out : 3D t & w-points depth
            &          pe3t    , pe3u    , pe3v   , pe3f   ,   &   !       vertical scale factors
            &          pe3w    , pe3uw   , pe3vw, k_bot, k_top             )     !           -      -      -
      !                                                     ! z-coordinate (3D arrays) from the 1D z-coord.
      IF ( ln_sco ) CALL zgr_sco( zht,   &                  ! in  : reference bathymetry
            &          pdept   , pdepw   ,                     &   ! out : 3D t & w-points depth
            &          pe3t    , pe3u    , pe3v   , pe3f   ,   &   !       vertical scale factors
            &          pe3w    , pe3uw   , pe3vw, k_bot, k_top            )     !           -      -      -

      !
   END SUBROUTINE usr_def_zgr


   SUBROUTINE zgr_z( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! 1D reference vertical coordinate
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_z  ***
      !!
      !! ** Purpose :   set the 1D depth of model levels and the resulting 
      !!              vertical scale factors.
      !!
      !! ** Method  :   1D z-coordinate system (use in all type of coordinate)
      !!       The depth of model levels is set from dep(k), an analytical function:
      !!                   w-level: depw_1d  = dep(k)
      !!                   t-level: dept_1d  = dep(k+0.5)
      !!       The scale factors are the discrete derivative of the depth:
      !!                   e3w_1d(jk) = dk[ dept_1d ] 
      !!                   e3t_1d(jk) = dk[ depw_1d ]
      !!           with at top and bottom :
      !!                   e3w_1d( 1 ) = 2 * ( dept_1d( 1 ) - depw_1d( 1 ) )
      !!                   e3t_1d(jpk) = 2 * ( dept_1d(jpk) - depw_1d(jpk) )
      !!       The depth are then re-computed from the sum of e3. This ensures 
      !!    that depths are identical when reading domain configuration file. 
      !!    Indeed, only e3. are saved in this file, depth are compute by a call
      !!    to the e3_to_depth subroutine.
      !!
      !!       Here the Madec & Imbard (1996) function is used.
      !!
      !! ** Action  : - pdept_1d, pdepw_1d : depth of T- and W-point (m)
      !!              - pe3t_1d , pe3w_1d  : scale factors at T- and W-levels (m)
      !!
      !! Reference : Marti, Madec & Delecluse, 1992, JGR, 97, No8, 12,763-12,766.
      !!             Madec and Imbard, 1996, Clim. Dyn.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pdept_1d, pdepw_1d   ! 1D grid-point depth        [m]
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pe3t_1d , pe3w_1d    ! 1D vertical scale factors  [m]
      !
      INTEGER  ::   jk       ! dummy loop indices
      REAL(wp) ::   zt, zw   ! local scalars
      !!----------------------------------------------------------------------
      !
      !
      IF(lwp) THEN            ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) '    zgr_z   : Reference vertical z-coordinates '
         WRITE(numout,*) '    ~~~~~~~'
      ENDIF

      !
      ! 1D Reference z-coordinate    (using Madec & Imbard 1996 function)
      ! -------------------------
      !
      DO jk = 1, jpk          ! depth at T and W-points
         zw = REAL( jk , wp ) - 1.0_wp
         zt = REAL( jk , wp ) - 0.5_wp
         pdepw_1d(jk) = zw * rn_bathy / ( REAL(jpk,wp) - 1.0_wp )
         pdept_1d(jk) = zt * rn_bathy / ( REAL(jpk,wp) - 1.0_wp )
      END DO
      !
      !                       ! e3t and e3w from depth
      CALL depth_to_e3( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d ) 
      !
      !                       ! recompute depths from SUM(e3)  <== needed
      CALL e3_to_depth( pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d ) 
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) '              Reference 1D z-coordinate depth and scale factors:'
         WRITE(numout, "(9x,' level  gdept_1d  gdepw_1d  e3t_1d   e3w_1d  ')" )
         WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, pdept_1d(jk), pdepw_1d(jk), pe3t_1d(jk), pe3w_1d(jk), jk = 1, jpk )
      ENDIF
      !
   END SUBROUTINE zgr_z


   SUBROUTINE zgr_msk_top_bot( k_top , k_bot )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE zgr_msk_top_bot  ***
      !!
      !! ** Purpose :   set the masked top and bottom ocean t-levels
      !!
      !! ** Method  :   GYRE case = closed flat box ocean without ocean cavities
      !!                   k_top = 1     except along north, south, east and west boundaries
      !!                   k_bot = jpk-1 except along north, south, east and west boundaries
      !!
      !! ** Action  : - k_top : first wet ocean level index
      !!              - k_bot : last  wet ocean level index
      !!----------------------------------------------------------------------
      INTEGER , DIMENSION(:,:), INTENT(out) ::   k_top , k_bot   ! first & last wet ocean level
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! 2D local workspace
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '    zgr_top_bot : defines the top and bottom wet ocean levels.'
      IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) '       SEAMOUNT case : terrain-following k_bot = jpkm1 for ocean points'
      !
      z2d(:,:) = REAL( jpkm1 , wp )                              ! flat bottom
      !
      k_bot(:,:) = NINT( z2d(:,:) )          ! =jpkm1 over the ocean point, =0 elsewhere
      !
      k_top(:,:) = MIN( 1 , k_bot(:,:) )     ! = 1    over the ocean point, =0 elsewhere
      !
   END SUBROUTINE zgr_msk_top_bot
   

   SUBROUTINE zgr_zco( pht, pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d,   &            ! in : 1D reference vertical coordinate
      &                pdept   , pdepw   ,                     &                 ! out: 3D t & w-points depth
      &                pe3t    , pe3u    , pe3v , pe3f,        &                 ! out: 3D vertical scale factors
      &                pe3w    , pe3uw   , pe3vw,              &                 !          -      -      -
      &                pk_bot  , pk_top                      )                   ! out: 2D Top and bottom level arrays
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zgr_zco  ***
      !!
      !! ** Purpose :   define the reference z-coordinate system
      !!
      !! ** Method  :   set 3D coord. arrays to reference 1D array 
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   pht                         ! 1D grid-point depth       [m]
      REAL(wp), DIMENSION(:)    , INTENT(in   ) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth       [m]
      REAL(wp), DIMENSION(:)    , INTENT(in   ) ::   pe3t_1d , pe3w_1d           ! 1D vertical scale factors [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! 3D grid-point depth       [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! 3D vertical scale factors [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         !    -       -      -
      INTEGER , DIMENSION(:,:)  , INTENT(inout) ::   pk_bot, pk_top              ! 2D Top and bottom level arrays
      !
      INTEGER  ::   jk
      !!----------------------------------------------------------------------
      !
      DO jk = 1, jpkm1
         WHERE( pdept_1d(jk) < pht(:,:) .AND. pht(:,:) <= pdept_1d(jk+1) )   pk_bot(:,:) = jk * pk_top(:,:)
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
      !
   END SUBROUTINE zgr_zco

   SUBROUTINE zgr_zps( pht, pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d,   &            ! in : 1D reference vertical coordinate
      &                pdept   , pdepw   ,                     &                 ! out: 3D t & w-points depth
      &                pe3t    , pe3u    , pe3v , pe3f,        &                 ! out: 3D vertical scale factors
      &                pe3w    , pe3uw   , pe3vw,              &                 !          -      -      -
      &                pk_bot  , pk_top                      )                   ! out: 2D Top and bottom level arrays
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zgr_zps  ***
      !!
      !! ** Purpose :   define the z-coordinate system
      !!
      !! ** Method  :   as per zco but with partial steps at lowest wet level
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   pht                         ! 1D grid-point depth       [m]
      REAL(wp), DIMENSION(:)    , INTENT(in   ) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth       [m]
      REAL(wp), DIMENSION(:)    , INTENT(in   ) ::   pe3t_1d , pe3w_1d           ! 1D vertical scale factors [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! 3D grid-point depth       [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! 3D vertical scale factors [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         !    -       -      -
      INTEGER , DIMENSION(:,:)  , INTENT(inout) ::   pk_bot, pk_top              ! 2D Top and bottom level arrays
      !
      INTEGER  ::   ji,jj,jk,ik
      REAL(wp) ::   ze3min
      !!----------------------------------------------------------------------
      !
      !
      ze3min = 0.1_wp * rn_dz
      IF(lwp) WRITE(numout,*) '   minimum thickness of the partial cells = 10 % of e3 = ', ze3min
      !
      !
      !                                !* bottom ocean compute from the depth of grid-points
      pk_bot(:,:) = jpkm1
      DO jk = jpkm1, 1, -1
         WHERE( pht(:,:) < pdepw_1d(jk) + ze3min )   pk_bot(:,:) = jk-1
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
      DO jj = 1, jpj                      ! bottom scale factors and depth at T- and W-points
         DO ji = 1, jpi
            ik = pk_bot(ji,jj)
               pdepw(ji,jj,ik+1) = MIN( pht(ji,jj) , pdepw_1d(ik+1) )
               pe3t (ji,jj,ik  ) = pdepw(ji,jj,ik+1) - pdepw(ji,jj,ik)
               pe3t (ji,jj,ik+1) = pe3t (ji,jj,ik  )
               !
               pdept(ji,jj,ik  ) = pdepw(ji,jj,ik  ) + pe3t (ji,jj,ik  ) * 0.5_wp
               pdept(ji,jj,ik+1) = pdepw(ji,jj,ik+1) + pe3t (ji,jj,ik+1) * 0.5_wp
               pe3w (ji,jj,ik+1) = pdept(ji,jj,ik+1) - pdept(ji,jj,ik)              ! = pe3t (ji,jj,ik  )
         END DO
      END DO

      ! Scale factors and depth at U-, V-, UW and VW-points
      DO jk = 1, jpk                        ! initialisation to z-scale factors
         e3u_0 (:,:,jk) = e3t_1d(jk)
         e3v_0 (:,:,jk) = e3t_1d(jk)
         e3uw_0(:,:,jk) = e3w_1d(jk)
         e3vw_0(:,:,jk) = e3w_1d(jk)
      END DO

      DO jk = 1,jpk                         ! Computed as the minimum of neighbooring scale factors
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               pe3u (ji,jj,jk) = MIN( pe3t(ji,jj,jk), pe3t(ji+1,jj,jk) )
               pe3v (ji,jj,jk) = MIN( pe3t(ji,jj,jk), pe3t(ji,jj+1,jk) )
               pe3uw(ji,jj,jk) = MIN( pe3w(ji,jj,jk), pe3w(ji+1,jj,jk) )
               pe3vw(ji,jj,jk) = MIN( pe3w(ji,jj,jk), pe3w(ji,jj+1,jk) )
            END DO
         END DO
      END DO
      CALL lbc_lnk('domzgr', e3u_0 , 'U', 1._wp )   ;   CALL lbc_lnk('domzgr', e3uw_0, 'U', 1._wp )   ! lateral boundary conditions
      CALL lbc_lnk('domzgr', e3v_0 , 'V', 1._wp )   ;   CALL lbc_lnk('domzgr', e3vw_0, 'V', 1._wp )
      ! Scale factor at F-point
      DO jk = 1, jpk                        ! initialisation to z-scale factors
         e3f_0(:,:,jk) = e3t_1d(jk)
      END DO
      DO jk = 1, jpk                        ! Computed as the minimum of neighbooring V-scale factors
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               e3f_0(ji,jj,jk) = MIN( e3v_0(ji,jj,jk), e3v_0(ji+1,jj,jk) )
            END DO
         END DO
      END DO
      CALL lbc_lnk('domzgr', e3f_0, 'F', 1._wp )       ! Lateral boundary conditions
      !      
      !
   END SUBROUTINE zgr_zps

   SUBROUTINE zgr_sco( pht,   &   ! in : reference bathymetry
      &                pdept   , pdepw   ,                     &   ! out: 3D t & w-points depth
      &                pe3t    , pe3u    , pe3v   , pe3f   ,   &   !      vertical scale factors
      &                pe3w    , pe3uw   , pe3vw,              &   !
      &                pk_bot  , pk_top                        )   !          -      -      -
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zgr_sco  ***
      !!
      !! ** Purpose :   define the z-coordinate system
      !!
      !! ** Method  :    
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   pht                         ! 1D grid-point depth       [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! grid-point depth          [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors    [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         !    -       -      -
      INTEGER , DIMENSION(:,:)  , INTENT(inout) ::   pk_bot, pk_top              !    -       -      -
      !
      INTEGER  ::   ji,jj,jk
      REAL(wp) ::   z1_jpkm1, zmax
      REAL(wp), DIMENSION(jpi, jpj) ::   zhu, zhv 
      REAL(wp), DIMENSION(jpk) ::   sigt_1d, sigw_1d 
      !!----------------------------------------------------------------------
      !
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
         zhu(ji,jj) = 0.5_wp * ( pht(ji,jj) + pht(ji+1,jj) )
         zhv(ji,jj) = 0.5_wp * ( pht(ji,jj) + pht(ji,jj+1) )
         END DO
      END DO
      CALL lbc_lnk( 'usrdef_zgr', zhu, 'U', 1. )
      CALL lbc_lnk( 'usrdef_zgr', zhv, 'V', 1. )
      z1_jpkm1 = 1._wp / REAL( jpkm1 , wp)
      DO jk = 1, jpk
         sigt_1d(jk) = ( REAL (jk-1, wp) + 0.5_wp ) / REAL ( jpkm1 )
         sigw_1d(jk) =   REAL (jk-1, wp)            / REAL ( jpkm1 )
         IF( lwp ) WRITE(numout, *) 'sigt_1d(jk), sigw_1d(jk)', jk, sigt_1d(jk), sigw_1d(jk)
      END DO
      DO jj = 1, jpj
         DO ji = 1, jpi
            DO jk = 1, jpk
               pdept(ji,jj,jk) = pht(ji,jj) * sigt_1d(jk)
               pdepw(ji,jj,jk) = pht(ji,jj) * sigw_1d(jk)
            END DO
            DO jk = 1, jpkm1
               pe3t (ji,jj,jk  ) = pdepw(ji,jj,jk+1) - pdepw(ji,jj,jk)
               pe3w (ji,jj,jk+1) = pdept(ji,jj,jk+1) - pdept(ji,jj,jk)
            END DO
            pe3t (ji,jj,jpk) = 2._wp * ( pdept(ji,jj,jpk) - pdepw(ji,jj,jpk) )
            pe3w (ji,jj,1  ) = 2._wp * ( pdept(ji,jj,1  ) - pdepw(ji,jj,1  ) )
         END DO
      END DO
      DO jk = 1, jpk
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               pe3u (ji,jj,jk) = 0.5_wp  * ( pe3t (ji,jj,jk)   + pe3t (ji+1,jj,  jk) )
               pe3v (ji,jj,jk) = 0.5_wp  * ( pe3w (ji,jj,jk)   + pe3t (ji,jj+1,  jk) )
               pe3uw(ji,jj,jk) = 0.5_wp  * ( pe3w (ji,jj,jk)   + pe3t (ji+1,jj,  jk) )
               pe3vw(ji,jj,jk) = 0.5_wp  * ( pe3w (ji,jj,jk)   + pe3t (ji,jj+1,  jk) )
               pe3f (ji,jj,jk) = 0.25 * ( pe3t (ji,jj,jk)   + pe3t (ji+1,jj,  jk) &
                                     &  + pe3t (ji,jj+1,jk) + pe3t (ji+1,jj+1,jk) ) 
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
      !
   END SUBROUTINE zgr_sco

  

   !!======================================================================
END MODULE usrdef_zgr
