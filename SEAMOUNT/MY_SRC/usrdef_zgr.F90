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
      &                    k_top  , k_bot                              )   ! top & bottom ocean level
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
      INTEGER                                 ::   ji, jj, jk                  ! dummy indices
      INTEGER                                 ::   ik                          ! local integers
      REAL(wp)                                ::   zlam_mid, zphi_mid          ! local scalar
      REAL(wp), DIMENSION(jpi,jpj)            ::   bathy                       ! model bathymetry
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
      zlam_mid = 0.5_wp * 1000._wp * rn_xdim
      zphi_mid = 0.5_wp * 1000._wp * rn_ydim
      bathy(:,:) = 0._wp
      DO jj = 1, jpj
         DO ji = 1, jpi
            bathy(ji,jj) = rn_bot_max - rn_smnt_H * EXP(                 &
               &           -(( 1000._wp * glamt(ji,jj) - zlam_mid)**2 + &
               &             ( 1000._wp * gphit(ji,jj) - zphi_mid)**2 ) & 
               &           / ( 1000._wp * rn_smnt_L)**2 )
         END DO
      END DO
      ! 
      ! ------------------------------------
      CALL zgr_z( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! Reference z-coordinate system
      !
      CALL zgr_msk_top_bot( k_top , k_bot )                 ! masked top and bottom ocean t-level indices
      !
      IF ( ln_zco ) CALL zgr_zco( bathy    ,                                    &  ! in : 2D bathymetry
            &                     pdept_1d , pdepw_1d, pe3t_1d , pe3w_1d,       &  !      1D ref. z-coord.
            &                     pdept    , pdepw   ,                          &  ! out: 3D T & W-points depth
            &                     pe3t     , pe3u    , pe3v    , pe3f   ,       &  !      vertical scale factors
            &                     pe3w     , pe3uw   , pe3vw   , k_bot  , k_top )  
      !
      IF ( ln_zps ) CALL zgr_zps( bathy    ,                                    &  ! in : 2D bathymetry
            &                     pdept_1d , pdepw_1d, pe3t_1d , pe3w_1d,       &  !      1D ref. z-coord.
            &                     pdept    , pdepw   ,                          &  ! out: 3D T & W-points depth
            &                     pe3t     , pe3u    , pe3v    , pe3f   ,       &  !      vertical scale factors
            &                     pe3w     , pe3uw   , pe3vw   , k_bot  , k_top )    
      !
      IF ( ln_sco ) CALL zgr_sco( bathy    ,                                    &  ! in : 2D bathymetry
            &                     pdept    , pdepw   ,                          &  ! out: 3D T & W-points depth
            &                     pe3t     , pe3u    , pe3v    , pe3f   ,       &  !      vertical scale factors
            &                     pe3w     , pe3uw   , pe3vw   , k_bot  , k_top )

      !
   END SUBROUTINE usr_def_zgr


   SUBROUTINE zgr_z( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )
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
         pdepw_1d(jk) = zw * rn_bot_max / ( REAL(jpk,wp) - 1.0_wp )
         pdept_1d(jk) = zt * rn_bot_max / ( REAL(jpk,wp) - 1.0_wp )
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
   

   SUBROUTINE zgr_zco( pht, pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d, &  
      &                pdept   , pdepw   ,                        &  
      &                pe3t    , pe3u    , pe3v , pe3f,           &
      &                pe3w    , pe3uw   , pe3vw,                 &
      &                pk_bot  , pk_top                           )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zgr_zco  ***
      !!
      !! ** Purpose :   define the reference z-coordinate system
      !!
      !! ** Method  :   set 3D coord. arrays to reference 1D array 
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   pht                         ! 2D bathymetry             [m]
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

   SUBROUTINE zgr_zps( pht     , pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d,   &   
      &                pdept   , pdepw   ,                               &        
      &                pe3t    , pe3u    , pe3v    , pe3f   ,            &        
      &                pe3w    , pe3uw   , pe3vw   ,                     &        
      &                pk_bot  , pk_top                                  )          
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zgr_zps  ***
      !!
      !! ** Purpose :   define the z-coordinate system
      !!
      !! ** Method  :   as per zco but with partial steps at lowest wet level
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   pht                         ! 2D bathymetry             [m]
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

      ! lateral boundary conditions
      CALL lbc_lnk('domzgr', e3u_0 , 'U', 1._wp ) 
      CALL lbc_lnk('domzgr', e3uw_0, 'U', 1._wp )
      CALL lbc_lnk('domzgr', e3v_0 , 'V', 1._wp )
      CALL lbc_lnk('domzgr', e3vw_0, 'V', 1._wp )

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

   SUBROUTINE zgr_sco( pht,                                    & 
      &                pdept   , pdepw   ,                     & 
      &                pe3t    , pe3u    , pe3v   , pe3f   ,   & 
      &                pe3w    , pe3uw   , pe3vw,              & 
      &                pk_bot  , pk_top                        )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zgr_sco  ***
      !!
      !! ** Purpose :   define the z-coordinate system
      !!
      !! ** Method  :    
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   pht                         ! 2D bathymetry          [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! grid-point depth       [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         !    -       -      -
      INTEGER , DIMENSION(:,:)  , INTENT(inout) ::   pk_bot, pk_top              !    -       -      -
      !
      INTEGER                                   ::   ji,jj,jk
      REAL(wp), DIMENSION(jpk)                  ::   sigt_1d, sigw_1d
      REAL(wp), DIMENSION(jpi,jpj)              ::   pht_w 
      !!----------------------------------------------------------------------
      !
      pht_w(:,:) = pht(:,:)
      ! Computing envelope bathymetry if using vanishing quasi-sigma levels
      IF( ln_vqs ) THEN 
        CALL s_vqs( pht, pht_w )
      ELSE
        pht_w(:,:) = pht(:,:)
      END IF

      ! Computing uniform sigma-coordinate 
      CALL sigma_coord(sigt_1d, sigw_1d)

      ! Computing SH94 stretched s-coordinate if requested
      IF( ln_s_sh94 ) CALL sh94_coord(sigt_1d, sigw_1d)

      ! Forcing to zero rn_hc in the case of 
      ! uniform sigma-coordinates
      IF( .NOT. ln_s_sh94 ) rn_hc = 0.0

      ! Computing depth of model levels
      ! and vertical scale factors
      DO jj = 1, jpj
         DO ji = 1, jpi
            ! Depth of model levels
            DO jk = 1, jpk
               pdept(ji,jj,jk) = rn_hc * sigt_1d(jk) + sigt_1d(jk) * (pht_w(ji,jj) - rn_hc)
               pdepw(ji,jj,jk) = rn_hc * sigw_1d(jk) + sigw_1d(jk) * (pht_W(ji,jj) - rn_hc)
            END DO
            ! Vertical scale factors as finite differences
            DO jk = 1, jpkm1
               pe3t (ji,jj,jk  ) = pdepw(ji,jj,jk+1) - pdepw(ji,jj,jk)
               pe3w (ji,jj,jk+1) = pdept(ji,jj,jk+1) - pdept(ji,jj,jk)
            END DO
            pe3t (ji,jj,jpk) = 2._wp * ( pdept(ji,jj,jpk) - pdepw(ji,jj,jpk) )
            pe3w (ji,jj,1  ) = 2._wp * ( pdept(ji,jj,1  ) - pdepw(ji,jj,1  ) )
         END DO
      END DO
 
      ! Adjust bottom ocean in case of VQS
      IF( ln_vqs ) THEN
        DO jj = 1, jpj
           DO ji = 1, jpi 
              DO jk = jpkm1, 1, -1
                 IF ( pht_w(ji,jj) < pdept(ji,jj,jk) ) THEN   
                    pk_bot(ji,jj) = jk-1
                    CYCLE
                 ENDIF
              END DO
           END DO
        END DO
      ENDIF

      ! Compute vertical scale factors for all grids
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

   SUBROUTINE s_vqs( pht, zenv )
      !!------------------------------------------------------------------------
      !!                  ***  ROUTINE s_vqs  ***
      !!
      !! ** Purpose :   compute the envelope bathymetry that will be used to
      !!                to define vanishing quasi-sigma (VQS) levels
      !!
      !! ** Method  :   Direct iterative method of Martinho and Batteen (2006).
      !!                The algorithm ensures that
      !!
      !!                              H_ij - H_n
      !!                              ---------- < rmax
      !!                              H_ij + H_n
      !!  
      !!                where H_ij is the depth at point (i,j) and H_n is the
      !!                neighbouring depth in the east, west, south or north 
      !!                direction.
      !! 
      !! Reference:     Martinho & Batteen, Oce. Mod. 13(2):166-175, 2006.
      !!------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   pht            ! 2D bathymetry [m]
      REAL(wp), DIMENSION(:,:)  , INTENT(  out) ::   zenv           ! 2D envelope   [m]
      !
      INTEGER                                   ::   ji,jj,jk,jl
      INTEGER                                   ::   iip1, ijp1
      REAL(wp)                                  ::   zrmax, zrfact
      REAL(wp), DIMENSION(jpi, jpj)             ::   ztmpi1, ztmpi2
      REAL(wp), DIMENSION(jpi, jpj)             ::   ztmpj1, ztmpj2
      REAL(wp), DIMENSION(jpi, jpj)             ::   zri, zrj
      !!----------------------------------------------------------------------
      !
      zenv(:,:) = pht(:,:)
      
      jl = 0
      zrmax = 1._wp
      !      
      ! set scaling factor used in reducing vertical gradients
      zrfact = ( 1._wp - rn_rmax ) / ( 1._wp + rn_rmax )
      !
      ! initialise temporary evelope depth arrays
      ztmpi1(:,:) = zenv(:,:)
      ztmpi2(:,:) = zenv(:,:)
      ztmpj1(:,:) = zenv(:,:)
      ztmpj2(:,:) = zenv(:,:)
      !
      ! initialise temporary r-value arrays
      zri(:,:) = 1._wp
      zrj(:,:) = 1._wp
      !
      DO WHILE( jl <= 10000 .AND. ( zrmax - rn_rmax ) > 1.e-8_wp ) !  Iterative loop  !
         !                                                         !
         !                                                         ================
         !                                                         !
         jl = jl + 1
         zrmax = 0._wp
         ! we set zrmax from previous r-values (zri and zrj) first
         ! if set after current r-value calculation (as previously)
         ! we could exit DO WHILE prematurely before checking r-value
         ! of current zenv
         DO jj = 1, jpj
            DO ji = 1, jpi
               zrmax = MAX( zrmax, ABS(zri(ji,jj)), ABS(zrj(ji,jj)) )
            END DO
         END DO
         zri(:,:) = 0._wp
         zrj(:,:) = 0._wp
         DO jj = 1, jpj
            DO ji = 1, jpi
               iip1 = MIN( ji+1, jpi )      ! force zri = 0 on last line (ji=ncli+1 to jpi)
               ijp1 = MIN( jj+1, jpj )      ! force zrj = 0 on last raw  (jj=nclj+1 to jpj)
               IF( (zenv(ji,jj) > 0._wp) .AND. (zenv(iip1,jj) > 0._wp)) THEN
                  zri(ji,jj) = ( zenv(iip1,jj  ) - zenv(ji,jj) ) / &
                    &          ( zenv(iip1,jj  ) + zenv(ji,jj) )
               END IF
               IF( (zenv(ji,jj) > 0._wp) .AND. (zenv(ji,ijp1) > 0._wp)) THEN
                  zrj(ji,jj) = ( zenv(ji,ijp1) - zenv(ji,jj) ) / &
                    &          ( zenv(ji,ijp1) + zenv(ji,jj) )
               END IF
               IF( zri(ji,jj) >  rn_rmax )   ztmpi1(ji  ,jj  ) = zenv(iip1,jj  ) * zrfact
               IF( zri(ji,jj) < -rn_rmax )   ztmpi2(iip1,jj  ) = zenv(ji  ,jj  ) * zrfact
               IF( zrj(ji,jj) >  rn_rmax )   ztmpj1(ji  ,jj  ) = zenv(ji  ,ijp1) * zrfact
               IF( zrj(ji,jj) < -rn_rmax )   ztmpj2(ji  ,ijp1) = zenv(ji  ,jj  ) * zrfact
            END DO
         END DO
         IF(lwp)WRITE(numout,*) 's_vqs :   iter= ',jl, ' rmax= ', zrmax
         !
         DO jj = 1, jpj
            DO ji = 1, jpi
               zenv(ji,jj) = MAX(zenv(ji,jj), ztmpi1(ji,jj), ztmpi2(ji,jj), ztmpj1(ji,jj), ztmpj2(ji,jj) )
            END DO
         END DO
         ! Apply lateral boundary condition
         CALL lbc_lnk( 'usrdef_zgr', zenv, 'T', 1. )
         ! CAUTION: keep the value when the lbc field is zero
         !CALL lbc_lnk( zenv, 'T', 1._wp, 'no0' )         
         !                                                  ! ================ !
      END DO                                                !     End loop     !               
      !
   END SUBROUTINE s_vqs

   SUBROUTINE sigma_coord( sigT, sigW )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE sigma  ***
      !!       
      !! ** Purpose :   provide the analytical function for sigma-coordinate
      !!                (not stretched s-coordinate).
      !!          
      !! ** Method  :   the function provide the non-dimensional position of
      !!                T and W points (i.e. between 0 and 1).
      !!
      !!----------------------------------------------------------------------
      REAL,   DIMENSION(:),   INTENT (  out) ::   sigT, sigW ! sigma coordinate 
                                                             ! at T and W points
      !
      INTEGER                                ::   jk
      !!----------------------------------------------------------------------

      DO jk = 1, jpk
         sigT(jk) = ( REAL (jk-1, wp) + 0.5_wp ) / REAL ( jpkm1 )
         sigW(jk) =   REAL (jk-1, wp)            / REAL ( jpkm1 )
         IF( lwp ) WRITE(numout, *) 'sigt_1d(jk), sigw_1d(jk)', jk, sigT(jk), sigW(jk)
      END DO

   END SUBROUTINE sigma_coord

   SUBROUTINE sh94_coord( sigT, sigW )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE sh94_coord ***
      !!
      !! ** Purpose :   provide the Song and Haidvogel 1994 analytical 
      !!                stretching function for s-coordinate.
      !!
      !!----------------------------------------------------------------------
      REAL,   DIMENSION(:),   INTENT (inout) ::   sigT, sigW ! uniform sigma-coord 
      INTEGER                                ::   jk
      REAL                                   ::   sT, sW     ! work variables
      !!----------------------------------------------------------------------                                          
      IF ( rn_theta > 0 ) THEN
         DO jk = 1, jpk 
            sT = sigT(jk)
            sW = sigW(jk)     
            sigT(jk) = (1._wp - rn_bb) * SINH(rn_theta * sT) / SINH(rn_theta) + rn_bb * &
              &        ( ( TANH(rn_theta * (sT + 0.5_wp)) - TANH(0.5_wp * rn_theta) ) / &
              &        (2._wp * TANH(0.5_wp * rn_theta) ) )
            sigW(jk) = (1._wp - rn_bb) * SINH(rn_theta * sW) / SINH(rn_theta) + rn_bb * &
              &        ( ( TANH(rn_theta * (sW + 0.5_wp)) - TANH(0.5_wp * rn_theta) ) / &
              &        (2._wp * TANH(0.5_wp * rn_theta) ) )
         END DO
      END IF

   END SUBROUTINE sh94_coord
   !!======================================================================
END MODULE usrdef_zgr
