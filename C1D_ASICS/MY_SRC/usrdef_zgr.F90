MODULE usrdef_zgr
   !!======================================================================
   !!                       ***  MODULE  usrdef_zgr  ***
   !!
   !!                       ===  C1D_ASICS configuration  ===
   !!
   !! User defined : vertical coordinate system of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-06  (R. Bourdalle-Badie)  Original code
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
   USE usrdef_nam     ! User defined : namelist variables
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_zgr        ! called by domzgr.F90

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
      REAL(wp) ::   zt, zw            ! local scalars
      REAL(wp) ::   zsur, za0, za1, zkth, zacr        ! Values for the Madec & Imbard (1996) function
      REAL(wp) ::   za2, zkth2, zacr2                 ! Values for optional double tanh function set from parameters
      REAL(wp), DIMENSION(jpi,jpj) ::   zht, zhu, z2d ! 2D workspace

      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr : C1D configuration (zps-coordinate closed box ocean without cavities)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      ! type of vertical coordinate
      ! ---------------------------
      ld_zco    = .FALSE.         ! C1D case:  z-coordinate without ocean cavities
      ld_zps    = .TRUE.
      ld_sco    = .FALSE.
      ld_isfcav = .FALSE.
      !
      ! Build the vertical coordinate system
      ! ------------------------------------
      !
      ! Set parameters of z(k) function
      ! -------------------------------
      zsur =   -3958.95137127683
      za0  =    103.953009600000
      za1  =    2.41595126900000
      zkth =    15.3510137000000
      zacr =    7.00000000000000
      za2  =    100.760928500000
      zkth2=    48.0298937200000
      zacr2=    13.0000000000000
      !
      IF(lwp) THEN            ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) '     zgr_z75L   : Reference vertical z-coordinates '
         WRITE(numout,*) '     ~~~~~~~'
         WRITE(numout,*) '       C1D case : L75 function with the following coefficients :'
         WRITE(numout,*) '                 zsur = ', zsur
         WRITE(numout,*) '                 za0  = ', za0
         WRITE(numout,*) '                 za1  = ', za1
         WRITE(numout,*) '                 zkth = ', zkth
         WRITE(numout,*) '                 zacr = ', zacr
         WRITE(numout,*) '                 za2  = ', za2
         WRITE(numout,*) '                 zkth2= ', zkth2
         WRITE(numout,*) '                 zacr2= ', zacr2
      ENDIF

      !                       !==  UNmasked meter bathymetry  ==!
      !
      zht(:,:) = rn_bathy
      !
      DO jk = 1, jpk          ! depth at T and W-points
         zw = REAL( jk , wp )
         zt = REAL( jk , wp ) + 0.5_wp
         pdepw_1d(jk) = ( zsur + za0 * zw + za1 * zacr * LOG ( COSH( (zw-zkth ) / zacr  ) )    &
                  &                    + za2 * zacr2* LOG ( COSH( (zw-zkth2) / zacr2 ) )  )
         pdept_1d(jk) = ( zsur + za0 * zt + za1 * zacr * LOG ( COSH( (zt-zkth ) / zacr  ) )    &
                  &                    + za2 * zacr2* LOG ( COSH( (zt-zkth2) / zacr2 ) )  )
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
      !                       !==  top masked level bathymetry  ==!  (all coordinates)
      !
      k_top(:,:) = 1 
      !                                   !* bottom ocean compute from the depth of grid-points
      k_bot(:,:) = jpkm1
      DO jk = jpkm1, 1, -1
        ze3min = 0.1_wp * pe3t_1d (jk)
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
      ! bottom scale factors and depth at T- and W-points
      DO_2D( 1, 1, 1, 1 )
         ik = k_bot(ji,jj)
         pdepw(ji,jj,ik+1) = MIN( zht(ji,jj) , pdepw_1d(ik+1) )
         pe3t (ji,jj,ik  ) = pdepw(ji,jj,ik+1) - pdepw(ji,jj,ik)
         pe3t (ji,jj,ik+1) = pe3t (ji,jj,ik  ) 
         !
         pdept(ji,jj,ik  ) = pdepw(ji,jj,ik  ) + pe3t (ji,jj,ik  ) * 0.5_wp
         pdept(ji,jj,ik+1) = pdepw(ji,jj,ik+1) + pe3t (ji,jj,ik+1) * 0.5_wp
         pe3w (ji,jj,ik+1) = pdept(ji,jj,ik+1) - pdept(ji,jj,ik)              ! = pe3t (ji,jj,ik  )
      END_2D        
      !                                   ! bottom scale factors and depth at  U-, V-, UW and VW-points
      !                                   ! usually Computed as the minimum of neighbooring scale factors
      pe3u (:,:,:) = pe3t(:,:,:)          ! HERE C1D configuration : 
      pe3v (:,:,:) = pe3t(:,:,:)          !    e3 increases with k-index 
      pe3f (:,:,:) = pe3t(:,:,:)          !    so e3 minimum of (i,i+1) points is (i) point
      pe3uw(:,:,:) = pe3w(:,:,:)          !    in j-direction e3v=e3t and e3f=e3v
      pe3vw(:,:,:) = pe3w(:,:,:)          !    ==>>  no need of lbc_lnk calls
      !      
      !
   END SUBROUTINE usr_def_zgr

   !!======================================================================
END MODULE usrdef_zgr
