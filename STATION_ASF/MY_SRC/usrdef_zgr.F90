MODULE usrdef_zgr
   !!======================================================================
   !!                       ***  MODULE  usrdef_zgr  ***
   !!
   !!                       ===  STATION_ASF case  ===
   !!
   !! User defined : vertical coordinate system of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-06  (G. Madec)  Original code
   !!            4.x  ! 2019-10  (L. Brodeau) Station ASF
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_zgr   : user defined vertical coordinate system
   !!      zgr_z      : reference 1D z-coordinate 
   !!      zgr_top_bot: ocean top and bottom level indices
   !!      zgr_zco    : 3D verticl coordinate in pure z-coordinate case
   !!---------------------------------------------------------------------
   USE oce            ! ocean variables
   USE usrdef_nam     ! User defined : namelist variables
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_zgr        ! called by domzgr.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_zgr.F90 12377 2020-02-12 14:39:06Z acc $
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
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr : STATION_ASF configuration, setting first level properties.'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      !
      ! type of vertical coordinate
      ! ---------------------------
      ld_zco    = .TRUE.         ! z-coordinate without ocean cavities
      ld_zps    = .FALSE.
      ld_sco    = .FALSE.
      ld_isfcav = .FALSE.

      !! 1st level (the only one that matters)
      pdept_1d(1) = rn_dept1 ! depth (m) at which the SST is taken/measured == depth of first T point!
      pdepw_1d(1) = 0._wp
      pe3t_1d(1)  = 2._wp*rn_dept1
      pe3w_1d(1)  = rn_dept1 ! LB???

      pdept(:,:,1) = rn_dept1
      pdepw(:,:,1) = 0._wp
      pe3t(:,:,1) = 2._wp*rn_dept1
      pe3u(:,:,1) = 2._wp*rn_dept1
      pe3v(:,:,1) = 2._wp*rn_dept1
      pe3f(:,:,1) = 2._wp*rn_dept1
      pe3w(:,:,1)  = rn_dept1  ! LB???
      pe3uw(:,:,1) = rn_dept1  ! LB???
      pe3vw(:,:,1) = rn_dept1  ! LB???
      
      !! 2nd level, technically useless (only for the sake of code stability)
      pdept_1d(2) = 3._wp*rn_dept1
      pdepw_1d(2) = 2._wp*rn_dept1
      pe3t_1d(2)  = 2._wp*rn_dept1
      pe3w_1d(2)  = 2._wp*rn_dept1

      pdept(:,:,2) = 3._wp*rn_dept1
      pdepw(:,:,2) = 2._wp*rn_dept1
      pe3t(:,:,2) = 2._wp*rn_dept1
      pe3u(:,:,2) = 2._wp*rn_dept1
      pe3v(:,:,2) = 2._wp*rn_dept1
      pe3f(:,:,2) = 2._wp*rn_dept1
      pe3w(:,:,2)  = 2._wp*rn_dept1
      pe3uw(:,:,2) = 2._wp*rn_dept1
      pe3vw(:,:,2) = 2._wp*rn_dept1

      k_top = 1
      k_bot = 1

   END SUBROUTINE usr_def_zgr
   !!======================================================================
END MODULE usrdef_zgr
