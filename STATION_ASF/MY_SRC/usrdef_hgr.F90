MODULE usrdef_hgr
   !!======================================================================
   !!                     ***  MODULE usrdef_hgr   ***
   !!
   !!                     ===  STATION_ASF configuration  ===
   !!
   !! User defined :   mesh and Coriolis parameter of a user configuration
   !!======================================================================
   !! History :  4.0 ! 2016-03  (S. Flavoni)
   !! History :  4.0 ! 2019-10  (L. Brodeau)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_hgr   : initialize the horizontal mesh
   !!----------------------------------------------------------------------
   USE dom_oce  , ONLY: nimpp, njmpp       ! ocean space and time domain
   USE c1d      ,  ONLY: rn_lon1d, rn_lat1d ! ocean lon/lat define by namelist
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE usrdef_nam     !
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_hgr   ! called in domhgr.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_hgr.F90 13216 2020-07-02 09:25:49Z rblod $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_hgr( plamt , plamu , plamv  , plamf  ,   &   ! geographic position (required)
      &                    pphit , pphiu , pphiv  , pphif  ,   &   !
      &                    kff   , pff_f , pff_t  ,            &   ! Coriolis parameter  (if domain not on the sphere)
      &                    pe1t  , pe1u  , pe1v   , pe1f   ,   &   ! scale factors       (required)
      &                    pe2t  , pe2u  , pe2v   , pe2f   ,   &   !
      &                    ke1e2u_v      , pe1e2u , pe1e2v     )   ! u- & v-surfaces (if gridsize reduction is used in strait(s))
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE usr_def_hgr  ***
      !!
      !! ** Purpose :   user defined mesh and Coriolis parameter
      !!
      !! ** Method  :   set all intent(out) argument to a proper value
      !!
      !!                Here STATION_ASF configuration :
      !!          Rectangular 3x3 domain
      !!          - Located at 150E-50N
      !!          - a constant horizontal resolution
      !!
      !! ** Action  : - define longitude & latitude of t-, u-, v- and f-points (in degrees)
      !!              - define coriolis parameter at f-point if the domain in not on the sphere (on beta-plane)
      !!              - define i- & j-scale factors at t-, u-, v- and f-points (in meters)
      !!              - define u- & v-surfaces (if gridsize reduction is used in some straits) (in m2)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   plamt, plamu, plamv, plamf   ! longitude outputs                     [degrees]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pphit, pphiu, pphiv, pphif   ! latitude outputs                      [degrees]
      INTEGER                 , INTENT(out) ::   kff                          ! =1 Coriolis parameter computed here, =0 otherwise
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pff_f, pff_t                 ! Coriolis factor at f-point                [1/s]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1t, pe1u, pe1v, pe1f       ! i-scale factors                             [m]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe2t, pe2u, pe2v, pe2f       ! j-scale factors                             [m]
      INTEGER                 , INTENT(out) ::   ke1e2u_v                     ! =1 u- & v-surfaces computed here, =0 otherwise
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1e2u, pe1e2v               ! u- & v-surfaces (if reduction in strait)   [m2]
      !
      REAL(wp) :: zf0   ! local scalar
      !!-------------------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_hgr : STATION_ASF configuration'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   uniform grid spacing WITH Coriolis force'
      IF(lwp) WRITE(numout,*) '              coordinates of station =>', rn_lon1d, rn_lat1d
      !
      !                       !==  grid point position  ==!   (Special case of 1D config: all lon/lat are the same)
      !                       ! longitude
      plamt(:,:) = rn_lon1d
      plamu(:,:) = rn_lon1d
      plamv(:,:) = rn_lon1d
      plamf(:,:) = rn_lon1d
      !                       ! latitude
      pphit(:,:) = rn_lat1d
      pphiu(:,:) = rn_lat1d
      pphiv(:,:) = rn_lat1d
      pphif(:,:) = rn_lat1d

      !                       !== Horizontal scale factors ==! (in meters)
      !                                         ! constant grid spacing
      pe1t(:,:) = 100.  ;   pe2t(:,:) = 100.
      pe1u(:,:) = 100.  ;   pe2u(:,:) = 100.
      pe1v(:,:) = 100.  ;   pe2v(:,:) = 100.
      pe1f(:,:) = 100.  ;   pe2f(:,:) = 100.
      !
      !                                         ! NO reduction of grid size in some straits
      ke1e2u_v = 0                              !    ==>> u_ & v_surfaces will be computed in dom_ghr routine
      pe1e2u(:,:) = 0._wp                       !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v(:,:) = 0._wp                       !             require an initialization of INTENT(out) arguments
      !
      !
      !                       !==  Coriolis parameter  ==!
      zf0   = 2._wp * omega * SIN( rad * rn_lat1d )
      pff_f(:,:) = zf0
      pff_t(:,:) = zf0
      kff = 1                 !  indicate to skip computing Coriolis parameter afterward
      !
   END SUBROUTINE usr_def_hgr

   !!======================================================================
END MODULE usrdef_hgr
