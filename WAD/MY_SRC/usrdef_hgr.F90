MODULE usrdef_hgr
   !!======================================================================
   !!                       ***  MODULE usrdef_hgr   ***
   !!
   !!                  ===  WAD_TEST_CASES configuration  ===
   !!
   !! User defined :   mesh and Coriolis parameter of a user configuration
   !!======================================================================
   !! History :  4.0 ! 2016-03  (S. Flavoni) 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_hgr    : initialize the horizontal mesh for WAD_TEST_CASES configuration
   !!----------------------------------------------------------------------
   USE dom_oce
   USE par_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE usrdef_nam, ONLY: rn_dx   ! horizontal resolution in meters
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_hgr   ! called by domhgr.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_hgr.F90 14223 2020-12-19 10:22:45Z smasson $ 
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
      !!                WAD_TEST_CASES configuration : uniform grid spacing (rn_dx)
      !!                without Coriolis force (f=0)
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
      INTEGER  ::   ji, jj   ! dummy loop indices
      REAL(wp) ::   zfact      ! local scalars
      !!-------------------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_hgr : WAD_TEST_CASES configuration basin'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   uniform grid spacing WITHOUT Coriolis force (f=0)'
      !
      !                       !==  grid point position  ==!   (in kilometers)
      zfact = rn_dx * 1.e-3         ! conversion in km
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         !                       ! longitude   (west coast at lon=0°)
         plamt(ji,jj) = zfact * (  - 0.5 + REAL( mig0(ji)-1 , wp )  )  
         plamu(ji,jj) = zfact * (          REAL( mig0(ji)-1 , wp )  )
         plamv(ji,jj) = plamt(ji,jj)
         plamf(ji,jj) = plamu(ji,jj)
         !                       ! latitude   (south coast at lat= 0°)
         pphit(ji,jj) = zfact * (  - 0.5 + REAL( mjg0(jj)-1 , wp )  )
         pphiu(ji,jj) = pphit(ji,jj)
         pphiv(ji,jj) = zfact * (          REAL( mjg0(jj)-1 , wp )  )
         pphif(ji,jj) = pphiv(ji,jj)
      END_2D
      !
      !                       !==  Horizontal scale factors  ==!   (in meters) 
      pe1t(:,:) = rn_dx   ;   pe2t(:,:) = rn_dx
      pe1u(:,:) = rn_dx   ;   pe2u(:,:) = rn_dx
      pe1v(:,:) = rn_dx   ;   pe2v(:,:) = rn_dx
      pe1f(:,:) = rn_dx   ;   pe2f(:,:) = rn_dx
      !
      !                             ! NO reduction of grid size in some straits 
      ke1e2u_v = 0                  !    ==>> u_ & v_surfaces will be computed in dom_ghr routine
      pe1e2u(:,:) = 0._wp           !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v(:,:) = 0._wp           !             require an initialization of INTENT(out) arguments
      !
      !
      !                       !==  Coriolis parameter  ==!
      kff = 1                       !  indicate not to compute Coriolis parameter afterward
      !
      pff_f(:,:) = 0._wp            ! here No earth rotation: f=0
      pff_t(:,:) = 0._wp
      !
   END SUBROUTINE usr_def_hgr

   !!======================================================================
END MODULE usrdef_hgr
