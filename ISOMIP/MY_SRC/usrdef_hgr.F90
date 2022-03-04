MODULE usrdef_hgr
   !!======================================================================
   !!                       ***  MODULE usrdef_hgr   ***
   !!
   !!                  ===  LOCK_EXCHANGE configuration  ===
   !!
   !! User defined :   mesh and Coriolis parameter of a user configuration
   !!======================================================================
   !! History :  NEMO  ! 2016-08  (S. Flavoni, G. Madec)   Original code
   !!                  ! 2017-02  (P. Mathiot, S. Flavoni) Adapt code to ISOMIP case
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_hgr    : initialize the horizontal mesh for ISOMIP configuration
   !!----------------------------------------------------------------------
   USE dom_oce
   USE par_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE usrdef_nam, ONLY: rn_e1deg, rn_e2deg   ! horizontal resolution in meters
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
      !!                ISOMIP configuration
      !!
      !! ** Action  : - define longitude & latitude of t-, u-, v- and f-points (in degrees) 
      !!              - define i- & j-scale factors at t-, u-, v- and f-points (in meters)
      !!              - define u- & v-surfaces (in m2)
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
      REAL(wp) ::   zfact, zti, zui, zvi, zfi, ztj, zuj, zvj, zfj      ! local scalars
      !!-------------------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'usr_def_hgr : ISOMIP configuration'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*)
         WRITE(numout,*) '   ===>>  geographical mesh on the sphere with regular grid-spacing'
         WRITE(numout,*) '          given by rn_e1deg and rn_e2deg'
      ENDIF
      !
      !                       !==  grid point position  ==!   (in degrees)
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         !                       ! longitude   (west coast at lon=0°)
         plamt(ji,jj) = rn_e1deg * (  - 0.5 + REAL( mig0(ji)-1 , wp )  )  
         plamu(ji,jj) = rn_e1deg * (          REAL( mig0(ji)-1 , wp )  )
         plamv(ji,jj) = plamt(ji,jj)
         plamf(ji,jj) = plamu(ji,jj)
         !                       ! latitude   (south coast at lat=-80°)
         pphit(ji,jj) = rn_e2deg * (  - 0.5 + REAL( mjg0(jj)-1 , wp )  ) - 80._wp
         pphiu(ji,jj) = pphit(ji,jj)
         pphiv(ji,jj) = rn_e2deg * (          REAL( mjg0(jj)-1 , wp )  ) - 80._wp
         pphif(ji,jj) = pphiv(ji,jj)
      END_2D
      !
      !                       !==  Horizontal scale factors  ==!   (in meters)
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         !                       ! e1   (zonal)
         pe1t(ji,jj) = ra * rad * COS( rad * pphit(ji,jj) ) * rn_e1deg
         pe1u(ji,jj) = ra * rad * COS( rad * pphiu(ji,jj) ) * rn_e1deg
         pe1v(ji,jj) = ra * rad * COS( rad * pphiv(ji,jj) ) * rn_e1deg
         pe1f(ji,jj) = ra * rad * COS( rad * pphif(ji,jj) ) * rn_e1deg
         !                       ! e2   (meridional)
         pe2t(ji,jj) = ra * rad * rn_e2deg
         pe2u(ji,jj) = ra * rad * rn_e2deg
         pe2v(ji,jj) = ra * rad * rn_e2deg
         pe2f(ji,jj) = ra * rad * rn_e2deg
      END_2D
      !                             ! NO reduction of grid size in some straits 
      ke1e2u_v    = 0               !    ==>> u_ & v_surfaces will be computed in dom_ghr routine
      pe1e2u(:,:) = 0._wp           !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v(:,:) = 0._wp           !             require an initialization of INTENT(out) arguments
      !
      !
      !                       !==  Coriolis parameter  ==!
      kff = 0                       ! Coriolis parameter calculated on the sphere
      pff_f(:,:) = 0._wp            ! CAUTION: set to zero to avoid error with some compilers that
      pff_t(:,:) = 0._wp            !             require an initialization of INTENT(out) arguments
      !
   END SUBROUTINE usr_def_hgr

   !!======================================================================
END MODULE usrdef_hgr
