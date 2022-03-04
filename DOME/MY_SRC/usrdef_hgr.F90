MODULE usrdef_hgr
   !!======================================================================
   !!                       ***  MODULE  usrdef_hgr  ***
   !!
   !!                      ===  DOME configuration  ===
   !!
   !! User defined :   mesh and Coriolis parameter of a user configuration
   !!======================================================================
   !! History :  NEMO  ! 2017-11  (J. Chanut)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_hgr    : initialize the horizontal mesh for DOME configuration
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE par_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE usrdef_nam, ONLY: rn_dx, rn_dy, rn_f0   ! horizontal resolution in meters
   !                                             and Coriolis freq. 
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   REAL(wp) :: roffsetx, roffsety ! Offset in km to first f-point

   PUBLIC   usr_def_hgr   ! called by domhgr.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_hgr.F90 13295 2020-07-10 18:24:21Z acc $ 
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
      !!                DOME configuration : beta-plance with uniform grid spacing (rn_dx)
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
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) ::   zti, ztj   ! local scalars
      !!-------------------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_hgr : DOME configuration bassin'
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '          f-plane with regular grid-spacing'
      IF(lwp) WRITE(numout,*) '          given by rn_dx and rn_dy' 
      !
      !                          
      ! Position coordinates (in kilometers)
      !                          ==========
      ! Offsets in km of the first south west f-point: 
      roffsetx = -1700._wp
      roffsety =  -800._wp 
#if defined key_agrif
      IF( .NOT.Agrif_Root() ) THEN
         ! deduce offset from parent:
         roffsetx = Agrif_Parent(roffsetx) &
              & + (-(nbghostcells_x_w - 1) + (Agrif_Parent(nbghostcells_x_w) &
              & + Agrif_Ix()-2)*Agrif_Rhox()) * 1.e-3 * rn_dx
         roffsety = Agrif_Parent(roffsety) &
              & + (-(nbghostcells_y_s - 1) + (Agrif_Parent(nbghostcells_y_s) &
              & + Agrif_Iy()-2)*Agrif_Rhoy()) * 1.e-3 * rn_dy
      ENDIF
#endif
         
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         zti = REAL( mig0(ji) - 1, wp )   ! start at i=0 in the global grid without halos
         ztj = REAL( mjg0(jj) - 1, wp )   ! start at j=0 in the global grid without halos
         
         plamt(ji,jj) = roffsetx + rn_dx * 1.e-3 * ( zti - 0.5_wp )
         plamu(ji,jj) = roffsetx + rn_dx * 1.e-3 *   zti 
         plamv(ji,jj) = plamt(ji,jj) 
         plamf(ji,jj) = plamu(ji,jj) 
         
         pphit(ji,jj) = roffsety + rn_dy * 1.e-3 * ( ztj - 0.5_wp )
         pphiv(ji,jj) = roffsety + rn_dy * 1.e-3 *   ztj
         pphiu(ji,jj) = pphit(ji,jj) 
         pphif(ji,jj) = pphiv(ji,jj) 
      END_2D
      !     
      ! Horizontal scale factors (in meters)
      !                              ======
      pe1t(:,:) = rn_dx  ;   pe2t(:,:) = rn_dy 
      pe1u(:,:) = rn_dx  ;   pe2u(:,:) = rn_dy 
      pe1v(:,:) = rn_dx  ;   pe2v(:,:) = rn_dy 
      pe1f(:,:) = rn_dx  ;   pe2f(:,:) = rn_dy 

      !                             ! NO reduction of grid size in some straits 
      ke1e2u_v = 0                  !    ==>> u_ & v_surfaces will be computed in dom_hgr routine
      pe1e2u(:,:) = 0._wp           !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v(:,:) = 0._wp           !             require an initialization of INTENT(out) arguments
      !
      !
      !                       !==  Coriolis parameter  ==!
      kff = 1                       !  indicate not to compute Coriolis parameter afterward
      !
      pff_f(:,:) = rn_f0
      pff_t(:,:) = rn_f0
      !
   END SUBROUTINE usr_def_hgr

   !!======================================================================
END MODULE usrdef_hgr
