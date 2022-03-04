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
   USE dom_oce  ,  ONLY: nimpp, njmpp       ! ocean space and time domain
   USE par_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE usrdef_nam, ONLY: rn_dx, rn_dy
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
   !! $Id: usrdef_hgr.F90 12740 2020-04-12 09:03:06Z smasson $ 
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
      REAL(wp) ::   zfact, zlam0, zphi0, zti, zui, zvi, zfi, ztj, zuj, zvj, zfj      ! local scalars
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
      zlam0 = 0.0
      zphi0 = 0.0
      DO_2D (1,1,1,1)
         zti = FLOAT( ji - 1 + nimpp - 1 )          ;  ztj = FLOAT( jj - 1 + njmpp - 1 )
         zui = FLOAT( ji - 1 + nimpp - 1 ) + 0.5_wp ;  zvj = FLOAT( jj - 1 + njmpp - 1 ) + 0.5_wp
         
         plamt(ji,jj) = zlam0 + rn_dx * zti
         plamu(ji,jj) = zlam0 + rn_dx * zui
         plamv(ji,jj) = plamt(ji,jj) 
         plamf(ji,jj) = plamu(ji,jj) 
         
         pphit(ji,jj) = zphi0 + rn_dy * ztj
         pphiv(ji,jj) = zphi0 + rn_dy * zvj
         pphiu(ji,jj) = pphit(ji,jj) 
         pphif(ji,jj) = pphiv(ji,jj) 
      END_2D
      !
      !                       !==  Horizontal scale factors  ==!   (in meters)
      DO_2D (1,1,1,1)
         !                       ! e1   (zonal)
         pe1t(ji,jj) = rn_dx
         pe1u(ji,jj) = rn_dx
         pe1v(ji,jj) = rn_dx
         pe1f(ji,jj) = rn_dx
         !                       ! e2   (meridional)
         pe2t(ji,jj) = rn_dy
         pe2u(ji,jj) = rn_dy
         pe2v(ji,jj) = rn_dy
         pe2f(ji,jj) = rn_dy
      END_2D
      !                             ! NO reduction of grid size in some straits 
      ke1e2u_v    = 0               !    ==>> u_ & v_surfaces will be computed in dom_ghr routine
      pe1e2u(:,:) = 0._wp           !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v(:,:) = 0._wp           !             require an initialization of INTENT(out) arguments
      !
      !
      !                       !==  Coriolis parameter  ==!
      kff = 1                       ! Coriolis parameter set to 0
      pff_f(:,:) = 0._wp            ! CAUTION: set to zero to avoid error with some compilers that
      pff_t(:,:) = 0._wp            !             require an initialization of INTENT(out) arguments
      !
   END SUBROUTINE usr_def_hgr

   !!======================================================================
END MODULE usrdef_hgr
