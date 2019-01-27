MODULE usrdef_hgr
   !!======================================================================
   !!                       ***  MODULE  usrdef_hgr  ***
   !!
   !!                      ===  BENCH configuration  ===
   !!
   !! User defined :   mesh and Coriolis parameter of a user configuration
   !!======================================================================
   !! History :  NEMO  ! 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_hgr    : initialize the horizontal mesh for BENCH configuration
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE par_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_hgr   ! called by domhgr.F90

   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0, NEMO Consortium (2016)
   !! $Id$ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
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
      !! ** Purpose :   square box mesh mesh and Coriolis parameter
      !!
      !! ** Method  :   set all intent(out) argument to a proper value
      !!                BENCH configuration : beta-plance with uniform grid spacing (zres)
      !!
      !! ** Action  : - define longitude & latitude of t-, u-, v- and f-points (in grid points) 
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
      REAL(wp) ::   zres, zf0
      REAL(wp) ::   zti, zui, ztj, zvj   ! local scalars
      !!-------------------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_hgr : BENCH configuration bassin'
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '          Beta-plane with regular grid-spacing'
      IF(lwp) WRITE(numout,*) '          given by rn_dx and rn_dy' 
      !
      !                          
      ! Position coordinates (in grid points)
      !                          ==========         
      DO jj = 1, jpj
         DO ji = 1, jpi
            
            zti = REAL( ji - 1 + nimpp - 1, wp )          ;  ztj = REAL( jj - 1 + njmpp - 1, wp )
            zui = REAL( ji - 1 + nimpp - 1, wp ) + 0.5_wp ;  zvj = REAL( jj - 1 + njmpp - 1, wp ) + 0.5_wp

            plamt(ji,jj) = zti
            plamu(ji,jj) = zui
            plamv(ji,jj) = zti
            plamf(ji,jj) = zui
   
            pphit(ji,jj) = ztj
            pphiv(ji,jj) = zvj
            pphiu(ji,jj) = ztj
            pphif(ji,jj) = zvj
            
         END DO
      END DO
      !     
      ! Horizontal scale factors (in meters)
      !                              ======
      zres = 1.e+5   !  100km
      pe1t(:,:) = zres  ;   pe2t(:,:) = zres 
      pe1u(:,:) = zres  ;   pe2u(:,:) = zres
      pe1v(:,:) = zres  ;   pe2v(:,:) = zres
      pe1f(:,:) = zres  ;   pe2f(:,:) = zres 

      !                             ! NO reduction of grid size in some straits 
      ke1e2u_v = 0                  !    ==>> u_ & v_surfaces will be computed in dom_hgr routine
      pe1e2u(:,:) = 0._wp           !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v(:,:) = 0._wp           !             require an initialization of INTENT(out) arguments
      !
      !
      !                       !==  Coriolis parameter  ==!
      kff = 1                       !  indicate not to compute Coriolis parameter afterward
      !
      zf0   = 2._wp * omega * SIN( rad * 45 )   ! constant coriolis factor corresponding to 45°N
      pff_f(:,:) = zf0
      pff_t(:,:) = zf0
      !
   END SUBROUTINE usr_def_hgr

   !!======================================================================
END MODULE usrdef_hgr
