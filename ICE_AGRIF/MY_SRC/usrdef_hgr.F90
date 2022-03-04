MODULE usrdef_hgr
   !!======================================================================
   !!                       ***  MODULE  usrdef_hgr  ***
   !!
   !!                      ===  ICE_AGRIF configuration  ===
   !!
   !! User defined :   mesh and Coriolis parameter of a user configuration
   !!======================================================================
   !! History :  NEMO  ! 2016-08  (S. Flavoni, G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_hgr    : initialize the horizontal mesh for ICE_AGRIF configuration
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE par_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE usrdef_nam, ONLY: rn_dx, rn_dy, ln_corio, rn_ppgphi0   ! horizontal resolution in meters
   !                                                            coriolis and reference latitude
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
   !! $Id: usrdef_hgr.F90 15119 2021-07-13 14:43:22Z jchanut $ 
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
      !!                ICE_AGRIF configuration : uniform grid spacing (rn_dx)
      !!                without Coriolis force (f=0)
      !!
      !! ** Action  : - define longitude & latitude of t-, u-, v- and f-points (in kms) 
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
      REAL(wp) ::   zbeta, zf0
      REAL(wp) ::   zti, ztj   ! local scalars
      !!-------------------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_hgr : ICE_AGRIF configuration bassin'
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '          f-plane with irregular grid-spacing (+- 10%)'
      IF(lwp) WRITE(numout,*) '          the max is given by rn_dx and rn_dy' 
      !
      !
      ! Position coordinates (in kilometers)
      !                          ==========
      ! Offset is given at first f-point, i.e. at (i,j) = (nn_hls+1, nn_hls+1)
      ! Here we assume the grid is centred around a T-point at the middle of
      ! of the domain (hence domain size is odd)
      roffsetx = (-REAL(Ni0glo-1, wp) + 1._wp) * 0.5 * 1.e-3 * rn_dx
      roffsety = (-REAL(Nj0glo-1, wp) + 1._wp) * 0.5 * 1.e-3 * rn_dy
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
         zti = REAL( mig0(ji)-1, wp )  ! start at i=0 in the global grid without halos
         ztj = REAL( mjg0(jj)-1, wp )  ! start at j=0 in the global grid without halos

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
!! ==> EITHER 1) variable scale factors
!! clem: This can be used with a 1proc simulation but I think it breaks repro when >1procs are used      
!!         DO jj = 1, jpj
!!            DO ji = 1, jpi
!!               !!pe1t(ji,jj) = rn_dx * EXP( -0.8/REAL(jpiglo**2) * (mi0(ji)-REAL(jpiglo+1)*0.5)**2 )  ! gaussian shape
!!               !!pe2t(ji,jj) = rn_dy * EXP( -0.8/REAL(jpjglo**2) * (mj0(jj)-REAL(jpjglo+1)*0.5)**2 )  ! gaussian shape
!!               pe1t(ji,jj) = rn_dx * ( 1. -0.1 * ABS(REAL(mi0(ji))-REAL(jpiglo+1)*0.5) / (1.-REAL(jpiglo+1)*0.5) ) ! linear shape
!!               pe2t(ji,jj) = rn_dy * ( 1. -0.1 * ABS(REAL(mj0(jj))-REAL(jpjglo+1)*0.5) / (1.-REAL(jpjglo+1)*0.5) ) ! linear shape
!!            END DO
!!         END DO
!!#if defined key_agrif 
!!         IF( .NOT. Agrif_Root() ) THEN ! only works if the zoom is positioned at the center of the parent grid
!!            DO jj = 1, jpj
!!               DO ji = 1, jpi
!!                  pe1t(ji,jj) = rn_dx * ( 1. -0.1 * ABS(REAL(mi0(ji))-REAL(jpiglo+1)*0.5) / (1.-REAL(jpiglo+1)*0.5)  &
!!                     &                            * REAL(jpiglo) / REAL(Agrif_Parent(jpiglo) * Agrif_Rhox()) )       ! factor to match parent grid
!!                  pe2t(ji,jj) = rn_dy * ( 1. -0.1 * ABS(REAL(mj0(jj))-REAL(jpjglo+1)*0.5) / (1.-REAL(jpjglo+1)*0.5)  &
!!                     &                            * REAL(jpjglo) / REAL(Agrif_Parent(jpjglo) * Agrif_Rhoy()) )       ! factor to match parent grid
!!               END DO
!!            END DO
!!         ENDIF
!!#endif
!! ==> OR 2) constant scale factors
         pe1t(:,:) = rn_dx
         pe2t(:,:) = rn_dy
!! ==> END
         
      pe1u(:,:) = pe1t(:,:)      ;      pe2u(:,:) = pe2t(:,:)
      pe1v(:,:) = pe1t(:,:)      ;      pe2v(:,:) = pe2t(:,:)
      pe1f(:,:) = pe1t(:,:)      ;      pe2f(:,:) = pe2t(:,:)

      !                             ! NO reduction of grid size in some straits 
      ke1e2u_v = 0                  !    ==>> u_ & v_surfaces will be computed in dom_ghr routine
      pe1e2u(:,:) = 0._wp           !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v(:,:) = 0._wp           !             require an initialization of INTENT(out) arguments
      !
      !
      !                       !==  Coriolis parameter  ==!
      kff = 1                       !  indicate not to compute Coriolis parameter afterward
      !
      IF( ln_corio ) THEN
         zbeta = 2._wp * omega * COS( rad * rn_ppgphi0 ) / ra
         zf0   = 2._wp * omega * SIN( rad * rn_ppgphi0 )
         pff_f(:,:) = zf0 + zbeta * pphif(:,:) * 1.e+3
         pff_t(:,:) = zf0 + zbeta * pphit(:,:) * 1.e+3
      ELSE
         pff_f(:,:) = 0.
         pff_t(:,:) = 0.
      ENDIF
      !
   END SUBROUTINE usr_def_hgr

   !!======================================================================
END MODULE usrdef_hgr
