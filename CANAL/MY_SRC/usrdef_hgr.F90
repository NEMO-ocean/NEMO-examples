MODULE usrdef_hgr
   !!======================================================================
   !!                       ***  MODULE  usrdef_hgr  ***
   !!
   !!                      ===  CANAL configuration  ===
   !!
   !! User defined :   mesh and Coriolis parameter of a user configuration
   !!======================================================================
   !! History :  NEMO  ! 2017-11  (J. Chanut)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_hgr    : initialize the horizontal mesh for CANAL configuration
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE par_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE usrdef_nam, ONLY: rn_dx, rn_dy, rn_0xratio, rn_0yratio, rn_ppgphi0, nn_fcase
   !                                                  and reference latitude
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_hgr   ! called by domhgr.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_hgr.F90 14433 2021-02-11 08:06:49Z smasson $ 
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
      !!                CANAL configuration : beta-plance with uniform grid spacing (rn_dx)
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
      INTEGER  ::   ii0, ij0   ! dummy loop indices
      REAL(wp) ::   zbeta, zf0
      REAL(wp) ::   zti, ztj   ! local scalars
      !!-------------------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_hgr : CANAL configuration bassin'
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '          Beta-plane with regular grid-spacing'
      IF(lwp) WRITE(numout,*) '          given by rn_dx and rn_dy' 
      !
      !                          
      ! Position coordinates (in kilometers)
      !                          ==========
      ii0 = NINT( REAL(Ni0glo, wp) * rn_0xratio )
      ij0 = NINT( REAL(Nj0glo, wp) * rn_0yratio )

#if defined key_agrif
      ! ! let lower left longitude and latitude from parent
      IF (.NOT.Agrif_root()) THEN
          to be coded...
      ENDIF 
#endif
         
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )         
         zti = REAL( mig0(ji)-ii0, wp )   ! =0 at i=ii0 in the global grid without halos
         ztj = REAL( mjg0(jj)-ij0, wp )   ! =0 at i=ij0 in the global grid without halos
         
         plamt(ji,jj) = rn_dx *   zti
         plamu(ji,jj) = rn_dx * ( zti + 0.5_wp ) 
         plamv(ji,jj) = plamt(ji,jj) 
         plamf(ji,jj) = plamu(ji,jj) 
         
         pphit(ji,jj) = rn_dy *   ztj
         pphiv(ji,jj) = rn_dy * ( ztj + 0.5_wp ) 
         pphiu(ji,jj) = pphit(ji,jj) 
         pphif(ji,jj) = pphiv(ji,jj) 
      END_2D
      !     
      ! Horizontal scale factors (in meters)
      !                              ======
      pe1t(:,:) = rn_dx * 1.e+3  ;   pe2t(:,:) = rn_dy * 1.e+3 
      pe1u(:,:) = rn_dx * 1.e+3  ;   pe2u(:,:) = rn_dy * 1.e+3
      pe1v(:,:) = rn_dx * 1.e+3  ;   pe2v(:,:) = rn_dy * 1.e+3
      pe1f(:,:) = rn_dx * 1.e+3  ;   pe2f(:,:) = rn_dy * 1.e+3 

      !                             ! NO reduction of grid size in some straits 
      ke1e2u_v = 0                  !    ==>> u_ & v_surfaces will be computed in dom_hgr routine
      pe1e2u(:,:) = 0._wp           !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v(:,:) = 0._wp           !             require an initialization of INTENT(out) arguments
      !
      !
      !                       !==  Coriolis parameter  ==!
      kff = 1                       !  indicate not to compute Coriolis parameter afterward
      !
      SELECT CASE(nn_fcase)
      CASE(0)
         zf0   = 2._wp * omega * SIN( rad * rn_ppgphi0 )
         pff_f(:,:) = zf0
         pff_t(:,:) = zf0
      CASE(1)
         zbeta = 2._wp * omega * COS( rad * rn_ppgphi0 ) / ra
         zf0   = 2._wp * omega * SIN( rad * rn_ppgphi0 )
         pff_f(:,:) = zf0 + zbeta * pphif(:,:) * 1.e+3
         pff_t(:,:) = zf0 + zbeta * pphit(:,:) * 1.e+3
      CASE(2)
         pff_f(:,:) = 2._wp * omega * SIN( rad * ( rn_ppgphi0 + pphif(:,:)/110. ) )
         pff_t(:,:) = 2._wp * omega * SIN( rad * ( rn_ppgphi0 + pphit(:,:)/110. ) )
      END SELECT
      !
   END SUBROUTINE usr_def_hgr

   !!======================================================================
END MODULE usrdef_hgr
