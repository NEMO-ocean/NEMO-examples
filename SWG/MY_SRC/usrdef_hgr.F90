MODULE usrdef_hgr
   !!======================================================================
   !!                     ***  MODULE usrdef_hgr   ***
   !!
   !!                     ===  SWG configuration  ===
   !!
   !! User defined :   mesh and Coriolis parameter of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-03  (S. Flavoni)
   !!             -   ! 2020-03  (A. Nasser) Shallow Water Eq. configuration
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_hgr   : initialize the horizontal mesh 
   !!----------------------------------------------------------------------
   USE dom_oce  , ONLY: nimpp, njmpp       ! ocean space and time domain
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
   !! $Id: usrdef_hgr.F90 10069 2018-08-28 14:12:24Z nicolasmartin $ 
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
      !!                Here SWG configuration :
      !!          Rectangular mid-latitude domain 
      !!          - with axes rotated by 45 degrees
      !!          - a constant horizontal resolution of 106 km 
      !!          - on a beta-plane
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
      INTEGER  ::   ji, jj               ! dummy loop indices
      REAL(wp) ::   zlam1, zlam0, zcos_theta, zim1 , zjm1 , ze1  , ze1deg ! local scalars
      REAL(wp) ::   zphi1, zphi0, zsin_theta, zim05, zjm05, znorme        !   -      -
      REAL(wp) ::   zgl, zbl      !   -      -

      !!-------------------------------------------------------------------------------
      !
      !     !==  beta-plane with regular grid-spacing and rotated domain ==!  (SWG configuration)
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_hgr : SWG configuration (beta-plane with rotated regular grid-spacing)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      !
      !                       !==  grid point position  ==!
      !
      ze1 =  rn_dx / REAL(nn_SWG, wp)                    ! [m] gridspacing used
      zgl =  rn_domsiz + 2._wp * REAL(nn_gc, wp) * ze1   ! [m] length of the square with ghostcells
      ! fit the best square around the square + ghost cells 
      zbl = zgl * ( COS( rn_theta * rad ) + SIN( rn_theta * rad ) )   ! length side bigger domain [m]                                 

      ! unrotated (0deg) 
      !zcos_theta = 1._wp
      !zsin_theta = 0._wp
      ! rotated case (45deg)
      !zcos_theta =   1._wp / SQRT( 2._wp ) 
      !zsin_theta =   1._wp / SQRT( 2._wp )  
      ! rotation angle
      zcos_theta = COS( rn_theta * rad)
      zsin_theta = SIN( rn_theta * rad)   

      ! exact origin in meters
      zlam1 =  zbl * COS((rn_theta + 45 )* rad ) / SQRT( 2._wp )  - rn_domsiz/2._wp 
      zphi1 =  zbl * SIN((rn_theta + 45 )* rad ) / SQRT( 2._wp )  - rn_domsiz/2._wp 
      ! origin put in the true corner of a cell so there will be no cropping 
      ! of the edge cells
      zlam0 = REAL( anint( zlam1 / ze1 ), wp ) * ze1
      zphi0 = REAl( anint( zphi1 / ze1 ), wp ) * ze1

      IF(lwp) WRITE(numout,*) '                  origin position    zlam0   = ', zlam0/1000,   ' km'
      IF(lwp) WRITE(numout,*) '                  origin position    zphi0   = ', zphi0/1000,   ' km'

      ! O1M = OM x rotation_theta + OO1
      ! zim1, zim05, zjm1, zjm05 fit for 2 ghost cells on each side
      DO jj = 1, jpj 
         DO ji = 1, jpi 
            zim1 = REAL( ji + nimpp - nn_hls )   ;   zim05 = REAL( ji + nimpp - nn_hls ) - 0.5
            zjm1 = REAL( jj + njmpp - nn_hls )   ;   zjm05 = REAL( jj + njmpp - nn_hls ) - 0.5
            !   
            !glamt(i,j) position (meters) at T-point 
            !gphit(i,j) position (meters) at T-point  
            plamt(ji,jj) =   zim05 * ze1 * zcos_theta - zjm05 * ze1 * zsin_theta - zlam0
            pphit(ji,jj) = + zim05 * ze1 * zsin_theta + zjm05 * ze1 * zcos_theta - zphi0
            !   
            !glamu(i,j) position (meters) at U-point
            !gphiu(i,j) position (meters) at U-point
            plamu(ji,jj) =   zim1  * ze1 * zcos_theta - zjm05 * ze1 * zsin_theta - zlam0
            pphiu(ji,jj) = + zim1  * ze1 * zsin_theta + zjm05 * ze1 * zcos_theta - zphi0
            !   
            !glamv(i,j) position (meters) at V-point
            !gphiv(i,j) position (meters) at V-point
            plamv(ji,jj) =   zim05 * ze1 * zcos_theta - zjm1  * ze1 * zsin_theta - zlam0
            pphiv(ji,jj) = + zim05 * ze1 * zsin_theta + zjm1  * ze1 * zcos_theta - zphi0
            !
            !glamf(i,j) position (meters) at F-point
            !gphif(i,j) position (meters) at F-point 
            plamf(ji,jj) =   zim1  * ze1 * zcos_theta - zjm1  * ze1 * zsin_theta - zlam0
            pphif(ji,jj) = + zim1  * ze1 * zsin_theta + zjm1  * ze1 * zcos_theta - zphi0
         END DO
      END DO
      !
      !                       !== Horizontal scale factors ==! (in meters)
      !                     
      !                                         ! constant grid spacing
      pe1t(:,:) =  ze1     ;      pe2t(:,:) = ze1
      pe1u(:,:) =  ze1     ;      pe2u(:,:) = ze1
      pe1v(:,:) =  ze1     ;      pe2v(:,:) = ze1
      pe1f(:,:) =  ze1     ;      pe2f(:,:) = ze1
      !
      !                                         ! NO reduction of grid size in some straits 
      ke1e2u_v = 0                              !    ==>> u_ & v_surfaces will be computed in dom_ghr routine
      pe1e2u(:,:) = 0._wp                       !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v(:,:) = 0._wp                       !             require an initialization of INTENT(out) arguments
      !
      !
      !                       !==  Coriolis parameter  ==!
      kff = 1                                            !  indicate not to compute ff afterward
      !
      pff_f(:,:) =  REAL( rn_f0, wp ) + REAL( rn_beta, wp ) * ABS( pphif(:,:) ) ! f = f0 +beta* y
      pff_t(:,:) =  REAL( rn_f0, wp ) + REAL( rn_beta, wp ) * ABS( pphit(:,:) ) ! f = f0 +beta* y
      ! 
      IF(lwp) WRITE(numout,*) '                           beta-plane used. f0   = ', rn_f0 ,  ' 1/s'
      IF(lwp) WRITE(numout,*) '                           beta-plane used. beta = ', rn_beta, ' 1/(s.m)'
      !
   END SUBROUTINE usr_def_hgr

   !!======================================================================
END MODULE usrdef_hgr
