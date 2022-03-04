MODULE bdyvol
   !!======================================================================
   !!                       ***  MODULE  bdyvol  ***
   !! Ocean dynamic :  Volume constraint when unstructured boundary 
   !!                  and filtered free surface are used
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2006-01  (J. Chanut) Bug correction
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!            4.0  !  2019-01  (P. Mathiot) adapted to time splitting     
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers 
   USE bdy_oce        ! ocean open boundary conditions
   USE sbc_oce        ! ocean surface boundary conditions
   USE isf_oce, ONLY : fwfisf_cav, fwfisf_par  ! ice shelf
   USE dom_oce        ! ocean space and time domain 
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! for mppsum
   USE lib_fortran    ! Fortran routines library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_vol2d    ! called by dynspg_ts

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: bdyvol.F90 12489 2020-02-28 15:55:11Z davestorkey $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_vol2d( kt, kc, pua2d, pva2d, phu, phv )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE bdyvol  ***
      !!
      !! ** Purpose :   This routine controls the volume of the system. 
      !!      A correction velocity is calculated to correct the total transport 
      !!      through the unstructured OBC. 
      !!
      !! ** Method  :   The correction velocity (zubtpecor here) is defined calculating
      !!      the total transport through all open boundaries (trans_bdy) minus
      !!      the cumulate E-P flux (z_cflxemp) divided by the total lateral 
      !!      surface (bdysurftot) of the unstructured boundary. 
      !!         zubtpecor = [trans_bdy - z_cflxemp ]*(1./bdysurftot)
      !!      with z_cflxemp => sum of (Evaporation minus Precipitation)
      !!                       over all the domain in m3/s at each time step.
      !!      z_cflxemp < 0 when precipitation dominate
      !!      z_cflxemp > 0 when evaporation dominate
      !!
      !!      There are 2 options (user's desiderata): 
      !!         1/ The volume changes according to E-P, this is the default
      !!            option. In this case the cumulate E-P flux are setting to
      !!            zero (z_cflxemp=0) to calculate the correction velocity. So
      !!            it will only balance the flux through open boundaries.
      !!            (set nn_volctl to 0 in tne namelist for this option)
      !!         2/ The volume is constant even with E-P flux. In this case
      !!            the correction velocity must balance both the flux 
      !!            through open boundaries and the ones through the free
      !!            surface. 
      !!            (set nn_volctl to 1 in tne namelist for this option)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, kc   ! ocean time-step index, cycle time-step
      !
      INTEGER  ::   ji, jj, jk, jb, jgrd
      INTEGER  ::   ib_bdy, ii, ij
      REAL(wp) ::   zubtpecor, ztranst
      REAL(wp), SAVE :: z_cflxemp                                  ! cumulated emp flux
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pua2d, pva2d  ! Barotropic velocities
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: phu, phv         ! Ocean depth at U- and V-points
      TYPE(OBC_INDEX), POINTER :: idx
      !!-----------------------------------------------------------------------------
      !
      ! Calculate the cumulate surface Flux z_cflxemp (m3/s) over all the domain
      ! -----------------------------------------------------------------------
      IF ( kc == 1 ) z_cflxemp = glob_sum( 'bdyvol', ( emp(:,:) - rnf(:,:) + fwfisf_cav(:,:) + fwfisf_par(:,:) ) * bdytmask(:,:) * e1e2t(:,:)  ) / rho0

      ! Compute bdy surface each cycle if non linear free surface
      ! ---------------------------------------------------------
      IF ( .NOT. ln_linssh ) THEN
         ! compute area each time step
         bdysurftot = bdy_segs_surf( phu, phv )
      ELSE
         ! compute area only the first time
         IF ( ( kt == nit000 ) .AND. ( kc == 1 ) ) bdysurftot = bdy_segs_surf( phu, phv )
      END IF

      ! Transport through the unstructured open boundary
      ! ------------------------------------------------
      zubtpecor = 0._wp
      DO ib_bdy = 1, nb_bdy
         idx => idx_bdy(ib_bdy)
         !
         jgrd = 2                               ! cumulate u component contribution first 
         DO jb = 1, idx%nblenrim(jgrd)
            ii = idx%nbi(jb,jgrd)
            ij = idx%nbj(jb,jgrd)
            IF( ii == 1 .OR. ii == jpi .OR. ij == 1 .OR. ij == jpj )  CYCLE   ! sum : else halo couted twice
            zubtpecor = zubtpecor + idx%flagu(jb,jgrd) * pua2d(ii,ij) * e2u(ii,ij) * phu(ii,ij) * tmask_i(ii,ij) * tmask_i(ii+1,ij)
         END DO
         jgrd = 3                               ! then add v component contribution
         DO jb = 1, idx%nblenrim(jgrd)
            ii = idx%nbi(jb,jgrd)
            ij = idx%nbj(jb,jgrd)
            IF( ii == 1 .OR. ii == jpi .OR. ij == 1 .OR. ij == jpj )  CYCLE   ! sum : else halo couted twice
            zubtpecor = zubtpecor + idx%flagv(jb,jgrd) * pva2d(ii,ij) * e1v(ii,ij) * phv(ii,ij) * tmask_i(ii,ij) * tmask_i(ii,ij+1)
         END DO
         !
      END DO
      IF( lk_mpp )   CALL mpp_sum( 'bdyvol', zubtpecor )   ! sum over the global domain

      ! The normal velocity correction
      ! ------------------------------
      IF( nn_volctl==1 ) THEN   ;   zubtpecor = ( zubtpecor - z_cflxemp ) / bdysurftot  ! maybe should be apply only once at the end
      ELSE                      ;   zubtpecor =   zubtpecor               / bdysurftot
      END IF

      ! Correction of the total velocity on the unstructured boundary to respect the mass flux conservation
      ! -------------------------------------------------------------
!      DO ib_bdy = 1, nb_bdy
!jc      
      DO ib_bdy = 2, nb_bdy
         idx => idx_bdy(ib_bdy)
         !
         jgrd = 2                               ! correct u component
         DO jb = 1, idx%nblen(jgrd)
               ii = idx%nbi(jb,jgrd)
               ij = idx%nbj(jb,jgrd)
               IF( ii == 1 .OR. ii == jpi .OR. ij == 1 .OR. ij == jpj )  CYCLE   ! to remove ?
!               pua2d(ii,ij) = pua2d(ii,ij) - idx%flagu(jb,jgrd) * zubtpecor * tmask_i(ii,ij) * tmask_i(ii+1,ij)
               pua2d(ii,ij) = pua2d(ii,ij) - zubtpecor * &
                             & tmask_i(ii+1,ij) * tmask_i(ii,ij) 
         END DO
         jgrd = 3                              ! correct v component
         DO jb = 1, idx%nblenrim(jgrd)
               ii = idx%nbi(jb,jgrd)
               ij = idx%nbj(jb,jgrd)
               IF( ii == 1 .OR. ii == jpi .OR. ij == 1 .OR. ij == jpj )  CYCLE   ! to remove ?
               pva2d(ii,ij) = pva2d(ii,ij) - idx%flagv(jb,jgrd) * zubtpecor * tmask_i(ii,ij) * tmask_i(ii,ij+1)
         END DO
         !
      END DO
      ! 
      ! Check the cumulated transport through unstructured OBC once barotropic velocities corrected
      ! ------------------------------------------------------
      IF( MOD( kt, MAX(nn_write,1) ) == 0 .AND. ( kc == 1 ) ) THEN
         !
         ! compute residual transport across boundary
         ztranst = 0._wp
         DO ib_bdy = 1, nb_bdy
            idx => idx_bdy(ib_bdy)
            !
            jgrd = 2                               ! correct u component
            DO jb = 1, idx%nblenrim(jgrd)
                  ii = idx%nbi(jb,jgrd)
                  ij = idx%nbj(jb,jgrd)
                  IF( ii == 1 .OR. ii == jpi .OR. ij == 1 .OR. ij == jpj )  CYCLE
                  ztranst = ztranst + idx%flagu(jb,jgrd) * pua2d(ii,ij) * e2u(ii,ij) * phu(ii,ij) * tmask_i(ii,ij) * tmask_i(ii+1,ij)
            END DO
            jgrd = 3                              ! correct v component
            DO jb = 1, idx%nblenrim(jgrd)
                  ii = idx%nbi(jb,jgrd)
                  ij = idx%nbj(jb,jgrd)
                  IF( ii == 1 .OR. ii == jpi .OR. ij == 1 .OR. ij == jpj )  CYCLE
                  ztranst = ztranst + idx%flagv(jb,jgrd) * pva2d(ii,ij) * e1v(ii,ij) * phv(ii,ij) * tmask_i(ii,ij) * tmask_i(ii,ij+1)
            END DO
            !
         END DO
         IF( lk_mpp )   CALL mpp_sum('bdyvol', ztranst )   ! sum over the global domain


         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'bdy_vol : time step :', kt
         IF(lwp) WRITE(numout,*)'~~~~~~~ '
         IF(lwp) WRITE(numout,*)'          cumulate flux EMP             =', z_cflxemp  , ' (m3/s)'
         IF(lwp) WRITE(numout,*)'          total lateral surface of OBC  =', bdysurftot, '(m2)'
         IF(lwp) WRITE(numout,*)'          correction velocity zubtpecor =', zubtpecor , '(m/s)'
         IF(lwp) WRITE(numout,*)'          cumulated transport ztranst   =', ztranst   , '(m3/s)'
      END IF 
      !
   END SUBROUTINE bdy_vol2d
   !
   REAL(wp) FUNCTION bdy_segs_surf(phu, phv)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE bdy_ctl_seg  ***
      !!
      !! ** Purpose :   Compute total lateral surface for volume correction
      !!
      !!----------------------------------------------------------------------

      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: phu, phv ! water column thickness at U and V points
      INTEGER            ::  igrd, ib_bdy, ib                ! loop indexes
      INTEGER , POINTER  ::  nbi, nbj                        ! short cuts
      REAL(wp), POINTER  ::  zflagu, zflagv                  !    -   -

      ! Compute total lateral surface for volume correction:
      ! ----------------------------------------------------
      bdy_segs_surf = 0._wp
      igrd = 2      ! Lateral surface at U-points
!      DO ib_bdy = 1, nb_bdy
!jc      
      DO ib_bdy = 2, nb_bdy

         DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)
            nbi => idx_bdy(ib_bdy)%nbi(ib,igrd)
            nbj => idx_bdy(ib_bdy)%nbj(ib,igrd)
            IF( nbi == 1 .OR. nbi == jpi .OR. nbj == 1 .OR. nbj == jpj )  CYCLE
            zflagu => idx_bdy(ib_bdy)%flagu(ib,igrd)
            bdy_segs_surf = bdy_segs_surf + phu(nbi, nbj)                              &
               &                            * e2u(nbi, nbj) * ABS( zflagu )            &
               &                            * tmask_i(nbi, nbj) * tmask_i(nbi+1, nbj)
         END DO
      END DO

      igrd=3 ! Add lateral surface at V-points
!      DO ib_bdy = 1, nb_bdy
!jc      
      DO ib_bdy = 2, nb_bdy

         DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)
            nbi => idx_bdy(ib_bdy)%nbi(ib,igrd)
            nbj => idx_bdy(ib_bdy)%nbj(ib,igrd)
            IF( nbi == 1 .OR. nbi == jpi .OR. nbj == 1 .OR. nbj == jpj )  CYCLE
            zflagv => idx_bdy(ib_bdy)%flagv(ib,igrd)
            bdy_segs_surf = bdy_segs_surf + phv(nbi, nbj)                              &
               &                            * e1v(nbi, nbj) * ABS( zflagv )            &
               &                            * tmask_i(nbi, nbj) * tmask_i(nbi, nbj+1)
         END DO
      END DO
      !
      ! redirect the time to bdyvol as this variable is only used by bdyvol
      IF( lk_mpp )   CALL mpp_sum( 'bdyvol', bdy_segs_surf ) ! sum over the global domain
      !
   END FUNCTION bdy_segs_surf
   !!======================================================================
END MODULE bdyvol
