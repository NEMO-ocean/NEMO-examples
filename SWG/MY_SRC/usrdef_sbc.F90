MODULE usrdef_sbc
   !!======================================================================
   !!                     ***  MODULE  usrdef_sbc  ***
   !!
   !!                     ===  SWG configuration  ===
   !!
   !! User defined :   surface forcing of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-03  (S. Flavoni, G. Madec)  user defined interface
   !!             -   ! 2020-03  (A. Nasser) Shallow Water Eq. configuration
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usrdef_sbc    : user defined surface bounday conditions in GYRE case
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE phycst         ! physical constants
   USE usrdef_nam
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran    ! Fortran library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usrdef_sbc_oce       ! routine called in sbcmod module
   PUBLIC   usrdef_sbc_ice_tau   ! routine called by icestp.F90 for ice dynamics
   PUBLIC   usrdef_sbc_ice_flx   ! routine called by icestp.F90 for ice thermo
      
   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_sbc.F90 10425 2018-12-19 21:54:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usrdef_sbc_oce( kt, Kbb )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE usrdef_sbc_oce  ***
      !!              
      !! ** Purpose :   provide at each time-step the GYRE surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   analytical seasonal cycle for GYRE configuration.
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - set the ocean surface boundary condition, i.e.   
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, sfx
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      INTEGER, INTENT(in) ::   Kbb  ! ocean time index
      !!
      INTEGER  ::   ji, jj                ! dummy loop indices
      REAL(wp) ::   ztauu, ztauv          ! wind intensity projeted
      REAL(wp) ::   zrhoa  = 1.22         ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3       ! drag coefficient
      REAL(wp) ::   ztx, zty, zmod, zcoef ! temporary variables
      !!---------------------------------------------------------------------

      ! ---------------------------- !
      !  heat and freshwater fluxes  !   (no fluxes)
      ! ---------------------------- !
      
      emp (:,:) = 0._wp
      sfx (:,:) = 0._wp
      qns (:,:) = 0._wp
      qsr (:,:) = 0._wp

      ! ---------------------------- !
      !       momentum fluxes        !
      ! ---------------------------- !
      ! rotated case (45deg)
      ! ztauu = 0.2_wp / SQRT( 2._wp )    ! N.m-2
      ! ztauv = 0.2_wp / SQRT( 2._wp)     ! N.m-2
      ! non rotated
      !ztauu = 0.2_wp                    ! N.m-2
      !ztauv = 0._wp                     ! N.m-2
      ! general case
      ztauu =   REAL( rn_tau, wp ) * COS( rn_theta * rad )   ! N.m-2
      ztauv = - REAL( rn_tau, wp ) * SIN( rn_theta * rad )   ! N.m-2
      
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         ! length of the domain : 2000km x 2000km 
         utau(ji,jj) = - ztauu * COS( rpi * gphiu(ji,jj) / 2000000._wp)
         vtau(ji,jj) = - ztauv * COS( rpi * gphiv(ji,jj) / 2000000._wp)
      END_2D
      
      ! module of wind stress and wind speed at T-point
      zcoef = 1. / ( zrhoa * zcdrag ) 
      DO_2D( 0, 0, 0, 0 )
         ztx = utau(ji-1,jj  ) + utau(ji,jj) 
         zty = vtau(ji  ,jj-1) + vtau(ji,jj) 
         zmod = 0.5 * SQRT( ztx * ztx + zty * zty )
         taum(ji,jj) = zmod
         wndm(ji,jj) = SQRT( zmod * zcoef )
      END_2D
      CALL lbc_lnk( 'usrdef_sbc', taum(:,:), 'T', 1. , wndm(:,:), 'T', 1. )
      !
   END SUBROUTINE usrdef_sbc_oce


   SUBROUTINE usrdef_sbc_ice_tau( kt )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
   END SUBROUTINE usrdef_sbc_ice_tau


   SUBROUTINE usrdef_sbc_ice_flx( kt, phs, phi )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phs    ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phi    ! ice thickness
   END SUBROUTINE usrdef_sbc_ice_flx

   !!======================================================================
END MODULE usrdef_sbc
