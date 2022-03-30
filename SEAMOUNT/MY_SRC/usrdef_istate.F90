MODULE usrdef_istate
   !!======================================================================
   !!                   ***  MODULE  usrdef_istate   ***
   !!
   !!                     ===  GYRE configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-03  (S. Flavoni) Original code
   !!                 ! 2020-11  (S. Techene, G. Madec) separate tsuv from ssh
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE dom_oce , ONLY : mi0, mig, mjg, glamt, gphit, ht_0
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE eosbn2, ONLY: rn_a0, rho0 
   USE usrdef_nam
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate       ! called in istate.F90
   PUBLIC   usr_def_istate_ssh   ! called by domqco.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_istate.F90 14053 2020-12-03 13:48:38Z techene $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

      REAL(wp) :: T0, dtem, delta  ! used to set the initial potential temperature profile

CONTAINS

   !!----------------------------------------------------------------------

   REAL FUNCTION theta_initial( pdept )
      REAL(wp), INTENT(in) ::   pdept       
      theta_initial = T0 + dtem*EXP( -pdept/delta )
   END FUNCTION theta_initial

   !!----------------------------------------------------------------------

   REAL FUNCTION ocean_depth( p_lam, p_phi )
      REAL(wp), INTENT(in) ::   p_lam      !  "longitude": depends on ji        
      REAL(wp), INTENT(in) ::   p_phi      !  "latitude":  depends on jj        

      ocean_depth = rn_bathy - rn_seamountheight * EXP( &
               &  - ( (1000._wp*p_lam - rn_lam_mid)**2 + (1000._wp*p_phi - rn_phi_mid)**2 ) / rn_l**2 ) 

   END FUNCTION ocean_depth

   !!----------------------------------------------------------------------
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here GYRE configuration example : (double gyre with rotated domain)
      !!
      !! ** Method  : - set temprature field
      !!              - set salinity   field
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
      !
      INTEGER                                              :: ji, jj, jk  ! dummy loop indices
      REAL(wp)                                             :: T0
      REAL(wp)                                             :: drho
      REAL(wp), DIMENSION(jpi,jpj)                         :: rho
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : SEAMOUNT_TEST_CASE configuration:'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   Ocean at rest with analytical initial stratification'
      IF(lwp) WRITE(numout,*) ''
      !
      pu  (:,:,:) = 0._wp        ! ocean at rest
      pv  (:,:,:) = 0._wp
      !
      SELECT CASE (nn_ini_cond)
         CASE (0)
            IF(lwp) WRITE(numout,*) '                 Shchepetkin & McWilliams (2003) initial density profile'
            IF(lwp) WRITE(numout,*) '                 and linear EOS only function of temperature.'
            !
            drho = 3._wp
            dtem = drho / rn_a0 
            delta = 500._wp
            T0 = 10._wp
         CASE (1)
            IF(lwp) WRITE(numout,*) '                 Ezer, Arango and Shchepetkin (2002) initial temperature profile'
            IF(lwp) WRITE(numout,*) '                 and non-linear TEOS10.'
            dtem  = 15._wp
            delta = 1000._wp
            T0 = 5._wp
      END SELECT
      !

      IF (ln_init_pt_val) THEN
         DO_3D (0, 0, 0, 0, 1, jpk-1 ) 
            pts(ji,jj,jk,jp_tem) = ptmask(ji,jj,jk) * theta_initial( pdept(ji,jj,jk) )  
         END_3D

      ELSE ! grid point mean values 

         CALL usr_def_ist_3d ( ptmask, pts (:,:,:, jp_tem) ) 

         CALL lbc_lnk ( 'usr_def_istate', pts(:,:,:,jp_tem), 'T', 1._wp ) 
	 
      END IF


      S0 = 35._wp
      pts(:,:,:,jp_sal) = 35._wp * ptmask(:,:,:)

      !IF(lwp) WRITE(numout,*) '                     *) Estimated Burger number              rn_Snum  = ', rn_Snum
      ! Estimating Burger Number for initial density profile:
      ! S = (N * H) / (f * L) = SQRT(g * H * drho / rho_ref) / (f * L)
      !rn_Snum = SQRT(grav * rn_bot_max * rn_drho / 1000._wp) / (rn_fplane * rn_smnt_L * 1000._wp)
      !!----------------------------------------------------------------------
      !
   END SUBROUTINE usr_def_istate

   
   SUBROUTINE usr_def_istate_ssh( ptmask, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate_ssh  ***
      !! 
      !! ** Purpose :   Initialization of ssh
      !!
      !! ** Method  :   Set ssh as null, ptmask is required for test cases
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask   [m]
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height   [m]
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate_ssh : GYRE configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~   Ocean at rest, ssh is zero'
      !
      ! Sea level:
      pssh(:,:) = 0._wp
      !
   END SUBROUTINE usr_def_istate_ssh

   SUBROUTINE usr_def_ist_3d( ptmask, ptheta )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_ist_3d  ***
      !! 
      !! ** Purpose :   Initialization of the tracer field using grid cell mean values for the SEAMOUNT test case 
      !!
      !! ** Method  : - set grid cell means values of the (potential) temperature field
      !!                The bathymetry and the tracer profile are specifie in usr_def_ist_1d (called by usr_def_ist_2d) 
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   ptheta   ! cell mean values of potential temperature
      !

      REAL(wp), DIMENSION(A2D(nn_hls)) :: zt_v_m, zt_v, zt_v_p   ! vol integral of pt at levels jk-1/2, jk and jk+1/2 
      REAL(wp), DIMENSION(A2D(nn_hls)) :: zvol_m, zvol, zvol_p   ! vol integral at levels jk-1/2, jk and jk+1/2

      INTEGER :: ji, jj, jk  ! dummy loop indices
      INTEGER :: jkk         ! half level index  
      
      jkk = 1 
      CALL usr_def_ist_2d ( jkk, zvol_m, zt_v_m) 
      
      DO jk  = 1, jpk-1  
         jkk = 2*jk 
         CALL usr_def_ist_2d ( jkk, zvol, zt_v) 
         jkk = jkk + 1 
         CALL usr_def_ist_2d ( jkk, zvol_p, zt_v_p) 

         DO_2D (0, 0, 0, 0) 
            ptheta(ji,jj,jk) = ptmask(ji,jj,jk) *   ( zt_v_m(ji,jj) + 4.*zt_v(ji,jj) + zt_v_p(ji,jj) ) & 
	    &                                     / ( zvol_m(ji,jj) + 4.*zvol(ji,jj) + zvol_p(ji,jj) )  
         END_2D 

         DO_2D (0, 0, 0, 0) 
            zt_v_m(ji,jj) = zt_v_p(ji,jj) 
            zvol_m(ji,jj) = zvol_p(ji,jj)
         END_2D 

      END DO ! jk 

      RETURN 

   END SUBROUTINE usr_def_ist_3d

!!----------------------------------------------------------------------------

   SUBROUTINE usr_def_ist_2d( kkk, pvol, pt_v)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_ist_2d  ***
      !! 
      !! ** Purpose :   Calculate grid cell volume integrals of potential temperature for the SEAMOUNT test case; 
      !!
      !! ** Method  : - Calculate integrals in the jj direction; call usr_def_ist_1d to calculate integrals in the ji direction
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kkk   ! vertical "half" index jkk = 2*jk or 2*jk+1
      REAL(wp), DIMENSION(A2D(nn_hls))     , INTENT(  out) ::   pvol  ! volume integral of the (partial) grid cell
      REAL(wp), DIMENSION(A2D(nn_hls))     , INTENT(  out) ::   pt_v  ! volume integral of theta for the same partial grid cell
      !

      REAL(wp), DIMENSION(A2D(nn_hls)) :: z_dz_t, z_dz_u, z_dz_v, z_dz_f   ! depth at the t, u, v and f points in the cell  
      REAL(wp), DIMENSION(A2D(nn_hls)) :: zt_t,   zt_u,   zt_v,   zt_f     ! depth times tracer value at the t, u, v and f points in the cell

      INTEGER  :: ji, jj             ! dummy loop indices

      REAL(wp) z_d_lam, z_d_phi      ! half grid spacing in ji and jj directions  
      REAL(wp) zeta_kk               ! sigma value at this jkk index level 
      REAL(wp) zpk                   ! real version of jpk 
      REAL(wp) z_lamt, z_phit        ! "longitude" and latitude at tracer points
      REAL(wp) z_lamu, z_phiv        ! "longitude" and latitude at u and v points respectively 
      REAL(wp) z_lam_tpi, z_phi_tpj  ! "longitude" at t(ji+1,jj) ; "latitude" at t(ji,jj+1)
      REAL(wp) z_dep_t, z_dep_u, z_dep_v, z_dep_f ! depth at this sigma level at t, u, v and f points 
      REAL(wp) z_dep_upj             ! depth at u(ji,jj+1) 
      
      REAL(wp) z_dz_a, z_dz_b        ! contribution to cell volume (depth) from cell "edge" and "corner" points 
      REAL(wp) zt_a,   zt_b          ! contribution to cell volume integrated potential temperature from corresponding points 
      
      zpk = REAL(jpk, wp) 

      z_d_lam = 0.5 * ( glamt(2,1) - glamt(1,1) )   ! grid spacing assumed uniform;  glamt depends on ji 
      z_d_phi = 0.5 * ( gphit(1,2) - gphit(1,1) )   ! grid spacing assumed uniform;  gphit depends on jj  

! calculate the sigma coordinate for this call ; this version assumes uniform spacing of sigma
      zeta_kk = REAL( kkk - 1, wp )  / ( 2._wp * REAL(jpk) )     

      DO_2D (1, 0, 1, 0)      !   need to calculate for lower halo pts  

         z_lamt = glamt(ji,jj)
         z_phit = gphit(ji,jj)
         z_dep_t = zeta_kk * ocean_depth( z_lamt, z_phit )

         IF ( ln_init_curved ) THEN  ! calculate depths allowing curvature inside the cells 
            z_lamu = z_lamt + z_d_lam 
            z_phiv = z_phit + z_d_phi 
            z_dep_u = zeta_kk * ocean_depth( z_lamu, z_phit )
            z_dep_v = zeta_kk * ocean_depth( z_lamt, z_phiv )
            z_dep_f = zeta_kk * ocean_depth( z_lamu, z_phiv )

         ELSE    ! interpolate depths at T points to U, V and F points (following the model's calculations)  
            z_lam_tpi = z_lamt + 2._wp * z_d_lam    ! next T point in ji   
            z_phi_tpj = z_phit + 2._wp * z_d_phi    ! next T point in jj 
            z_dep_u   = 0.5_wp * ( z_dep_t + zeta_kk * ocean_depth( z_lam_tpi, z_phit ) )  
            z_dep_v   = 0.5_wp * ( z_dep_t + zeta_kk * ocean_depth( z_lamt, z_phi_tpj ) )
            z_dep_upj = 0.5_wp * zeta_kk *( ocean_depth( z_lamt,    z_phi_tpj ) +  &
	    &                               ocean_depth( z_lam_tpi, z_phi_tpj ) )
	    z_dep_f   = 0.5_wp * ( z_dep_u + z_dep_upj ) 

	 END IF 

! again assumes uniform spacing of sigma  
         z_dz_t(ji,jj) = z_dep_t / zpk 
         z_dz_u(ji,jj) = z_dep_u / zpk 
         z_dz_v(ji,jj) = z_dep_v / zpk 
         z_dz_f(ji,jj) = z_dep_f / zpk 

         zt_t(ji,jj) =  z_dz_t(ji,jj) * theta_initial( z_dep_t )
         zt_u(ji,jj) =  z_dz_u(ji,jj) * theta_initial( z_dep_u )
         zt_v(ji,jj) =  z_dz_v(ji,jj) * theta_initial( z_dep_v )
         zt_f(ji,jj) =  z_dz_f(ji,jj) * theta_initial( z_dep_f )

      END_2D 

      DO_2D (0, 0, 0, 0)      !   
         z_dz_a = 4._wp * ( z_dz_u(ji,jj) + z_dz_u(ji-1,jj) + z_dz_v(ji,jj) + z_dz_v(ji,jj-1) )
         zt_a   = 4._wp * (   zt_u(ji,jj) +   zt_u(ji-1,jj) +   zt_v(ji,jj) +   zt_v(ji,jj-1) ) 

         z_dz_b = z_dz_f(ji,jj) + z_dz_f(ji-1,jj) + z_dz_f(ji,jj-1) + z_dz_f(ji-1,jj-1) 
         zt_b   =   zt_f(ji,jj) +   zt_f(ji-1,jj) +   zt_f(ji,jj-1) +   zt_f(ji-1,jj-1) 
      
	 pvol(ji,jj) = 16._wp*z_dz_t(ji,jj) + z_dz_a + z_dz_b 
         pt_v(ji,jj) = 16._wp*  zt_t(ji,jj) +   zt_a +   zt_b    

      END_2D  

      RETURN 

   END SUBROUTINE usr_def_ist_2d   
   
   !!======================================================================

   !!======================================================================
END MODULE usrdef_istate
