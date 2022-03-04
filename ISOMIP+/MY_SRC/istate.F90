MODULE istate
   !!======================================================================
   !!                     ***  MODULE  istate  ***
   !! Ocean state   :  initial state setting
   !!=====================================================================
   !! History :  OPA  !  1989-12  (P. Andrich)  Original code
   !!            5.0  !  1991-11  (G. Madec)  rewritting
   !!            6.0  !  1996-01  (G. Madec)  terrain following coordinates
   !!            8.0  !  2001-09  (M. Levy, M. Ben Jelloul)  istate_eel
   !!            8.0  !  2001-09  (M. Levy, M. Ben Jelloul)  istate_uvg
   !!   NEMO     1.0  !  2003-08  (G. Madec, C. Talandier)  F90: Free form, modules + EEL R5
   !!             -   !  2004-05  (A. Koch-Larrouy)  istate_gyre 
   !!            2.0  !  2006-07  (S. Masson)  distributed restart using iom
   !!            3.3  !  2010-10  (C. Ethe) merge TRC-TRA
   !!            3.4  !  2011-04  (G. Madec) Merge of dtatem and dtasal & suppression of tb,tn/sb,sn 
   !!            3.7  !  2016-04  (S. Flavoni) introduce user defined initial state 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   istate_init   : initial state setting
   !!   istate_uvg    : initial velocity in geostropic balance
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers 
   USE dom_oce        ! ocean space and time domain 
   USE daymod         ! calendar
   USE dtatsd         ! data temperature and salinity   (dta_tsd routine)
   USE dtauvd         ! data: U & V current             (dta_uvd routine)
   USE domvvl          ! varying vertical mesh
   USE wet_dry         ! wetting and drying (needed for wad_istate)
   USE usrdef_istate   ! User defined initial state
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O library
   USE lib_mpp         ! MPP library
   USE restart         ! restart

#if defined key_agrif
   USE agrif_oce       ! initial state interpolation
   USE agrif_oce_interp
#endif   

   IMPLICIT NONE
   PRIVATE

   PUBLIC   istate_init   ! routine called by nemogcm.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: istate.F90 11423 2019-08-08 14:02:49Z mathiot $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE istate_init( Kbb, Kmm, Kaa )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE istate_init  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracer fields.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::  Kbb, Kmm, Kaa   ! ocean time level indices
      !
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zgdept     ! 3D table for qco substitute
!!gm see comment further down
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   zuvd    ! U & V data workspace
!!gm end
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'istate_init : Initialization of the dynamics and tracers'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'

      CALL dta_tsd_init                 ! Initialisation of T & S input data
      IF( ln_c1d) CALL dta_uvd_init     ! Initialisation of U & V input data (c1d only)

      rhd  (:,:,:      ) = 0._wp   ;   rhop (:,:,:  ) = 0._wp      ! set one for all to 0 at level jpk
      rn2b (:,:,:      ) = 0._wp   ;   rn2  (:,:,:  ) = 0._wp      ! set one for all to 0 at levels 1 and jpk
      ts   (:,:,:,:,Kaa) = 0._wp                                   ! set one for all to 0 at level jpk
      rab_b(:,:,:,:    ) = 0._wp   ;   rab_n(:,:,:,:) = 0._wp      ! set one for all to 0 at level jpk
#if defined key_agrif
      uu   (:,:,:  ,Kaa) = 0._wp   ! used in agrif_oce_sponge at initialization
      vv   (:,:,:  ,Kaa) = 0._wp   ! used in agrif_oce_sponge at initialization    
#endif

#if defined key_agrif
      IF ( (.NOT.Agrif_root()).AND.ln_init_chfrpar ) THEN
         numror = 0                           ! define numror = 0 -> no restart file to read
         ln_1st_euler = .true.                ! Set time-step indicator at nit000 (euler forward)
         CALL day_init 
         CALL agrif_istate_oce( Kbb, Kmm, Kaa )   ! Interp from parent
         !
         ts (:,:,:,:,Kmm) = ts (:,:,:,:,Kbb)
         uu (:,:,:  ,Kmm) = uu (:,:,:  ,Kbb)
         vv (:,:,:  ,Kmm) = vv (:,:,:  ,Kbb)
      ELSE
#endif
         IF( ln_rstart ) THEN                    ! Restart from a file
            !                                    ! -------------------
            CALL rst_read( Kbb, Kmm )            ! Read the restart file
            CALL day_init                        ! model calendar (using both namelist and restart infos)
            !
         ELSE                                    ! Start from rest
            !                                    ! ---------------
            numror = 0                           ! define numror = 0 -> no restart file to read
            l_1st_euler = .true.                 ! Set time-step indicator at nit000 (euler forward)
            CALL day_init                        ! model calendar (using both namelist and restart infos)
            !                                    ! Initialization of ocean to zero
            !
            IF( ln_tsd_init ) THEN               
               CALL dta_tsd( nit000, 'ini', ts(:,:,:,:,Kbb) )       ! read 3D T and S data at nit000
            ENDIF
            !
            IF( ln_uvd_init .AND. ln_c1d ) THEN               
               CALL dta_uvd( nit000, Kbb, uu(:,:,:,Kbb), vv(:,:,:,Kbb) )   ! read 3D U and V data at nit000
            ELSE
               uu  (:,:,:,Kbb) = 0._wp               ! set the ocean at rest
               vv  (:,:,:,Kbb) = 0._wp  
            ENDIF
               !
               !
            IF( .NOT. ln_tsd_init .AND. .NOT. ln_uvd_init ) THEN
               DO jk = 1, jpk
                  zgdept(:,:,jk) = gdept(:,:,jk,Kbb)
               END DO
               CALL usr_def_istate( zgdept, tmask, ts(:,:,:,:,Kbb), uu(:,:,:,Kbb), vv(:,:,:,Kbb) )         
            ENDIF
            ts  (:,:,:,:,Kmm) = ts (:,:,:,:,Kbb)       ! set now values from to before ones
            uu    (:,:,:,Kmm) = uu   (:,:,:,Kbb)
            vv    (:,:,:,Kmm) = vv   (:,:,:,Kbb)

         ENDIF 
#if defined key_agrif
      ENDIF
#endif
      ! 
      ! Initialize "now" and "before" barotropic velocities:
      ! Do it whatever the free surface method, these arrays being eventually used
      !
      uu_b(:,:,Kmm) = 0._wp   ;   vv_b(:,:,Kmm) = 0._wp
      uu_b(:,:,Kbb) = 0._wp   ;   vv_b(:,:,Kbb) = 0._wp
      !
!!gm  the use of umsak & vmask is not necessary below as uu(:,:,:,Kmm), vv(:,:,:,Kmm), uu(:,:,:,Kbb), vv(:,:,:,Kbb) are always masked
      DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpkm1 )
         uu_b(ji,jj,Kmm) = uu_b(ji,jj,Kmm) + e3u(ji,jj,jk,Kmm) * uu(ji,jj,jk,Kmm) * umask(ji,jj,jk)
         vv_b(ji,jj,Kmm) = vv_b(ji,jj,Kmm) + e3v(ji,jj,jk,Kmm) * vv(ji,jj,jk,Kmm) * vmask(ji,jj,jk)
         !
         uu_b(ji,jj,Kbb) = uu_b(ji,jj,Kbb) + e3u(ji,jj,jk,Kbb) * uu(ji,jj,jk,Kbb) * umask(ji,jj,jk)
         vv_b(ji,jj,Kbb) = vv_b(ji,jj,Kbb) + e3v(ji,jj,jk,Kbb) * vv(ji,jj,jk,Kbb) * vmask(ji,jj,jk)
      END_3D
      !
      uu_b(:,:,Kmm) = uu_b(:,:,Kmm) * r1_hu(:,:,Kmm)
      vv_b(:,:,Kmm) = vv_b(:,:,Kmm) * r1_hv(:,:,Kmm)
      !
      uu_b(:,:,Kbb) = uu_b(:,:,Kbb) * r1_hu(:,:,Kbb)
      vv_b(:,:,Kbb) = vv_b(:,:,Kbb) * r1_hv(:,:,Kbb)
      !
   END SUBROUTINE istate_init

   !!======================================================================
END MODULE istate
