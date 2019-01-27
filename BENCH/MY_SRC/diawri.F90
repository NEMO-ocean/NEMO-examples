MODULE diawri
   !!======================================================================
   !!                     ***  MODULE  diawri  ***
   !! Ocean diagnostics :  write ocean output files
   !!=====================================================================
   !! History :  OPA  ! 1991-03  (M.-A. Foujols)  Original code
   !!            4.0  ! 1991-11  (G. Madec)
   !!                 ! 1992-06  (M. Imbard)  correction restart file
   !!                 ! 1992-07  (M. Imbard)  split into diawri and rstwri
   !!                 ! 1993-03  (M. Imbard)  suppress writibm
   !!                 ! 1998-01  (C. Levy)  NETCDF format using ioipsl INTERFACE
   !!                 ! 1999-02  (E. Guilyardi)  name of netCDF files + variables
   !!            8.2  ! 2000-06  (M. Imbard)  Original code (diabort.F)
   !!   NEMO     1.0  ! 2002-06  (A.Bozec, E. Durand)  Original code (diainit.F)
   !!             -   ! 2002-09  (G. Madec)  F90: Free form and module
   !!             -   ! 2002-12  (G. Madec)  merge of diabort and diainit, F90
   !!                 ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.2  ! 2008-11  (B. Lemaire) creation from old diawri
   !!            3.7  ! 2014-01  (G. Madec) remove eddy induced velocity from no-IOM output
   !!                 !                     change name of output variables in dia_wri_state
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_wri       : create the standart output files
   !!   dia_wri_state : create an output NetCDF file for a single instantaeous ocean state and forcing fields
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers 
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE dianam         ! build name of file (routine)
   USE diahth         ! thermocline diagnostics
   USE dynadv   , ONLY: ln_dynadv_vec
   USE icb_oce        ! Icebergs
   USE icbdia         ! Iceberg budgets
   USE ldftra         ! lateral physics: eddy diffusivity coef.
   USE ldfdyn         ! lateral physics: eddy viscosity   coef.
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE sbc_ice        ! Surface boundary condition: ice fields
   USE sbcssr         ! restoring term toward SST/SSS climatology
   USE sbcwave        ! wave parameters
   USE wet_dry        ! wetting and drying
   USE zdf_oce        ! ocean vertical physics
   USE zdfdrg         ! ocean vertical physics: top/bottom friction
   USE zdfmxl         ! mixed layer
   !
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager ! I/O manager
   USE diatmb         ! Top,middle,bottom output
   USE dia25h         ! 25h Mean output
   USE iom            ! 
   USE ioipsl         ! 

#if defined key_si3
   USE ice
   USE icewri 
#endif
   USE lib_mpp         ! MPP library
   USE timing          ! preformance summary
   USE diurnal_bulk    ! diurnal warm layer
   USE cool_skin       ! Cool skin

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_wri                 ! routines called by step.F90
   PUBLIC   dia_wri_state
   PUBLIC   dia_wri_alloc           ! Called by nemogcm module

   INTEGER ::   nid_T, nz_T, nh_T, ndim_T, ndim_hT   ! grid_T file
   INTEGER ::          nb_T              , ndim_bT   ! grid_T file
   INTEGER ::   nid_U, nz_U, nh_U, ndim_U, ndim_hU   ! grid_U file
   INTEGER ::   nid_V, nz_V, nh_V, ndim_V, ndim_hV   ! grid_V file
   INTEGER ::   nid_W, nz_W, nh_W                    ! grid_W file
   INTEGER ::   ndex(1)                              ! ???
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_hT, ndex_hU, ndex_hV
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_T, ndex_U, ndex_V
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_bT

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: diawri.F90 9598 2018-05-15 22:47:16Z nicolasmartin $
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

#if defined key_iomput
   !!----------------------------------------------------------------------
   !!   'key_iomput'                                        use IOM library
   !!----------------------------------------------------------------------
   INTEGER FUNCTION dia_wri_alloc()
      !
      dia_wri_alloc = 0
      !
   END FUNCTION dia_wri_alloc

   
   SUBROUTINE dia_wri( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_wri  ***
      !!                   
      !! ** Purpose :   Standard output of opa: dynamics and tracer fields 
      !!      NETCDF format is used by default 
      !!
      !! ** Method  :  use iom_put
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      INTEGER ::   ji, jj, jk       ! dummy loop indices
      INTEGER ::   ikbot            ! local integer
      REAL(wp)::   zztmp , zztmpx   ! local scalar
      REAL(wp)::   zztmp2, zztmpy   !   -      -
      REAL(wp), DIMENSION(jpi,jpj)     ::   z2d   ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   z3d   ! 3D workspace
      !!----------------------------------------------------------------------
      ! 
      IF( ln_timing )   CALL timing_start('dia_wri')
      ! 
      ! Output the initial state and forcings
      IF( ninist == 1 ) THEN                       
         CALL dia_wri_state( 'output.init' )
         ninist = 0
      ENDIF

      ! Output of initial vertical scale factor
      CALL iom_put("e3t_0", e3t_0(:,:,:) )
      CALL iom_put("e3u_0", e3u_0(:,:,:) )
      CALL iom_put("e3v_0", e3v_0(:,:,:) )
      !
      CALL iom_put( "e3t" , e3t_n(:,:,:) )
      CALL iom_put( "e3u" , e3u_n(:,:,:) )
      CALL iom_put( "e3v" , e3v_n(:,:,:) )
      CALL iom_put( "e3w" , e3w_n(:,:,:) )
      IF( iom_use("e3tdef") )   &
         CALL iom_put( "e3tdef"  , ( ( e3t_n(:,:,:) - e3t_0(:,:,:) ) / e3t_0(:,:,:) * 100 * tmask(:,:,:) ) ** 2 )

      IF( ll_wd ) THEN
         CALL iom_put( "ssh" , (sshn+ssh_ref)*tmask(:,:,1) )   ! sea surface height (brought back to the reference used for wetting and drying)
      ELSE
         CALL iom_put( "ssh" , sshn )              ! sea surface height
      ENDIF

      IF( iom_use("wetdep") )   &                  ! wet depth
         CALL iom_put( "wetdep" , ht_0(:,:) + sshn(:,:) )
      
      CALL iom_put( "toce", tsn(:,:,:,jp_tem) )    ! 3D temperature
      CALL iom_put(  "sst", tsn(:,:,1,jp_tem) )    ! surface temperature
      IF ( iom_use("sbt") ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               ikbot = mbkt(ji,jj)
               z2d(ji,jj) = tsn(ji,jj,ikbot,jp_tem)
            END DO
         END DO
         CALL iom_put( "sbt", z2d )                ! bottom temperature
      ENDIF
      
      CALL iom_put( "soce", tsn(:,:,:,jp_sal) )    ! 3D salinity
      CALL iom_put(  "sss", tsn(:,:,1,jp_sal) )    ! surface salinity
      IF ( iom_use("sbs") ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               ikbot = mbkt(ji,jj)
               z2d(ji,jj) = tsn(ji,jj,ikbot,jp_sal)
            END DO
         END DO
         CALL iom_put( "sbs", z2d )                ! bottom salinity
      ENDIF

      IF ( iom_use("taubot") ) THEN                ! bottom stress
         zztmp = rau0 * 0.25
         z2d(:,:) = 0._wp
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zztmp2 = (  ( rCdU_bot(ji+1,jj)+rCdU_bot(ji  ,jj) ) * un(ji  ,jj,mbku(ji  ,jj))  )**2   &
                  &   + (  ( rCdU_bot(ji  ,jj)+rCdU_bot(ji-1,jj) ) * un(ji-1,jj,mbku(ji-1,jj))  )**2   &
                  &   + (  ( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj  ) ) * vn(ji,jj  ,mbkv(ji,jj  ))  )**2   &
                  &   + (  ( rCdU_bot(ji,jj  )+rCdU_bot(ji,jj-1) ) * vn(ji,jj-1,mbkv(ji,jj-1))  )**2
               z2d(ji,jj) = zztmp * SQRT( zztmp2 ) * tmask(ji,jj,1) 
               !
            END DO
         END DO
         CALL lbc_lnk( 'diawri', z2d, 'T', 1. )
         CALL iom_put( "taubot", z2d )           
      ENDIF
         
      CALL iom_put( "uoce", un(:,:,:) )            ! 3D i-current
      CALL iom_put(  "ssu", un(:,:,1) )            ! surface i-current
      IF ( iom_use("sbu") ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               ikbot = mbku(ji,jj)
               z2d(ji,jj) = un(ji,jj,ikbot)
            END DO
         END DO
         CALL iom_put( "sbu", z2d )                ! bottom i-current
      ENDIF
      
      CALL iom_put( "voce", vn(:,:,:) )            ! 3D j-current
      CALL iom_put(  "ssv", vn(:,:,1) )            ! surface j-current
      IF ( iom_use("sbv") ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               ikbot = mbkv(ji,jj)
               z2d(ji,jj) = vn(ji,jj,ikbot)
            END DO
         END DO
         CALL iom_put( "sbv", z2d )                ! bottom j-current
      ENDIF

      CALL iom_put( "woce", wn )                   ! vertical velocity
      IF( iom_use('w_masstr') .OR. iom_use('w_masstr2') ) THEN   ! vertical mass transport & its square value
         ! Caution: in the VVL case, it only correponds to the baroclinic mass transport.
         z2d(:,:) = rau0 * e1e2t(:,:)
         DO jk = 1, jpk
            z3d(:,:,jk) = wn(:,:,jk) * z2d(:,:)
         END DO
         CALL iom_put( "w_masstr" , z3d )  
         IF( iom_use('w_masstr2') )   CALL iom_put( "w_masstr2", z3d(:,:,:) * z3d(:,:,:) )
      ENDIF

      CALL iom_put( "avt" , avt )                  ! T vert. eddy diff. coef.
      CALL iom_put( "avs" , avs )                  ! S vert. eddy diff. coef.
      CALL iom_put( "avm" , avm )                  ! T vert. eddy visc. coef.

      IF( iom_use('logavt') )   CALL iom_put( "logavt", LOG( MAX( 1.e-20_wp, avt(:,:,:) ) ) )
      IF( iom_use('logavs') )   CALL iom_put( "logavs", LOG( MAX( 1.e-20_wp, avs(:,:,:) ) ) )

      IF ( iom_use("sstgrad") .OR. iom_use("sstgrad2") ) THEN
         DO jj = 2, jpjm1                                    ! sst gradient
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zztmp  = tsn(ji,jj,1,jp_tem)
               zztmpx = ( tsn(ji+1,jj,1,jp_tem) - zztmp ) * r1_e1u(ji,jj) + ( zztmp - tsn(ji-1,jj  ,1,jp_tem) ) * r1_e1u(ji-1,jj)
               zztmpy = ( tsn(ji,jj+1,1,jp_tem) - zztmp ) * r1_e2v(ji,jj) + ( zztmp - tsn(ji  ,jj-1,1,jp_tem) ) * r1_e2v(ji,jj-1)
               z2d(ji,jj) = 0.25 * ( zztmpx * zztmpx + zztmpy * zztmpy )   &
                  &              * umask(ji,jj,1) * umask(ji-1,jj,1) * vmask(ji,jj,1) * umask(ji,jj-1,1)
            END DO
         END DO
         CALL lbc_lnk( 'diawri', z2d, 'T', 1. )
         CALL iom_put( "sstgrad2",  z2d )          ! square of module of sst gradient
         z2d(:,:) = SQRT( z2d(:,:) )
         CALL iom_put( "sstgrad" ,  z2d )          ! module of sst gradient
      ENDIF
         
      ! heat and salt contents
      IF( iom_use("heatc") ) THEN
         z2d(:,:)  = 0._wp 
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  z2d(ji,jj) = z2d(ji,jj) + e3t_n(ji,jj,jk) * tsn(ji,jj,jk,jp_tem) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL iom_put( "heatc", rau0_rcp * z2d )   ! vertically integrated heat content (J/m2)
      ENDIF

      IF( iom_use("saltc") ) THEN
         z2d(:,:)  = 0._wp 
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  z2d(ji,jj) = z2d(ji,jj) + e3t_n(ji,jj,jk) * tsn(ji,jj,jk,jp_sal) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL iom_put( "saltc", rau0 * z2d )          ! vertically integrated salt content (PSU*kg/m2)
      ENDIF
      !
      IF ( iom_use("eken") ) THEN
         z3d(:,:,jpk) = 0._wp 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zztmp  = 0.25_wp * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
                  z3d(ji,jj,jk) = zztmp * (  un(ji-1,jj,jk)**2 * e2u(ji-1,jj) * e3u_n(ji-1,jj,jk)   &
                     &                     + un(ji  ,jj,jk)**2 * e2u(ji  ,jj) * e3u_n(ji  ,jj,jk)   &
                     &                     + vn(ji,jj-1,jk)**2 * e1v(ji,jj-1) * e3v_n(ji,jj-1,jk)   &
                     &                     + vn(ji,jj  ,jk)**2 * e1v(ji,jj  ) * e3v_n(ji,jj  ,jk)   )
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'diawri', z3d, 'T', 1. )
         CALL iom_put( "eken", z3d )                 ! kinetic energy
      ENDIF
      !
      CALL iom_put( "hdiv", hdivn )                  ! Horizontal divergence
      !
      IF( iom_use("u_masstr") .OR. iom_use("u_masstr_vint") .OR. iom_use("u_heattr") .OR. iom_use("u_salttr") ) THEN
         z3d(:,:,jpk) = 0.e0
         z2d(:,:) = 0.e0
         DO jk = 1, jpkm1
            z3d(:,:,jk) = rau0 * un(:,:,jk) * e2u(:,:) * e3u_n(:,:,jk) * umask(:,:,jk)
            z2d(:,:) = z2d(:,:) + z3d(:,:,jk)
         END DO
         CALL iom_put( "u_masstr"     , z3d )         ! mass transport in i-direction
         CALL iom_put( "u_masstr_vint", z2d )         ! mass transport in i-direction vertical sum
      ENDIF
      
      IF( iom_use("u_heattr") ) THEN
         z2d(:,:) = 0._wp 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z2d(ji,jj) = z2d(ji,jj) + z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_tem) + tsn(ji+1,jj,jk,jp_tem) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'diawri', z2d, 'U', -1. )
         CALL iom_put( "u_heattr", 0.5*rcp * z2d )    ! heat transport in i-direction
      ENDIF

      IF( iom_use("u_salttr") ) THEN
         z2d(:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z2d(ji,jj) = z2d(ji,jj) + z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_sal) + tsn(ji+1,jj,jk,jp_sal) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'diawri', z2d, 'U', -1. )
         CALL iom_put( "u_salttr", 0.5 * z2d )        ! heat transport in i-direction
      ENDIF

      
      IF( iom_use("v_masstr") .OR. iom_use("v_heattr") .OR. iom_use("v_salttr") ) THEN
         z3d(:,:,jpk) = 0.e0
         DO jk = 1, jpkm1
            z3d(:,:,jk) = rau0 * vn(:,:,jk) * e1v(:,:) * e3v_n(:,:,jk) * vmask(:,:,jk)
         END DO
         CALL iom_put( "v_masstr", z3d )              ! mass transport in j-direction
      ENDIF
      
      IF( iom_use("v_heattr") ) THEN
         z2d(:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z2d(ji,jj) = z2d(ji,jj) + z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_tem) + tsn(ji,jj+1,jk,jp_tem) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'diawri', z2d, 'V', -1. )
         CALL iom_put( "v_heattr", 0.5*rcp * z2d )    !  heat transport in j-direction
      ENDIF

      IF( iom_use("v_salttr") ) THEN
         z2d(:,:) = 0._wp 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z2d(ji,jj) = z2d(ji,jj) + z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_sal) + tsn(ji,jj+1,jk,jp_sal) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'diawri', z2d, 'V', -1. )
         CALL iom_put( "v_salttr", 0.5 * z2d )        !  heat transport in j-direction
      ENDIF

      IF( iom_use("tosmint") ) THEN
         z2d(:,:) = 0._wp
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z2d(ji,jj) = z2d(ji,jj) + e3t_n(ji,jj,jk) *  tsn(ji,jj,jk,jp_tem)
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'diawri', z2d, 'T', -1. )
         CALL iom_put( "tosmint", rau0 * z2d )        ! Vertical integral of temperature
      ENDIF
      IF( iom_use("somint") ) THEN
         z2d(:,:)=0._wp
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z2d(ji,jj) = z2d(ji,jj) + e3t_n(ji,jj,jk) * tsn(ji,jj,jk,jp_sal)
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'diawri', z2d, 'T', -1. )
         CALL iom_put( "somint", rau0 * z2d )         ! Vertical integral of salinity
      ENDIF

      CALL iom_put( "bn2", rn2 )                      ! Brunt-Vaisala buoyancy frequency (N^2)
      !

      IF (ln_diatmb)   CALL dia_tmb                   ! tmb values 
          
      IF (ln_dia25h)   CALL dia_25h( kt )             ! 25h averaging

      IF( ln_timing )   CALL timing_stop('dia_wri')
      !
   END SUBROUTINE dia_wri

#else
   !!----------------------------------------------------------------------
   !!   Default option                                  use IOIPSL  library
   !!----------------------------------------------------------------------

   INTEGER FUNCTION dia_wri_alloc()
      !
      dia_wri_alloc = 0
      !
   END FUNCTION dia_wri_alloc

   SUBROUTINE dia_wri( kt )
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      IF( ninist == 1 ) THEN     !==  Output the initial state and forcings  ==!
         CALL dia_wri_state( 'output.init' )
         ninist = 0
      ENDIF

   END SUBROUTINE dia_wri

#endif

   SUBROUTINE dia_wri_state( cdfile_name )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dia_wri_state  ***
      !!        
      !! ** Purpose :   create a NetCDF file named cdfile_name which contains 
      !!      the instantaneous ocean state and forcing fields.
      !!        Used to find errors in the initial state or save the last
      !!      ocean state in case of abnormal end of a simulation
      !!
      !! ** Method  :   NetCDF files using ioipsl
      !!      File 'output.init.nc'  is created if ninist = 1 (namelist)
      !!      File 'output.abort.nc' is created in case of abnormal job end
      !!----------------------------------------------------------------------
      CHARACTER (len=* ), INTENT( in ) ::   cdfile_name      ! name of the file created
      !!
      INTEGER :: inum
      !!----------------------------------------------------------------------
      ! 
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dia_wri_state : single instantaneous ocean state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~   and forcing fields file created '
      IF(lwp) WRITE(numout,*) '                and named :', cdfile_name, '...nc'

#if defined key_si3
     CALL iom_open( TRIM(cdfile_name), inum, ldwrt = .TRUE., kdlev = jpl )
#else
     CALL iom_open( TRIM(cdfile_name), inum, ldwrt = .TRUE. )
#endif

      CALL iom_rstput( 0, 0, inum, 'votemper', tsn(:,:,:,jp_tem) )    ! now temperature
      CALL iom_rstput( 0, 0, inum, 'vosaline', tsn(:,:,:,jp_sal) )    ! now salinity
      CALL iom_rstput( 0, 0, inum, 'sossheig', sshn              )    ! sea surface height
      CALL iom_rstput( 0, 0, inum, 'vozocrtx', un                )    ! now i-velocity
      CALL iom_rstput( 0, 0, inum, 'vomecrty', vn                )    ! now j-velocity
      CALL iom_rstput( 0, 0, inum, 'vovecrtz', wn                )    ! now k-velocity
      IF( ALLOCATED(ahtu) ) THEN
         CALL iom_rstput( 0, 0, inum,  'ahtu', ahtu              )    ! aht at u-point
         CALL iom_rstput( 0, 0, inum,  'ahtv', ahtv              )    ! aht at v-point
      ENDIF
      IF( ALLOCATED(ahmt) ) THEN
         CALL iom_rstput( 0, 0, inum,  'ahmt', ahmt              )    ! ahmt at u-point
         CALL iom_rstput( 0, 0, inum,  'ahmf', ahmf              )    ! ahmf at v-point
      ENDIF
      CALL iom_rstput( 0, 0, inum, 'sowaflup', emp - rnf         )    ! freshwater budget
      CALL iom_rstput( 0, 0, inum, 'sohefldo', qsr + qns         )    ! total heat flux
      CALL iom_rstput( 0, 0, inum, 'soshfldo', qsr               )    ! solar heat flux
      CALL iom_rstput( 0, 0, inum, 'soicecov', fr_i              )    ! ice fraction
      CALL iom_rstput( 0, 0, inum, 'sozotaux', utau              )    ! i-wind stress
      CALL iom_rstput( 0, 0, inum, 'sometauy', vtau              )    ! j-wind stress
      IF(  .NOT.ln_linssh  ) THEN             
         CALL iom_rstput( 0, 0, inum, 'vovvldep', gdept_n        )    !  T-cell depth 
         CALL iom_rstput( 0, 0, inum, 'vovvle3t', e3t_n          )    !  T-cell thickness  
      END IF
      IF( ln_wave .AND. ln_sdw ) THEN
         CALL iom_rstput( 0, 0, inum, 'sdzocrtx', usd            )    ! now StokesDrift i-velocity
         CALL iom_rstput( 0, 0, inum, 'sdmecrty', vsd            )    ! now StokesDrift j-velocity
         CALL iom_rstput( 0, 0, inum, 'sdvecrtz', wsd            )    ! now StokesDrift k-velocity
      ENDIF
 
#if defined key_si3
      IF( nn_ice == 2 ) THEN   ! condition needed in case agrif + ice-model but no-ice in child grid
         CALL ice_wri_state( inum )
      ENDIF
#endif
      !
      CALL iom_close( inum )
      ! 
   END SUBROUTINE dia_wri_state

   !!======================================================================
END MODULE diawri
