MODULE isf_oce
   !!======================================================================
   !!                       ***  MODULE  isf_oce  ***
   !! Ice shelves :  ice shelves variables defined in memory
   !!======================================================================
   !! History :  3.2  !  2011-02  (C.Harris  ) Original code isf cav
   !!            X.X  !  2006-02  (C. Wang   ) Original code bg03
   !!            3.4  !  2013-03  (P. Mathiot) Merging + parametrization
   !!            4.1  !  2019-09  (P. Mathiot) Split param/explicit ice shelf and re-organisation
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   isf          : define and allocate ice shelf variables
   !!----------------------------------------------------------------------

   USE par_oce       , ONLY: jpi, jpj, jpk
   USE in_out_manager, ONLY: wp, jpts ! I/O manager
   USE lib_mpp       , ONLY: ctl_stop, mpp_sum      ! MPP library
   USE fldread        ! read input fields

   IMPLICIT NONE

   PRIVATE

   PUBLIC   isf_alloc, isf_alloc_par, isf_alloc_cav, isf_alloc_cpl, isf_dealloc_cpl
   !
   !-------------------------------------------------------
   ! 0 :              namelist parameter
   !-------------------------------------------------------
   !
   ! 0.1 -------- ice shelf cavity parameter --------------
   CHARACTER(LEN=256), PUBLIC :: cn_isfdir
   LOGICAL           , PUBLIC :: ln_isf
   LOGICAL           , PUBLIC :: ln_isfdebug
   !
   ! 0.2 -------- ice shelf cavity opened namelist parameter -------------
   LOGICAL           , PUBLIC :: ln_isfcav_mlt   !: logical for the use of ice shelf parametrisation
   REAL(wp)          , PUBLIC :: rn_gammat0      !: temperature exchange coeficient    []
   REAL(wp)          , PUBLIC :: rn_gammas0      !: salinity    exchange coeficient    []
   REAL(wp)          , PUBLIC :: rn_vtide        !: tidal background velocity (can be different to what is used in the 
   REAL(wp)          , PUBLIC :: rn_htbl         !: Losch top boundary layer thickness [m]
   REAL(wp)          , PUBLIC :: rn_isfload_T    !: 
   REAL(wp)          , PUBLIC :: rn_isfload_S    !: 
   CHARACTER(LEN=256), PUBLIC :: cn_gammablk     !: gamma formulation
   CHARACTER(LEN=256), PUBLIC :: cn_isfcav_mlt   !: melt formulation (cavity/param)
   CHARACTER(LEN=256), PUBLIC :: cn_isfload      !: ice shelf load computation method
   TYPE(FLD_N)       , PUBLIC :: sn_isfcav_fwf   !: information about the isf melting file to be read
   !
   ! 0.3 -------- ice shelf cavity parametrised namelist parameter -------------
   LOGICAL           , PUBLIC :: ln_isfpar_mlt      !: logical for the computation of melt inside the cavity
   REAL(wp)          , PUBLIC :: rn_isfpar_bg03_gt0 !: temperature exchange coeficient [m/s]
   CHARACTER(LEN=256), PUBLIC :: cn_isfpar_mlt      !: melt formulation (cavity/param)
   TYPE(FLD_N)       , PUBLIC :: sn_isfpar_fwf      !: information about the isf melting file to be read
   TYPE(FLD_N)       , PUBLIC :: sn_isfpar_zmax     !: information about the grounding line depth file to be read
   TYPE(FLD_N)       , PUBLIC :: sn_isfpar_zmin     !: information about the calving   line depth file to be read
   TYPE(FLD_N)       , PUBLIC :: sn_isfpar_Leff     !: information about the effective length     file to be read
   !
   ! 0.4 -------- coupling namelist parameter -------------
   LOGICAL, PUBLIC :: ln_isfcpl      !:
   LOGICAL, PUBLIC :: ln_isfcpl_cons !:
   INTEGER, PUBLIC :: nn_drown       !:
   !
   !-------------------------------------------------------
   ! 1 :              ice shelf parameter
   !-------------------------------------------------------
   !
   REAL(wp), PARAMETER, PUBLIC :: rLfusisf = 0.334e6_wp    !: latent heat of fusion of ice shelf     [J/kg]
   REAL(wp), PARAMETER, PUBLIC :: rcpisf = 2000.0_wp       !: specific heat of ice shelf             [J/kg/K]
   REAL(wp), PARAMETER, PUBLIC :: rkappa = 0.0_wp          !: ISOMIP+ no heat diffusivity through the ice-shelf [m2/s]
   REAL(wp), PARAMETER, PUBLIC :: rhoisf = 920.0_wp        !: volumic mass of ice shelf              [kg/m3]
   REAL(wp), PARAMETER, PUBLIC :: rtsurf = -20.0           !: surface temperature                    [C]
   !
   !-------------------------------------------------------
   ! 2 :              ice shelf global variables
   !-------------------------------------------------------
   !
   ! 2.1 -------- ice shelf cavity parameter --------------
   LOGICAL , PUBLIC            :: l_isfoasis = .FALSE.
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)    ::   risfload                    !: ice shelf load
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)    ::   fwfisf_oasis
   !
   ! 2.2 -------- ice shelf cavity melt namelist parameter -------------
   INTEGER  , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: mskisf_cav                    !:
   INTEGER  , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: misfkt_cav   , misfkb_cav     !: 
   REAL(wp) , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: rhisf_tbl_cav, rfrac_tbl_cav  !: 
   REAL(wp) , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fwfisf_cav   , fwfisf_cav_b   !: before and now net fwf from the ice shelf        [kg/m2/s]
   REAL(wp) , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: risf_cav_tsc , risf_cav_tsc_b !: before and now T & S isf contents [K.m/s & PSU.m/s]  
   TYPE(FLD), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)     :: sf_isfcav_fwf                 !:
   !
   REAL(wp) , PUBLIC                                      :: risf_lamb1, risf_lamb2, risf_lamb3  ! freezing point linearization coeficient
   !
   ! 2.3 -------- ice shelf param. melt namelist parameter -------------
   INTEGER  , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: mskisf_par                    !:
   INTEGER  , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: misfkt_par   , misfkb_par     !:
   REAL(wp) , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: rhisf_tbl_par, rfrac_tbl_par  !: 
   REAL(wp) , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fwfisf_par   , fwfisf_par_b   !: before and now net fwf from the ice shelf        [kg/m2/s]
   REAL(wp) , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: risf_par_tsc , risf_par_tsc_b !: before and now T & S isf contents [K.m/s & PSU.m/s]  
   TYPE(FLD), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)     :: sf_isfpar_fwf                 !:
   !
   REAL(wp) , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: rhisf0_tbl_par                !: thickness of tbl (initial value)  [m]
   REAL(wp) , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: risfLeff                      !:
   !
   ! 2.4 -------- coupling namelist parameter -------------
   INTEGER , PUBLIC                                        ::   nstp_iscpl   !:
   REAL(wp), PUBLIC                                        ::   rdt_iscpl    !: 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   risfcpl_ssh, risfcpl_cons_ssh, risfcpl_cons_ssh_b               !:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   risfcpl_vol, risfcpl_cons_vol, risfcpl_cons_vol_b  !:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   risfcpl_tsc, risfcpl_cons_tsc, risfcpl_cons_tsc_b  !:
   !
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcisf.F90 10536 2019-01-16 19:21:09Z mathiot $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE isf_alloc_par()
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isf_alloc_par  ***
      !!
      !! ** Purpose : 
      !!
      !! ** Method  : 
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ierr, ialloc
      !!----------------------------------------------------------------------
      ierr = 0       ! set to zero if no array to be allocated
      !
      ALLOCATE(risfLeff(jpi,jpj), STAT=ialloc)
      ierr = ierr + ialloc
      !
      ALLOCATE(misfkt_par(jpi,jpj), misfkb_par(jpi,jpj), STAT=ialloc )
      ierr = ierr + ialloc
      !
      ALLOCATE( rfrac_tbl_par(jpi,jpj), STAT=ialloc)
      ierr = ierr + ialloc
      !
      ALLOCATE( rhisf_tbl_par(jpi,jpj), rhisf0_tbl_par(jpi,jpj), STAT=ialloc)
      ierr = ierr + ialloc
      !
      ALLOCATE( mskisf_par(jpi,jpj), STAT=ialloc)
      ierr = ierr + ialloc
      !
      CALL mpp_sum ( 'isf', ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'isf: failed to allocate arrays.' )
      !
   END SUBROUTINE isf_alloc_par

   
   SUBROUTINE isf_alloc_cav()
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isf_alloc_cav  ***
      !!
      !! ** Purpose : 
      !!
      !! ** Method  : 
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ierr, ialloc
      !!----------------------------------------------------------------------
      ierr = 0       ! set to zero if no array to be allocated
      !
      ALLOCATE(misfkt_cav(jpi,jpj), misfkb_cav(jpi,jpj), STAT=ialloc )
      ierr = ierr + ialloc
      !
      ALLOCATE( rfrac_tbl_cav(jpi,jpj), STAT=ialloc)
      ierr = ierr + ialloc
      !
      ALLOCATE( rhisf_tbl_cav(jpi,jpj), STAT=ialloc)
      ierr = ierr + ialloc
      !
      CALL mpp_sum ( 'isf', ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'isf: failed to allocate arrays.' )
      !
   END SUBROUTINE isf_alloc_cav

   
   SUBROUTINE isf_alloc_cpl()
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isf_alloc_cpl  ***
      !!
      !! ** Purpose : allocate array use for the ice sheet coupling
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ierr, ialloc
      !!----------------------------------------------------------------------
      ierr = 0
      !
      ALLOCATE( risfcpl_ssh(jpi,jpj) , risfcpl_tsc(jpi,jpj,jpk,jpts) , risfcpl_vol(jpi,jpj,jpk) , STAT=ialloc )
      ierr = ierr + ialloc
      !
      risfcpl_tsc(:,:,:,:) = 0._wp ; risfcpl_vol(:,:,:) = 0._wp ; risfcpl_ssh(:,:) = 0._wp

      IF ( ln_isfcpl_cons ) THEN
         ALLOCATE( risfcpl_cons_tsc(jpi,jpj,jpk,jpts) , risfcpl_cons_vol(jpi,jpj,jpk) , risfcpl_cons_ssh(jpi,jpj) , STAT=ialloc )
         ierr = ierr + ialloc
         !
         risfcpl_cons_tsc(:,:,:,:) = 0._wp ; risfcpl_cons_vol(:,:,:) = 0._wp ; risfcpl_cons_ssh(:,:) = 0._wp
         !
      END IF
      !
      CALL mpp_sum ( 'isf', ierr )
      IF( ierr /= 0 )   CALL ctl_stop('STOP','isfcpl: failed to allocate arrays.')
      !
   END SUBROUTINE isf_alloc_cpl

   
   SUBROUTINE isf_dealloc_cpl()
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isf_dealloc_cpl  ***
      !!
      !! ** Purpose : de-allocate useless public 3d array used for ice sheet coupling
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ierr, ialloc
      !!----------------------------------------------------------------------
      ierr = 0
      !
      DEALLOCATE( risfcpl_ssh , risfcpl_tsc , risfcpl_vol , STAT=ialloc )
      ierr = ierr + ialloc
      !
      CALL mpp_sum ( 'isf', ierr )
      IF( ierr /= 0 )   CALL ctl_stop('STOP','isfcpl: failed to deallocate arrays.')
      !
   END SUBROUTINE isf_dealloc_cpl

   
   SUBROUTINE isf_alloc()
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isf_alloc  ***
      !!
      !! ** Purpose : allocate array used for the ice shelf cavity (cav and par)
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ierr, ialloc
      !!----------------------------------------------------------------------
      !
      ierr = 0       ! set to zero if no array to be allocated
      !
      ALLOCATE( fwfisf_par  (jpi,jpj) , fwfisf_par_b(jpi,jpj) ,     &
         &      fwfisf_cav  (jpi,jpj) , fwfisf_cav_b(jpi,jpj) ,     &
         &      fwfisf_oasis(jpi,jpj)                         , STAT=ialloc )
      ierr = ierr + ialloc
      !
      ALLOCATE( risf_par_tsc(jpi,jpj,jpts) , risf_par_tsc_b(jpi,jpj,jpts) , STAT=ialloc )
      ierr = ierr + ialloc
      !
      ALLOCATE( risf_cav_tsc(jpi,jpj,jpts) , risf_cav_tsc_b(jpi,jpj,jpts) , STAT=ialloc )
      ierr = ierr + ialloc
      !
      ALLOCATE( risfload(jpi,jpj) , STAT=ialloc )
      ierr = ierr + ialloc
      !
      ALLOCATE( mskisf_cav(jpi,jpj) , STAT=ialloc )
      ierr = ierr + ialloc
      !
      CALL mpp_sum ( 'isf', ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'isf: failed to allocate arrays.' )
      !
      ! initalisation of fwf and tsc array to 0
      risfload    (:,:)   = 0._wp
      fwfisf_oasis(:,:)   = 0._wp
      fwfisf_par  (:,:)   = 0._wp   ;   fwfisf_par_b  (:,:)   = 0._wp
      fwfisf_cav  (:,:)   = 0._wp   ;   fwfisf_cav_b  (:,:)   = 0._wp
      risf_cav_tsc(:,:,:) = 0._wp   ;   risf_cav_tsc_b(:,:,:) = 0._wp
      risf_par_tsc(:,:,:) = 0._wp   ;   risf_par_tsc_b(:,:,:) = 0._wp
      !
   END SUBROUTINE isf_alloc
   
   !!======================================================================
END MODULE isf_oce
