MODULE bdyini
   !!======================================================================
   !!                       ***  MODULE  bdyini  ***
   !! Unstructured open boundaries : initialisation
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2007-01  (D. Storkey) Update to use IOM module
   !!             -   !  2007-01  (D. Storkey) Tidal forcing
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (E.O'Dea) updates for Shelf configurations
   !!            3.3  !  2010-09  (D.Storkey) add ice boundary conditions
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!            3.4  !  2012     (J. Chanut) straight open boundary case update
   !!            3.5  !  2012     (S. Mocavero, I. Epicoco) optimization of BDY communications
   !!            3.7  !  2016     (T. Lovato) Remove bdy macro, call here init for dta and tides
   !!----------------------------------------------------------------------
   !!   bdy_init      : Initialization of unstructured open boundaries
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain
   USE bdy_oce        ! unstructured open boundary conditions
   USE bdydta         ! open boundary cond. setting   (bdy_dta_init routine)
   USE bdytides       ! open boundary cond. setting   (bdytide_init routine)
   USE sbctide        ! Tidal forcing or not
   USE phycst   , ONLY: rday
   !
   USE in_out_manager ! I/O units
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! for mpp_sum  
   USE iom            ! I/O

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_init   ! routine called in nemo_init

   INTEGER, PARAMETER ::   jp_nseg = 100   ! 
   INTEGER, PARAMETER ::   nrimmax =  20   ! maximum rimwidth in structured
                                               ! open boundary data files
   ! Straight open boundary segment parameters:
   INTEGER  ::   nbdysege, nbdysegw, nbdysegn, nbdysegs 
   INTEGER, DIMENSION(jp_nseg) ::   jpieob, jpjedt, jpjeft, npckge   !
   INTEGER, DIMENSION(jp_nseg) ::   jpiwob, jpjwdt, jpjwft, npckgw   !
   INTEGER, DIMENSION(jp_nseg) ::   jpjnob, jpindt, jpinft, npckgn   !
   INTEGER, DIMENSION(jp_nseg) ::   jpjsob, jpisdt, jpisft, npckgs   !
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: bdyini.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE bdy_init  ***
      !!
      !! ** Purpose :   Initialization of the dynamics and tracer fields with
      !!              unstructured open boundaries.
      !!
      !! ** Method  :   Read initialization arrays (mask, indices) to identify
      !!              an unstructured open boundary
      !!
      !! ** Input   :  bdy_init.nc, input file for unstructured open boundaries
      !!----------------------------------------------------------------------
      NAMELIST/nambdy/ ln_bdy, nb_bdy, ln_coords_file, cn_coords_file,         &
         &             ln_mask_file, cn_mask_file, cn_dyn2d, nn_dyn2d_dta,     &
         &             cn_dyn3d, nn_dyn3d_dta, cn_tra, nn_tra_dta,             &
         &             ln_tra_dmp, ln_dyn3d_dmp, rn_time_dmp, rn_time_dmp_out, &
         &             cn_ice, nn_ice_dta,                                     &
         &             rn_ice_tem, rn_ice_sal, rn_ice_age,                     &
         &             ln_vol, nn_volctl, nn_rimwidth, nb_jpk_bdy
         !
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------

      ! ------------------------
      ! Read namelist parameters
      ! ------------------------
      REWIND( numnam_ref )              ! Namelist nambdy in reference namelist :Unstructured open boundaries
      READ  ( numnam_ref, nambdy, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nambdy in reference namelist', lwp )
      REWIND( numnam_cfg )              ! Namelist nambdy in configuration namelist :Unstructured open boundaries
      READ  ( numnam_cfg, nambdy, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nambdy in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nambdy )

      IF( .NOT. Agrif_Root() ) ln_bdy = .FALSE.   ! forced for Agrif children
      
      ! -----------------------------------------
      ! unstructured open boundaries use control
      ! -----------------------------------------
      IF ( ln_bdy ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'bdy_init : initialization of open boundaries'
         IF(lwp) WRITE(numout,*) '~~~~~~~~'
         !
         ! Open boundaries definition (arrays and masks)
         CALL bdy_segs
         !
         ! Open boundaries initialisation of external data arrays
         CALL bdy_dta_init
         !
         ! Open boundaries initialisation of tidal harmonic forcing
         IF( ln_tide ) CALL bdytide_init
         !
      ELSE
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'bdy_init : open boundaries not used (ln_bdy = F)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~'
         !
      ENDIF
      !
   END SUBROUTINE bdy_init


   SUBROUTINE bdy_segs
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE bdy_init  ***
      !!         
      !! ** Purpose :   Definition of unstructured open boundaries.
      !!
      !! ** Method  :   Read initialization arrays (mask, indices) to identify 
      !!              an unstructured open boundary
      !!
      !! ** Input   :  bdy_init.nc, input file for unstructured open boundaries
      !!----------------------------------------------------------------------      
      INTEGER  ::   ib_bdy, ii, ij, ik, igrd, ib, ir, iseg ! dummy loop indices
      INTEGER  ::   icount, icountr, ibr_max, ilen1, ibm1  ! local integers
      INTEGER  ::   iwe, ies, iso, ino, inum, id_dummy     !   -       -
      INTEGER  ::   igrd_start, igrd_end, jpbdta           !   -       -
      INTEGER  ::   jpbdtau, jpbdtas                       !   -       -
      INTEGER  ::   ib_bdy1, ib_bdy2, ib1, ib2             !   -       -
      INTEGER  ::   i_offset, j_offset                     !   -       -
      INTEGER , POINTER  ::  nbi, nbj, nbr                 ! short cuts
      REAL(wp), POINTER  ::  flagu, flagv                  !    -   -
      REAL(wp), POINTER, DIMENSION(:,:)       ::   pmask    ! pointer to 2D mask fields
      REAL(wp) ::   zefl, zwfl, znfl, zsfl                 ! local scalars
      INTEGER, DIMENSION (2)                  ::   kdimsz
      INTEGER, DIMENSION(jpbgrd,jp_bdy)       ::   nblendta         ! Length of index arrays 
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:)  ::   nbidta, nbjdta   ! Index arrays: i and j indices of bdy dta
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:)  ::   nbrdta           ! Discrete distance from rim points
      CHARACTER(LEN=1),DIMENSION(jpbgrd)      ::   cgrid
      INTEGER :: com_east, com_west, com_south, com_north          ! Flags for boundaries sending
      INTEGER :: com_east_b, com_west_b, com_south_b, com_north_b  ! Flags for boundaries receiving
      INTEGER :: iw_b(4), ie_b(4), is_b(4), in_b(4)                ! Arrays for neighbours coordinates
      REAL(wp), TARGET, DIMENSION(jpi,jpj) ::   zfmask  ! temporary fmask array excluding coastal boundary condition (shlat)
      !!
      CHARACTER(LEN=1)                     ::   ctypebdy   !     -        - 
      INTEGER                              ::   nbdyind, nbdybeg, nbdyend
      !!
      NAMELIST/nambdy_index/ ctypebdy, nbdyind, nbdybeg, nbdyend
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------
      !
      cgrid = (/'t','u','v'/)

      ! -----------------------------------------
      ! Check and write out namelist parameters
      ! -----------------------------------------
      IF( jperio /= 0 )   CALL ctl_stop( 'bdy_segs: Cyclic or symmetric,',   &
         &                               ' and general open boundary condition are not compatible' )

      IF( nb_bdy == 0 ) THEN 
        IF(lwp) WRITE(numout,*) 'nb_bdy = 0, NO OPEN BOUNDARIES APPLIED.'
      ELSE
        IF(lwp) WRITE(numout,*) 'Number of open boundary sets : ', nb_bdy
      ENDIF

      DO ib_bdy = 1,nb_bdy
        IF(lwp) WRITE(numout,*) ' ' 
        IF(lwp) WRITE(numout,*) '------ Open boundary data set ',ib_bdy,'------' 

        IF( ln_coords_file(ib_bdy) ) THEN
           IF(lwp) WRITE(numout,*) 'Boundary definition read from file '//TRIM(cn_coords_file(ib_bdy))
        ELSE
           IF(lwp) WRITE(numout,*) 'Boundary defined in namelist.'
        ENDIF
        IF(lwp) WRITE(numout,*)

        IF(lwp) WRITE(numout,*) 'Boundary conditions for barotropic solution:  '
        SELECT CASE( cn_dyn2d(ib_bdy) )                  
          CASE( 'none' )         
             IF(lwp) WRITE(numout,*) '      no open boundary condition'        
             dta_bdy(ib_bdy)%ll_ssh = .false.
             dta_bdy(ib_bdy)%ll_u2d = .false.
             dta_bdy(ib_bdy)%ll_v2d = .false.
          CASE( 'frs' )          
             IF(lwp) WRITE(numout,*) '      Flow Relaxation Scheme'
             dta_bdy(ib_bdy)%ll_ssh = .false.
             dta_bdy(ib_bdy)%ll_u2d = .true.
             dta_bdy(ib_bdy)%ll_v2d = .true.
          CASE( 'flather' )      
             IF(lwp) WRITE(numout,*) '      Flather radiation condition'
             dta_bdy(ib_bdy)%ll_ssh = .true.
             dta_bdy(ib_bdy)%ll_u2d = .true.
             dta_bdy(ib_bdy)%ll_v2d = .true.
          CASE( 'orlanski' )     
             IF(lwp) WRITE(numout,*) '      Orlanski (fully oblique) radiation condition with adaptive nudging'
             dta_bdy(ib_bdy)%ll_ssh = .false.
             dta_bdy(ib_bdy)%ll_u2d = .true.
             dta_bdy(ib_bdy)%ll_v2d = .true.
          CASE( 'orlanski_npo' ) 
             IF(lwp) WRITE(numout,*) '      Orlanski (NPO) radiation condition with adaptive nudging'
             dta_bdy(ib_bdy)%ll_ssh = .false.
             dta_bdy(ib_bdy)%ll_u2d = .true.
             dta_bdy(ib_bdy)%ll_v2d = .true.
          CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for cn_dyn2d' )
        END SELECT
        IF( cn_dyn2d(ib_bdy) /= 'none' ) THEN
           SELECT CASE( nn_dyn2d_dta(ib_bdy) )                   ! 
              CASE( 0 )      ;   IF(lwp) WRITE(numout,*) '      initial state used for bdy data'        
              CASE( 1 )      ;   IF(lwp) WRITE(numout,*) '      boundary data taken from file'
              CASE( 2 )      ;   IF(lwp) WRITE(numout,*) '      tidal harmonic forcing taken from file'
              CASE( 3 )      ;   IF(lwp) WRITE(numout,*) '      boundary data AND tidal harmonic forcing taken from files'
              CASE DEFAULT   ;   CALL ctl_stop( 'nn_dyn2d_dta must be between 0 and 3' )
           END SELECT
           IF (( nn_dyn2d_dta(ib_bdy) .ge. 2 ).AND.(.NOT.ln_tide)) THEN
             CALL ctl_stop( 'You must activate with ln_tide to add tidal forcing at open boundaries' )
           ENDIF
        ENDIF
        IF(lwp) WRITE(numout,*)

        IF(lwp) WRITE(numout,*) 'Boundary conditions for baroclinic velocities:  '
        SELECT CASE( cn_dyn3d(ib_bdy) )                  
          CASE('none')
             IF(lwp) WRITE(numout,*) '      no open boundary condition'        
             dta_bdy(ib_bdy)%ll_u3d = .false.
             dta_bdy(ib_bdy)%ll_v3d = .false.
          CASE('frs')       
             IF(lwp) WRITE(numout,*) '      Flow Relaxation Scheme'
             dta_bdy(ib_bdy)%ll_u3d = .true.
             dta_bdy(ib_bdy)%ll_v3d = .true.
          CASE('specified')
             IF(lwp) WRITE(numout,*) '      Specified value'
             dta_bdy(ib_bdy)%ll_u3d = .true.
             dta_bdy(ib_bdy)%ll_v3d = .true.
          CASE('neumann')
             IF(lwp) WRITE(numout,*) '      Neumann conditions'
             dta_bdy(ib_bdy)%ll_u3d = .false.
             dta_bdy(ib_bdy)%ll_v3d = .false.
          CASE('zerograd')
             IF(lwp) WRITE(numout,*) '      Zero gradient for baroclinic velocities'
             dta_bdy(ib_bdy)%ll_u3d = .false.
             dta_bdy(ib_bdy)%ll_v3d = .false.
          CASE('zero')
             IF(lwp) WRITE(numout,*) '      Zero baroclinic velocities (runoff case)'
             dta_bdy(ib_bdy)%ll_u3d = .false.
             dta_bdy(ib_bdy)%ll_v3d = .false.
          CASE('orlanski')
             IF(lwp) WRITE(numout,*) '      Orlanski (fully oblique) radiation condition with adaptive nudging'
             dta_bdy(ib_bdy)%ll_u3d = .true.
             dta_bdy(ib_bdy)%ll_v3d = .true.
          CASE('orlanski_npo')
             IF(lwp) WRITE(numout,*) '      Orlanski (NPO) radiation condition with adaptive nudging'
             dta_bdy(ib_bdy)%ll_u3d = .true.
             dta_bdy(ib_bdy)%ll_v3d = .true.
          CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for cn_dyn3d' )
        END SELECT
        IF( cn_dyn3d(ib_bdy) /= 'none' ) THEN
           SELECT CASE( nn_dyn3d_dta(ib_bdy) )                   ! 
              CASE( 0 )      ;   IF(lwp) WRITE(numout,*) '      initial state used for bdy data'        
              CASE( 1 )      ;   IF(lwp) WRITE(numout,*) '      boundary data taken from file'
              CASE DEFAULT   ;   CALL ctl_stop( 'nn_dyn3d_dta must be 0 or 1' )
           END SELECT
        ENDIF

        IF ( ln_dyn3d_dmp(ib_bdy) ) THEN
           IF ( cn_dyn3d(ib_bdy) == 'none' ) THEN
              IF(lwp) WRITE(numout,*) 'No open boundary condition for baroclinic velocities: ln_dyn3d_dmp is set to .false.'
              ln_dyn3d_dmp(ib_bdy)=.false.
           ELSEIF ( cn_dyn3d(ib_bdy) == 'frs' ) THEN
              CALL ctl_stop( 'Use FRS OR relaxation' )
           ELSE
              IF(lwp) WRITE(numout,*) '      + baroclinic velocities relaxation zone'
              IF(lwp) WRITE(numout,*) '      Damping time scale: ',rn_time_dmp(ib_bdy),' days'
              IF((lwp).AND.rn_time_dmp(ib_bdy)<0) CALL ctl_stop( 'Time scale must be positive' )
              dta_bdy(ib_bdy)%ll_u3d = .true.
              dta_bdy(ib_bdy)%ll_v3d = .true.
           ENDIF
        ELSE
           IF(lwp) WRITE(numout,*) '      NO relaxation on baroclinic velocities'
        ENDIF
        IF(lwp) WRITE(numout,*)

        IF(lwp) WRITE(numout,*) 'Boundary conditions for temperature and salinity:  '
        SELECT CASE( cn_tra(ib_bdy) )                  
          CASE('none')
             IF(lwp) WRITE(numout,*) '      no open boundary condition'        
             dta_bdy(ib_bdy)%ll_tem = .false.
             dta_bdy(ib_bdy)%ll_sal = .false.
          CASE('frs')
             IF(lwp) WRITE(numout,*) '      Flow Relaxation Scheme'
             dta_bdy(ib_bdy)%ll_tem = .true.
             dta_bdy(ib_bdy)%ll_sal = .true.
          CASE('specified')
             IF(lwp) WRITE(numout,*) '      Specified value'
             dta_bdy(ib_bdy)%ll_tem = .true.
             dta_bdy(ib_bdy)%ll_sal = .true.
          CASE('neumann')
             IF(lwp) WRITE(numout,*) '      Neumann conditions'
             dta_bdy(ib_bdy)%ll_tem = .false.
             dta_bdy(ib_bdy)%ll_sal = .false.
          CASE('runoff')
             IF(lwp) WRITE(numout,*) '      Runoff conditions : Neumann for T and specified to 0.1 for salinity'
             dta_bdy(ib_bdy)%ll_tem = .false.
             dta_bdy(ib_bdy)%ll_sal = .false.
          CASE('orlanski')
             IF(lwp) WRITE(numout,*) '      Orlanski (fully oblique) radiation condition with adaptive nudging'
             dta_bdy(ib_bdy)%ll_tem = .true.
             dta_bdy(ib_bdy)%ll_sal = .true.
          CASE('orlanski_npo')
             IF(lwp) WRITE(numout,*) '      Orlanski (NPO) radiation condition with adaptive nudging'
             dta_bdy(ib_bdy)%ll_tem = .true.
             dta_bdy(ib_bdy)%ll_sal = .true.
          CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for cn_tra' )
        END SELECT
        IF( cn_tra(ib_bdy) /= 'none' ) THEN
           SELECT CASE( nn_tra_dta(ib_bdy) )                   ! 
              CASE( 0 )      ;   IF(lwp) WRITE(numout,*) '      initial state used for bdy data'        
              CASE( 1 )      ;   IF(lwp) WRITE(numout,*) '      boundary data taken from file'
              CASE DEFAULT   ;   CALL ctl_stop( 'nn_tra_dta must be 0 or 1' )
           END SELECT
        ENDIF

        IF ( ln_tra_dmp(ib_bdy) ) THEN
           IF ( cn_tra(ib_bdy) == 'none' ) THEN
              IF(lwp) WRITE(numout,*) 'No open boundary condition for tracers: ln_tra_dmp is set to .false.'
              ln_tra_dmp(ib_bdy)=.false.
           ELSEIF ( cn_tra(ib_bdy) == 'frs' ) THEN
              CALL ctl_stop( 'Use FRS OR relaxation' )
           ELSE
              IF(lwp) WRITE(numout,*) '      + T/S relaxation zone'
              IF(lwp) WRITE(numout,*) '      Damping time scale: ',rn_time_dmp(ib_bdy),' days'
              IF(lwp) WRITE(numout,*) '      Outflow damping time scale: ',rn_time_dmp_out(ib_bdy),' days'
              IF((lwp).AND.rn_time_dmp(ib_bdy)<0) CALL ctl_stop( 'Time scale must be positive' )
              dta_bdy(ib_bdy)%ll_tem = .true.
              dta_bdy(ib_bdy)%ll_sal = .true.
           ENDIF
        ELSE
           IF(lwp) WRITE(numout,*) '      NO T/S relaxation'
        ENDIF
        IF(lwp) WRITE(numout,*)

#if defined key_si3
         IF(lwp) WRITE(numout,*) 'Boundary conditions for sea ice:  '
         SELECT CASE( cn_ice(ib_bdy) )                  
         CASE('none')
             IF(lwp) WRITE(numout,*) '      no open boundary condition'        
             dta_bdy(ib_bdy)%ll_a_i = .false.
             dta_bdy(ib_bdy)%ll_h_i = .false.
             dta_bdy(ib_bdy)%ll_h_s = .false.
         CASE('frs')
             IF(lwp) WRITE(numout,*) '      Flow Relaxation Scheme'
             dta_bdy(ib_bdy)%ll_a_i = .true.
             dta_bdy(ib_bdy)%ll_h_i = .true.
             dta_bdy(ib_bdy)%ll_h_s = .true.
         CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for cn_ice' )
         END SELECT
        IF( cn_ice(ib_bdy) /= 'none' ) THEN 
           SELECT CASE( nn_ice_dta(ib_bdy) )                   ! 
              CASE( 0 )      ;   IF(lwp) WRITE(numout,*) '      initial state used for bdy data'        
              CASE( 1 )      ;   IF(lwp) WRITE(numout,*) '      boundary data taken from file'
              CASE DEFAULT   ;   CALL ctl_stop( 'nn_ice_dta must be 0 or 1' )
           END SELECT
        ENDIF
        IF(lwp) WRITE(numout,*)
        IF(lwp) WRITE(numout,*) '      tem of bdy sea-ice = ', rn_ice_tem(ib_bdy)         
        IF(lwp) WRITE(numout,*) '      sal of bdy sea-ice = ', rn_ice_sal(ib_bdy)         
        IF(lwp) WRITE(numout,*) '      age of bdy sea-ice = ', rn_ice_age(ib_bdy)         
#endif

        IF(lwp) WRITE(numout,*) '      Width of relaxation zone = ', nn_rimwidth(ib_bdy)
        IF(lwp) WRITE(numout,*)
         !
      END DO

     IF( nb_bdy > 0 ) THEN
        IF( ln_vol ) THEN                     ! check volume conservation (nn_volctl value)
          IF(lwp) WRITE(numout,*) 'Volume correction applied at open boundaries'
          IF(lwp) WRITE(numout,*)
          SELECT CASE ( nn_volctl )
            CASE( 1 )      ;   IF(lwp) WRITE(numout,*) '      The total volume will be constant'
            CASE( 0 )      ;   IF(lwp) WRITE(numout,*) '      The total volume will vary according to the surface E-P flux'
            CASE DEFAULT   ;   CALL ctl_stop( 'nn_volctl must be 0 or 1' )
          END SELECT
          IF(lwp) WRITE(numout,*)
        ELSE
          IF(lwp) WRITE(numout,*) 'No volume correction applied at open boundaries'
          IF(lwp) WRITE(numout,*)
        ENDIF
        IF( nb_jpk_bdy > 0 ) THEN
           IF(lwp) WRITE(numout,*) '*** open boundary will be interpolate in the vertical onto the native grid ***'
        ELSE
           IF(lwp) WRITE(numout,*) '*** open boundary will be read straight onto the native grid without vertical interpolation ***'
        ENDIF
     ENDIF

      ! -------------------------------------------------
      ! Initialise indices arrays for open boundaries
      ! -------------------------------------------------

      ! Work out global dimensions of boundary data
      ! ---------------------------------------------
      REWIND( numnam_cfg )     

      nblendta(:,:) = 0
      nbdysege = 0
      nbdysegw = 0
      nbdysegn = 0
      nbdysegs = 0
      icount   = 0 ! count user defined segments
      ! Dimensions below are used to allocate arrays to read external data
      jpbdtas = 1 ! Maximum size of boundary data (structured case)
      jpbdtau = 1 ! Maximum size of boundary data (unstructured case)

      DO ib_bdy = 1, nb_bdy

         IF( .NOT. ln_coords_file(ib_bdy) ) THEN ! Work out size of global arrays from namelist parameters
 
            icount = icount + 1
            ! No REWIND here because may need to read more than one nambdy_index namelist.
            ! Read only namelist_cfg to avoid unseccessfull overwrite 
            ! keep full control of the configuration namelist
            READ  ( numnam_cfg, nambdy_index, IOSTAT = ios, ERR = 904 )
904         IF( ios /= 0 )   CALL ctl_nam ( ios , 'nambdy_index in configuration namelist', lwp )
            IF(lwm) WRITE ( numond, nambdy_index )

            SELECT CASE ( TRIM(ctypebdy) )
              CASE( 'N' )
                 IF( nbdyind == -1 ) THEN  ! Automatic boundary definition: if nbdysegX = -1
                    nbdyind  = jpjglo - 2  ! set boundary to whole side of model domain.
                    nbdybeg  = 2
                    nbdyend  = jpiglo - 1
                 ENDIF
                 nbdysegn = nbdysegn + 1
                 npckgn(nbdysegn) = ib_bdy ! Save bdy package number
                 jpjnob(nbdysegn) = nbdyind
                 jpindt(nbdysegn) = nbdybeg
                 jpinft(nbdysegn) = nbdyend
                 !
              CASE( 'S' )
                 IF( nbdyind == -1 ) THEN  ! Automatic boundary definition: if nbdysegX = -1
                    nbdyind  = 2           ! set boundary to whole side of model domain.
                    nbdybeg  = 2
                    nbdyend  = jpiglo - 1
                 ENDIF
                 nbdysegs = nbdysegs + 1
                 npckgs(nbdysegs) = ib_bdy ! Save bdy package number
                 jpjsob(nbdysegs) = nbdyind
                 jpisdt(nbdysegs) = nbdybeg
                 jpisft(nbdysegs) = nbdyend
                 !
              CASE( 'E' )
                 IF( nbdyind == -1 ) THEN  ! Automatic boundary definition: if nbdysegX = -1
                    nbdyind  = jpiglo - 2  ! set boundary to whole side of model domain.
                    nbdybeg  = 2
                    nbdyend  = jpjglo - 1
                 ENDIF
                 nbdysege = nbdysege + 1 
                 npckge(nbdysege) = ib_bdy ! Save bdy package number
                 jpieob(nbdysege) = nbdyind
                 jpjedt(nbdysege) = nbdybeg
                 jpjeft(nbdysege) = nbdyend
                 !
              CASE( 'W' )
                 IF( nbdyind == -1 ) THEN  ! Automatic boundary definition: if nbdysegX = -1
                    nbdyind  = 2           ! set boundary to whole side of model domain.
                    nbdybeg  = 2
                    nbdyend  = jpjglo - 1
                 ENDIF
                 nbdysegw = nbdysegw + 1
                 npckgw(nbdysegw) = ib_bdy ! Save bdy package number
                 jpiwob(nbdysegw) = nbdyind
                 jpjwdt(nbdysegw) = nbdybeg
                 jpjwft(nbdysegw) = nbdyend
                 !
              CASE DEFAULT   ;   CALL ctl_stop( 'ctypebdy must be N, S, E or W' )
            END SELECT

            ! For simplicity we assume that in case of straight bdy, arrays have the same length
            ! (even if it is true that last tangential velocity points
            ! are useless). This simplifies a little bit boundary data format (and agrees with format
            ! used so far in obc package)

            nblendta(1:jpbgrd,ib_bdy) =  (nbdyend - nbdybeg + 1) * nn_rimwidth(ib_bdy)
            jpbdtas = MAX(jpbdtas, (nbdyend - nbdybeg + 1))
            IF (lwp.and.(nn_rimwidth(ib_bdy)>nrimmax)) &
            & CALL ctl_stop( 'rimwidth must be lower than nrimmax' )

         ELSE            ! Read size of arrays in boundary coordinates file.
            CALL iom_open( cn_coords_file(ib_bdy), inum )
            DO igrd = 1, jpbgrd
               id_dummy = iom_varid( inum, 'nbi'//cgrid(igrd), kdimsz=kdimsz )  
               nblendta(igrd,ib_bdy) = MAXVAL(kdimsz)
               jpbdtau = MAX(jpbdtau, MAXVAL(kdimsz))
            END DO
            CALL iom_close( inum )
            !
         ENDIF 
         !
      END DO ! ib_bdy

      IF (nb_bdy>0) THEN
         jpbdta = MAXVAL(nblendta(1:jpbgrd,1:nb_bdy))

         ! Allocate arrays
         !---------------
         ALLOCATE( nbidta(jpbdta, jpbgrd, nb_bdy), nbjdta(jpbdta, jpbgrd, nb_bdy),    &
            &      nbrdta(jpbdta, jpbgrd, nb_bdy) )

         IF( nb_jpk_bdy>0 ) THEN
            ALLOCATE( dta_global(jpbdtau, 1, nb_jpk_bdy) )
            ALLOCATE( dta_global_z(jpbdtau, 1, nb_jpk_bdy) )
            ALLOCATE( dta_global_dz(jpbdtau, 1, nb_jpk_bdy) )
         ELSE
            ALLOCATE( dta_global(jpbdtau, 1, jpk) )
            ALLOCATE( dta_global_z(jpbdtau, 1, jpk) ) ! needed ?? TODO
            ALLOCATE( dta_global_dz(jpbdtau, 1, jpk) )! needed ?? TODO
         ENDIF

         IF ( icount>0 ) THEN
            IF( nb_jpk_bdy>0 ) THEN
               ALLOCATE( dta_global2(jpbdtas, nrimmax, nb_jpk_bdy) )
               ALLOCATE( dta_global2_z(jpbdtas, nrimmax, nb_jpk_bdy) )
               ALLOCATE( dta_global2_dz(jpbdtas, nrimmax, nb_jpk_bdy) )
            ELSE
               ALLOCATE( dta_global2(jpbdtas, nrimmax, jpk) )
               ALLOCATE( dta_global2_z(jpbdtas, nrimmax, jpk) ) ! needed ?? TODO
               ALLOCATE( dta_global2_dz(jpbdtas, nrimmax, jpk) )! needed ?? TODO  
            ENDIF
         ENDIF
         ! 
      ENDIF

      ! Now look for crossings in user (namelist) defined open boundary segments:
      !--------------------------------------------------------------------------
      IF( icount>0 )   CALL bdy_ctl_seg

      ! Calculate global boundary index arrays or read in from file
      !------------------------------------------------------------               
      ! 1. Read global index arrays from boundary coordinates file.
      DO ib_bdy = 1, nb_bdy
         !
         IF( ln_coords_file(ib_bdy) ) THEN
            !
            CALL iom_open( cn_coords_file(ib_bdy), inum )
            DO igrd = 1, jpbgrd
               CALL iom_get( inum, jpdom_unknown, 'nbi'//cgrid(igrd), dta_global(1:nblendta(igrd,ib_bdy),:,1) )
               DO ii = 1,nblendta(igrd,ib_bdy)
                  nbidta(ii,igrd,ib_bdy) = INT( dta_global(ii,1,1) )
               END DO
               CALL iom_get( inum, jpdom_unknown, 'nbj'//cgrid(igrd), dta_global(1:nblendta(igrd,ib_bdy),:,1) )
               DO ii = 1,nblendta(igrd,ib_bdy)
                  nbjdta(ii,igrd,ib_bdy) = INT( dta_global(ii,1,1) )
               END DO
               CALL iom_get( inum, jpdom_unknown, 'nbr'//cgrid(igrd), dta_global(1:nblendta(igrd,ib_bdy),:,1) )
               DO ii = 1,nblendta(igrd,ib_bdy)
                  nbrdta(ii,igrd,ib_bdy) = INT( dta_global(ii,1,1) )
               END DO
               !
               ibr_max = MAXVAL( nbrdta(:,igrd,ib_bdy) )
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) ' Maximum rimwidth in file is ', ibr_max
               IF(lwp) WRITE(numout,*) ' nn_rimwidth from namelist is ', nn_rimwidth(ib_bdy)
               IF (ibr_max < nn_rimwidth(ib_bdy))   &
                     CALL ctl_stop( 'nn_rimwidth is larger than maximum rimwidth in file',cn_coords_file(ib_bdy) )
            END DO
            CALL iom_close( inum )
            !
         ENDIF 
         !
      END DO      
    
      ! 2. Now fill indices corresponding to straight open boundary arrays:
      ! East
      !-----
      DO iseg = 1, nbdysege
         ib_bdy = npckge(iseg)
         !
         ! ------------ T points -------------
         igrd=1
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ij = jpjedt(iseg), jpjeft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = jpieob(iseg) + 2 - ir
               nbjdta(icount, igrd, ib_bdy) = ij
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
         !
         ! ------------ U points -------------
         igrd=2
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ij = jpjedt(iseg), jpjeft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = jpieob(iseg) + 1 - ir
               nbjdta(icount, igrd, ib_bdy) = ij
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
         !
         ! ------------ V points -------------
         igrd=3
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
!            DO ij = jpjedt(iseg), jpjeft(iseg) - 1
            DO ij = jpjedt(iseg), jpjeft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = jpieob(iseg) + 2 - ir
               nbjdta(icount, igrd, ib_bdy) = ij
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
            nbidta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
            nbjdta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
         ENDDO
      ENDDO
      !
      ! West
      !-----
      DO iseg = 1, nbdysegw
         ib_bdy = npckgw(iseg)
         !
         ! ------------ T points -------------
         igrd=1
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ij = jpjwdt(iseg), jpjwft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = jpiwob(iseg) + ir - 1
               nbjdta(icount, igrd, ib_bdy) = ij
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
         !
         ! ------------ U points -------------
         igrd=2
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ij = jpjwdt(iseg), jpjwft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = jpiwob(iseg) + ir - 1
               nbjdta(icount, igrd, ib_bdy) = ij
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
         !
         ! ------------ V points -------------
         igrd=3
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
!            DO ij = jpjwdt(iseg), jpjwft(iseg) - 1
            DO ij = jpjwdt(iseg), jpjwft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = jpiwob(iseg) + ir - 1
               nbjdta(icount, igrd, ib_bdy) = ij
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
            nbidta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
            nbjdta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
         ENDDO
      ENDDO
      !
      ! North
      !-----
      DO iseg = 1, nbdysegn
         ib_bdy = npckgn(iseg)
         !
         ! ------------ T points -------------
         igrd=1
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ii = jpindt(iseg), jpinft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = ii
               nbjdta(icount, igrd, ib_bdy) = jpjnob(iseg) + 2 - ir 
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
         !
         ! ------------ U points -------------
         igrd=2
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
!            DO ii = jpindt(iseg), jpinft(iseg) - 1
            DO ii = jpindt(iseg), jpinft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = ii
               nbjdta(icount, igrd, ib_bdy) = jpjnob(iseg) + 2 - ir
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
            nbidta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
            nbjdta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
         ENDDO
         !
         ! ------------ V points -------------
         igrd=3
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ii = jpindt(iseg), jpinft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = ii
               nbjdta(icount, igrd, ib_bdy) = jpjnob(iseg) + 1 - ir
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
      ENDDO
      !
      ! South
      !-----
      DO iseg = 1, nbdysegs
         ib_bdy = npckgs(iseg)
         !
         ! ------------ T points -------------
         igrd=1
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ii = jpisdt(iseg), jpisft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = ii
               nbjdta(icount, igrd, ib_bdy) = jpjsob(iseg) + ir - 1
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
         !
         ! ------------ U points -------------
         igrd=2
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
!            DO ii = jpisdt(iseg), jpisft(iseg) - 1
            DO ii = jpisdt(iseg), jpisft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = ii
               nbjdta(icount, igrd, ib_bdy) = jpjsob(iseg) + ir - 1
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
            nbidta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
            nbjdta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
         ENDDO
         !
         ! ------------ V points -------------
         igrd=3
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ii = jpisdt(iseg), jpisft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = ii
               nbjdta(icount, igrd, ib_bdy) = jpjsob(iseg) + ir - 1
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
      ENDDO

      !  Deal with duplicated points
      !-----------------------------
      ! We assign negative indices to duplicated points (to remove them from bdy points to be updated)
      ! if their distance to the bdy is greater than the other
      ! If their distance are the same, just keep only one to avoid updating a point twice
      DO igrd = 1, jpbgrd
         DO ib_bdy1 = 1, nb_bdy
            DO ib_bdy2 = 1, nb_bdy
               IF (ib_bdy1/=ib_bdy2) THEN
                  DO ib1 = 1, nblendta(igrd,ib_bdy1)
                     DO ib2 = 1, nblendta(igrd,ib_bdy2)
                        IF ((nbidta(ib1, igrd, ib_bdy1)==nbidta(ib2, igrd, ib_bdy2)).AND. &
                        &   (nbjdta(ib1, igrd, ib_bdy1)==nbjdta(ib2, igrd, ib_bdy2))) THEN
!                           IF ((lwp).AND.(igrd==1)) WRITE(numout,*) ' found coincident point ji, jj:', & 
!                                                       &              nbidta(ib1, igrd, ib_bdy1),      & 
!                                                       &              nbjdta(ib2, igrd, ib_bdy2)
                           ! keep only points with the lowest distance to boundary:
                           IF (nbrdta(ib1, igrd, ib_bdy1)<nbrdta(ib2, igrd, ib_bdy2)) THEN
                             nbidta(ib2, igrd, ib_bdy2) =-ib_bdy2
                             nbjdta(ib2, igrd, ib_bdy2) =-ib_bdy2
                           ELSEIF (nbrdta(ib1, igrd, ib_bdy1)>nbrdta(ib2, igrd, ib_bdy2)) THEN
                             nbidta(ib1, igrd, ib_bdy1) =-ib_bdy1
                             nbjdta(ib1, igrd, ib_bdy1) =-ib_bdy1
                           ! Arbitrary choice if distances are the same:
                           ELSE
                             nbidta(ib1, igrd, ib_bdy1) =-ib_bdy1
                             nbjdta(ib1, igrd, ib_bdy1) =-ib_bdy1
                           ENDIF
                        END IF
                     END DO
                  END DO
               ENDIF
            END DO
         END DO
      END DO

      ! Work out dimensions of boundary data on each processor
      ! ------------------------------------------------------

      ! Rather assume that boundary data indices are given on global domain
      ! TO BE DISCUSSED ?
!      iw = mig(1) + 1            ! if monotasking and no zoom, iw=2
!      ie = mig(1) + nlci-1 - 1   ! if monotasking and no zoom, ie=jpim1
!      is = mjg(1) + 1            ! if monotasking and no zoom, is=2
!      in = mjg(1) + nlcj-1 - 1   ! if monotasking and no zoom, in=jpjm1      
      iwe = mig(1) - 1 + 2         ! if monotasking and no zoom, iw=2
      ies = mig(1) + nlci-1 - 1  ! if monotasking and no zoom, ie=jpim1
      iso = mjg(1) - 1 + 2         ! if monotasking and no zoom, is=2
      ino = mjg(1) + nlcj-1 - 1  ! if monotasking and no zoom, in=jpjm1

      ALLOCATE( nbondi_bdy(nb_bdy))
      ALLOCATE( nbondj_bdy(nb_bdy))
      nbondi_bdy(:)=2
      nbondj_bdy(:)=2
      ALLOCATE( nbondi_bdy_b(nb_bdy))
      ALLOCATE( nbondj_bdy_b(nb_bdy))
      nbondi_bdy_b(:)=2
      nbondj_bdy_b(:)=2

      ! Work out dimensions of boundary data on each neighbour process
      IF(nbondi == 0) THEN
         iw_b(1) = 1 + nimppt(nowe+1)
         ie_b(1) = 1 + nimppt(nowe+1)+nlcit(nowe+1)-3
         is_b(1) = 1 + njmppt(nowe+1)
         in_b(1) = 1 + njmppt(nowe+1)+nlcjt(nowe+1)-3

         iw_b(2) = 1 + nimppt(noea+1)
         ie_b(2) = 1 + nimppt(noea+1)+nlcit(noea+1)-3
         is_b(2) = 1 + njmppt(noea+1)
         in_b(2) = 1 + njmppt(noea+1)+nlcjt(noea+1)-3
      ELSEIF(nbondi == 1) THEN
         iw_b(1) = 1 + nimppt(nowe+1)
         ie_b(1) = 1 + nimppt(nowe+1)+nlcit(nowe+1)-3
         is_b(1) = 1 + njmppt(nowe+1)
         in_b(1) = 1 + njmppt(nowe+1)+nlcjt(nowe+1)-3
      ELSEIF(nbondi == -1) THEN
         iw_b(2) = 1 + nimppt(noea+1)
         ie_b(2) = 1 + nimppt(noea+1)+nlcit(noea+1)-3
         is_b(2) = 1 + njmppt(noea+1)
         in_b(2) = 1 + njmppt(noea+1)+nlcjt(noea+1)-3
      ENDIF

      IF(nbondj == 0) THEN
         iw_b(3) = 1 + nimppt(noso+1)
         ie_b(3) = 1 + nimppt(noso+1)+nlcit(noso+1)-3
         is_b(3) = 1 + njmppt(noso+1)
         in_b(3) = 1 + njmppt(noso+1)+nlcjt(noso+1)-3

         iw_b(4) = 1 + nimppt(nono+1)
         ie_b(4) = 1 + nimppt(nono+1)+nlcit(nono+1)-3
         is_b(4) = 1 + njmppt(nono+1)
         in_b(4) = 1 + njmppt(nono+1)+nlcjt(nono+1)-3
      ELSEIF(nbondj == 1) THEN
         iw_b(3) = 1 + nimppt(noso+1)
         ie_b(3) = 1 + nimppt(noso+1)+nlcit(noso+1)-3
         is_b(3) = 1 + njmppt(noso+1)
         in_b(3) = 1 + njmppt(noso+1)+nlcjt(noso+1)-3
      ELSEIF(nbondj == -1) THEN
         iw_b(4) = 1 + nimppt(nono+1)
         ie_b(4) = 1 + nimppt(nono+1)+nlcit(nono+1)-3
         is_b(4) = 1 + njmppt(nono+1)
         in_b(4) = 1 + njmppt(nono+1)+nlcjt(nono+1)-3
      ENDIF

      DO ib_bdy = 1, nb_bdy
         DO igrd = 1, jpbgrd
            icount  = 0
            icountr = 0
            idx_bdy(ib_bdy)%nblen(igrd)    = 0
            idx_bdy(ib_bdy)%nblenrim(igrd) = 0
            DO ib = 1, nblendta(igrd,ib_bdy)
               ! check that data is in correct order in file
               ibm1 = MAX(1,ib-1)
               IF(lwp) THEN         ! Since all procs read global data only need to do this check on one proc...
                  IF( nbrdta(ib,igrd,ib_bdy) < nbrdta(ibm1,igrd,ib_bdy) ) THEN
                     CALL ctl_stop('bdy_segs : ERROR : boundary data in file must be defined ', &
                          &        ' in order of distance from edge nbr A utility for re-ordering ', &
                          &        ' boundary coordinates and data files exists in the TOOLS/OBC directory')
                  ENDIF    
               ENDIF
               ! check if point is in local domain
               IF(  nbidta(ib,igrd,ib_bdy) >= iwe .AND. nbidta(ib,igrd,ib_bdy) <= ies .AND.   &
                  & nbjdta(ib,igrd,ib_bdy) >= iso .AND. nbjdta(ib,igrd,ib_bdy) <= ino      ) THEN
                  !
                  icount = icount  + 1
                  !
                  IF( nbrdta(ib,igrd,ib_bdy) == 1 )   icountr = icountr+1
               ENDIF
            END DO
            idx_bdy(ib_bdy)%nblenrim(igrd) = icountr !: length of rim boundary data on each proc
            idx_bdy(ib_bdy)%nblen   (igrd) = icount  !: length of boundary data on each proc        
         END DO  ! igrd

         ! Allocate index arrays for this boundary set
         !--------------------------------------------
         ilen1 = MAXVAL( idx_bdy(ib_bdy)%nblen(:) )
         ALLOCATE( idx_bdy(ib_bdy)%nbi   (ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%nbj   (ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%nbr   (ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%nbd   (ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%nbdout(ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%nbmap (ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%nbw   (ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%flagu (ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%flagv (ilen1,jpbgrd) )

         ! Dispatch mapping indices and discrete distances on each processor
         ! -----------------------------------------------------------------

         com_east  = 0
         com_west  = 0
         com_south = 0
         com_north = 0

         com_east_b  = 0
         com_west_b  = 0
         com_south_b = 0
         com_north_b = 0

         DO igrd = 1, jpbgrd
            icount  = 0
            ! Loop on rimwidth to ensure outermost points come first in the local arrays.
            DO ir=1, nn_rimwidth(ib_bdy)
               DO ib = 1, nblendta(igrd,ib_bdy)
                  ! check if point is in local domain and equals ir
                  IF(  nbidta(ib,igrd,ib_bdy) >= iwe .AND. nbidta(ib,igrd,ib_bdy) <= ies .AND.   &
                     & nbjdta(ib,igrd,ib_bdy) >= iso .AND. nbjdta(ib,igrd,ib_bdy) <= ino .AND.   &
                     & nbrdta(ib,igrd,ib_bdy) == ir  ) THEN
                     !
                     icount = icount  + 1

                     ! Rather assume that boundary data indices are given on global domain
                     ! TO BE DISCUSSED ?
!                     idx_bdy(ib_bdy)%nbi(icount,igrd)   = nbidta(ib,igrd,ib_bdy)- mig(1)+1
!                     idx_bdy(ib_bdy)%nbj(icount,igrd)   = nbjdta(ib,igrd,ib_bdy)- mjg(1)+1
                     idx_bdy(ib_bdy)%nbi(icount,igrd)   = nbidta(ib,igrd,ib_bdy)- mig(1)+1
                     idx_bdy(ib_bdy)%nbj(icount,igrd)   = nbjdta(ib,igrd,ib_bdy)- mjg(1)+1
                     ! check if point has to be sent
                     ii = idx_bdy(ib_bdy)%nbi(icount,igrd)
                     ij = idx_bdy(ib_bdy)%nbj(icount,igrd)
                     if((com_east .ne. 1) .and. (ii == (nlci-1)) .and. (nbondi .le. 0)) then
                        com_east = 1
                     elseif((com_west .ne. 1) .and. (ii == 2) .and. (nbondi .ge. 0) .and. (nbondi .ne. 2)) then
                        com_west = 1
                     endif 
                     if((com_south .ne. 1) .and. (ij == 2) .and. (nbondj .ge. 0) .and. (nbondj .ne. 2)) then
                        com_south = 1
                     elseif((com_north .ne. 1) .and. (ij == (nlcj-1)) .and. (nbondj .le. 0)) then
                        com_north = 1
                     endif 
                     idx_bdy(ib_bdy)%nbr(icount,igrd)   = nbrdta(ib,igrd,ib_bdy)
                     idx_bdy(ib_bdy)%nbmap(icount,igrd) = ib
                  ENDIF
                  ! check if point has to be received from a neighbour
                  IF(nbondi == 0) THEN
                     IF( nbidta(ib,igrd,ib_bdy) >= iw_b(1) .AND. nbidta(ib,igrd,ib_bdy) <= ie_b(1) .AND.   &
                       & nbjdta(ib,igrd,ib_bdy) >= is_b(1) .AND. nbjdta(ib,igrd,ib_bdy) <= in_b(1) .AND.   &
                       & nbrdta(ib,igrd,ib_bdy) == ir  ) THEN
                       ii = nbidta(ib,igrd,ib_bdy)- iw_b(1)+2
                       if((com_west_b .ne. 1) .and. (ii == (nlcit(nowe+1)-1))) then
                          ij = nbjdta(ib,igrd,ib_bdy) - is_b(1)+2
                          if((ij == 2) .and. (nbondj == 0 .or. nbondj == 1)) then
                            com_south = 1
                          elseif((ij == nlcjt(nowe+1)-1) .and. (nbondj == 0 .or. nbondj == -1)) then
                            com_north = 1
                          endif
                          com_west_b = 1
                       endif 
                     ENDIF
                     IF( nbidta(ib,igrd,ib_bdy) >= iw_b(2) .AND. nbidta(ib,igrd,ib_bdy) <= ie_b(2) .AND.   &
                       & nbjdta(ib,igrd,ib_bdy) >= is_b(2) .AND. nbjdta(ib,igrd,ib_bdy) <= in_b(2) .AND.   &
                       & nbrdta(ib,igrd,ib_bdy) == ir  ) THEN
                       ii = nbidta(ib,igrd,ib_bdy)- iw_b(2)+2
                       if((com_east_b .ne. 1) .and. (ii == 2)) then
                          ij = nbjdta(ib,igrd,ib_bdy) - is_b(2)+2
                          if((ij == 2) .and. (nbondj == 0 .or. nbondj == 1)) then
                            com_south = 1
                          elseif((ij == nlcjt(noea+1)-1) .and. (nbondj == 0 .or. nbondj == -1)) then
                            com_north = 1
                          endif
                          com_east_b = 1
                       endif 
                     ENDIF
                  ELSEIF(nbondi == 1) THEN
                     IF( nbidta(ib,igrd,ib_bdy) >= iw_b(1) .AND. nbidta(ib,igrd,ib_bdy) <= ie_b(1) .AND.   &
                       & nbjdta(ib,igrd,ib_bdy) >= is_b(1) .AND. nbjdta(ib,igrd,ib_bdy) <= in_b(1) .AND.   &
                       & nbrdta(ib,igrd,ib_bdy) == ir  ) THEN
                       ii = nbidta(ib,igrd,ib_bdy)- iw_b(1)+2
                       if((com_west_b .ne. 1) .and. (ii == (nlcit(nowe+1)-1))) then
                          ij = nbjdta(ib,igrd,ib_bdy) - is_b(1)+2
                          if((ij == 2) .and. (nbondj == 0 .or. nbondj == 1)) then
                            com_south = 1
                          elseif((ij == nlcjt(nowe+1)-1) .and. (nbondj == 0 .or. nbondj == -1)) then
                            com_north = 1
                          endif
                          com_west_b = 1
                       endif 
                     ENDIF
                  ELSEIF(nbondi == -1) THEN
                     IF( nbidta(ib,igrd,ib_bdy) >= iw_b(2) .AND. nbidta(ib,igrd,ib_bdy) <= ie_b(2) .AND.   &
                       & nbjdta(ib,igrd,ib_bdy) >= is_b(2) .AND. nbjdta(ib,igrd,ib_bdy) <= in_b(2) .AND.   &
                       & nbrdta(ib,igrd,ib_bdy) == ir  ) THEN
                       ii = nbidta(ib,igrd,ib_bdy)- iw_b(2)+2
                       if((com_east_b .ne. 1) .and. (ii == 2)) then
                          ij = nbjdta(ib,igrd,ib_bdy) - is_b(2)+2
                          if((ij == 2) .and. (nbondj == 0 .or. nbondj == 1)) then
                            com_south = 1
                          elseif((ij == nlcjt(noea+1)-1) .and. (nbondj == 0 .or. nbondj == -1)) then
                            com_north = 1
                          endif
                          com_east_b = 1
                       endif 
                     ENDIF
                  ENDIF
                  IF(nbondj == 0) THEN
                     IF(com_north_b .ne. 1 .AND. (nbidta(ib,igrd,ib_bdy) == iw_b(4)-1  &
                       & .OR. nbidta(ib,igrd,ib_bdy) == ie_b(4)+1) .AND. &
                       & nbjdta(ib,igrd,ib_bdy) == is_b(4) .AND. nbrdta(ib,igrd,ib_bdy) == ir) THEN
                       com_north_b = 1 
                     ENDIF
                     IF(com_south_b .ne. 1 .AND. (nbidta(ib,igrd,ib_bdy) == iw_b(3)-1  &
                       &.OR. nbidta(ib,igrd,ib_bdy) == ie_b(3)+1) .AND. &
                       & nbjdta(ib,igrd,ib_bdy) == in_b(3) .AND. nbrdta(ib,igrd,ib_bdy) == ir) THEN
                       com_south_b = 1 
                     ENDIF
                     IF( nbidta(ib,igrd,ib_bdy) >= iw_b(3) .AND. nbidta(ib,igrd,ib_bdy) <= ie_b(3) .AND.   &
                       & nbjdta(ib,igrd,ib_bdy) >= is_b(3) .AND. nbjdta(ib,igrd,ib_bdy) <= in_b(3) .AND.   &
                       & nbrdta(ib,igrd,ib_bdy) == ir  ) THEN
                       ij = nbjdta(ib,igrd,ib_bdy)- is_b(3)+2
                       if((com_south_b .ne. 1) .and. (ij == (nlcjt(noso+1)-1))) then
                          com_south_b = 1
                       endif 
                     ENDIF
                     IF( nbidta(ib,igrd,ib_bdy) >= iw_b(4) .AND. nbidta(ib,igrd,ib_bdy) <= ie_b(4) .AND.   &
                       & nbjdta(ib,igrd,ib_bdy) >= is_b(4) .AND. nbjdta(ib,igrd,ib_bdy) <= in_b(4) .AND.   &
                       & nbrdta(ib,igrd,ib_bdy) == ir  ) THEN
                       ij = nbjdta(ib,igrd,ib_bdy)- is_b(4)+2
                       if((com_north_b .ne. 1) .and. (ij == 2)) then
                          com_north_b = 1
                       endif 
                     ENDIF
                  ELSEIF(nbondj == 1) THEN
                     IF( com_south_b .ne. 1 .AND. (nbidta(ib,igrd,ib_bdy) == iw_b(3)-1 .OR. &
                       & nbidta(ib,igrd,ib_bdy) == ie_b(3)+1) .AND. &
                       & nbjdta(ib,igrd,ib_bdy) == in_b(3) .AND. nbrdta(ib,igrd,ib_bdy) == ir) THEN
                       com_south_b = 1 
                     ENDIF
                     IF( nbidta(ib,igrd,ib_bdy) >= iw_b(3) .AND. nbidta(ib,igrd,ib_bdy) <= ie_b(3) .AND.   &
                       & nbjdta(ib,igrd,ib_bdy) >= is_b(3) .AND. nbjdta(ib,igrd,ib_bdy) <= in_b(3) .AND.   &
                       & nbrdta(ib,igrd,ib_bdy) == ir  ) THEN
                       ij = nbjdta(ib,igrd,ib_bdy)- is_b(3)+2
                       if((com_south_b .ne. 1) .and. (ij == (nlcjt(noso+1)-1))) then
                          com_south_b = 1
                       endif 
                     ENDIF
                  ELSEIF(nbondj == -1) THEN
                     IF(com_north_b .ne. 1 .AND. (nbidta(ib,igrd,ib_bdy) == iw_b(4)-1  &
                       & .OR. nbidta(ib,igrd,ib_bdy) == ie_b(4)+1) .AND. &
                       & nbjdta(ib,igrd,ib_bdy) == is_b(4) .AND. nbrdta(ib,igrd,ib_bdy) == ir) THEN
                       com_north_b = 1 
                     ENDIF
                     IF( nbidta(ib,igrd,ib_bdy) >= iw_b(4) .AND. nbidta(ib,igrd,ib_bdy) <= ie_b(4) .AND.   &
                       & nbjdta(ib,igrd,ib_bdy) >= is_b(4) .AND. nbjdta(ib,igrd,ib_bdy) <= in_b(4) .AND.   &
                       & nbrdta(ib,igrd,ib_bdy) == ir  ) THEN
                       ij = nbjdta(ib,igrd,ib_bdy)- is_b(4)+2
                       if((com_north_b .ne. 1) .and. (ij == 2)) then
                          com_north_b = 1
                       endif 
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO 

         ! definition of the i- and j- direction local boundaries arrays used for sending the boundaries
         IF(     (com_east  == 1) .and. (com_west  == 1) ) THEN   ;   nbondi_bdy(ib_bdy) =  0
         ELSEIF( (com_east  == 1) .and. (com_west  == 0) ) THEN   ;   nbondi_bdy(ib_bdy) = -1
         ELSEIF( (com_east  == 0) .and. (com_west  == 1) ) THEN   ;   nbondi_bdy(ib_bdy) =  1
         ENDIF
         IF(     (com_north == 1) .and. (com_south == 1) ) THEN   ;   nbondj_bdy(ib_bdy) =  0
         ELSEIF( (com_north == 1) .and. (com_south == 0) ) THEN   ;   nbondj_bdy(ib_bdy) = -1
         ELSEIF( (com_north == 0) .and. (com_south == 1) ) THEN   ;   nbondj_bdy(ib_bdy) =  1
         ENDIF

         ! definition of the i- and j- direction local boundaries arrays used for receiving the boundaries
         IF(     (com_east_b  == 1) .and. (com_west_b  == 1) ) THEN   ;   nbondi_bdy_b(ib_bdy) =  0
         ELSEIF( (com_east_b  == 1) .and. (com_west_b  == 0) ) THEN   ;   nbondi_bdy_b(ib_bdy) = -1
         ELSEIF( (com_east_b  == 0) .and. (com_west_b  == 1) ) THEN   ;   nbondi_bdy_b(ib_bdy) =  1
         ENDIF
         IF(     (com_north_b == 1) .and. (com_south_b == 1) ) THEN   ;   nbondj_bdy_b(ib_bdy) =  0
         ELSEIF( (com_north_b == 1) .and. (com_south_b == 0) ) THEN   ;   nbondj_bdy_b(ib_bdy) = -1
         ELSEIF( (com_north_b == 0) .and. (com_south_b == 1) ) THEN   ;   nbondj_bdy_b(ib_bdy) =  1
         ENDIF

         ! Compute rim weights for FRS scheme
         ! ----------------------------------
         DO igrd = 1, jpbgrd
            DO ib = 1, idx_bdy(ib_bdy)%nblen(igrd)
               nbr => idx_bdy(ib_bdy)%nbr(ib,igrd)
               idx_bdy(ib_bdy)%nbw(ib,igrd) = 1.- TANH( REAL( nbr - 1 ) *0.5 )      ! tanh formulation
!               idx_bdy(ib_bdy)%nbw(ib,igrd) = (REAL(nn_rimwidth(ib_bdy)+1-nbr)/REAL(nn_rimwidth(ib_bdy)))**2.  ! quadratic
!               idx_bdy(ib_bdy)%nbw(ib,igrd) =  REAL(nn_rimwidth(ib_bdy)+1-nbr)/REAL(nn_rimwidth(ib_bdy))       ! linear
            END DO
         END DO 

         ! Compute damping coefficients
         ! ----------------------------
         DO igrd = 1, jpbgrd
            DO ib = 1, idx_bdy(ib_bdy)%nblen(igrd)
               nbr => idx_bdy(ib_bdy)%nbr(ib,igrd)
               idx_bdy(ib_bdy)%nbd(ib,igrd) = 1. / ( rn_time_dmp(ib_bdy) * rday ) & 
               & *(REAL(nn_rimwidth(ib_bdy)+1-nbr)/REAL(nn_rimwidth(ib_bdy)))**2.   ! quadratic
               idx_bdy(ib_bdy)%nbdout(ib,igrd) = 1. / ( rn_time_dmp_out(ib_bdy) * rday ) & 
               & *(REAL(nn_rimwidth(ib_bdy)+1-nbr)/REAL(nn_rimwidth(ib_bdy)))**2.   ! quadratic
            END DO
         END DO 

      END DO

      ! ------------------------------------------------------
      ! Initialise masks and find normal/tangential directions
      ! ------------------------------------------------------

      ! Read global 2D mask at T-points: bdytmask
      ! -----------------------------------------
      ! bdytmask = 1  on the computational domain AND on open boundaries
      !          = 0  elsewhere   
 
      bdytmask(:,:) = ssmask(:,:)

      ! Derive mask on U and V grid from mask on T grid

      bdyumask(:,:) = 0._wp
      bdyvmask(:,:) = 0._wp
      DO ij = 1, jpjm1
         DO ii = 1, jpim1
            bdyumask(ii,ij) = bdytmask(ii,ij) * bdytmask(ii+1, ij )
            bdyvmask(ii,ij) = bdytmask(ii,ij) * bdytmask(ii  ,ij+1)  
         END DO
      END DO
      CALL lbc_lnk_multi( 'bdyini', bdyumask, 'U', 1. , bdyvmask, 'V', 1. )   ! Lateral boundary cond. 

      ! bdy masks are now set to zero on boundary points:
      !
      igrd = 1
      DO ib_bdy = 1, nb_bdy
        DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)      
          bdytmask(idx_bdy(ib_bdy)%nbi(ib,igrd), idx_bdy(ib_bdy)%nbj(ib,igrd)) = 0._wp
        END DO
      END DO
      !
      igrd = 2
      DO ib_bdy = 1, nb_bdy
        DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)
          bdyumask(idx_bdy(ib_bdy)%nbi(ib,igrd), idx_bdy(ib_bdy)%nbj(ib,igrd)) = 0._wp
        END DO
      END DO
      !
      igrd = 3
      DO ib_bdy = 1, nb_bdy
        DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)
          bdyvmask(idx_bdy(ib_bdy)%nbi(ib,igrd), idx_bdy(ib_bdy)%nbj(ib,igrd)) = 0._wp
        END DO
      END DO

      ! For the flagu/flagv calculation below we require a version of fmask without
      ! the land boundary condition (shlat) included:
      zfmask(:,:) = 0
      DO ij = 2, jpjm1
         DO ii = 2, jpim1
            zfmask(ii,ij) = tmask(ii,ij  ,1) * tmask(ii+1,ij  ,1)   &
           &              * tmask(ii,ij+1,1) * tmask(ii+1,ij+1,1)
         END DO      
      END DO

      ! Lateral boundary conditions
      CALL lbc_lnk( 'bdyini', zfmask, 'F', 1. ) 
      CALL lbc_lnk_multi( 'bdyini', bdyumask, 'U', 1. , bdyvmask, 'V', 1., bdytmask, 'T', 1. )
      DO ib_bdy = 1, nb_bdy       ! Indices and directions of rim velocity components

         idx_bdy(ib_bdy)%flagu(:,:) = 0._wp
         idx_bdy(ib_bdy)%flagv(:,:) = 0._wp
         icount = 0 

         ! Calculate relationship of U direction to the local orientation of the boundary
         ! flagu = -1 : u component is normal to the dynamical boundary and its direction is outward
         ! flagu =  0 : u is tangential
         ! flagu =  1 : u is normal to the boundary and is direction is inward
  
         DO igrd = 1, jpbgrd 
            SELECT CASE( igrd )
               CASE( 1 )   ;   pmask => umask   (:,:,1)   ;   i_offset = 0
               CASE( 2 )   ;   pmask => bdytmask(:,:)     ;   i_offset = 1
               CASE( 3 )   ;   pmask => zfmask  (:,:)     ;   i_offset = 0
            END SELECT 
            icount = 0
            DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)  
               nbi => idx_bdy(ib_bdy)%nbi(ib,igrd)
               nbj => idx_bdy(ib_bdy)%nbj(ib,igrd)
               zefl = pmask(nbi+i_offset-1,nbj)
               zwfl = pmask(nbi+i_offset,nbj)
               ! This error check only works if you are using the bdyXmask arrays
               IF( i_offset == 1 .and. zefl + zwfl == 2 ) THEN
                  icount = icount + 1
                  IF(lwp) WRITE(numout,*) 'Problem with igrd = ',igrd,' at (global) nbi, nbj : ',mig(nbi),mjg(nbj)
               ELSE
                  idx_bdy(ib_bdy)%flagu(ib,igrd) = -zefl + zwfl
               ENDIF
            END DO
            IF( icount /= 0 ) THEN
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) ' E R R O R : Some ',cgrid(igrd),' grid points,',   &
                  ' are not boundary points (flagu calculation). Check nbi, nbj, indices for boundary set ',ib_bdy
               IF(lwp) WRITE(numout,*) ' ========== '
               IF(lwp) WRITE(numout,*)
               nstop = nstop + 1
            ENDIF 
         END DO

         ! Calculate relationship of V direction to the local orientation of the boundary
         ! flagv = -1 : v component is normal to the dynamical boundary but its direction is outward
         ! flagv =  0 : v is tangential
         ! flagv =  1 : v is normal to the boundary and is direction is inward

         DO igrd = 1, jpbgrd 
            SELECT CASE( igrd )
               CASE( 1 )   ;   pmask => vmask (:,:,1)   ;   j_offset = 0
               CASE( 2 )   ;   pmask => zfmask(:,:)     ;   j_offset = 0
               CASE( 3 )   ;   pmask => bdytmask        ;   j_offset = 1
            END SELECT 
            icount = 0
            DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)  
               nbi => idx_bdy(ib_bdy)%nbi(ib,igrd)
               nbj => idx_bdy(ib_bdy)%nbj(ib,igrd)
               znfl = pmask(nbi,nbj+j_offset-1)
               zsfl = pmask(nbi,nbj+j_offset  )
               ! This error check only works if you are using the bdyXmask arrays
               IF( j_offset == 1 .and. znfl + zsfl == 2 ) THEN
                  IF(lwp) WRITE(numout,*) 'Problem with igrd = ',igrd,' at (global) nbi, nbj : ',mig(nbi),mjg(nbj)
                  icount = icount + 1
               ELSE
                  idx_bdy(ib_bdy)%flagv(ib,igrd) = -znfl + zsfl
               END IF
            END DO
            IF( icount /= 0 ) THEN
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) ' E R R O R : Some ',cgrid(igrd),' grid points,',   &
                  ' are not boundary points (flagv calculation). Check nbi, nbj, indices for boundary set ',ib_bdy
               IF(lwp) WRITE(numout,*) ' ========== '
               IF(lwp) WRITE(numout,*)
               nstop = nstop + 1
            ENDIF 
         END DO
         !
      END DO

      ! Compute total lateral surface for volume correction:
      ! ----------------------------------------------------
      ! JC: this must be done at each time step with non-linear free surface
      bdysurftot = 0._wp 
      IF( ln_vol ) THEN  
         igrd = 2      ! Lateral surface at U-points
         DO ib_bdy = 1, nb_bdy
            DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)
               nbi => idx_bdy(ib_bdy)%nbi(ib,igrd)
               nbj => idx_bdy(ib_bdy)%nbj(ib,igrd)
               flagu => idx_bdy(ib_bdy)%flagu(ib,igrd)
               bdysurftot = bdysurftot + hu_n   (nbi  , nbj)                           &
                  &                    * e2u    (nbi  , nbj) * ABS( flagu ) &
                  &                    * tmask_i(nbi  , nbj)                           &
                  &                    * tmask_i(nbi+1, nbj)                   
            END DO
         END DO

         igrd=3 ! Add lateral surface at V-points
         DO ib_bdy = 1, nb_bdy
            DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)
               nbi => idx_bdy(ib_bdy)%nbi(ib,igrd)
               nbj => idx_bdy(ib_bdy)%nbj(ib,igrd)
               flagv => idx_bdy(ib_bdy)%flagv(ib,igrd)
               bdysurftot = bdysurftot + hv_n   (nbi, nbj  )                           &
                  &                    * e1v    (nbi, nbj  ) * ABS( flagv ) &
                  &                    * tmask_i(nbi, nbj  )                           &
                  &                    * tmask_i(nbi, nbj+1)
            END DO
         END DO
         !
         CALL mpp_sum( 'bdyini', bdysurftot )      ! sum over the global domain
      END IF   
      !
      ! Tidy up
      !--------
      IF( nb_bdy>0 )   DEALLOCATE( nbidta, nbjdta, nbrdta )
      !
   END SUBROUTINE bdy_segs


   SUBROUTINE bdy_ctl_seg
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE bdy_ctl_seg  ***
      !!
      !! ** Purpose :   Check straight open boundary segments location
      !!
      !! ** Method  :   - Look for open boundary corners
      !!                - Check that segments start or end on land 
      !!----------------------------------------------------------------------
      INTEGER  ::   ib, ib1, ib2, ji ,jj, itest  
      INTEGER, DIMENSION(jp_nseg,2) :: icorne, icornw, icornn, icorns  
      REAL(wp), DIMENSION(2) ::   ztestmask
      !!----------------------------------------------------------------------
      !
      IF (lwp) WRITE(numout,*) ' '
      IF (lwp) WRITE(numout,*) 'bdy_ctl_seg: Check analytical segments'
      IF (lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      !
      IF(lwp) WRITE(numout,*) 'Number of east  segments     : ', nbdysege
      IF(lwp) WRITE(numout,*) 'Number of west  segments     : ', nbdysegw
      IF(lwp) WRITE(numout,*) 'Number of north segments     : ', nbdysegn
      IF(lwp) WRITE(numout,*) 'Number of south segments     : ', nbdysegs
      ! 1. Check bounds
      !----------------
      DO ib = 1, nbdysegn
         IF (lwp) WRITE(numout,*) '**check north seg bounds pckg: ', npckgn(ib)
         IF ((jpjnob(ib).ge.jpjglo-1).or.& 
            &(jpjnob(ib).le.1))        CALL ctl_stop( 'nbdyind out of domain' )
         IF (jpindt(ib).ge.jpinft(ib)) CALL ctl_stop( 'Bdy start index is greater than end index' )
         IF (jpindt(ib).le.1     )     CALL ctl_stop( 'Start index out of domain' )
         IF (jpinft(ib).ge.jpiglo)     CALL ctl_stop( 'End index out of domain' )
      END DO
      !
      DO ib = 1, nbdysegs
         IF (lwp) WRITE(numout,*) '**check south seg bounds pckg: ', npckgs(ib)
         IF ((jpjsob(ib).ge.jpjglo-1).or.& 
            &(jpjsob(ib).le.1))        CALL ctl_stop( 'nbdyind out of domain' )
         IF (jpisdt(ib).ge.jpisft(ib)) CALL ctl_stop( 'Bdy start index is greater than end index' )
         IF (jpisdt(ib).le.1     )     CALL ctl_stop( 'Start index out of domain' )
         IF (jpisft(ib).ge.jpiglo)     CALL ctl_stop( 'End index out of domain' )
      END DO
      !
      DO ib = 1, nbdysege
         IF (lwp) WRITE(numout,*) '**check east  seg bounds pckg: ', npckge(ib)
         IF ((jpieob(ib).ge.jpiglo-1).or.& 
            &(jpieob(ib).le.1))        CALL ctl_stop( 'nbdyind out of domain' )
         IF (jpjedt(ib).ge.jpjeft(ib)) CALL ctl_stop( 'Bdy start index is greater than end index' )
         IF (jpjedt(ib).lt.1     )     CALL ctl_stop( 'Start index out of domain' )!ACC
         IF (jpjeft(ib).gt.jpjglo)     CALL ctl_stop( 'End index out of domain' )!ACC
      END DO
      !
      DO ib = 1, nbdysegw
         IF (lwp) WRITE(numout,*) '**check west  seg bounds pckg: ', npckgw(ib)
         IF ((jpiwob(ib).ge.jpiglo-1).or.& 
            &(jpiwob(ib).le.1))        CALL ctl_stop( 'nbdyind out of domain' )
         IF (jpjwdt(ib).ge.jpjwft(ib)) CALL ctl_stop( 'Bdy start index is greater than end index' )
         IF (jpjwdt(ib).lt.1     )     CALL ctl_stop( 'Start index out of domain' ) !ACC
         IF (jpjwft(ib).gt.jpjglo)     CALL ctl_stop( 'End index out of domain' ) !ACC
      ENDDO
      !
      !      
      ! 2. Look for segment crossings
      !------------------------------ 
      IF (lwp) WRITE(numout,*) '**Look for segments corners  :'
      !
      itest = 0 ! corner number
      !
      ! flag to detect if start or end of open boundary belongs to a corner
      ! if not (=0), it must be on land.
      ! if a corner is detected, save bdy package number for further tests
      icorne(:,:)=0. ; icornw(:,:)=0. ; icornn(:,:)=0. ; icorns(:,:)=0.
      ! South/West crossings
      IF ((nbdysegw > 0).AND.(nbdysegs > 0)) THEN
         DO ib1 = 1, nbdysegw        
            DO ib2 = 1, nbdysegs
               IF (( jpisdt(ib2)<=jpiwob(ib1)).AND. &
                &  ( jpisft(ib2)>=jpiwob(ib1)).AND. &
                &  ( jpjwdt(ib1)<=jpjsob(ib2)).AND. &
                &  ( jpjwft(ib1)>=jpjsob(ib2))) THEN
                  IF ((jpjwdt(ib1)==jpjsob(ib2)).AND.(jpisdt(ib2)==jpiwob(ib1))) THEN 
                     ! We have a possible South-West corner                      
!                     WRITE(numout,*) ' Found a South-West corner at (i,j): ', jpisdt(ib2), jpjwdt(ib1) 
!                     WRITE(numout,*) ' between segments: ', npckgw(ib1), npckgs(ib2)
                     icornw(ib1,1) = npckgs(ib2)
                     icorns(ib2,1) = npckgw(ib1)
                  ELSEIF ((jpisft(ib2)==jpiwob(ib1)).AND.(jpjwft(ib1)==jpjsob(ib2))) THEN
                     IF(lwp) WRITE(numout,*)
                     IF(lwp) WRITE(numout,*) ' E R R O R : Found an acute open boundary corner at point (i,j)= ', &
                     &                                     jpisft(ib2), jpjwft(ib1)
                     IF(lwp) WRITE(numout,*) ' ==========  Not allowed yet'
                     IF(lwp) WRITE(numout,*) '             Crossing problem with West segment: ',npckgw(ib1), & 
                     &                                                    ' and South segment: ',npckgs(ib2)
                     IF(lwp) WRITE(numout,*)
                     nstop = nstop + 1
                  ELSE
                     IF(lwp) WRITE(numout,*)
                     IF(lwp) WRITE(numout,*) ' E R R O R : Check South and West Open boundary indices'
                     IF(lwp) WRITE(numout,*) ' ==========  Crossing problem with West segment: ',npckgw(ib1) , &
                     &                                                    ' and South segment: ',npckgs(ib2)
                     IF(lwp) WRITE(numout,*)
                     nstop = nstop+1
                  END IF
               END IF
            END DO
         END DO
      END IF
      !
      ! South/East crossings
      IF ((nbdysege > 0).AND.(nbdysegs > 0)) THEN
         DO ib1 = 1, nbdysege
            DO ib2 = 1, nbdysegs
               IF (( jpisdt(ib2)<=jpieob(ib1)+1).AND. &
                &  ( jpisft(ib2)>=jpieob(ib1)+1).AND. &
                &  ( jpjedt(ib1)<=jpjsob(ib2)  ).AND. &
                &  ( jpjeft(ib1)>=jpjsob(ib2)  )) THEN
                  IF ((jpjedt(ib1)==jpjsob(ib2)).AND.(jpisft(ib2)==jpieob(ib1)+1)) THEN
                     ! We have a possible South-East corner 
!                     WRITE(numout,*) ' Found a South-East corner at (i,j): ', jpisft(ib2), jpjedt(ib1) 
!                     WRITE(numout,*) ' between segments: ', npckge(ib1), npckgs(ib2)
                     icorne(ib1,1) = npckgs(ib2)
                     icorns(ib2,2) = npckge(ib1)
                  ELSEIF ((jpjeft(ib1)==jpjsob(ib2)).AND.(jpisdt(ib2)==jpieob(ib1)+1)) THEN
                     IF(lwp) WRITE(numout,*)
                     IF(lwp) WRITE(numout,*) ' E R R O R : Found an acute open boundary corner at point (i,j)= ', &
                     &                                     jpisdt(ib2), jpjeft(ib1)
                     IF(lwp) WRITE(numout,*) ' ==========  Not allowed yet'
                     IF(lwp) WRITE(numout,*) '             Crossing problem with East segment: ',npckge(ib1), &
                     &                                                    ' and South segment: ',npckgs(ib2)
                     IF(lwp) WRITE(numout,*)
                     nstop = nstop + 1
                  ELSE
                     IF(lwp) WRITE(numout,*)
                     IF(lwp) WRITE(numout,*) ' E R R O R : Check South and East Open boundary indices'
                     IF(lwp) WRITE(numout,*) ' ==========  Crossing problem with East segment: ',npckge(ib1), &
                     &                                                    ' and South segment: ',npckgs(ib2)
                     IF(lwp) WRITE(numout,*)
                     nstop = nstop + 1
                  END IF
               END IF
            END DO
         END DO
      END IF
      !
      ! North/West crossings
      IF ((nbdysegn > 0).AND.(nbdysegw > 0)) THEN
         DO ib1 = 1, nbdysegw        
            DO ib2 = 1, nbdysegn
               IF (( jpindt(ib2)<=jpiwob(ib1)  ).AND. &
                &  ( jpinft(ib2)>=jpiwob(ib1)  ).AND. &
                &  ( jpjwdt(ib1)<=jpjnob(ib2)+1).AND. &
                &  ( jpjwft(ib1)>=jpjnob(ib2)+1)) THEN
                  IF ((jpjwft(ib1)==jpjnob(ib2)+1).AND.(jpindt(ib2)==jpiwob(ib1))) THEN
                     ! We have a possible North-West corner 
!                     WRITE(numout,*) ' Found a North-West corner at (i,j): ', jpindt(ib2), jpjwft(ib1) 
!                     WRITE(numout,*) ' between segments: ', npckgw(ib1), npckgn(ib2)
                     icornw(ib1,2) = npckgn(ib2)
                     icornn(ib2,1) = npckgw(ib1)
                  ELSEIF ((jpjwdt(ib1)==jpjnob(ib2)+1).AND.(jpinft(ib2)==jpiwob(ib1))) THEN
                     IF(lwp) WRITE(numout,*)
                     IF(lwp) WRITE(numout,*) ' E R R O R : Found an acute open boundary corner at point (i,j)= ', &
                     &                                     jpinft(ib2), jpjwdt(ib1)
                     IF(lwp) WRITE(numout,*) ' ==========  Not allowed yet'
                     IF(lwp) WRITE(numout,*) '             Crossing problem with West segment: ',npckgw(ib1), &
                     &                                                    ' and North segment: ',npckgn(ib2)
                     IF(lwp) WRITE(numout,*)
                     nstop = nstop + 1
                  ELSE
                     IF(lwp) WRITE(numout,*)
                     IF(lwp) WRITE(numout,*) ' E R R O R : Check North and West Open boundary indices'
                     IF(lwp) WRITE(numout,*) ' ==========  Crossing problem with West segment: ',npckgw(ib1), &
                     &                                                    ' and North segment: ',npckgn(ib2)
                     IF(lwp) WRITE(numout,*)
                     nstop = nstop + 1
                  END IF
               END IF
            END DO
         END DO
      END IF
      !
      ! North/East crossings
      IF ((nbdysegn > 0).AND.(nbdysege > 0)) THEN
         DO ib1 = 1, nbdysege        
            DO ib2 = 1, nbdysegn
               IF (( jpindt(ib2)<=jpieob(ib1)+1).AND. &
                &  ( jpinft(ib2)>=jpieob(ib1)+1).AND. &
                &  ( jpjedt(ib1)<=jpjnob(ib2)+1).AND. &
                &  ( jpjeft(ib1)>=jpjnob(ib2)+1)) THEN
                  IF ((jpjeft(ib1)==jpjnob(ib2)+1).AND.(jpinft(ib2)==jpieob(ib1)+1)) THEN
                     ! We have a possible North-East corner 
!                     WRITE(numout,*) ' Found a North-East corner at (i,j): ', jpinft(ib2), jpjeft(ib1)
!                     WRITE(numout,*) ' between segments: ', npckge(ib1), npckgn(ib2)
                     icorne(ib1,2) = npckgn(ib2)
                     icornn(ib2,2) = npckge(ib1)
                  ELSEIF ((jpjedt(ib1)==jpjnob(ib2)+1).AND.(jpindt(ib2)==jpieob(ib1)+1)) THEN
                     IF(lwp) WRITE(numout,*)
                     IF(lwp) WRITE(numout,*) ' E R R O R : Found an acute open boundary corner at point (i,j)= ', &
                     &                                     jpindt(ib2), jpjedt(ib1)
                     IF(lwp) WRITE(numout,*) ' ==========  Not allowed yet'
                     IF(lwp) WRITE(numout,*) '             Crossing problem with East segment: ',npckge(ib1), &
                     &                                                    ' and North segment: ',npckgn(ib2)
                     IF(lwp) WRITE(numout,*)
                     nstop = nstop + 1
                  ELSE
                     IF(lwp) WRITE(numout,*)
                     IF(lwp) WRITE(numout,*) ' E R R O R : Check North and East Open boundary indices'
                     IF(lwp) WRITE(numout,*) ' ==========  Crossing problem with East segment: ',npckge(ib1), &
                     &                                                    ' and North segment: ',npckgn(ib2)
                     IF(lwp) WRITE(numout,*)
                     nstop = nstop + 1
                  END IF
               END IF
            END DO
         END DO
      END IF
      !
      ! 3. Check if segment extremities are on land
      !-------------------------------------------- 
      !
      ! West segments
      DO ib = 1, nbdysegw
         ! get mask at boundary extremities:
         ztestmask(1:2)=0.
         DO ji = 1, jpi
            DO jj = 1, jpj             
              IF (((ji + nimpp - 1) == jpiwob(ib)).AND. & 
               &  ((jj + njmpp - 1) == jpjwdt(ib))) ztestmask(1)=tmask(ji,jj,1)
              IF (((ji + nimpp - 1) == jpiwob(ib)).AND. & 
               &  ((jj + njmpp - 1) == jpjwft(ib))) ztestmask(2)=tmask(ji,jj,1)  
            END DO
         END DO
         CALL mpp_sum( 'bdyini', ztestmask, 2 )   ! sum over the global domain

         IF (ztestmask(1)==1) THEN 
            IF (icornw(ib,1)==0) THEN
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) ' E R R O R : Open boundary segment ', npckgw(ib)
               IF(lwp) WRITE(numout,*) ' ==========  does not start on land or on a corner'                                                  
               IF(lwp) WRITE(numout,*)
               nstop = nstop + 1
            ELSE
               ! This is a corner
               IF(lwp) WRITE(numout,*) 'Found a South-West corner at (i,j): ', jpiwob(ib), jpjwdt(ib)
               CALL bdy_ctl_corn(npckgw(ib), icornw(ib,1))
               itest=itest+1
            ENDIF
         ENDIF
         IF (ztestmask(2)==1) THEN
            IF (icornw(ib,2)==0) THEN
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) ' E R R O R : Open boundary segment ', npckgw(ib)
               IF(lwp) WRITE(numout,*) ' ==========  does not end on land or on a corner'                                                  
               IF(lwp) WRITE(numout,*)
               nstop = nstop + 1
            ELSE
               ! This is a corner
               IF(lwp) WRITE(numout,*) 'Found a North-West corner at (i,j): ', jpiwob(ib), jpjwft(ib)
               CALL bdy_ctl_corn(npckgw(ib), icornw(ib,2))
               itest=itest+1
            ENDIF
         ENDIF
      END DO
      !
      ! East segments
      DO ib = 1, nbdysege
         ! get mask at boundary extremities:
         ztestmask(1:2)=0.
         DO ji = 1, jpi
            DO jj = 1, jpj             
              IF (((ji + nimpp - 1) == jpieob(ib)+1).AND. & 
               &  ((jj + njmpp - 1) == jpjedt(ib))) ztestmask(1)=tmask(ji,jj,1)
              IF (((ji + nimpp - 1) == jpieob(ib)+1).AND. & 
               &  ((jj + njmpp - 1) == jpjeft(ib))) ztestmask(2)=tmask(ji,jj,1)  
            END DO
         END DO
         CALL mpp_sum( 'bdyini', ztestmask, 2 )   ! sum over the global domain

         IF (ztestmask(1)==1) THEN
            IF (icorne(ib,1)==0) THEN
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) ' E R R O R : Open boundary segment ', npckge(ib)
               IF(lwp) WRITE(numout,*) ' ==========  does not start on land or on a corner'                                                  
               IF(lwp) WRITE(numout,*)
               nstop = nstop + 1 
            ELSE
               ! This is a corner
               IF(lwp) WRITE(numout,*) 'Found a South-East corner at (i,j): ', jpieob(ib)+1, jpjedt(ib)
               CALL bdy_ctl_corn(npckge(ib), icorne(ib,1))
               itest=itest+1
            ENDIF
         ENDIF
         IF (ztestmask(2)==1) THEN
            IF (icorne(ib,2)==0) THEN
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) ' E R R O R : Open boundary segment ', npckge(ib)
               IF(lwp) WRITE(numout,*) ' ==========  does not end on land or on a corner'                                                  
               IF(lwp) WRITE(numout,*)
               nstop = nstop + 1
            ELSE
               ! This is a corner
               IF(lwp) WRITE(numout,*) 'Found a North-East corner at (i,j): ', jpieob(ib)+1, jpjeft(ib)
               CALL bdy_ctl_corn(npckge(ib), icorne(ib,2))
               itest=itest+1
            ENDIF
         ENDIF
      END DO
      !
      ! South segments
      DO ib = 1, nbdysegs
         ! get mask at boundary extremities:
         ztestmask(1:2)=0.
         DO ji = 1, jpi
            DO jj = 1, jpj             
              IF (((jj + njmpp - 1) == jpjsob(ib)).AND. & 
               &  ((ji + nimpp - 1) == jpisdt(ib))) ztestmask(1)=tmask(ji,jj,1)
              IF (((jj + njmpp - 1) == jpjsob(ib)).AND. & 
               &  ((ji + nimpp - 1) == jpisft(ib))) ztestmask(2)=tmask(ji,jj,1)  
            END DO
         END DO
         CALL mpp_sum( 'bdyini', ztestmask, 2 )   ! sum over the global domain

         IF ((ztestmask(1)==1).AND.(icorns(ib,1)==0)) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) ' E R R O R : Open boundary segment ', npckgs(ib)
            IF(lwp) WRITE(numout,*) ' ==========  does not start on land or on a corner'                                                  
            IF(lwp) WRITE(numout,*)
            nstop = nstop + 1
         ENDIF
         IF ((ztestmask(2)==1).AND.(icorns(ib,2)==0)) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) ' E R R O R : Open boundary segment ', npckgs(ib)
            IF(lwp) WRITE(numout,*) ' ==========  does not end on land or on a corner'                                                  
            IF(lwp) WRITE(numout,*)
            nstop = nstop + 1
         ENDIF
      END DO
      !
      ! North segments
      DO ib = 1, nbdysegn
         ! get mask at boundary extremities:
         ztestmask(1:2)=0.
         DO ji = 1, jpi
            DO jj = 1, jpj             
              IF (((jj + njmpp - 1) == jpjnob(ib)+1).AND. & 
               &  ((ji + nimpp - 1) == jpindt(ib))) ztestmask(1)=tmask(ji,jj,1)
              IF (((jj + njmpp - 1) == jpjnob(ib)+1).AND. & 
               &  ((ji + nimpp - 1) == jpinft(ib))) ztestmask(2)=tmask(ji,jj,1)  
            END DO
         END DO
         CALL mpp_sum( 'bdyini', ztestmask, 2 )   ! sum over the global domain

         IF ((ztestmask(1)==1).AND.(icornn(ib,1)==0)) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) ' E R R O R : Open boundary segment ', npckgn(ib)
            IF(lwp) WRITE(numout,*) ' ==========  does not start on land'                                                  
            IF(lwp) WRITE(numout,*)
            nstop = nstop + 1
         ENDIF
         IF ((ztestmask(2)==1).AND.(icornn(ib,2)==0)) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) ' E R R O R : Open boundary segment ', npckgn(ib)
            IF(lwp) WRITE(numout,*) ' ==========  does not end on land'                                                  
            IF(lwp) WRITE(numout,*)
            nstop = nstop + 1
         ENDIF
      END DO
      !
      IF ((itest==0).AND.(lwp)) WRITE(numout,*) 'NO open boundary corner found'
      !
      ! Other tests TBD: 
      ! segments completly on land
      ! optimized open boundary array length according to landmask
      ! Nudging layers that overlap with interior domain
      !
   END SUBROUTINE bdy_ctl_seg


   SUBROUTINE bdy_ctl_corn( ib1, ib2 )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE bdy_ctl_corn  ***
      !!
      !! ** Purpose :   Check numerical schemes consistency between
      !!                segments having a common corner
      !!
      !! ** Method  :   
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)  ::   ib1, ib2
      INTEGER :: itest
      !!----------------------------------------------------------------------
      itest = 0

      IF( cn_dyn2d(ib1) /= cn_dyn2d(ib2) )   itest = itest + 1
      IF( cn_dyn3d(ib1) /= cn_dyn3d(ib2) )   itest = itest + 1
      IF( cn_tra  (ib1) /= cn_tra  (ib2) )   itest = itest + 1
      !
      IF( nn_dyn2d_dta(ib1) /= nn_dyn2d_dta(ib2) )   itest = itest + 1
      IF( nn_dyn3d_dta(ib1) /= nn_dyn3d_dta(ib2) )   itest = itest + 1
      IF( nn_tra_dta  (ib1) /= nn_tra_dta  (ib2) )   itest = itest + 1
      !
      IF( nn_rimwidth(ib1) /= nn_rimwidth(ib2) )   itest = itest + 1   
      !
      IF( itest>0 ) THEN
         IF(lwp) WRITE(numout,*) ' E R R O R : Segments ', ib1, 'and ', ib2
         IF(lwp) WRITE(numout,*) ' ==========  have different open bdy schemes'                                                  
         IF(lwp) WRITE(numout,*)
         nstop = nstop + 1
      ENDIF
      !
   END SUBROUTINE bdy_ctl_corn

   !!=================================================================================
END MODULE bdyini
