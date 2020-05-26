MODULE sbcssm
   !!======================================================================
   !!                       ***  MODULE  sbcssm  ***
   !! Off-line : interpolation of the physical fields
   !!======================================================================
   !! History :  3.4  ! 2012-03 (S. Alderson)  original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_ssm_init  : initialization, namelist read, and SAVEs control
   !!   sbc_ssm       : Interpolation of the fields
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE c1d            ! 1D configuration: lk_c1d
   USE dom_oce        ! ocean domain: variables
   USE sbc_oce        ! surface module: variables
   USE phycst         ! physical constants
   USE eosbn2         ! equation of state - Brunt Vaisala frequency
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE lib_mpp        ! distributed memory computing library
   USE prtctl         ! print control
   USE fldread        ! read input fields
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_ssm_init   ! called by sbc_init
   PUBLIC   sbc_ssm        ! called by sbc

   CHARACTER(len=100) ::   cn_dir        ! Root directory for location of ssm files
   LOGICAL            ::   ln_3d_uve     ! specify whether input velocity data is 3D
   LOGICAL            ::   ln_read_frq   ! specify whether we must read frq or not

   LOGICAL            ::   l_sasread     ! Ice intilisation: =T read a file ; =F anaytical initilaistion
   LOGICAL            ::   l_initdone = .false.
   INTEGER     ::   nfld_3d
   INTEGER     ::   nfld_2d

   INTEGER     ::   jf_tem         ! index of temperature
   INTEGER     ::   jf_sal         ! index of salinity
   INTEGER     ::   jf_usp         ! index of u velocity component
   INTEGER     ::   jf_vsp         ! index of v velocity component
   INTEGER     ::   jf_ssh         ! index of sea surface height
   INTEGER     ::   jf_e3t         ! index of first T level thickness
   INTEGER     ::   jf_frq         ! index of fraction of qsr absorbed in the 1st T level

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_ssm_3d  ! structure of input fields (file information, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_ssm_2d  ! structure of input fields (file information, fields read)

   !!----------------------------------------------------------------------
   !! NEMO/SAS 4.0 , NEMO Consortium (2018)
   !! $Id: sbcssm.F90 12615 2020-03-26 15:18:49Z laurent $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_ssm( kt, Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_ssm  ***
      !!
      !! ** Purpose :  Prepares dynamics and physics fields from a NEMO run
      !!               for an off-line simulation using surface processes only
      !!
      !! ** Method : calculates the position of data
      !!             - interpolates data if needed
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm   ! ocean time level indices
      ! (not needed for SAS but needed to keep a consistent interface in sbcmod.F90)
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) ::   ztinta     ! ratio applied to after  records when doing time interpolation
      REAL(wp) ::   ztintb     ! ratio applied to before records when doing time interpolation
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start( 'sbc_ssm')

      IF ( l_sasread ) THEN
         IF( nfld_3d > 0 ) CALL fld_read( kt, 1, sf_ssm_3d )      !==   read data at kt time step   ==!
         IF( nfld_2d > 0 ) CALL fld_read( kt, 1, sf_ssm_2d )      !==   read data at kt time step   ==!
         !
         IF( ln_3d_uve ) THEN
            IF( .NOT. ln_linssh ) THEN
               e3t_m(:,:) = sf_ssm_3d(jf_e3t)%fnow(:,:,1) * tmask(:,:,1) ! vertical scale factor
            ELSE
               e3t_m(:,:) = e3t_0(:,:,1)                                 ! vertical scale factor
            ENDIF
            ssu_m(:,:) = sf_ssm_3d(jf_usp)%fnow(:,:,1) * umask(:,:,1)    ! u-velocity
            ssv_m(:,:) = sf_ssm_3d(jf_vsp)%fnow(:,:,1) * vmask(:,:,1)    ! v-velocity
         ELSE
            IF( .NOT. ln_linssh ) THEN
               e3t_m(:,:) = sf_ssm_2d(jf_e3t)%fnow(:,:,1) * tmask(:,:,1) ! vertical scale factor
            ELSE
               e3t_m(:,:) = e3t_0(:,:,1)                                 ! vertical scale factor
            ENDIF
            ssu_m(:,:) = sf_ssm_2d(jf_usp)%fnow(:,:,1) * umask(:,:,1)    ! u-velocity
            ssv_m(:,:) = sf_ssm_2d(jf_vsp)%fnow(:,:,1) * vmask(:,:,1)    ! v-velocity
         ENDIF
         !
         sst_m(:,:) = sf_ssm_2d(jf_tem)%fnow(:,:,1) * tmask(:,:,1)    ! temperature
         sss_m(:,:) = sf_ssm_2d(jf_sal)%fnow(:,:,1) * tmask(:,:,1)    ! salinity
         ssh_m(:,:) = sf_ssm_2d(jf_ssh)%fnow(:,:,1) * tmask(:,:,1)    ! sea surface height
         IF( ln_read_frq ) THEN
            frq_m(:,:) = sf_ssm_2d(jf_frq)%fnow(:,:,1) * tmask(:,:,1) ! solar penetration
         ELSE
            frq_m(:,:) = 1._wp
         ENDIF
      ELSE
         sss_m(:,:) = 35._wp                             ! =35. to obtain a physical value for the freezing point
         CALL eos_fzp( sss_m(:,:), sst_m(:,:) )          ! sst_m is set at the freezing point
         ssu_m(:,:) = 0._wp
         ssv_m(:,:) = 0._wp
         ssh_m(:,:) = 0._wp
         IF( .NOT. ln_linssh ) e3t_m(:,:) = e3t_0(:,:,1) !clem: necessary at least for sas2D
         frq_m(:,:) = 1._wp                              !              - -
         ssh  (:,:,Kmm) = 0._wp                              !              - -
      ENDIF

      IF ( nn_ice == 1 ) THEN
         ts(:,:,1,jp_tem,Kmm) = sst_m(:,:)
         ts(:,:,1,jp_sal,Kmm) = sss_m(:,:)
         ts(:,:,1,jp_tem,Kbb) = sst_m(:,:)
         ts(:,:,1,jp_sal,Kbb) = sss_m(:,:)
      ENDIF
      uu (:,:,1,Kbb) = ssu_m(:,:)
      vv (:,:,1,Kbb) = ssv_m(:,:)

      IF(sn_cfctl%l_prtctl) THEN            ! print control
         CALL prt_ctl(tab2d_1=sst_m, clinfo1=' sst_m   - : ', mask1=tmask   )
         CALL prt_ctl(tab2d_1=sss_m, clinfo1=' sss_m   - : ', mask1=tmask   )
         CALL prt_ctl(tab2d_1=ssu_m, clinfo1=' ssu_m   - : ', mask1=umask   )
         CALL prt_ctl(tab2d_1=ssv_m, clinfo1=' ssv_m   - : ', mask1=vmask   )
         CALL prt_ctl(tab2d_1=ssh_m, clinfo1=' ssh_m   - : ', mask1=tmask   )
         IF( .NOT.ln_linssh )   CALL prt_ctl(tab2d_1=ssh_m, clinfo1=' e3t_m   - : ', mask1=tmask   )
         IF( ln_read_frq    )   CALL prt_ctl(tab2d_1=frq_m, clinfo1=' frq_m   - : ', mask1=tmask   )
      ENDIF
      !
      IF( l_initdone ) THEN          !   Mean value at each nn_fsbc time-step   !
         CALL iom_put( 'ssu_m', ssu_m )
         CALL iom_put( 'ssv_m', ssv_m )
         CALL iom_put( 'sst_m', sst_m )
         CALL iom_put( 'sss_m', sss_m )
         CALL iom_put( 'ssh_m', ssh_m )
         IF( .NOT.ln_linssh )   CALL iom_put( 'e3t_m', e3t_m )
         IF( ln_read_frq    )   CALL iom_put( 'frq_m', frq_m )
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop( 'sbc_ssm')
      !
   END SUBROUTINE sbc_ssm


   SUBROUTINE sbc_ssm_init( Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_ssm_init  ***
      !!
      !! ** Purpose :   Initialisation of sea surface mean data
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Kmm   ! ocean time level indices
      ! (not needed for SAS but needed to keep a consistent interface in sbcmod.F90)
      INTEGER  :: ierr, ierr0, ierr1, ierr2, ierr3   ! return error code
      INTEGER  :: ifpr                               ! dummy loop indice
      INTEGER  :: inum, idv, idimv, jpm              ! local integer
      INTEGER  ::   ios                              ! Local integer output status for namelist read
      !!
      CHARACTER(len=100)                     ::  cn_dir       ! Root directory for location of core files
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::  slf_3d       ! array of namelist information on the fields to read
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::  slf_2d       ! array of namelist information on the fields to read
      TYPE(FLD_N) ::   sn_tem, sn_sal                     ! information about the fields to be read
      TYPE(FLD_N) ::   sn_usp, sn_vsp
      TYPE(FLD_N) ::   sn_ssh, sn_e3t, sn_frq
      !!
      NAMELIST/namsbc_sas/ l_sasread, cn_dir, ln_3d_uve, ln_read_frq,   &
         &                 sn_tem, sn_sal, sn_usp, sn_vsp, sn_ssh, sn_e3t, sn_frq
      !!----------------------------------------------------------------------
      !
      IF( ln_rstart .AND. nn_components == jp_iam_sas )   RETURN
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_ssm_init : sea surface mean data initialisation '
         WRITE(numout,*) '~~~~~~~~~~~~ '
      ENDIF
      !
      READ  ( numnam_ref, namsbc_sas, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_sas in reference namelist' )
      READ  ( numnam_cfg, namsbc_sas, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_sas in configuration namelist' )
      IF(lwm) WRITE ( numond, namsbc_sas )
      !
      IF(lwp) THEN                              ! Control print
         WRITE(numout,*) '   Namelist namsbc_sas'
         WRITE(numout,*) '      Initialisation using an input file                                 l_sasread   = ', l_sasread
         WRITE(numout,*) '      Are we supplying a 3D u,v and e3 field                             ln_3d_uve   = ', ln_3d_uve
         WRITE(numout,*) '      Are we reading frq (fraction of qsr absorbed in the 1st T level)   ln_read_frq = ', ln_read_frq
      ENDIF
      !
      !! switch off stuff that isn't sensible with a standalone module
      !! note that we need sbc_ssm called first in sbc
      !
      IF( ln_apr_dyn ) THEN
         IF( lwp ) WRITE(numout,*) '         ==>>>   No atmospheric gradient needed with StandAlone Surface scheme'
         ln_apr_dyn = .FALSE.
      ENDIF
      IF( ln_rnf ) THEN
         IF( lwp ) WRITE(numout,*) '         ==>>>   No runoff needed with StandAlone Surface scheme'
         ln_rnf = .FALSE.
      ENDIF
      IF( ln_ssr ) THEN
         IF( lwp ) WRITE(numout,*) '         ==>>>   No surface relaxation needed with StandAlone Surface scheme'
         ln_ssr = .FALSE.
      ENDIF
      IF( nn_fwb > 0 ) THEN
         IF( lwp ) WRITE(numout,*) '         ==>>>   No freshwater budget adjustment needed with StandAlone Surface scheme'
         nn_fwb = 0
      ENDIF

      !
      IF( l_sasread ) THEN                       ! store namelist information in an array
         !
         !! following code is a bit messy, but distinguishes between when u,v are 3d arrays and
         !! when we have other 3d arrays that we need to read in
         !! so if a new field is added i.e. jf_new, just give it the next integer in sequence
         !! for the corresponding dimension (currently if ln_3d_uve is true, 4 for 2d and 3 for 3d,
         !! alternatively if ln_3d_uve is false, 6 for 2d and 1 for 3d), reset nfld_3d, nfld_2d,
         !! and the rest of the logic should still work
         !
         jf_tem = 1   ;   jf_ssh = 3   ! default 2D fields index
         jf_sal = 2   ;   jf_frq = 4   !
         !
         IF( ln_3d_uve ) THEN
            jf_usp = 1   ;   jf_vsp = 2   ;   jf_e3t = 3     ! define 3D fields index
            nfld_3d  = 2 + COUNT( (/.NOT.ln_linssh/) )       ! number of 3D fields to read
            nfld_2d  = 3 + COUNT( (/ln_read_frq/) )          ! number of 2D fields to read
         ELSE
            jf_usp = 4   ;   jf_e3t = 6                      ! update 2D fields index
            jf_vsp = 5   ;   jf_frq = 6 + COUNT( (/.NOT.ln_linssh/) )
            !
            nfld_3d  = 0                                     ! no 3D fields to read
            nfld_2d  = 5 + COUNT( (/.NOT.ln_linssh/) ) + COUNT( (/ln_read_frq/) )    ! number of 2D fields to read
         ENDIF
         !
         IF( nfld_3d > 0 ) THEN
            ALLOCATE( slf_3d(nfld_3d), STAT=ierr )         ! set slf structure
            IF( ierr > 0 ) THEN
               CALL ctl_stop( 'sbc_ssm_init: unable to allocate slf 3d structure' )   ;   RETURN
            ENDIF
            slf_3d(jf_usp) = sn_usp
            slf_3d(jf_vsp) = sn_vsp
            IF( .NOT.ln_linssh )   slf_3d(jf_e3t) = sn_e3t
         ENDIF
         !
         IF( nfld_2d > 0 ) THEN
            ALLOCATE( slf_2d(nfld_2d), STAT=ierr )         ! set slf structure
            IF( ierr > 0 ) THEN
               CALL ctl_stop( 'sbc_ssm_init: unable to allocate slf 2d structure' )   ;   RETURN
            ENDIF
            slf_2d(jf_tem) = sn_tem   ;   slf_2d(jf_sal) = sn_sal   ;   slf_2d(jf_ssh) = sn_ssh
            IF( ln_read_frq )   slf_2d(jf_frq) = sn_frq
            IF( .NOT. ln_3d_uve ) THEN
               slf_2d(jf_usp) = sn_usp ; slf_2d(jf_vsp) = sn_vsp
               IF( .NOT.ln_linssh )   slf_2d(jf_e3t) = sn_e3t
            ENDIF
         ENDIF
         !
         ierr1 = 0    ! default definition if slf_?d(ifpr)%ln_tint = .false.
         IF( nfld_3d > 0 ) THEN
            ALLOCATE( sf_ssm_3d(nfld_3d), STAT=ierr )         ! set sf structure
            IF( ierr > 0 ) THEN
               CALL ctl_stop( 'sbc_ssm_init: unable to allocate sf structure' )   ;   RETURN
            ENDIF
            DO ifpr = 1, nfld_3d
               ALLOCATE( sf_ssm_3d(ifpr)%fnow(jpi,jpj,jpk)    , STAT=ierr0 )
               IF( slf_3d(ifpr)%ln_tint )   ALLOCATE( sf_ssm_3d(ifpr)%fdta(jpi,jpj,jpk,2)  , STAT=ierr1 )
               IF( ierr0 + ierr1 > 0 ) THEN
                  CALL ctl_stop( 'sbc_ssm_init : unable to allocate sf_ssm_3d array structure' )   ;   RETURN
               ENDIF
            END DO
            !                                         ! fill sf with slf_i and control print
            CALL fld_fill( sf_ssm_3d, slf_3d, cn_dir, 'sbc_ssm_init', '3D Data in file', 'namsbc_ssm' )
         ENDIF
         !
         IF( nfld_2d > 0 ) THEN
            ALLOCATE( sf_ssm_2d(nfld_2d), STAT=ierr )         ! set sf structure
            IF( ierr > 0 ) THEN
               CALL ctl_stop( 'sbc_ssm_init: unable to allocate sf 2d structure' )   ;   RETURN
            ENDIF
            DO ifpr = 1, nfld_2d
               ALLOCATE( sf_ssm_2d(ifpr)%fnow(jpi,jpj,1)    , STAT=ierr0 )
               IF( slf_2d(ifpr)%ln_tint )   ALLOCATE( sf_ssm_2d(ifpr)%fdta(jpi,jpj,1,2)  , STAT=ierr1 )
               IF( ierr0 + ierr1 > 0 ) THEN
                  CALL ctl_stop( 'sbc_ssm_init : unable to allocate sf_ssm_2d array structure' )   ;   RETURN
               ENDIF
            END DO
            !
            CALL fld_fill( sf_ssm_2d, slf_2d, cn_dir, 'sbc_ssm_init', '2D Data in file', 'namsbc_ssm' )
         ENDIF
         !
         IF( nfld_3d > 0 )   DEALLOCATE( slf_3d, STAT=ierr )
         IF( nfld_2d > 0 )   DEALLOCATE( slf_2d, STAT=ierr )
         !
      ENDIF
      !
      CALL sbc_ssm( nit000, Kbb, Kmm )   ! need to define ss?_m arrays used in iceistate
      l_initdone = .TRUE.
      !
   END SUBROUTINE sbc_ssm_init

   !!======================================================================
END MODULE sbcssm
