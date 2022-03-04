MODULE dtatsd
   !!======================================================================
   !!                     ***  MODULE  dtatsd  ***
   !! Ocean data  :  read ocean Temperature & Salinity Data from gridded data
   !!======================================================================
   !! History :  OPA  ! 1991-03  ()  Original code
   !!             -   ! 1992-07  (M. Imbard)
   !!            8.0  ! 1999-10  (M.A. Foujols, M. Imbard)  NetCDF FORMAT
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            3.3  ! 2010-10  (C. Bricaud, S. Masson)  use of fldread
   !!            3.4  ! 2010-11  (G. Madec, C. Ethe) Merge of dtatem and dtasal + remove CPP keys
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dta_tsd      : read and time interpolated ocean Temperature & Salinity Data
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE phycst          ! physical constants
   USE dom_oce         ! ocean space and time domain
   USE domtile
   USE fldread         ! read input fields
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dta_tsd_init   ! called by opa.F90
   PUBLIC   dta_tsd        ! called by istate.F90 and tradmp.90

   !                                  !!* namtsd  namelist : Temperature & Salinity Data *
   LOGICAL , PUBLIC ::   ln_tsd_init   !: T & S data flag
   LOGICAL , PUBLIC ::   ln_tsd_dmp    !: internal damping toward input data flag

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_tsdini ! structure of input SST (file informations, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_tsddmp ! structure of input SST (file informations, fields read)

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dtatsd.F90 10213 2018-10-23 14:40:09Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dta_tsd_init( ld_tradmp )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_tsd_init  ***
      !!
      !! ** Purpose :   initialisation of T & S input data
      !!
      !! ** Method  : - Read namtsd namelist
      !!              - allocates T & S data structure
      !!----------------------------------------------------------------------
      LOGICAL, INTENT(in), OPTIONAL ::   ld_tradmp   ! force the initialization when tradp is used
      !
      INTEGER ::   ios, ierr0, ierr1, ierr2, ierr3   ! local integers
      !!
      CHARACTER(len=100)            ::   cn_dir          ! Root directory for location of ssr files
      TYPE(FLD_N), DIMENSION( jpts) ::   slf_i           ! array of namelist informations on the fields to read
      TYPE(FLD_N)                   ::   sn_tem, sn_sal
      TYPE(FLD_N)                   ::   sn_dmpt, sn_dmps
      !!
      NAMELIST/namtsd/   ln_tsd_init, ln_tsd_dmp, cn_dir, sn_tem, sn_sal, sn_dmpt, sn_dmps
      !!----------------------------------------------------------------------
      !
      !  Initialisation
      ierr0 = 0  ;  ierr1 = 0  ;  ierr2 = 0  ;  ierr3 = 0
      !
      READ  ( numnam_ref, namtsd, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtsd in reference namelist' )
      READ  ( numnam_cfg, namtsd, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtsd in configuration namelist' )
      IF(lwm) WRITE ( numond, namtsd )

      IF( PRESENT( ld_tradmp ) )   ln_tsd_dmp = .TRUE.     ! forces the initialization when tradmp is used

      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dta_tsd_init : Temperature & Salinity data '
         WRITE(numout,*) '~~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namtsd'
         WRITE(numout,*) '      Initialisation of ocean T & S with T &S input data   ln_tsd_init = ', ln_tsd_init
         WRITE(numout,*) '      damping of ocean T & S toward T &S input data        ln_tsd_dmp  = ', ln_tsd_dmp
         WRITE(numout,*)
         IF( .NOT.ln_tsd_init .AND. .NOT.ln_tsd_dmp ) THEN
            WRITE(numout,*)
            WRITE(numout,*) '   ===>>   T & S data not used'
         ENDIF
      ENDIF
      !
      IF( ln_rstart .AND. ln_tsd_init ) THEN
         CALL ctl_warn( 'dta_tsd_init: ocean restart and T & S data intialisation, ',   &
            &           'we keep the restart T & S values and set ln_tsd_init to FALSE' )
         ln_tsd_init = .FALSE.
      ENDIF
      !
      !                             ! allocate the arrays (if necessary)
      IF( ln_tsd_init ) THEN
         !
         ALLOCATE( sf_tsdini(jpts), STAT=ierr0 )
         IF( ierr0 > 0 ) THEN
            CALL ctl_stop( 'dta_tsd_init: unable to allocate sf_tsd structure' )   ;   RETURN
         ENDIF
         !
                                ALLOCATE( sf_tsdini(jp_tem)%fnow(jpi,jpj,jpk)   , STAT=ierr0 )
         IF( sn_tem%ln_tint )   ALLOCATE( sf_tsdini(jp_tem)%fdta(jpi,jpj,jpk,2) , STAT=ierr1 )
                                ALLOCATE( sf_tsdini(jp_sal)%fnow(jpi,jpj,jpk)   , STAT=ierr2 )
         IF( sn_sal%ln_tint )   ALLOCATE( sf_tsdini(jp_sal)%fdta(jpi,jpj,jpk,2) , STAT=ierr3 )
         !
         IF( ierr0 + ierr1 + ierr2 + ierr3 > 0 ) THEN
            CALL ctl_stop( 'dta_tsd : unable to allocate T & S data arrays' )   ;   RETURN
         ENDIF
         !                         ! fill sf_tsd with sn_tem & sn_sal and control print
         slf_i(jp_tem) = sn_tem   ;   slf_i(jp_sal) = sn_sal
         CALL fld_fill( sf_tsdini, slf_i, cn_dir, 'dta_tsd', 'Temperature & Salinity data', 'namtsd', no_print )
         !
      END IF

      IF( ln_tsd_dmp ) THEN
         !
         ALLOCATE( sf_tsddmp(jpts), STAT=ierr0 )
         IF( ierr0 > 0 ) THEN
            CALL ctl_stop( 'dta_tsd_init: unable to allocate sf_tsddmp structure' )   ;   RETURN
         ENDIF
         !
         ! dmp file
                                 ALLOCATE( sf_tsddmp(jp_tem)%fnow(jpi,jpj,jpk)   , STAT=ierr0 )
         IF( sn_dmpt%ln_tint )   ALLOCATE( sf_tsddmp(jp_tem)%fdta(jpi,jpj,jpk,2) , STAT=ierr1 )
                                 ALLOCATE( sf_tsddmp(jp_sal)%fnow(jpi,jpj,jpk)   , STAT=ierr2 )
         IF( sn_dmps%ln_tint )   ALLOCATE( sf_tsddmp(jp_sal)%fdta(jpi,jpj,jpk,2) , STAT=ierr3 )
         !
         IF( ierr0 + ierr1 + ierr2 + ierr3 > 0 ) THEN
            CALL ctl_stop( 'dta_tsd : unable to allocate T & S dmp data arrays' )   ;   RETURN
         ENDIF
         !
         !                         ! fill sf_tsd with sn_tem & sn_sal and control print
         slf_i(jp_tem) = sn_dmpt   ;   slf_i(jp_sal) = sn_dmps
         CALL fld_fill( sf_tsddmp, slf_i, cn_dir, 'dta_tsd', 'Temperature & Salinity dmp data', 'namtsd', no_print )
         !
      ENDIF
      !
   END SUBROUTINE dta_tsd_init


   SUBROUTINE dta_tsd( kt, cddta, ptsd )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_tsd  ***
      !!
      !! ** Purpose :   provides T and S data at kt
      !!
      !! ** Method  : - call fldread routine
      !!              - ORCA_R2: add some hand made alteration to read data
      !!              - 'key_orca_lev10' interpolates on 10 times more levels
      !!              - s- or mixed z-s coordinate: vertical interpolation on model mesh
      !!              - ln_tsd_dmp=F: deallocates the T-S data structure
      !!                as T-S data are no are used
      !!
      !! ** Action  :   ptsd   T-S data on medl mesh and interpolated at time-step kt
      !!----------------------------------------------------------------------
      INTEGER                          , INTENT(in   ) ::   kt     ! ocean time-step
      CHARACTER(LEN=3)                 , INTENT(in   ) ::   cddta  ! dmp or ini
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,jpts), INTENT(  out) ::   ptsd   ! T & S data
      !
      INTEGER ::   ji, jj, jk, jl, jkk   ! dummy loop indicies
      INTEGER ::   ik, il0, il1, ii0, ii1, ij0, ij1   ! local integers
      REAL(wp)::   zl, zi                             ! local scalars
      REAL(wp), DIMENSION(jpk) ::  ztp, zsp   ! 1D workspace
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                                         ! Do only for the full domain
         IF( ln_tile ) CALL dom_tile_stop( ldhold=.TRUE. )             ! Use full domain

         SELECT CASE(cddta)
         CASE('ini')
            CALL fld_read( kt, 1, sf_tsdini ) !==   read T & S data at kt time step   ==!
         CASE('dmp')
            CALL fld_read( kt, 1, sf_tsddmp ) !==   read T & S data at kt time step   ==!
         CASE DEFAULT
            CALL ctl_stop('STOP', 'dta_tsd: cddta case unknown')
         END SELECT

         IF( ln_tile ) CALL dom_tile_start( ldhold=.TRUE. )            ! Revert to tile domain
      ENDIF
      !
      SELECT CASE(cddta)
      CASE('ini')
         DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpk )
            ptsd(ji,jj,jk,jp_tem) = sf_tsdini(jp_tem)%fnow(ji,jj,jk)    ! NO mask
            ptsd(ji,jj,jk,jp_sal) = sf_tsdini(jp_sal)%fnow(ji,jj,jk)
         END_3D
      CASE('dmp')
         DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpk )
            ptsd(ji,jj,jk,jp_tem) = sf_tsddmp(jp_tem)%fnow(ji,jj,jk)    ! NO mask
            ptsd(ji,jj,jk,jp_sal) = sf_tsddmp(jp_sal)%fnow(ji,jj,jk)
         END_3D
      CASE DEFAULT
         CALL ctl_stop('STOP', 'dta_tsd: cddta case unknown')
      END SELECT
      !
      IF( ln_sco ) THEN                   !==   s- or mixed s-zps-coordinate   ==!
         !
         IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
            IF( kt == nit000 .AND. lwp )THEN
               WRITE(numout,*)
               WRITE(numout,*) 'dta_tsd: interpolates T & S data onto the s- or mixed s-z-coordinate mesh'
            ENDIF
         ENDIF
         !
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )                  ! vertical interpolation of T & S
            DO jk = 1, jpk                        ! determines the intepolated T-S profiles at each (i,j) points
               zl = gdept_0(ji,jj,jk)
               IF(     zl < gdept_1d(1  ) ) THEN          ! above the first level of data
                  ztp(jk) =  ptsd(ji,jj,1    ,jp_tem)
                  zsp(jk) =  ptsd(ji,jj,1    ,jp_sal)
               ELSEIF( zl > gdept_1d(jpk) ) THEN          ! below the last level of data
                  ztp(jk) =  ptsd(ji,jj,jpkm1,jp_tem)
                  zsp(jk) =  ptsd(ji,jj,jpkm1,jp_sal)
               ELSE                                      ! inbetween : vertical interpolation between jkk & jkk+1
                  DO jkk = 1, jpkm1                                  ! when  gdept(jkk) < zl < gdept(jkk+1)
                     IF( (zl-gdept_1d(jkk)) * (zl-gdept_1d(jkk+1)) <= 0._wp ) THEN
                        zi = ( zl - gdept_1d(jkk) ) / (gdept_1d(jkk+1)-gdept_1d(jkk))
                        ztp(jk) = ptsd(ji,jj,jkk,jp_tem) + ( ptsd(ji,jj,jkk+1,jp_tem) - ptsd(ji,jj,jkk,jp_tem) ) * zi
                        zsp(jk) = ptsd(ji,jj,jkk,jp_sal) + ( ptsd(ji,jj,jkk+1,jp_sal) - ptsd(ji,jj,jkk,jp_sal) ) * zi
                     ENDIF
                  END DO
               ENDIF
            END DO
            DO jk = 1, jpkm1
               ptsd(ji,jj,jk,jp_tem) = ztp(jk) * tmask(ji,jj,jk)     ! mask required for mixed zps-s-coord
               ptsd(ji,jj,jk,jp_sal) = zsp(jk) * tmask(ji,jj,jk)
            END DO
            ptsd(ji,jj,jpk,jp_tem) = 0._wp
            ptsd(ji,jj,jpk,jp_sal) = 0._wp
         END_2D
         !
      ELSE                                !==   z- or zps- coordinate   ==!
         !
         DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpk )
            ptsd(ji,jj,jk,jp_tem) = ptsd(ji,jj,jk,jp_tem) * tmask(ji,jj,jk)    ! Mask
            ptsd(ji,jj,jk,jp_sal) = ptsd(ji,jj,jk,jp_sal) * tmask(ji,jj,jk)
         END_3D
         !
         IF( ln_zps ) THEN                      ! zps-coordinate (partial steps) interpolation at the last ocean level
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               ik = mbkt(ji,jj)
               IF( ik > 1 ) THEN
                  zl = ( gdept_1d(ik) - gdept_0(ji,jj,ik) ) / ( gdept_1d(ik) - gdept_1d(ik-1) )
                  ptsd(ji,jj,ik,jp_tem) = (1.-zl) * ptsd(ji,jj,ik,jp_tem) + zl * ptsd(ji,jj,ik-1,jp_tem)
                  ptsd(ji,jj,ik,jp_sal) = (1.-zl) * ptsd(ji,jj,ik,jp_sal) + zl * ptsd(ji,jj,ik-1,jp_sal)
               ENDIF
               ik = mikt(ji,jj)
               IF( ik > 1 ) THEN
                  zl = ( gdept_0(ji,jj,ik) - gdept_1d(ik) ) / ( gdept_1d(ik+1) - gdept_1d(ik) )
                  ptsd(ji,jj,ik,jp_tem) = (1.-zl) * ptsd(ji,jj,ik,jp_tem) + zl * ptsd(ji,jj,ik+1,jp_tem)
                  ptsd(ji,jj,ik,jp_sal) = (1.-zl) * ptsd(ji,jj,ik,jp_sal) + zl * ptsd(ji,jj,ik+1,jp_sal)
               END IF
            END_2D
         ENDIF
         !
      ENDIF
      !
      SELECT CASE(cddta)
      CASE('ini') 
         !                        !==   deallocate T & S structure   ==! 
         !                                              (data used only for initialisation)
         IF(lwp) WRITE(numout,*) 'dta_tsd: deallocte T & S arrays as they are only use to initialize the run'
                                        DEALLOCATE( sf_tsdini(jp_tem)%fnow )     ! T arrays in the structure
         IF( sf_tsdini(jp_tem)%ln_tint )   DEALLOCATE( sf_tsdini(jp_tem)%fdta )
                                        DEALLOCATE( sf_tsdini(jp_sal)%fnow )     ! S arrays in the structure
         IF( sf_tsdini(jp_sal)%ln_tint )   DEALLOCATE( sf_tsdini(jp_sal)%fdta )
                                        DEALLOCATE( sf_tsdini              )     ! the structure itself
         !
      END SELECT
      !
   END SUBROUTINE dta_tsd

   !!======================================================================
END MODULE dtatsd
