MODULE stpctl
   !!======================================================================
   !!                       ***  MODULE  stpctl  ***
   !! Ocean run control :  gross check of the ocean time stepping
   !!======================================================================
   !! History :  OPA  ! 1991-03  (G. Madec) Original code
   !!            6.0  ! 1992-06  (M. Imbard)
   !!            8.0  ! 1997-06  (A.M. Treguier)
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            2.0  ! 2009-07  (G. Madec)  Add statistic for time-spliting
   !!            3.7  ! 2016-09  (G. Madec)  Remove solver
   !!            4.0  ! 2017-04  (G. Madec)  regroup global communications
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   stp_ctl      : Control the run
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE c1d             ! 1D vertical configuration
   USE diawri          ! Standard run outputs       (dia_wri_state routine)
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! distributed memory computing
   USE zdf_oce ,  ONLY : ln_zad_Aimp       ! ocean vertical physics variables
   USE wet_dry,   ONLY : ll_wd, ssh_ref    ! reference depth for negative bathy

   USE netcdf          ! NetCDF library
   IMPLICIT NONE
   PRIVATE

   PUBLIC stp_ctl           ! routine called by step.F90

   INTEGER  ::   idrun, idtime, idssh, idu, ids1, ids2, idt1, idt2, idc1, idw1, istatus
   LOGICAL  ::   lsomeoce
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: stpctl.F90 10572 2019-01-24 15:37:13Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE stp_ctl( kt, kindic )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE stp_ctl  ***
      !!                     
      !! ** Purpose :   Control the run
      !!
      !! ** Method  : - Save the time step in numstp
      !!              - Print it each 50 time steps
      !!              - Stop the run IF problem encountered by setting indic=-3
      !!                Problems checked: |ssh| maximum larger than 10 m
      !!                                  |U|   maximum larger than 10 m/s 
      !!                                  negative sea surface salinity
      !!
      !! ** Actions :   "time.step" file = last ocean time-step
      !!                "run.stat"  file = run statistics
      !!                nstop indicator sheared among all local domain (lk_mpp=T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER, INTENT(inout) ::   kindic   ! error indicator
      !!
      INTEGER                ::   ji, jj, jk          ! dummy loop indices
      INTEGER, DIMENSION(2)  ::   ih                  ! min/max loc indices
      INTEGER, DIMENSION(3)  ::   iu, is1, is2        ! min/max loc indices
      REAL(wp)               ::   zzz                 ! local real 
      REAL(wp), DIMENSION(9) ::   zmax
      LOGICAL                ::   ll_wrtstp, ll_colruns, ll_wrtruns
      CHARACTER(len=20) :: clname
      !!----------------------------------------------------------------------
      !
      ll_wrtstp  = ( MOD( kt, sn_cfctl%ptimincr ) == 0 ) .OR. ( kt == nitend )
      ll_colruns = ll_wrtstp .AND. ( ln_ctl .OR. sn_cfctl%l_runstat )
      ll_wrtruns = ll_colruns .AND. lwm
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'stp_ctl : time-stepping control'
         WRITE(numout,*) '~~~~~~~'
         !                                ! open time.step file
         IF( lwm ) CALL ctl_opn( numstp, 'time.step', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
         !                                ! open run.stat file(s) at start whatever
         !                                ! the value of sn_cfctl%ptimincr
         IF( lwm .AND. ( ln_ctl .OR. sn_cfctl%l_runstat ) ) THEN
            CALL ctl_opn( numrun, 'run.stat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
            clname = 'run.stat.nc'
            IF( .NOT. Agrif_Root() )   clname = TRIM(Agrif_CFixed())//"_"//TRIM(clname)
            istatus = NF90_CREATE( TRIM(clname), NF90_CLOBBER, idrun )
            istatus = NF90_DEF_DIM( idrun, 'time', NF90_UNLIMITED, idtime )
            istatus = NF90_DEF_VAR( idrun, 'abs_ssh_max', NF90_DOUBLE, (/ idtime /), idssh )
            istatus = NF90_DEF_VAR( idrun,   'abs_u_max', NF90_DOUBLE, (/ idtime /), idu   )
            istatus = NF90_DEF_VAR( idrun,       's_min', NF90_DOUBLE, (/ idtime /), ids1  )
            istatus = NF90_DEF_VAR( idrun,       's_max', NF90_DOUBLE, (/ idtime /), ids2  )
            istatus = NF90_DEF_VAR( idrun,       't_min', NF90_DOUBLE, (/ idtime /), idt1  )
            istatus = NF90_DEF_VAR( idrun,       't_max', NF90_DOUBLE, (/ idtime /), idt2  )
            IF( ln_zad_Aimp ) THEN
               istatus = NF90_DEF_VAR( idrun,   'abs_wi_max', NF90_DOUBLE, (/ idtime /), idw1  )
               istatus = NF90_DEF_VAR( idrun,       'Cu_max', NF90_DOUBLE, (/ idtime /), idc1  )
            ENDIF
            istatus = NF90_ENDDEF(idrun)
            zmax(8:9) = 0._wp    ! initialise to zero in case ln_zad_Aimp option is not in use
         ENDIF
      ENDIF
      IF( kt == nit000 )   lsomeoce = COUNT( ssmask(:,:) == 1._wp ) > 0
      !
      IF(lwm .AND. ll_wrtstp) THEN        !==  current time step  ==!   ("time.step" file)
         WRITE ( numstp, '(1x, i8)' )   kt
         REWIND( numstp )
      ENDIF
      !
      !                                   !==  test of extrema  ==!
      IF( ll_wd ) THEN
         zmax(1) = MAXVAL(  ABS( sshn(:,:) + ssh_ref*tmask(:,:,1) )  )        ! ssh max 
      ELSE
         zmax(1) = MAXVAL(  ABS( sshn(:,:) )  )                               ! ssh max
      ENDIF
      zmax(2) = MAXVAL(  ABS( un(:,:,:) )  )                                  ! velocity max (zonal only)
      zmax(3) = MAXVAL( -tsn(:,:,:,jp_sal) , mask = tmask(:,:,:) == 1._wp )   ! minus salinity max
      zmax(4) = MAXVAL(  tsn(:,:,:,jp_sal) , mask = tmask(:,:,:) == 1._wp )   !       salinity max
      zmax(5) = MAXVAL( -tsn(:,:,:,jp_tem) , mask = tmask(:,:,:) == 1._wp )   ! minus temperature max
      zmax(6) = MAXVAL(  tsn(:,:,:,jp_tem) , mask = tmask(:,:,:) == 1._wp )   !       temperature max
      zmax(7) = REAL( nstop , wp )                                            ! stop indicator
      IF( ln_zad_Aimp ) THEN
         zmax(8) = MAXVAL(  ABS( wi(:,:,:) ) , mask = wmask(:,:,:) == 1._wp ) ! implicit vertical vel. max
         zmax(9) = MAXVAL(   Cu_adv(:,:,:)   , mask = tmask(:,:,:) == 1._wp ) !       cell Courant no. max
      ENDIF
      !
      IF( ll_colruns ) THEN
         CALL mpp_max( "stpctl", zmax )          ! max over the global domain
         nstop = NINT( zmax(7) )                 ! nstop indicator sheared among all local domains
      ENDIF
      !                                   !==  run statistics  ==!   ("run.stat" files)
      IF( ll_wrtruns ) THEN
         WRITE(numrun,9500) kt, zmax(1), zmax(2), -zmax(3), zmax(4)
         istatus = NF90_PUT_VAR( idrun, idssh, (/ zmax(1)/), (/kt/), (/1/) )
         istatus = NF90_PUT_VAR( idrun,   idu, (/ zmax(2)/), (/kt/), (/1/) )
         istatus = NF90_PUT_VAR( idrun,  ids1, (/-zmax(3)/), (/kt/), (/1/) )
         istatus = NF90_PUT_VAR( idrun,  ids2, (/ zmax(4)/), (/kt/), (/1/) )
         istatus = NF90_PUT_VAR( idrun,  idt1, (/-zmax(5)/), (/kt/), (/1/) )
         istatus = NF90_PUT_VAR( idrun,  idt2, (/ zmax(6)/), (/kt/), (/1/) )
         IF( ln_zad_Aimp ) THEN
            istatus = NF90_PUT_VAR( idrun,  idw1, (/ zmax(8)/), (/kt/), (/1/) )
            istatus = NF90_PUT_VAR( idrun,  idc1, (/ zmax(9)/), (/kt/), (/1/) )
         ENDIF
         IF( MOD( kt , 100 ) == 0 ) istatus = NF90_SYNC(idrun)
         IF( kt == nitend         ) istatus = NF90_CLOSE(idrun)
      END IF
      !                                   !==  error handling  ==!
      IF( ( ln_ctl .OR. lsomeoce ) .AND. (   &             ! have use mpp_max (because ln_ctl=.T.) or contains some ocean points
         &  zmax(1) >   50._wp .OR.   &                    ! too large sea surface height ( > 50 m )
         &  zmax(2) >   20._wp .OR.   &                    ! too large velocity ( > 20 m/s)
!!$         &  zmax(3) >=   0._wp .OR.   &                    ! negative or zero sea surface salinity
!!$         &  zmax(4) >= 100._wp .OR.   &                    ! too large sea surface salinity ( > 100 )
!!$         &  zmax(4) <    0._wp .OR.   &                    ! too large sea surface salinity (keep this line for sea-ice)
         &  ISNAN( zmax(1) + zmax(2) + zmax(3) ) ) ) THEN   ! NaN encounter in the tests
         IF( lk_mpp .AND. ln_ctl ) THEN
            CALL mpp_maxloc( 'stpctl', ABS(sshn)        , ssmask(:,:)  , zzz, ih  )
            CALL mpp_maxloc( 'stpctl', ABS(un)          , umask (:,:,:), zzz, iu  )
            CALL mpp_minloc( 'stpctl', tsn(:,:,:,jp_sal), tmask (:,:,:), zzz, is1 )
            CALL mpp_maxloc( 'stpctl', tsn(:,:,:,jp_sal), tmask (:,:,:), zzz, is2 )
         ELSE
            ih(:)  = MAXLOC( ABS( sshn(:,:)   )                              ) + (/ nimpp - 1, njmpp - 1    /)
            iu(:)  = MAXLOC( ABS( un  (:,:,:) )                              ) + (/ nimpp - 1, njmpp - 1, 0 /)
            is1(:) = MINLOC( tsn(:,:,:,jp_sal), mask = tmask(:,:,:) == 1._wp ) + (/ nimpp - 1, njmpp - 1, 0 /)
            is2(:) = MAXLOC( tsn(:,:,:,jp_sal), mask = tmask(:,:,:) == 1._wp ) + (/ nimpp - 1, njmpp - 1, 0 /)
         ENDIF
         
         WRITE(ctmp1,*) ' stp_ctl: |ssh| > 50 m  or  |U| > 20 m/s  or  NaN encounter in the tests'
         WRITE(ctmp2,9100) kt,   zmax(1), ih(1) , ih(2)
         WRITE(ctmp3,9200) kt,   zmax(2), iu(1) , iu(2) , iu(3)
         WRITE(ctmp4,9300) kt, - zmax(3), is1(1), is1(2), is1(3)
         WRITE(ctmp5,9400) kt,   zmax(4), is2(1), is2(2), is2(3)
         WRITE(ctmp6,*) '      ===> output of last computed fields in output.abort.nc file'
         
         CALL dia_wri_state( 'output.abort' )     ! create an output.abort file
         
         IF( .NOT. ln_ctl ) THEN
            WRITE(ctmp8,*) 'E R R O R message from sub-domain: ', narea
            CALL ctl_stop( 'STOP', ctmp1, ' ', ctmp8, ' ', ctmp2, ctmp3, ctmp4, ctmp5, ctmp6 )
         ELSE
            CALL ctl_stop( ctmp1, ' ', ctmp2, ctmp3, ctmp4, ctmp5, ' ', ctmp6, ' ' )
         ENDIF

         kindic = -3
         !
      ENDIF
      !
9100  FORMAT (' kt=',i8,'   |ssh| max: ',1pg11.4,', at  i j  : ',2i5)
9200  FORMAT (' kt=',i8,'   |U|   max: ',1pg11.4,', at  i j k: ',3i5)
9300  FORMAT (' kt=',i8,'   S     min: ',1pg11.4,', at  i j k: ',3i5)
9400  FORMAT (' kt=',i8,'   S     max: ',1pg11.4,', at  i j k: ',3i5)
9500  FORMAT(' it :', i8, '    |ssh|_max: ', D23.16, ' |U|_max: ', D23.16,' S_min: ', D23.16,' S_max: ', D23.16)
      !
   END SUBROUTINE stp_ctl

   !!======================================================================
END MODULE stpctl
