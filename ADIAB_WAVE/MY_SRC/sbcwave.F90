MODULE sbcwave
   !!======================================================================
   !!                       ***  MODULE  sbcwave  ***
   !! Wave module
   !!======================================================================
   !! History :  3.3  !  2011-09  (M. Adani)  Original code: Drag Coefficient
   !!         :  3.4  !  2012-10  (M. Adani)  Stokes Drift
   !!            3.6  !  2014-09  (E. Clementi,P. Oddo) New Stokes Drift Computation
   !!             -   !  2016-12  (G. Madec, E. Clementi) update Stoke drift computation
   !!                                                    + add sbc_wave_ini routine
   !!            4.2  !  2020-12  (G. Madec, E. Clementi) updates, new Stoke drift computation
   !!                                                    according to Couvelard et al.,2019
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_stokes    : calculate 3D Stokes-drift velocities
   !!   sbc_wave      : wave data from wave model: forced (netcdf files) or coupled mode
   !!   sbc_wave_init : initialisation fo surface waves
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE oce            ! ocean variables
   USE dom_oce        ! ocean domain variables
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE bdy_oce        ! open boundary condition variables
   USE domvvl         ! domain: variable volume layers
   USE usrdef_nam , ONLY: ln_STOKES_ADIAB   
   !
   USE iom            ! I/O manager library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE fldread        ! read input fields
   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_stokes      ! routine called in sbccpl
   PUBLIC   sbc_wave        ! routine called in sbcmod
   PUBLIC   sbc_wave_init   ! routine called in sbcmod

   ! Variables checking if the wave parameters are coupled (if not, they are read from file)
   LOGICAL, PUBLIC ::   cpl_hsig          = .FALSE.
   LOGICAL, PUBLIC ::   cpl_phioc         = .FALSE.
   LOGICAL, PUBLIC ::   cpl_sdrftx        = .FALSE.
   LOGICAL, PUBLIC ::   cpl_sdrfty        = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wper          = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wnum          = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wstrf         = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wdrag         = .FALSE.
   LOGICAL, PUBLIC ::   cpl_charn         = .FALSE.
   LOGICAL, PUBLIC ::   cpl_taw           = .FALSE.
   LOGICAL, PUBLIC ::   cpl_bhd           = .FALSE.
   LOGICAL, PUBLIC ::   cpl_tusd          = .FALSE.
   LOGICAL, PUBLIC ::   cpl_tvsd          = .FALSE.

   INTEGER ::   jpfld    ! number of files to read for stokes drift
   INTEGER ::   jp_usd   ! index of stokes drift  (i-component) (m/s)    at T-point
   INTEGER ::   jp_vsd   ! index of stokes drift  (j-component) (m/s)    at T-point
   INTEGER ::   jp_hsw   ! index of significant wave hight      (m)      at T-point
   INTEGER ::   jp_wmp   ! index of mean wave period            (s)      at T-point

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_cd      ! structure of input fields (file informations, fields read) Drag Coefficient
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_sd      ! structure of input fields (file informations, fields read) Stokes Drift
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_wn      ! structure of input fields (file informations, fields read) wave number for Qiao
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_tauoc   ! structure of input fields (file informations, fields read) normalized wave stress into the ocean

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   cdn_wave        !: Neutral drag coefficient at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   hsw             !: Significant Wave Height at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   wmp             !: Wave Mean Period at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   wnum            !: Wave Number at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tauoc_wave      !: stress reduction factor  at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tsd2d           !: Surface Stokes Drift module at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   div_sd          !: barotropic stokes drift divergence
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   ut0sd, vt0sd    !: surface Stokes drift velocities at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   usd, vsd, wsd   !: Stokes drift velocities at u-, v- & w-points, resp.u
!
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   charn           !: charnock coefficient at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tawx            !: Net wave-supported stress, u
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tawy            !: Net wave-supported stress, v
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   twox            !: wave-ocean momentum flux, u
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   twoy            !: wave-ocean momentum flux, v
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tauoc_wavex     !: stress reduction factor  at, u component
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tauoc_wavey     !: stress reduction factor  at, v component
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   phioc           !: tke flux from wave model
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   KZN2            !: Kz*N2
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   bhd_wave        !: Bernoulli head. wave induce pression
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tusd, tvsd      !: Stokes drift transport
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   ZMX             !: Kz*N2
   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcwave.F90 14433 2021-02-11 08:06:49Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_stokes( Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_stokes  ***
      !!
      !! ** Purpose :   compute the 3d Stokes Drift according to Breivik et al.,
      !!                2014 (DOI: 10.1175/JPO-D-14-0020.1)
      !!
      !! ** Method  : - Calculate the horizontal Stokes drift velocity (Breivik et al. 2014)
      !!              - Calculate its horizontal divergence
      !!              - Calculate the vertical Stokes drift velocity
      !!              - Calculate the barotropic Stokes drift divergence
      !!
      !! ** action  : - tsd2d         : module of the surface Stokes drift velocity
      !!              - usd, vsd, wsd : 3 components of the Stokes drift velocity
      !!              - div_sd        : barotropic Stokes drift divergence
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: Kmm ! ocean time level index
      INTEGER  ::   jj, ji, jk   ! dummy loop argument
      INTEGER  ::   ik           ! local integer
      REAL(wp) ::  ztransp, zfac, ztemp, zsp0, zsqrt, zbreiv16_w
      REAL(wp) ::  zdep_u, zdep_v, zkh_u, zkh_v, zda_u, zda_v, sdtrp
      REAL(wp) ::  zkd_u, zkd_v, zdep_ue, zdep_ve  ! variable for Stokes drift in SHALLOW/INTER WATER
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE ::   zk_t, zk_u, zk_v, zu0_sd, zv0_sd ! 2D workspace
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ze3divh, zInt_w                  ! 3D workspace
      !!---------------------------------------------------------------------
      !
      ALLOCATE( ze3divh(jpi,jpj,jpkm1) ) ! jpkm1 -> avoid lbc_lnk on jpk that is not defined
      ALLOCATE( zInt_w(jpi,jpj,jpk) )
      ALLOCATE( zk_t(jpi,jpj), zk_u(jpi,jpj), zk_v(jpi,jpj), zu0_sd(jpi,jpj), zv0_sd(jpi,jpj) )
      zk_t    (:,:) = 0._wp
      zk_u    (:,:) = 0._wp
      zk_v    (:,:) = 0._wp
      zu0_sd  (:,:) = 0._wp
      zv0_sd  (:,:) = 0._wp
      ze3divh (:,:,:) = 0._wp

      !
      ! select parameterization for the calculation of vertical Stokes drift
      ! exp. wave number at t-point
      IF( ln_breivikFV_2016 ) THEN
      ! Assumptions :  ut0sd and vt0sd are surface Stokes drift at T-points
      !                sdtrp is the norm of Stokes transport
      !
         zfac = 0.166666666667_wp
         DO_2D( 1, 1, 1, 1 ) ! In the deep-water limit we have ke = ||ust0||/( 6 * ||transport|| )
            zsp0          = SQRT( ut0sd(ji,jj)*ut0sd(ji,jj) + vt0sd(ji,jj)*vt0sd(ji,jj) ) !<-- norm of Surface Stokes drift
            tsd2d(ji,jj)  = zsp0
            IF( cpl_tusd .AND. cpl_tvsd ) THEN  !stokes transport is provided in coupled mode
               sdtrp      = SQRT( tusd(ji,jj)*tusd(ji,jj) + tvsd(ji,jj)*tvsd(ji,jj) )  !<-- norm of Surface Stokes drift transport
            ELSE
               ! Stokes drift transport estimated from Hs and Tmean
               sdtrp      = 2.0_wp * rpi / 16.0_wp *                             &
                   &        hsw(ji,jj)*hsw(ji,jj) / MAX( wmp(ji,jj), 0.0000001_wp )
            ENDIF
            zk_t (ji,jj)  = zfac * zsp0 / MAX ( sdtrp, 0.0000001_wp ) !<-- ke = ||ust0||/( 6 * ||transport|| )
         END_2D
      !# define zInt_w ze3divh
         DO_3D( 1, 1, 1, 1, 1, jpk ) ! Compute the primitive of Breivik 2016 function at W-points
            zfac             = - 2._wp * zk_t (ji,jj) * gdepw(ji,jj,jk,Kmm)  !<-- zfac should be negative definite
            ztemp            = EXP ( zfac )
            zsqrt            = SQRT( -zfac )
            zbreiv16_w       = ztemp - SQRT(rpi)*zsqrt*ERFC(zsqrt) !Eq. 16 Breivik 2016
            zInt_w(ji,jj,jk) = ztemp - 4._wp * zk_t (ji,jj) * gdepw(ji,jj,jk,Kmm) * zbreiv16_w
         END_3D
!
         DO jk = 1, jpkm1
            zfac = 0.166666666667_wp
            DO_2D( 1, 1, 1, 1 ) !++ Compute the FV Breivik 2016 function at T-points
               zsp0          = zfac / MAX(zk_t (ji,jj),0.0000001_wp)
               ztemp         = zInt_w(ji,jj,jk) - zInt_w(ji,jj,jk+1)
               zu0_sd(ji,jj) = ut0sd(ji,jj) * zsp0 * ztemp * tmask(ji,jj,jk)
               zv0_sd(ji,jj) = vt0sd(ji,jj) * zsp0 * ztemp * tmask(ji,jj,jk)
            END_2D
            DO_2D( 1, 0, 1, 0 ) ! ++ Interpolate at U/V points
               zfac          =  1.0_wp / e3u(ji  ,jj,jk,Kmm)
               usd(ji,jj,jk) =  0.5_wp * zfac * ( zu0_sd(ji,jj)+zu0_sd(ji+1,jj) ) * umask(ji,jj,jk)
               zfac          =  1.0_wp / e3v(ji  ,jj,jk,Kmm)
               vsd(ji,jj,jk) =  0.5_wp * zfac * ( zv0_sd(ji,jj)+zv0_sd(ji,jj+1) ) * vmask(ji,jj,jk)
            END_2D
         ENDDO
      !# undef zInt_w
      !
      ELSE IF (ln_STOKES_ADIAB) THEN
         DO_2D( 1, 0, 1, 0 )
       !  velocity at u- & v-points
           zk_u(ji,jj) = 0.5_wp * ( wnum(ji,jj) + wnum(ji+1,jj) )
           zk_v(ji,jj) = 0.5_wp * ( wnum(ji,jj) + wnum(ji,jj+1) )
            !
           zu0_sd(ji,jj) = 0.5_wp * ( ut0sd(ji,jj) + ut0sd(ji+1,jj) )
           zv0_sd(ji,jj) = 0.5_wp * ( vt0sd(ji,jj) + vt0sd(ji,jj+1) )
         END_2D

      !                       !==  horizontal Stokes Drift 3D velocity  ==!
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
            zdep_u = 0.5_wp * ( gdept(ji,jj,jk,Kmm) + gdept(ji+1,jj,jk,Kmm) )
            zdep_v = 0.5_wp * ( gdept(ji,jj,jk,Kmm) + gdept(ji,jj+1,jk,Kmm) )
            zdep_ue= 0.5_wp * ( gdept(ji,jj,jpkm1,Kmm) + gdept(ji+1,jj,jpkm1,Kmm) )
            zdep_ve= 0.5_wp * ( gdept(ji,jj,jpkm1,Kmm) + gdept(ji+1,jj,jpkm1,Kmm) )
            !
            zkh_u = zk_u(ji,jj) * zdep_u     ! k * depth
            zkh_v = zk_v(ji,jj) * zdep_v
            zkd_u = zk_u(ji,jj) * zdep_ue
            zkd_v = zk_v(ji,jj) * zdep_ve
                                ! Depth attenuation
            zda_u = COSH(-2.0_wp*zkh_u+2.0_wp*zkd_u)/COSH(2.0_wp*zkd_u)
            zda_v = COSH(-2.0_wp*zkh_v+2.0_wp*zkd_v)/COSH(2.0_wp*zkd_v)
            !
            usd(ji,jj,jk) = zda_u * zu0_sd(ji,jj) * umask(ji,jj,jk)
            vsd(ji,jj,jk) = zda_v * zv0_sd(ji,jj) * vmask(ji,jj,jk)
         END_3D
      ELSE
         zfac = 2.0_wp * rpi / 16.0_wp
         DO_2D( 1, 1, 1, 1 )
            ! Stokes drift velocity estimated from Hs and Tmean
            ztransp = zfac * hsw(ji,jj)*hsw(ji,jj) / MAX( wmp(ji,jj), 0.0000001_wp )
            ! Stokes surface speed
            tsd2d(ji,jj) = SQRT( ut0sd(ji,jj)*ut0sd(ji,jj) + vt0sd(ji,jj)*vt0sd(ji,jj))
            ! Wavenumber scale
            zk_t(ji,jj) = ABS( tsd2d(ji,jj) ) / MAX( ABS( 5.97_wp*ztransp ), 0.0000001_wp )
         END_2D
         DO_2D( 1, 0, 1, 0 )          ! exp. wave number & Stokes drift velocity at u- & v-points
            zk_u(ji,jj) = 0.5_wp * ( zk_t(ji,jj) + zk_t(ji+1,jj) )
            zk_v(ji,jj) = 0.5_wp * ( zk_t(ji,jj) + zk_t(ji,jj+1) )
            !
            zu0_sd(ji,jj) = 0.5_wp * ( ut0sd(ji,jj) + ut0sd(ji+1,jj) )
            zv0_sd(ji,jj) = 0.5_wp * ( vt0sd(ji,jj) + vt0sd(ji,jj+1) )
         END_2D

      !                       !==  horizontal Stokes Drift 3D velocity  ==!

         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
            zdep_u = 0.5_wp * ( gdept(ji,jj,jk,Kmm) + gdept(ji+1,jj,jk,Kmm) )
            zdep_v = 0.5_wp * ( gdept(ji,jj,jk,Kmm) + gdept(ji,jj+1,jk,Kmm) )
            !
            zkh_u = zk_u(ji,jj) * zdep_u     ! k * depth
            zkh_v = zk_v(ji,jj) * zdep_v
            !                                ! Depth attenuation
            zda_u = EXP( -2.0_wp*zkh_u ) / ( 1.0_wp + 8.0_wp*zkh_u )
            zda_v = EXP( -2.0_wp*zkh_v ) / ( 1.0_wp + 8.0_wp*zkh_v )
            !
            usd(ji,jj,jk) = zda_u * zu0_sd(ji,jj) * umask(ji,jj,jk)
            vsd(ji,jj,jk) = zda_v * zv0_sd(ji,jj) * vmask(ji,jj,jk)
         END_3D
      ENDIF

      CALL lbc_lnk( 'sbcwave', usd, 'U', -1.0_wp, vsd, 'V', -1.0_wp )

      !
      !                       !==  vertical Stokes Drift 3D velocity  ==!
      !
      DO_3D( 0, 1, 0, 1, 1, jpkm1 )    ! Horizontal e3*divergence
         ze3divh(ji,jj,jk) = (  e2u(ji  ,jj) * e3u(ji  ,jj,jk,Kmm) * usd(ji  ,jj,jk)    &
            &                 - e2u(ji-1,jj) * e3u(ji-1,jj,jk,Kmm) * usd(ji-1,jj,jk)    &
            &                 + e1v(ji,jj  ) * e3v(ji,jj  ,jk,Kmm) * vsd(ji,jj  ,jk)    &
            &                 - e1v(ji,jj-1) * e3v(ji,jj-1,jk,Kmm) * vsd(ji,jj-1,jk)  ) &
            &                * r1_e1e2t(ji,jj)
      END_3D
      !
      CALL lbc_lnk( 'sbcwave', ze3divh, 'T', 1.0_wp )
      !
      IF( ln_linssh ) THEN   ;   ik = 1   ! none zero velocity through the sea surface
      ELSE                   ;   ik = 2   ! w=0 at the surface (set one for all in sbc_wave_init)
      ENDIF
      DO jk = jpkm1, ik, -1          ! integrate from the bottom the hor. divergence (NB: at k=jpk w is always zero)
         wsd(:,:,jk) = wsd(:,:,jk+1) - ze3divh(:,:,jk)
      END DO
      !
      IF( ln_bdy ) THEN
         DO jk = 1, jpkm1
            wsd(:,:,jk) = wsd(:,:,jk) * bdytmask(:,:)
         END DO
      ENDIF
      !                       !==  Horizontal divergence of barotropic Stokes transport  ==!
      div_sd(:,:) = 0._wp
      DO jk = 1, jpkm1                                 !
        div_sd(:,:) = div_sd(:,:) + ze3divh(:,:,jk)
      END DO
      !
      CALL iom_put( "ustokes",  usd  )
      CALL iom_put( "vstokes",  vsd  )
      CALL iom_put( "wstokes",  wsd  )
!      !
      DEALLOCATE( ze3divh, zInt_w )
      DEALLOCATE( zk_t, zk_u, zk_v, zu0_sd, zv0_sd )
      !
   END SUBROUTINE sbc_stokes
!
!
   SUBROUTINE sbc_wave( kt, Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_wave  ***
      !!
      !! ** Purpose :   read wave parameters from wave model in netcdf files
      !!                or from a coupled wave mdoel
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt   ! ocean time step
      INTEGER, INTENT(in   ) ::   Kmm  ! ocean time index
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_wave : update the read waves fields'
         WRITE(numout,*) '~~~~~~~~ '
      ENDIF
      !
      IF( ln_cdgw .AND. .NOT. cpl_wdrag ) THEN     !==  Neutral drag coefficient  ==!
         CALL fld_read( kt, nn_fsbc, sf_cd )             ! read from external forcing
         cdn_wave(:,:) = sf_cd(1)%fnow(:,:,1) * tmask(:,:,1)
      ENDIF

      IF( ln_tauoc .AND. .NOT. cpl_wstrf ) THEN    !==  Wave induced stress  ==!
         CALL fld_read( kt, nn_fsbc, sf_tauoc )          ! read stress reduction factor due to wave from external forcing
         tauoc_wave(:,:) = sf_tauoc(1)%fnow(:,:,1) * tmask(:,:,1)
      ELSEIF ( ln_taw .AND. cpl_taw ) THEN
         IF (kt < 1) THEN ! The first fields gave by OASIS have very high erroneous values ....
            twox(:,:)=0._wp
            twoy(:,:)=0._wp
            tawx(:,:)=0._wp
            tawy(:,:)=0._wp
            tauoc_wavex(:,:) = 1._wp
            tauoc_wavey(:,:) = 1._wp
         ELSE
            tauoc_wavex(:,:) = abs(twox(:,:)/tawx(:,:))
            tauoc_wavey(:,:) = abs(twoy(:,:)/tawy(:,:))
         ENDIF
      ENDIF

      ! Read also wave number if needed, so that it is available in coupling routines
      IF( ln_STOKES_ADIAB .AND. .NOT. cpl_wnum ) THEN     !==wavenumber==!
         CALL fld_read( kt, nn_fsbc, sf_wn )             ! read wave parameters from external forcing
         wnum(:,:) = sf_wn(1)%fnow(:,:,1) * tmask(:,:,1)
      ENDIF

      IF ( ln_phioc .and. cpl_phioc .and.  kt == nit000 ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_wave : PHIOC from wave model'
         WRITE(numout,*) '~~~~~~~~ '
      ENDIF

      IF( ln_sdw .AND. .NOT. cpl_sdrftx)  THEN       !==  Computation of the 3d Stokes Drift  ==!
         !
         IF( jpfld > 0 ) THEN                            ! Read from file only if the field is not coupled
            CALL fld_read( kt, nn_fsbc, sf_sd )          ! read wave parameters from external forcing
            !                                            ! NB: test case mode, not read as jpfld=0
            IF( jp_hsw > 0 )   hsw  (:,:) = sf_sd(jp_hsw)%fnow(:,:,1) * tmask(:,:,1)  ! significant wave height
            IF( jp_wmp > 0 )   wmp  (:,:) = sf_sd(jp_wmp)%fnow(:,:,1) * tmask(:,:,1)  ! wave mean period
            IF( jp_usd > 0 )   ut0sd(:,:) = sf_sd(jp_usd)%fnow(:,:,1) * tmask(:,:,1)  ! 2D zonal Stokes Drift at T point
            IF( jp_vsd > 0 )   vt0sd(:,:) = sf_sd(jp_vsd)%fnow(:,:,1) * tmask(:,:,1)  ! 2D meridional Stokes Drift at T point
         ENDIF
         !
         IF( jpfld == 4 .OR. ln_wave_test )   &
            &      CALL sbc_stokes( Kmm )                 ! Calculate only if all required fields are read
            !                                            ! or in wave test case
         !  !                                            ! In coupled case the call is done after (in sbc_cpl)
      ENDIF
         !
   END SUBROUTINE sbc_wave


   SUBROUTINE sbc_wave_init
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_wave_init  ***
      !!
      !! ** Purpose :   Initialisation fo surface waves
      !!
      !! ** Method  : - Read namelist namsbc_wave
      !!              - create the structure used to read required wave fields
      !!                (its size depends on namelist options)
      !! ** action
      !!---------------------------------------------------------------------
      INTEGER ::   ierror, ios   ! local integer
      INTEGER ::   ifpr
      !!
      CHARACTER(len=100)     ::  cn_dir                            ! Root directory for location of drag coefficient files
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::   slf_i            ! array of namelist informations on the fields to read
      TYPE(FLD_N)            ::  sn_cdg, sn_usd, sn_vsd,  &
                             &   sn_hsw, sn_wmp, sn_wnum, sn_tauoc    ! informations about the fields to be read
      !
      NAMELIST/namsbc_wave/ cn_dir, sn_cdg, sn_usd, sn_vsd, sn_hsw, sn_wmp, sn_wnum, sn_tauoc,   &
         &                  ln_cdgw, ln_sdw, ln_tauoc, ln_stcor, ln_charn, ln_taw, ln_phioc,     &
         &                  ln_wave_test, ln_bern_srfc, ln_breivikFV_2016, ln_vortex_force, ln_stshear
      !!---------------------------------------------------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_wave_init : surface waves in the system'
         WRITE(numout,*) '~~~~~~~~~~~~~ '
      ENDIF
      !
      READ  ( numnam_ref, namsbc_wave, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_wave in reference namelist')

      READ  ( numnam_cfg, namsbc_wave, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_wave in configuration namelist' )
      IF(lwm) WRITE ( numond, namsbc_wave )
      !
      IF(lwp) THEN
         WRITE(numout,*) '   Namelist namsbc_wave'
         WRITE(numout,*) '      Stokes drift                                  ln_sdw = ', ln_sdw
         WRITE(numout,*) '      Breivik 2016                       ln_breivikFV_2016 = ', ln_breivikFV_2016
         WRITE(numout,*) '      Stokes Coriolis & tracer advection terms    ln_stcor = ', ln_stcor
         WRITE(numout,*) '      Vortex Force                         ln_vortex_force = ', ln_vortex_force
         WRITE(numout,*) '      Bernouilli Head Pressure                ln_bern_srfc = ', ln_bern_srfc
         WRITE(numout,*) '      wave modified ocean stress                  ln_tauoc = ', ln_tauoc
         WRITE(numout,*) '      neutral drag coefficient (CORE bulk only)    ln_cdgw = ', ln_cdgw
         WRITE(numout,*) '      charnock coefficient                        ln_charn = ', ln_charn
         WRITE(numout,*) '      Stress modificated by wave                    ln_taw = ', ln_taw
         WRITE(numout,*) '      TKE flux from wave                          ln_phioc = ', ln_phioc
         WRITE(numout,*) '      Surface shear with Stokes drift           ln_stshear = ', ln_stshear
         WRITE(numout,*) '      Test with constant wave fields          ln_wave_test = ', ln_wave_test
      ENDIF

      !                                ! option check
      IF( .NOT.( ln_cdgw .OR. ln_sdw .OR. ln_tauoc .OR. ln_stcor .OR. ln_charn) )   &
         &     CALL ctl_warn( 'Ask for wave coupling but ln_cdgw=F, ln_sdw=F, ln_tauoc=F, ln_stcor=F')
      IF( ln_cdgw .AND. ln_blk )   &
         &     CALL ctl_stop( 'drag coefficient read from wave model NOT available yet with aerobulk package')
      IF( ln_stcor .AND. .NOT.ln_sdw )   &
         &     CALL ctl_stop( 'Stokes-Coriolis term calculated only if activated Stokes Drift ln_sdw=T')

      !                             !==  Allocate wave arrays  ==!
      ALLOCATE( ut0sd (jpi,jpj)    , vt0sd (jpi,jpj) )
      ALLOCATE( hsw   (jpi,jpj)    , wmp   (jpi,jpj) )
      ALLOCATE( wnum  (jpi,jpj) )
      ALLOCATE( tsd2d (jpi,jpj)    , div_sd(jpi,jpj)    , bhd_wave(jpi,jpj)     )
      ALLOCATE( usd   (jpi,jpj,jpk), vsd   (jpi,jpj,jpk), wsd     (jpi,jpj,jpk) )
      ALLOCATE( tusd  (jpi,jpj)    , tvsd  (jpi,jpj)    , ZMX     (jpi,jpj,jpk) )
      usd   (:,:,:) = 0._wp
      vsd   (:,:,:) = 0._wp
      wsd   (:,:,:) = 0._wp
      hsw     (:,:) = 0._wp
      wmp     (:,:) = 0._wp
      ut0sd   (:,:) = 0._wp
      vt0sd   (:,:) = 0._wp
      tusd    (:,:) = 0._wp
      tvsd    (:,:) = 0._wp
      bhd_wave(:,:) = 0._wp
      ZMX   (:,:,:) = 0._wp
!
      IF( ln_wave_test ) THEN       !==  Wave TEST case  ==!   set uniform waves fields
         jpfld    = 0                   ! No field read
         ln_cdgw  = .FALSE.             ! No neutral wave drag input
         ln_tauoc = .FALSE.             ! No wave induced drag reduction factor
         ut0sd(:,:) = 0.13_wp * tmask(:,:,1)   ! m/s
         vt0sd(:,:) = 0.00_wp                  ! m/s
         hsw  (:,:) = 2.80_wp                  ! meters
         wmp  (:,:) = 8.00_wp                  ! seconds
         !
      ELSE                          !==  create the structure associated with fields to be read  ==!
         IF( ln_cdgw ) THEN                       ! wave drag
            IF( .NOT. cpl_wdrag ) THEN
               ALLOCATE( sf_cd(1), STAT=ierror )               !* allocate and fill sf_wave with sn_cdg
               IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wave structure' )
               !
                                      ALLOCATE( sf_cd(1)%fnow(jpi,jpj,1)   )
               IF( sn_cdg%ln_tint )   ALLOCATE( sf_cd(1)%fdta(jpi,jpj,1,2) )
               CALL fld_fill( sf_cd, (/ sn_cdg /), cn_dir, 'sbc_wave_init', 'Wave module ', 'namsbc_wave' )
            ENDIF
            ALLOCATE( cdn_wave(jpi,jpj) )
            cdn_wave(:,:) = 0._wp
         ENDIF
         IF( ln_charn ) THEN                     ! wave drag
            IF( .NOT. cpl_charn ) THEN
               CALL ctl_stop( 'STOP', 'Charnock based wind stress can be used in coupled mode only' )
            ENDIF
            ALLOCATE( charn(jpi,jpj) )
            charn(:,:) = 0._wp
         ENDIF
         IF( ln_taw ) THEN                     ! wind stress
            IF( .NOT. cpl_taw ) THEN
               CALL ctl_stop( 'STOP', 'wind stress from wave model can be used in coupled mode only, use ln_cdgw instead' )
            ENDIF
            ALLOCATE( tawx(jpi,jpj) )
            ALLOCATE( tawy(jpi,jpj) )
            ALLOCATE( twox(jpi,jpj) )
            ALLOCATE( twoy(jpi,jpj) )
            ALLOCATE( tauoc_wavex(jpi,jpj) )
            ALLOCATE( tauoc_wavey(jpi,jpj) )
            tawx(:,:) = 0._wp
            tawy(:,:) = 0._wp
            twox(:,:) = 0._wp
            twoy(:,:) = 0._wp
            tauoc_wavex(:,:) = 1._wp
            tauoc_wavey(:,:) = 1._wp
         ENDIF

         IF( ln_phioc ) THEN                     ! TKE flux
            IF( .NOT. cpl_phioc ) THEN
                CALL ctl_stop( 'STOP', 'phioc can be used in coupled mode only' )
            ENDIF
            ALLOCATE( phioc(jpi,jpj) )
            phioc(:,:) = 0._wp
         ENDIF

         IF( ln_tauoc ) THEN                    ! normalized wave stress into the ocean
            IF( .NOT. cpl_wstrf ) THEN
               ALLOCATE( sf_tauoc(1), STAT=ierror )           !* allocate and fill sf_wave with sn_tauoc
               IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_tauoc structure' )
               !
                                       ALLOCATE( sf_tauoc(1)%fnow(jpi,jpj,1)   )
               IF( sn_tauoc%ln_tint )  ALLOCATE( sf_tauoc(1)%fdta(jpi,jpj,1,2) )
               CALL fld_fill( sf_tauoc, (/ sn_tauoc /), cn_dir, 'sbc_wave_init', 'Wave module', 'namsbc_wave' )
            ENDIF
            ALLOCATE( tauoc_wave(jpi,jpj) )
            tauoc_wave(:,:) = 0._wp
         ENDIF

         IF( ln_sdw ) THEN                      ! Stokes drift
            ! 1. Find out how many fields have to be read from file if not coupled
            jpfld=0
            jp_usd=0   ;   jp_vsd=0   ;   jp_hsw=0   ;   jp_wmp=0
            IF( .NOT. cpl_sdrftx ) THEN
               jpfld  = jpfld + 1
               jp_usd = jpfld
            ENDIF
            IF( .NOT. cpl_sdrfty ) THEN
               jpfld  = jpfld + 1
               jp_vsd = jpfld
            ENDIF
            IF( .NOT. cpl_hsig ) THEN
               jpfld  = jpfld + 1
               jp_hsw = jpfld
            ENDIF
            IF( .NOT. cpl_wper ) THEN
               jpfld  = jpfld + 1
               jp_wmp = jpfld
            ENDIF
            ! 2. Read from file only the non-coupled fields
            IF( jpfld > 0 ) THEN
               ALLOCATE( slf_i(jpfld) )
               IF( jp_usd > 0 )   slf_i(jp_usd) = sn_usd
               IF( jp_vsd > 0 )   slf_i(jp_vsd) = sn_vsd
               IF( jp_hsw > 0 )   slf_i(jp_hsw) = sn_hsw
               IF( jp_wmp > 0 )   slf_i(jp_wmp) = sn_wmp
               ALLOCATE( sf_sd(jpfld), STAT=ierror )   !* allocate and fill sf_sd with stokes drift
               IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wave structure' )
               !
               DO ifpr= 1, jpfld
                  ALLOCATE( sf_sd(ifpr)%fnow(jpi,jpj,1) )
                  IF( slf_i(ifpr)%ln_tint )   ALLOCATE( sf_sd(ifpr)%fdta(jpi,jpj,1,2) )
               END DO
               !
               CALL fld_fill( sf_sd, slf_i, cn_dir, 'sbc_wave_init', 'Wave module ', 'namsbc_wave' )
               sf_sd(jp_usd)%zsgn = -1._wp   ;  sf_sd(jp_vsd)%zsgn = -1._wp   ! vector field at T point: overwrite default definition of zsgn
            ENDIF
                       !
            ! 3. Wave number (only needed for Qiao parametrisation,
            ! ln_zdfqiao=T)
            IF( ln_STOKES_ADIAB ) THEN
             IF( .NOT. cpl_wnum ) THEN
                ALLOCATE( sf_wn(1), STAT=ierror )           !* allocate and fill sf_wave with sn_wnum
                IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wn structure' )
                                       ALLOCATE( sf_wn(1)%fnow(jpi,jpj,1)   )
                IF( sn_wnum%ln_tint )  ALLOCATE( sf_wn(1)%fdta(jpi,jpj,1,2) )
                CALL fld_fill( sf_wn, (/ sn_wnum /), cn_dir, 'sbc_wave_init', 'Wave module', 'namsbc_wave' )
             ENDIF
            ENDIF
            !
            !
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE sbc_wave_init

   !!======================================================================
END MODULE sbcwave
