MODULE usrdef_nam
   !!======================================================================
   !!                     ***  MODULE usrdef_nam   ***
   !!
   !!                     ===  SWG configuration  ===
   !!
   !! User defined : set the domain characteristics of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-03  (S. Flavoni, G. Madec)  Original code
   !!             -   ! 2020-03  (A. Nasser) Shallow Water Eq. configuration
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_nam   : read user defined namelist and set global domain size
   !!   usr_def_hgr   : initialize the horizontal mesh 
   !!----------------------------------------------------------------------
   USE dom_oce  , ONLY: nimpp, njmpp       ! ocean space and time domain
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_nam   ! called in nemogcm.F90 module
   !                              !!* namusr_def namelist *!!
   INTEGER , PUBLIC ::   nn_SWG             ! 1/nn_SWG = the resolution chosen in degrees and thus defining the horizontal domain size
   REAL(wp), PUBLIC ::   rn_theta           ! rotation angle (in degree) of the grid 
   INTEGER , PUBLIC ::   nn_gc              ! number of ghostcells
   REAL(wp), PUBLIC ::   rn_domsiz          ! size of the domain (default 2000km)  [m]
   REAL(wp), PUBLIC ::   rn_dx              ! gridspacing (default 100km)          [m]
   REAL(wp), PUBLIC ::   rn_tau             ! wind stress on the surface           [N/m2]
   REAL(wp), PUBLIC ::   rn_f0              !    f-plan coriolis parameter         [1/s] 
   REAL(wp), PUBLIC ::   rn_beta            ! beta-plan coriolis parameter         [1/m.s] 
   REAL(wp), PUBLIC ::   rn_modified_grav   ! modified gravity                     [m/s2] 
   REAL(wp), PUBLIC ::   rn_rfr             ! layer friction                       [1/s]
   !                              !!* temporary *!!
   INTEGER , PUBLIC ::   nn_dynldf_lap_typ            ! choose type of laplacian (ideally from namelist)
   !                                                  ! = 1   divrot    laplacian
   !                                                  ! = 2   symmetric laplacian (Griffies&Hallberg 2000)
   !                                                  ! = 3   symmetric laplacian (cartesian)
   LOGICAL , PUBLIC ::   ln_dynldf_lap_PM             ! if True - apply the P.Marchand boundary condition on the laplacian
   !                              !!* penalisation *!!
   REAL(wp), PUBLIC ::   rn_abp             ! alpha boundary parameter                                       [-]
   INTEGER , PUBLIC ::   nn_cnp             ! number of cell on which is smoothed the porosity (phi)         [-]
   REAL(wp), PUBLIC ::   rn_fsp             ! friction parameter 1/epsilon of the permeability               [1/s]
   !
   REAL(wp), PUBLIC ::   r1_abp             ! inverse alpha boundary parameter                            [-]
   !
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_nam.F90 11536 2019-09-11 13:54:18Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_nam( cd_cfg, kk_cfg, kpi, kpj, kpk, ldIperio, ldJperio, ldNFold, cdNFtype )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!                    
      !! ** Purpose :   read user defined namelist and define the domain size
      !!
      !! ** Method  :   read in namusr_def containing all the user specific namelist parameter
      !!
      !!                Here SWG configuration
      !!
      !! ** input   : - namusr_def namelist found in namelist_cfg
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(out) ::   cd_cfg               ! configuration name
      INTEGER         , INTENT(out) ::   kk_cfg               ! configuration resolution
      INTEGER         , INTENT(out) ::   kpi, kpj, kpk        ! global domain sizes
      LOGICAL         , INTENT(out) ::   ldIperio, ldJperio   ! i- and j- periodicity
      LOGICAL         , INTENT(out) ::   ldNFold              ! North pole folding
      CHARACTER(len=1), INTENT(out) ::   cdNFtype             ! Folding type: T or F
      !
      INTEGER  ::   ios             ! Local integer
      REAL(wp) ::   ze1, zgl, zbl   ! gridspacing, length of the biggest square
      !!
      NAMELIST/namusr_def/ nn_SWG, rn_theta, jpkglo,            &   ! 
         &                 nn_gc ,rn_domsiz, rn_dx,             &   ! domain parameters
         &                 rn_f0 ,rn_beta,                      &   ! coriolis parameter
         &                 rn_modified_grav, rn_rfr , rn_tau,   &   ! reduced gravity, friction, wind
         &                 nn_dynldf_lap_typ, ln_dynldf_lap_PM, &   ! temporary parameter
         &                 rn_abp, nn_cnp, rn_fsp                   ! penalisation parameters
      !!----------------------------------------------------------------------
      !
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist' )
      !
      IF(lwm)   WRITE( numond, namusr_def )
      !
      cd_cfg = 'SWG'               ! name & resolution (not used)
      
#if defined key_agrif
      IF (.NOT.Agrif_root()) nn_SWG = Agrif_parent(nn_SWG) * Agrif_irhox()
#endif

      kk_cfg = nn_SWG
       !
      ze1 =  rn_dx / REAL(nn_SWG, wp)                    ! [m] gridspacing used
      zgl =  rn_domsiz + 2._wp * REAL(nn_gc, wp) * ze1   ! [m] length of the square with ghostcells
      ! rotation
      zbl = zgl * ( COS( rn_theta * rad ) + SIN( rn_theta * rad ) )   ! length side bigger domain [m]  
      !
      kpi = ceiling(zbl / ze1 )    ! Global Domain size + ghost cells                
      kpj = ceiling(zbl / ze1 )    ! Global Domain size + ghost cells              
      !
      IF( rn_modified_grav /= 0._wp) grav = rn_modified_grav   ! update the gravity 
      !
      kpk = jpkglo
      !                             ! Set the lateral boundary condition of the global domain
      ldIperio = .FALSE.   ;   ldJperio = .FALSE.   ! SWG configuration : closed domain
      ldNFold  = .FALSE.   ;   cdNFtype = '-'
      !
# if defined key_bvp
      r1_abp = 1._wp / rn_abp
#endif
      !                             ! control print
      IF(lwp) THEN
         WRITE(numout,*) '   '
         WRITE(numout,*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'
         WRITE(numout,*) '~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namusr_def : SWG case'
         WRITE(numout,*) '                                   domain size       rn_domsiz   = ', rn_domsiz, 'm'
         WRITE(numout,*) '                                   gridspacing           rn_dx   = ', rn_dx, 'm'
         WRITE(numout,*) '      inverse resolution & implied domain size          nn_SWG   = ', nn_SWG
         WRITE(numout,*) '                           implied gridspacing           rn_dx   = ', rn_dx, 'm'
         WRITE(numout,*) '                          number of ghostcells           nn_gc   = ', nn_gc
         WRITE(numout,*) '   '
         WRITE(numout,*) '                    rotation angle chosen             rn_theta   = ', rn_theta, 'deg'
         WRITE(numout,*) '                    modified gravity          rn_modified_grav   = ', rn_modified_grav, 'm2/s'
         WRITE(numout,*) '      number of model levels                              jpkglo = ', kpk
         WRITE(numout,*) '   '
      ENDIF
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
