MODULE usrdef_nam
   !!======================================================================
   !!                       ***  MODULE  usrdef_nam  ***
   !!
   !!                     ===  ISOMIP configuration  ===
   !!
   !! User defined : set the domain characteristics of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2016-03  (S. Flavoni, G. Madec)   Original code
   !!                 ! 2017-02  (P. Mathiot, S. Flavoni) Adapt code to ISOMIP case
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_nam   : read user defined namelist and set global domain size
   !!   usr_def_hgr   : initialize the horizontal mesh 
   !!----------------------------------------------------------------------
   USE dom_oce  , ONLY: nimpp , njmpp            ! i- & j-indices of the local domain
   USE dom_oce  , ONLY: ln_zco, ln_zps, ln_sco   ! flag of type of coordinate
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_nam   ! called by nemogcm.F90

   !                                         !!* namusr_def namelist *!!
   REAL(wp), PUBLIC ::   rn_e1deg, rn_e2deg   !: horizontal resolution   [degrees]
   REAL(wp), PUBLIC ::   rn_e3                !: vertical   resolution         [m]
   
   REAL(wp), PARAMETER, PUBLIC ::   rbathy = 900._wp   !: depth of the seafloor   [m]

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_nam.F90 10074 2018-08-28 16:15:49Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_nam( ldtxt, ldnam, cd_cfg, kk_cfg, kpi, kpj, kpk, kperio )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!                    
      !! ** Purpose :   read user defined namelist and define the domain size
      !!
      !! ** Method  :   read in namusr_def containing all the user specific namelist parameter
      !!
      !!                Here ISOMIP configuration
      !!
      !! ** input   : - namusr_def namelist found in namelist_cfg
      !!----------------------------------------------------------------------
      CHARACTER(len=*), DIMENSION(:), INTENT(out) ::   ldtxt, ldnam    ! stored print information
      CHARACTER(len=*)              , INTENT(out) ::   cd_cfg          ! configuration name
      INTEGER                       , INTENT(out) ::   kk_cfg          ! configuration resolution
      INTEGER                       , INTENT(out) ::   kpi, kpj, kpk   ! global domain sizes 
      INTEGER                       , INTENT(out) ::   kperio          ! lateral global domain b.c. 
      !
      INTEGER ::   ios, ii   ! Local integer
      !!
      NAMELIST/namusr_def/ ln_zco, ln_zps, ln_sco, rn_e1deg, rn_e2deg, rn_e3
      !!----------------------------------------------------------------------
      !
      ii = 1
      !
      REWIND( numnam_cfg )          ! Namelist namusr_def (exist in namelist_cfg only)
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist', .TRUE. )
      !
      WRITE( ldnam(:), namusr_def )
      !
      cd_cfg = 'ISOMIP'           ! name & resolution (not used)
      kk_cfg = INT( rn_e3 )
      !
      ! Global Domain size:  ISOMIP domain is  15° x 10° x 900 m
      kpi = INT(  15.0  / rn_e1deg ) + 2     ! add 2 for t-point in the east  & west  coasts
      kpj = INT(  10.0  / rn_e2deg ) + 2     !     -        -           north & south   -
      kpk = INT( rbathy / rn_e3    ) + 1     ! add 1 for t-point in the seafloor
      !
      !                             ! control print
      WRITE(ldtxt(ii),*) '   '                                                                              ;   ii = ii + 1
      WRITE(ldtxt(ii),*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'       ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '~~~~~~~~~~~ '                                                                     ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Namelist namusr_def : ISOMIP test case'                                        ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      type of vertical coordinate : '                                             ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         z-coordinate flag                     ln_zco   = ', ln_zco               ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         z-partial-step coordinate flag        ln_zps   = ', ln_zps               ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         s-coordinate flag                     ln_sco   = ', ln_sco               ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      resolution'                                                                 ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         zonal      resolution                 rn_e1deg = ', rn_e1deg, ' degrees' ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         meridional resolution                 rn_e1deg = ', rn_e1deg, ' degrees' ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         vertical   resolution                 rn_e3    = ', rn_e3   , ' meters'  ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      ISOMIP domain = 15° x 10° x 900 m'                                          ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         resulting global domain size :        jpiglo   = ', kpi                  ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '                                               jpjglo   = ', kpj                  ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '                                               jpkglo   = ', kpk                  ;   ii = ii + 1
      !
      !
      !                             ! Set the lateral boundary condition of the global domain
      kperio = 0                    ! ISOMIP configuration : close basin
      !
      WRITE(ldtxt(ii),*) '   '                                                                          ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Lateral boundary condition of the global domain'                           ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      ISOMIP : closed basin                    jperio   = ', kperio             ;   ii = ii + 1
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
