MODULE usrdef_nam
   !!======================================================================
   !!                       ***  MODULE usrdef_nam  ***
   !!
   !!                  ===  LOCK_EXCHANGE configuration  ===
   !!
   !! User defined : set the domain characteristics of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2016-03  (S. Flavoni, G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_nam   : read user defined namelist and set global domain size
   !!   usr_def_hgr   : initialize the horizontal mesh 
   !!----------------------------------------------------------------------
   USE dom_oce  , ONLY: nimpp , njmpp            ! i- & j-indices of the local domain
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_nam   ! called by nemogcm.F90

   !                              !!* namusr_def namelist *!!
   REAL(wp), PUBLIC ::   rn_dx     ! resolution in meters defining the horizontal domain size
   REAL(wp), PUBLIC ::   rn_dz     ! resolution in meters defining the vertical   domain size

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
      !!                Here LOCK_EXCHANGE configuration
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
      NAMELIST/namusr_def/ rn_dx, rn_dz
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
      !
      cd_cfg = 'LOCK_EXCHANGE'      ! name & resolution (not used)
      kk_cfg = INT( rn_dx )
      !
      ! Global Domain size:  LOCK_EXCHANGE domain is 64 km x 3 grid-points x 20 m
      kpi = INT(  64.e3 / rn_dx ) + 2
      kpj = 3
      kpk = INT(  20.  / rn_dz ) + 1
      !
      !                             ! control print
      WRITE(ldtxt(ii),*) '   '                                                                          ;   ii = ii + 1
      WRITE(ldtxt(ii),*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'   ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '~~~~~~~~~~~ '                                                                 ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Namelist namusr_def : LOCK_EXCHANGE test case'                             ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      horizontal resolution                    rn_dx  = ', rn_dx, ' meters'   ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      vertical   resolution                    rn_dz  = ', rn_dz, ' meters'   ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      LOCK_EXCHANGE domain = 64 km  x  3 grid-points  x  20 m'                ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         resulting global domain size :        jpiglo = ', kpi                ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '                                               jpjglo = ', kpj                ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '                                               jpkglo = ', kpk                ;   ii = ii + 1
      !
      !                             ! Set the lateral boundary condition of the global domain
      kperio = 0                    ! LOCK_EXCHANGE configuration : closed domain
      !
      WRITE(ldtxt(ii),*) '   '                                                                          ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Lateral boundary condition of the global domain'                           ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      closed                                   jperio = ', kperio             ;   ii = ii + 1
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
