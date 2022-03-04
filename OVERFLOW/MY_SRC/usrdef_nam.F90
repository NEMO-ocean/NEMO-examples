MODULE usrdef_nam
   !!======================================================================
   !!                       ***  MODULE  usrdef_nam  ***
   !!
   !!                     ===  OVERFLOW configuration  ===
   !!
   !! User defined : set the domain characteristics of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2016-03  (S. Flavoni, G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_nam   : read user defined namelist and set global domain size
   !!   usr_def_hgr   : initialize the horizontal mesh 
   !!----------------------------------------------------------------------
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

   !                              !!* namusr_def namelist *!!
   REAL(wp), PUBLIC ::   rn_dx     ! resolution in meters defining the horizontal domain size
   REAL(wp), PUBLIC ::   rn_dz     ! resolution in meters defining the vertical   domain size

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_nam.F90 14433 2021-02-11 08:06:49Z smasson $ 
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
      !!                Here OVERFLOW configuration
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
      INTEGER ::   ios   ! Local integer
      !!
      NAMELIST/namusr_def/ ln_zco, ln_zps, ln_sco, rn_dx, rn_dz
      !!----------------------------------------------------------------------
      !
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist' )
      !
      IF(lwm)   WRITE( numond, namusr_def )
      !
      cd_cfg = 'OVERFLOW'           ! name & resolution (not used)
      kk_cfg = INT( rn_dx )
      !
      ! Global Domain size:  OVERFLOW domain is  200 km x 3 grid-points x 2000 m
      kpi = INT( 200.e3 / rn_dx ) + 2
      kpj = 3
      kpk = INT(  2000. / rn_dz ) + 1
      !
      !                             ! control print
      WRITE(numout,*) '   '
      WRITE(numout,*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'
      WRITE(numout,*) '~~~~~~~~~~~ '
      WRITE(numout,*) '   Namelist namusr_def : OVERFLOW test case'
      WRITE(numout,*) '      type of vertical coordinate : '
      WRITE(numout,*) '         z-coordinate flag                     ln_zco = ', ln_zco
      WRITE(numout,*) '         z-partial-step coordinate flag        ln_zps = ', ln_zps
      WRITE(numout,*) '         s-coordinate flag                     ln_sco = ', ln_sco
      WRITE(numout,*) '      horizontal resolution                    rn_dx  = ', rn_dx, ' meters'
      WRITE(numout,*) '      vertical   resolution                    rn_dz  = ', rn_dz, ' meters'
      WRITE(numout,*) '      OVERFLOW domain = 200 km x 3 grid-points x 2000 m'
      WRITE(numout,*) '         resulting global domain size :        Ni0glo = ', kpi
      WRITE(numout,*) '                                               Nj0glo = ', kpj
      WRITE(numout,*) '                                               jpkglo = ', kpk
      WRITE(numout,*) '   '
      !
      !                             ! Set the lateral boundary condition of the global domain
      ldIperio = .FALSE.   ;   ldJperio = .FALSE.   ! OVERFLOW configuration : closed domain
      ldNFold  = .FALSE.   ;   cdNFtype = '-'
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
