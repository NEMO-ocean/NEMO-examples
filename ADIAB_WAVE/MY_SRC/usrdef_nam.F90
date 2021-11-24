MODULE usrdef_nam
   !!======================================================================
   !!                     ***  MODULE usrdef_nam   ***
   !!
   !!                     ===  ADIAB_WAVE configuration  ===
   !!
   !! User defined : set the domain characteristics of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-03  (S. Flavoni, G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_nam   : read user defined namelist and set global domain size
   !!   usr_def_hgr   : initialize the horizontal mesh
   !!----------------------------------------------------------------------
   USE dom_oce
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_nam   ! called in nemogcm.F90 module

   !                              !!* namusr_def namelist *!!
   REAL(wp), PUBLIC ::   rn_dx     ! resolution in meters defining the horizontal domain size
   REAL(wp), PUBLIC ::   rn_dz     ! resolution in meters defining the vertical domain size
   REAL(wp), PUBLIC ::   rn_dy
   LOGICAL , PUBLIC ::   ln_STOKES_ADIAB         !Shallow/Inter water formula for the Stokes drift
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
      !!                Here ADIAB configuration
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
      NAMELIST/namusr_def/ ln_zco, ln_zps, ln_sco, rn_dx, rn_dy, rn_dz, ln_STOKES_ADIAB
      !!----------------------------------------------------------------------
      !
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist' )
      !
      IF(lwm)   WRITE( numond, namusr_def )
      !
      cd_cfg = 'ADIAB_WAVE'           ! name & resolution
      kk_cfg = INT( rn_dx )
      !
      ! Global Domain size:  ADIAB_WAVE domain is  780 km x 5 grid-points x 4-6 m
      kpi = INT( 780. / rn_dx )
      kpj = 5. !
      kpk = INT(  6. / rn_dz )!+1

      !                             ! Set the lateral boundary condition of the global domain
      ldIperio = .FALSE.    ;   ldJperio = .FALSE.   !  closed domain
      ldNFold  = .FALSE.   ;   cdNFtype = '-'
      !
      !                             ! control print
      IF(lwp) THEN
      WRITE(numout,*) '   '
      WRITE(numout,*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'
      WRITE(numout,*) '~~~~~~~~~~~ '
      WRITE(numout,*) '   Namelist namusr_def : ADIAB_WAVE test case'
      WRITE(numout,*) '      type of vertical coordinate : '
      WRITE(numout,*) '         z-coordinate flag                     ln_zco = ', ln_zco
      WRITE(numout,*) '         z-partial-step coordinate flag        ln_zps = ', ln_zps
      WRITE(numout,*) '         s-coordinate flag                     ln_sco = ', ln_sco
      WRITE(numout,*) '      horizontal resolution                    rn_dx  = ', rn_dx, ' meters'
      WRITE(numout,*) '      horizontal resolution                    rn_dy  = ', rn_dy, ' meters'
      WRITE(numout,*) '      vertical   resolution                    rn_dz  = ', rn_dz, ' meters'
      WRITE(numout,*) '      ADIAB_WAVE domain = 780 m x 5 grid-points x 4-6 m'
      WRITE(numout,*) '         resulting global domain size :        Ni0glo = ', kpi
      WRITE(numout,*) '                                               Nj0glo = ', kpj
      WRITE(numout,*) '                                               jpkglo = ', kpk
      WRITE(numout,*) ' Stokes drift for Shallow/Inter water (ln_STOKES_ADIAB) = ', ln_STOKES_ADIAB    
      !
      ! Set the lateral boundary condition of the
      ! ADIAB_WAVE configuration : E-W periodic basin
      !
      WRITE(numout,*) '   '
      WRITE(numout,*) '   Lateral boundary condition of the global domain'
      WRITE(numout,*) '      ADIAB_WAVE : closed basin                  Iperio =', ldIperio
      WRITE(numout,*) '      ADIAB_WAVE : closed basin                  Jperio =', ldJperio
      ENDIF
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
