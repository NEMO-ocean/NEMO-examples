MODULE usrdef_nam
   !!======================================================================
   !!                     ***  MODULE usrdef_nam   ***
   !!
   !!                   ===  SEAMOUNT configuration  ===
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
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_nam   ! called in nemogcm.F90 module

   !                              !!* namusr_def namelist *!!
   REAL(wp), PUBLIC ::   rn_dx     ! resolution in meters defining the horizontal domain size
   REAL(wp), PUBLIC ::   rn_dz     ! resolution in meters defining the vertical   domain size
   REAL(wp), PUBLIC ::   rn_length
   REAL(wp), PUBLIC ::   rn_width
   REAL(wp), PUBLIC ::   rn_drho     ! resolution in meters defining the horizontal domain size
   REAL(wp), PUBLIC ::   rn_initrho
   REAL(wp), PUBLIC ::   rn_s
   REAL(wp), PUBLIC ::   rn_bathy
   REAL(wp), PUBLIC ::   rn_seamountheight
   REAL(wp), PUBLIC ::   rn_l
   REAL(wp), PUBLIC ::   rn_f
   LOGICAL,  PUBLIC ::   ln_exp_init
   LOGICAL,  PUBLIC ::   ln_linear_init

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_nam.F90 14072 2020-12-04 07:48:38Z laurent $
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
      !!                Here SEAMOUNT configuration defined as per Beckmann
      !!                and Haidvogel (1993) with uniformly spaced sigma coordinates
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
      NAMELIST/namusr_def/ ln_zco, ln_zps, ln_sco, rn_length, rn_width, rn_dx, rn_dz, rn_initrho, rn_s, rn_bathy, rn_seamountheight, rn_l, rn_f, ln_exp_init, ln_linear_init
      !!----------------------------------------------------------------------
      !
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist' )
      !
      IF(lwm)   WRITE( numond, namusr_def )
      !
      !
      cd_cfg = 'Seamount'      ! name & resolution (not used)
      kk_cfg = 0

      !                             ! Set the lateral boundary condition of the global domain
      ldIperio = .TRUE.    ;   ldJperio = .true.   ! C1D configuration : 3x3 basin with cyclic Est-West and Norht-South condition
      ldNFold  = .FALSE.   ;   cdNFtype = '-'

      !
      ! Global Domain size:  SEAMOUNT_TEST_CASE domain is rn_length km x rn_width km x rn_bathy m
      kpi = INT(  1000._wp * rn_length / rn_dx ) + 1
      kpj = INT(  1000._wp * rn_width / rn_dx ) + 2
      kpk = INT(  rn_bathy  / rn_dz ) + 1
      ! Calculating the density difference from the given Burger Number in the namelist_cfg
      ! rn_drho =  rho_ref * depth * (S * f * L)^2 / g 
      rn_drho = 1000._wp * rn_bathy * ( rn_s * rn_f * rn_l / rn_bathy ) ** 2._wp / grav
      !
      !                             ! SEAMOUNT_TEST_CASE configuration : closed domain
      !                             ! control print
      IF(lwp) THEN
         WRITE(numout,*) '   '                                                                          
         WRITE(numout,*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'   
         WRITE(numout,*) '~~~~~~~~~~~ '                                                                 
         WRITE(numout,*) '   Namelist namusr_def : SEAMOUNT_TEST_CASE test case'                        
         WRITE(numout,*) '      horizontal resolution                    rn_dx   =  ', rn_dx, ' meters' 
         WRITE(numout,*) '      vertical   resolution                    rn_dz   =  ', rn_dz, ' meters' 
         WRITE(numout,*) '      SEAMOUNT_TEST_CASE domain'                                              
         WRITE(numout,*) '         resulting global domain size :        jpiglo  =  ', kpi              
         WRITE(numout,*) '                                               jpjglo  =  ', kpj              
         WRITE(numout,*) '                                               jpkglo  =  ', kpk              
         WRITE(numout,*) '   For Burger Number S = ', rn_s,            ' rn_drho =  ', rn_drho          
         WRITE(numout,*) '   '                                                                          
         WRITE(numout,*) '   Lateral boundary condition of the global domain'                           
         WRITE(numout,*) '      east-west                                   ldIperio = ', ldIperio     
         WRITE(numout,*) '      north-south                                 ldJperio = ', ldJperio   
         WRITE(numout,*) '   Initial condition:                          ln_exp_init =    ', ln_exp_init
         WRITE(numout,*) '                                               ln_linear_init = ', ln_linear_init
      ENDIF
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
