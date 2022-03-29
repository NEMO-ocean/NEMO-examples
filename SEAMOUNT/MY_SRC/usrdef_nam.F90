MODULE usrdef_nam
   !!======================================================================
   !!                     ***  MODULE usrdef_nam   ***
   !!
   !!                   ===  SEAMOUNT configuration  === test
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
   !! PHYSICAL DOMAIN
   REAL(wp), PUBLIC :: rn_xdim     ! x-dimension of the domain [km]
   REAL(wp), PUBLIC :: rn_ydim     ! y-dimension of the domain [km]
   REAL(wp), PUBLIC :: rn_bot_max  ! max ocean depth (>0) [m]
   REAL(wp), PUBLIC :: rn_smnt_H   ! seamount height (>0) [m]
   REAL(wp), PUBLIC :: rn_smnt_L   ! seamount width [km]
   REAL(wp), PUBLIC :: rn_fplane   ! Coriolis parameter for f-plane approximation   
   REAL(wp), PUBLIC :: rn_Snum     ! Burger number of the initial stratification
   !! NUMERICAL DISCRETIZATION
   REAL(wp), PUBLIC :: rn_dx       ! horizontal resolution [m]          
   REAL(wp), PUBLIC :: rn_dz       ! vertical resolution far from the seamount
                                   ! and assuming no stretching [m]
   !! VERTICAL COORDINATE
   ! Stretched s-levels (ln_sco = .true.)
   LOGICAL,  PUBLIC :: ln_s_sh94   ! TRUE:  s-levels using Song & Haidvogel 1994 (SH94)
                                   !        stretching function
                                   ! FALSE: uniform sigma-levels
   REAL(wp), PUBLIC :: rn_theta    ! SH94 surface control parameter (0<=theta<=20)
   REAL(wp), PUBLIC :: rn_bb       ! SH94 bottom control parameter (0<=bb<=1)
   REAL(wp), PUBLIC :: rn_hc       ! SH94 critical depth for transition to
                                   ! stretched coordinates [m]
   ! Paramaters for vqs-coordinate (ln_sco = .true.)
   LOGICAL,  PUBLIC :: ln_vqs      ! activating vanishing quasi-sigma levels (TRUE)
   REAL(wp), PUBLIC :: rn_rmax     ! maximum cut-off slope parameter value allowed
                                   ! if using vqs-levels (0<r_max<1)
   !! INITIAL CONDITION
   INTEGER,  PUBLIC :: nn_ini_cond ! 0: Shchepetkin and McWilliams (2002)
                                   ! 1: Ezer, Arango and Shchepetkin (2002)

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
      !!                Here SEAMOUNT configuration defined as in 
      !!                Shchepetkin & McWilliams (2003) or Ezer, Arango and Shchepetkin (2002)
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
      INTEGER                       ::   ios    ! Local integer
      INTEGER                       ::   ioptio ! Local string
      !!
      NAMELIST/namusr_def/ rn_xdim  , rn_ydim  , rn_bot_max , rn_smnt_H   , rn_smnt_L, & 
         &                 rn_fplane, rn_dx    , rn_dz      , ln_zco      , ln_zps   , &
         &                 ln_sco   , ln_s_sh94, rn_theta   , rn_bb       , rn_hc    , &
         &                 ln_vqs   , rn_rmax  , nn_ini_cond
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

      !                        ! Set the lateral boundary condition of the global domain
      SELECT CASE (nn_ini_cond)
         CASE (0)             ! Shchepetkin and McWilliams (2002):
            ldIperio = .TRUE. ! basin with cyclic Est-West
         CASE (1)             ! Ezer, Arango and Shchepetkin (2002):   
            ldIperio = .FALSE.! basin with closed boundaries
      END SELECT
      
      ldJperio = .FALSE.
      ldNFold  = .FALSE.
      cdNFtype = '-'

      !
      ! Global Domain size:
      kpi = INT(  1000._wp * rn_xdim / rn_dx ) + 1
      kpj = INT(  1000._wp * rn_ydim / rn_dx ) + 1
      kpk = INT(  rn_bot_max  / rn_dz ) + 1
      
      ! Check vertical coordinate options
      ioptio = 0           
      IF( ln_zco      )   ioptio = ioptio + 1
      IF( ln_zps      )   ioptio = ioptio + 1
      IF( ln_sco      )   ioptio = ioptio + 1
      IF( ioptio /= 1 )   CALL ctl_stop( ' none or several vertical coordinate options used' )
      
      !                             ! SEAMOUNT_TEST_CASE configuration : closed domain
      !                             ! control print
      IF(lwp) THEN
         WRITE(numout,*) '   '                                                                          
         WRITE(numout,*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'   
         WRITE(numout,*) '~~~~~~~~~~~ '                                                                 
         WRITE(numout,*) '   Namelist namusr_def : SEAMOUNT test case'                        
         WRITE(numout,*) '      *) horizontal resolution                    rn_dx     =  ', rn_dx, ' m' 
         WRITE(numout,*) '      *) vertical   resolution                    rn_dz     =  ', rn_dz, ' m'
         WRITE(numout,*) '         (far from the seamount and assuming no stretching)'
         WRITE(numout,*) '      *) vertical   coordinate:                '
         WRITE(numout,*) '           z-coordinate - full steps            ln_zco    = ', ln_zco
         WRITE(numout,*) '           z-coordinate - partial steps         ln_zps    = ', ln_zps
         WRITE(numout,*) '           s-coordinate                         ln_sco    = ', ln_sco
         WRITE(numout,*) '           stretched s-levels                   ln_s_sh94 = ', ln_s_sh94
         IF( ln_s_sh94 ) THEN
           WRITE(numout,*) '                                                rn_theta  = ', rn_theta
           WRITE(numout,*) '                                                rn_bb     = ', rn_bb
           WRITE(numout,*) '                                                rn_hc     = ', rn_hc
         ENDIF
         WRITE(numout,*) '           using vanishing quasi-sigma levels   ln_vqs    = ', ln_vqs
         IF( ln_vqs    ) THEN
           WRITE(numout,*) '           maximum slope parameter allowed      rn_rmax   = ', rn_rmax
         ENDIF
         WRITE(numout,*) '      *) initial condition:                    '
         IF( nn_ini_cond == 0 ) THEN
           WRITE(numout,*) '           Shchepetkin and McWilliams (2002)'
         ELSE IF( nn_ini_cond == 1 ) THEN
           WRITE(numout,*) '           Ezer, Arango and Shchepetkin (2002)'
         ENDIF         
         WRITE(numout,*) ''
         WRITE(numout,*) '   Global domain size:                   '
         WRITE(numout,*) '          1st dimension of glob. dom. i-dir.     jpiglo  =  ', kpi
         WRITE(numout,*) '          2nd dimension of glob. dom. j-dir.     jpjglo  =  ', kpj        
         WRITE(numout,*) '          3rd dimension of glob. dom. k-dir.     jpkglo  =  ', kpk              
         WRITE(numout,*) '   Lateral boundary condition of the global domain'                           
         WRITE(numout,*) '      east-west                                   ldIperio = ', ldIperio     
         WRITE(numout,*) '      north-south                                 ldJperio = ', ldJperio   
      ENDIF
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
