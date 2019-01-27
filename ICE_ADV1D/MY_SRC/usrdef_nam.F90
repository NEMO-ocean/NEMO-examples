MODULE usrdef_nam
   !!======================================================================
   !!                       ***  MODULE  usrdef_nam  ***
   !!
   !!                      ===  ICE_ADV1D configuration  ===
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

   !                               !!* namusr_def namelist *!!
   REAL(wp), PUBLIC ::   rn_dx      ! resolution in meters defining the horizontal domain size
   REAL(wp), PUBLIC ::   rn_dy      ! resolution in meters defining the horizontal domain size
   REAL(wp), PUBLIC ::   rn_ppgphi0 ! reference latitude for beta-plane 
   LOGICAL , PUBLIC ::   ln_corio   ! set coriolis at 0 (ln_corio=F) or not 

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
      !!                Here ICE_ADV1D configuration
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
      REAL(wp)::   zlx, zly  ! Local scalars
      !!
      NAMELIST/namusr_def/ rn_dx, rn_dy, ln_corio, rn_ppgphi0
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
      cd_cfg = 'ICE_ADV1D'           ! name & resolution (not used)
      kk_cfg = INT( rn_dx )
      !
      ! Global Domain size:  ICE_ADV1D domain is  480 m x 480 m x 10 m
      kpi = INT( 480.*0.5 / rn_dx ) -1
      kpj = INT( 480.*0.5 / rn_dy ) -1
      kpk = 1
      !
      zlx = kpi*rn_dx*1.e-3
      zly = kpj*rn_dy*1.e-3
      !                             ! control print
      WRITE(ldtxt(ii),*) '   '                                                                          ;   ii = ii + 1
      WRITE(ldtxt(ii),*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'   ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '~~~~~~~~~~~ '                                                                 ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Namelist namusr_def : ICE_ADV1D test case'                                 ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      horizontal resolution                    rn_dx  = ', rn_dx, ' meters'   ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      horizontal resolution                    rn_dy  = ', rn_dy, ' meters'   ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      ICE_ADV1D domain  '                                                     ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         LX [km]: ', zlx                                                      ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         LY [km]: ', zly                                                      ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         resulting global domain size :        jpiglo = ', kpi                ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '                                               jpjglo = ', kpj                ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '                                               jpkglo = ', kpk                ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         Coriolis:', ln_corio                                                 ;   ii = ii + 1
      !
      !                             ! Set the lateral boundary condition of the global domain
      kperio = 0                    ! ICE_ADV1D configuration : bi-periodic basin
      !
      WRITE(ldtxt(ii),*) '   '                                                                          ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Lateral boundary condition of the global domain'                           ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      ICE_ADV1D : closed basin                    jperio = ', kperio          ;   ii = ii + 1
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
