MODULE usrdef_nam
   !!======================================================================
   !!                       ***  MODULE  usrdef_nam  ***
   !!
   !!                      ===  ICE_ADV2D configuration  ===
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
   !! $Id: usrdef_nam.F90 10161 2018-10-01 09:48:55Z clem $ 
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
      !!                Here ICE_ADV2D configuration
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
#if defined key_agrif 
      ! Domain parameters are taken from parent:
      IF( .NOT. Agrif_Root() ) THEN
         rn_dx = Agrif_Parent(rn_dx)/Agrif_Rhox()
         rn_dy = Agrif_Parent(rn_dy)/Agrif_Rhoy()
         rn_ppgphi0 = Agrif_Parent(rn_ppgphi0)
      ENDIF
#endif
      !
      WRITE( ldnam(:), namusr_def )
      !
      cd_cfg = 'ICE_ADV2D'           ! name & resolution (not used)
      kk_cfg = INT( rn_dx )
      !
      ! Global Domain size:  ICE_ADV2D domain is  300 km x 300 Km x 10 m
      kpi = INT( 300.e3 / rn_dx ) -1
      kpj = INT( 300.e3 / rn_dy ) -1
#if defined key_agrif
      IF( .NOT. Agrif_Root() ) THEN
         kpi = nbcellsx + 2 + 2*nbghostcells
         kpj = nbcellsy + 2 + 2*nbghostcells
      ENDIF
#endif
      kpk = 1
      !
!!      zlx = (kpi-2)*rn_dx*1.e-3
!!      zly = (kpj-2)*rn_dy*1.e-3
      zlx = kpi*rn_dx*1.e-3
      zly = kpj*rn_dy*1.e-3
      !                             ! control print
      WRITE(ldtxt(ii),*) '   '                                                                          ;   ii = ii + 1
      WRITE(ldtxt(ii),*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'   ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '~~~~~~~~~~~ '                                                                 ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Namelist namusr_def : ICE_ADV2D test case'                                 ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      horizontal resolution                    rn_dx  = ', rn_dx, ' meters'   ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      horizontal resolution                    rn_dy  = ', rn_dy, ' meters'   ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      ICE_ADV2D domain = 300 km x 300Km x 1 grid-point '                      ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         LX [km]: ', zlx                                                      ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         LY [km]: ', zly                                                      ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         resulting global domain size :        jpiglo = ', kpi                ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '                                               jpjglo = ', kpj                ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '                                               jpkglo = ', kpk                ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '         Coriolis:', ln_corio                                                 ;   ii = ii + 1
      !
      !                             ! Set the lateral boundary condition of the global domain
      kperio = 7                    ! ICE_ADV2D configuration : bi-periodic basin
#if defined key_agrif
      IF( .NOT. Agrif_Root() ) THEN
      kperio = 0
      ENDIF
#endif
      !
      WRITE(ldtxt(ii),*) '   '                                                                          ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Lateral boundary condition of the global domain'                           ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      ICE_ADV2D : bi-periodic basin               jperio = ', kperio          ;   ii = ii + 1
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
