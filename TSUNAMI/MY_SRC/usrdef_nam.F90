MODULE usrdef_nam
   !!======================================================================
   !!                       ***  MODULE  usrdef_nam  ***
   !!
   !!                      ===  TSUNAMI configuration  ===
   !!
   !! User defined : set the domain characteristics of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2017-10  (J. Chanut)  Original code
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

   PUBLIC   usr_def_nam   ! called by nemogcm.F90

   !                              !!* namusr_def namelist *!!
   REAL(wp)         ::   rn_domszx  = 1800.  ! x horizontal size         [km]
   REAL(wp)         ::   rn_domszy  = 1800.  ! y horizontal size         [km]
   REAL(wp), PUBLIC ::   rn_domszz  = 5000.  ! z horizontal size          [m]
   REAL(wp), PUBLIC ::   rn_dx      =   30.  ! x horizontal resolution   [km]
   REAL(wp), PUBLIC ::   rn_dy      =   30.  ! y horizontal resolution   [km]
   REAL(wp), PUBLIC ::   rn_0xratio =    0.5 ! x domain ratio of the 0
   REAL(wp), PUBLIC ::   rn_0yratio =    0.5 ! x domain ratio of the 0
   INTEGER , PUBLIC ::   nn_fcase   =    1   ! F computation (0:f0, 1:Beta, 2:real)
   REAL(wp), PUBLIC ::   rn_ppgphi0 =   38.5 ! reference latitude for beta-plane 

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_nam.F90 13472 2020-09-16 13:05:19Z smasson $ 
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
      !!                Here TSUNAMI configuration
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
      INTEGER ::   ios      ! Local integer
      LOGICAL ::   ln_Iperio, ln_Jperio
      !!
      NAMELIST/namusr_def/  rn_domszx, rn_domszy, rn_domszz, rn_dx, rn_dy, rn_0xratio, rn_0yratio   &
         &                 , nn_fcase, rn_ppgphi0, ln_Iperio, ln_Jperio
      !!----------------------------------------------------------------------
      !
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist' )
      !
      IF(lwm)   WRITE( numond, namusr_def )
      !
#if defined key_agrif 
      ! Domain parameters are taken from parent:
      IF( .NOT. Agrif_Root() ) THEN
         rn_dx = Agrif_Parent(rn_dx)/Agrif_Rhox()
         rn_dy = Agrif_Parent(rn_dy)/Agrif_Rhoy()
         rn_ppgphi0 = Agrif_Parent(rn_ppgphi0)
      ENDIF
      rn_0xratio = 0.5   ! not really working I guess... 
      rn_0yratio = 0.5
#endif
      !
      cd_cfg = 'TSUNAMI'             ! name & resolution (not used)
      kk_cfg = INT( rn_dx )
      !
      IF( Agrif_Root() ) THEN        ! Global Domain size:  TSUNAMI global domain is  1800 km x 1800 Km x 5000 m
         kpi = NINT( rn_domszx / rn_dx ) + 1
         kpj = NINT( rn_domszy / rn_dy ) + 1
      ELSE                           ! Global Domain size: add nbghostcells + 1 "land" point on each side
         kpi  = nbcellsx + nbghostcells_x_w + nbghostcells_x_e + 2
         kpj  = nbcellsy + nbghostcells_y_s + nbghostcells_y_n + 2
      ENDIF
      kpk = 2
      !                              ! Set the lateral boundary condition of the global domain
      !
      ldIperio = ln_Iperio   ;   ldJperio = ln_Jperio
      ldNFold  =  .FALSE.    ;   cdNFtype = '-'
      !
      !                              ! control print
      IF(lwp) THEN
         WRITE(numout,*) '   '
         WRITE(numout,*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'
         WRITE(numout,*) '~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namusr_def : TSUNAMI test case'
         WRITE(numout,*) '      horizontal domain size-x          rn_domszx  = ', rn_domszx, ' km'
         WRITE(numout,*) '      horizontal domain size-y          rn_domszy  = ', rn_domszy, ' km'
         WRITE(numout,*) '      vertical   domain size-z          rn_domszz  = ', rn_domszz, '  m'
         WRITE(numout,*) '      horizontal x-resolution           rn_dx      = ',     rn_dx, ' km'
         WRITE(numout,*) '      horizontal y-resolution           rn_dy      = ',     rn_dy, ' km'
         WRITE(numout,*) '      x-domain ratio of the 0           rn_0xratio = ', rn_0xratio
         WRITE(numout,*) '      y-domain ratio of the 0           rn_0yratio = ', rn_0yratio
         WRITE(numout,*) '      F computation                     nn_fcase   = ',   nn_fcase
         WRITE(numout,*) '      Reference latitude                rn_ppgphi0 = ', rn_ppgphi0
         WRITE(numout,*) '   '
      ENDIF
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
