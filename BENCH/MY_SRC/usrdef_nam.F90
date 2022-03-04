MODULE usrdef_nam
   !!======================================================================
   !!                       ***  MODULE  usrdef_nam  ***
   !!
   !!                      ===  BENCH configuration  ===
   !!
   !! User defined : set the domain characteristics of a user configuration
   !!======================================================================
   !! History :  NEMO !
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_nam   : read user defined namelist and set global domain size
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! to get ctl_nam
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_nam   ! called by nemogcm.F90
   
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2016)
   !! $Id$ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
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
      !!                Here EW_CANAL configuration
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
      INTEGER ::   ios         ! Local integer
      !                              !!* namusr_def namelist *!!
      INTEGER ::   nn_isize    ! number of point in i-direction of global(local) domain if >0 (<0)  
      INTEGER ::   nn_jsize    ! number of point in j-direction of global(local) domain if >0 (<0)  
      INTEGER ::   nn_ksize    ! total number of point in k-direction
      !                              !!* nammpp namelist *!!
      INTEGER          ::   jpni, jpnj
      LOGICAL          ::   ln_listonly
      LOGICAL          ::   ln_Iperio, ln_Jperio
      LOGICAL          ::   ln_NFold
      character(len=1) ::   cn_NFtype
      !!
      NAMELIST/namusr_def/ nn_isize, nn_jsize, nn_ksize, ln_Iperio, ln_Jperio, ln_NFold, cn_NFtype
      NAMELIST/nammpp/ jpni, jpnj, nn_hls, ln_nnogather, ln_listonly, nn_comm
      !!----------------------------------------------------------------------     
      !
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 903 )
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist' )
      !
      IF(lwm)   WRITE( numond, namusr_def )      
      !
      cd_cfg = 'BENCH'             ! name & resolution (not used)
      kk_cfg = 0
      !
      IF( nn_isize < 0 .AND. nn_jsize < 0 ) THEN
      !
         READ  ( numnam_ref, nammpp, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 )   CALL ctl_nam ( ios , 'nammpp in reference namelist' )
         !
         READ  ( numnam_cfg, nammpp, IOSTAT = ios, ERR = 902 )
902      IF( ios >  0 )   CALL ctl_nam ( ios , 'nammpp in configuration namelist' )

         kpi = -nn_isize * jpni
         kpj = -nn_jsize * jpnj
      ELSE
         kpi = nn_isize
         kpj = nn_jsize
      ENDIF
      kpk = nn_ksize
      !
      ldIperio = ln_Iperio   ;   ldJperio = ln_Jperio
      ldNFold  = ln_NFold    ;   cdNFtype = cn_NFtype
      !
      !                             ! control print
      IF(lwp) THEN
         WRITE(numout,*) '   '
         WRITE(numout,*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'
         WRITE(numout,*) '~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namusr_def : BENCH test case'
         IF( nn_isize > 0 ) THEN
            WRITE(numout,*) '      global domain size-x            nn_isize = ',  nn_isize
         ELSE
            WRITE(numout,*) '                                          jpni = ', jpni
            WRITE(numout,*) '       local domain size-x           -nn_isize = ', -nn_isize
            WRITE(numout,*) '      global domain size-x                 kpi = ', kpi
         ENDIF
         IF( nn_jsize > 0 ) THEN
            WRITE(numout,*) '      global domain size-y            nn_jsize = ', nn_jsize
         ELSE
            WRITE(numout,*) '                                          jpnj = ', jpnj
            WRITE(numout,*) '       local domain size-y           -nn_jsize = ', -nn_jsize
            WRITE(numout,*) '      global domain size-y                 kpj = ', kpj
         ENDIF
         WRITE(numout,*) '      global domain size-z            nn_ksize = ', nn_ksize
         WRITE(numout,*) '   '
      ENDIF
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
