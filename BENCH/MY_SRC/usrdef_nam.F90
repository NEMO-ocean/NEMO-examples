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

   SUBROUTINE usr_def_nam( ldtxt, ldnam, cd_cfg, kk_cfg, kpi, kpj, kpk, kperio )
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
      CHARACTER(len=*), DIMENSION(:), INTENT(out) ::   ldtxt, ldnam    ! stored print information
      CHARACTER(len=*)              , INTENT(out) ::   cd_cfg          ! configuration name
      INTEGER                       , INTENT(out) ::   kk_cfg          ! configuration resolution
      INTEGER                       , INTENT(out) ::   kpi, kpj, kpk   ! global domain sizes 
      INTEGER                       , INTENT(out) ::   kperio          ! lateral global domain b.c. 
      !
      !
      INTEGER ::   ios, ii      ! Local integer
      !                              !!* namusr_def namelist *!!
      INTEGER ::   nn_isize    ! number of point in i-direction of global(local) domain if >0 (<0)  
      INTEGER ::   nn_jsize    ! number of point in j-direction of global(local) domain if >0 (<0)  
      INTEGER ::   nn_ksize    ! total number of point in k-direction
      INTEGER ::   nn_perio    ! periodicity
      !                              !!* nammpp namelist *!!
      CHARACTER(len=1) ::   cn_mpi_send
      INTEGER          ::   nn_buffer, jpni, jpnj
      LOGICAL          ::   ln_nnogather
      !!
      NAMELIST/namusr_def/ nn_isize, nn_jsize, nn_ksize, nn_perio
      NAMELIST/nammpp/ cn_mpi_send, nn_buffer, jpni, jpnj, ln_nnogather
      !!----------------------------------------------------------------------     
      !
      REWIND( numnam_cfg )          ! Namelist namusr_def (exist in namelist_cfg only)
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 903 )
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist', .TRUE. )
      WRITE( ldnam(:), namusr_def )      
      !
      cd_cfg = 'BENCH'             ! name & resolution (not used)
      kk_cfg = 0

      IF( nn_isize < 0 .AND. nn_jsize < 0 ) THEN
      !
         REWIND( numnam_ref )              ! Namelist nammpp in reference namelist: mpi variables
         READ  ( numnam_ref, nammpp, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 )   CALL ctl_nam ( ios , 'nammpp in reference namelist', lwp )
         !
         REWIND( numnam_cfg )              ! Namelist nammpp in configuration namelist: mpi variables
         READ  ( numnam_cfg, nammpp, IOSTAT = ios, ERR = 902 )
902      IF( ios >  0 )   CALL ctl_nam ( ios , 'nammpp in configuration namelist', lwp )

         kpi = ( -nn_isize - 2*nn_hls ) * jpni + 2*nn_hls
         kpj = ( -nn_jsize - 2*nn_hls ) * jpnj + 2*nn_hls
      ELSE
         kpi = nn_isize
         kpj = nn_jsize
      ENDIF
 
      kpk = nn_ksize
      kperio = nn_perio

      !                             ! control print
      ii = 1
      WRITE(ldtxt(ii),*) '   '                                                                          ;   ii = ii + 1
      WRITE(ldtxt(ii),*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'   ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '~~~~~~~~~~~ '                                                                 ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Namelist namusr_def : BENCH test case'                                     ;   ii = ii + 1
      IF( nn_isize > 0 ) THEN
         WRITE(ldtxt(ii),*) '      global domain size-x            nn_isize = ',  nn_isize              ;   ii = ii + 1
      ELSE
         WRITE(ldtxt(ii),*) '                                          jpni = ', jpni                   ;   ii = ii + 1
         WRITE(ldtxt(ii),*) '       local domain size-x           -nn_isize = ', -nn_isize              ;   ii = ii + 1
         WRITE(ldtxt(ii),*) '      global domain size-x                 kpi = ', kpi                    ;   ii = ii + 1
      ENDIF
      IF( nn_jsize > 0 ) THEN
         WRITE(ldtxt(ii),*) '      global domain size-y            nn_jsize = ', nn_jsize               ;   ii = ii + 1
      ELSE
         WRITE(ldtxt(ii),*) '                                          jpnj = ', jpnj                   ;   ii = ii + 1
         WRITE(ldtxt(ii),*) '       local domain size-y           -nn_jsize = ', -nn_jsize              ;   ii = ii + 1
         WRITE(ldtxt(ii),*) '      global domain size-y                 kpj = ', kpj                    ;   ii = ii + 1
      ENDIF
      WRITE(ldtxt(ii),*) '      global domain size-z            nn_ksize = ', nn_ksize                  ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      LBC of the global domain          kperio = ', kperio                    ;   ii = ii + 1
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
