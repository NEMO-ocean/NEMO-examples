MODULE step_c1d
   !!======================================================================
   !!                       ***  MODULE step_c1d  ***
   !! Time-stepping    : manager of the ocean, tracer and ice time stepping - c1d case
   !!======================================================================
   !! History :   2.0  !  2004-04  (C. Ethe)  adapted from step.F90 for C1D
   !!             3.0  !  2008-04  (G. Madec)  redo the adaptation to include SBC
   !!             4.1  !  2019-08  (A. Coward, D. Storkey) rewrite in preparation for new timestepping scheme
   !!             4.1  !  2019-12  (L. Brodeau) STATION_ASF test-case
   !!----------------------------------------------------------------------
#if defined key_c1d
   !!----------------------------------------------------------------------
   !!   'key_c1d'                                       1D Configuration
   !!----------------------------------------------------------------------
   !!   stp_c1d        : NEMO system time-stepping in c1d case
   !!----------------------------------------------------------------------
   USE step_oce        ! time stepping definition modules
   USE step, ONLY : Nbb, Nnn, Naa, Nrhs ! time level indices
   USE restart         ! restart

   IMPLICIT NONE
   PRIVATE

   PUBLIC stp_c1d   ! called by nemogcm.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: step_c1d.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE stp_c1d( kstp )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp_c1d  ***
      !!
      !! ** Purpose :  - Time stepping of SBC including sea ice (dynamic and thermodynamic eqs.)
      !!               - Time stepping of OPA (momentum and active tracer eqs.)
      !!               - Time stepping of TOP (passive tracer eqs.)
      !!
      !! ** Method  : -1- Update forcings and data
      !!              -2- Update vertical ocean physics
      !!              -3- Compute the t and s trends
      !!              -4- Update t and s
      !!              -5- Compute the momentum trends
      !!              -6- Update the horizontal velocity
      !!              -7- Compute the diagnostics variables (rd,N2, div,cur,w)
      !!              -8- Outputs and diagnostics
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kstp   ! ocean time-step index
      !
      INTEGER ::   jk       ! dummy loop indice
      INTEGER ::   indic    ! error indicator if < 0
      !! ---------------------------------------------------------------------

      indic = 0                ! reset to no error condition
      IF( kstp == nit000 )   CALL iom_init( "nemo")   ! iom_put initialization (must be done after nemo_init for AGRIF+XIOS+OASIS)
      IF( kstp /= nit000 )   CALL day( kstp )         ! Calendar (day was already called at nit000 in day_init)
      CALL iom_setkt( kstp - nit000 + 1, "nemo" )   ! say to iom that we are at time step kstp

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Update data, open boundaries, surface boundary condition (including sea-ice)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      CALL sbc    ( kstp, Nbb, Nnn )  ! Sea Boundary Condition (including sea-ice)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! diagnostics and outputs
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      CALL dia_wri( kstp, Nnn )  ! ocean model: outputs

      ! Swap time levels
      Nrhs = Nbb
      Nbb = Nnn
      Nnn = Naa
      Naa = Nrhs

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Control and restarts
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      CALL stp_ctl( kstp, Nnn )

      IF( kstp == nit000 )   CALL iom_close( numror )          ! close input  ocean restart file
      IF( lrst_oce       )   CALL rst_write( kstp, Nbb, Nnn )  ! write output ocean restart file
      !
#if defined key_iomput
      IF( kstp == nitend .OR. nstop > 0 )   CALL xios_context_finalize()   ! needed for XIOS
      !
#endif
   END SUBROUTINE stp_c1d

#else
   !!----------------------------------------------------------------------
   !!   Default key                                            NO 1D Config
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE stp_c1d ( kt )      ! dummy routine
      IMPLICIT NONE
      INTEGER, INTENT( in ) :: kt
      WRITE(*,*) 'stp_c1d: You should not have seen this print! error?', kt
   END SUBROUTINE stp_c1d
#endif

   !!======================================================================
END MODULE step_c1d
