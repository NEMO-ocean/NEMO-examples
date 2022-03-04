MODULE isfstp
   !!======================================================================
   !!                       ***  MODULE  isfstp  ***
   !! Ice Shelves :  compute iceshelf load, melt and heat flux
   !!======================================================================
   !! History :  3.2  !  2011-02  (C.Harris  ) Original code isf cav
   !!            X.X  !  2006-02  (C. Wang   ) Original code bg03
   !!            3.4  !  2013-03  (P. Mathiot) Merging + parametrization
   !!            4.1  !  2019-09  (P. Mathiot) Split param/explicit ice shelf and re-organisation
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   isfstp       : compute iceshelf melt and heat flux
   !!----------------------------------------------------------------------
   USE isf_oce                                      ! isf variables
   USE isfload, ONLY: isf_load                      ! ice shelf load
   USE isftbl , ONLY: isf_tbl_lvl                   ! ice shelf boundary layer
   USE isfpar , ONLY: isf_par, isf_par_init         ! ice shelf parametrisation
   USE isfcav , ONLY: isf_cav, isf_cav_init         ! ice shelf cavity
   USE isfcpl , ONLY: isfcpl_rst_write, isfcpl_init ! isf variables

   USE dom_oce        ! ocean space and time domain
   USE oce      , ONLY: ssh                           ! sea surface height
   USE domvvl,  ONLY: ln_vvl_zstar                      ! zstar logical
   USE zdfdrg,  ONLY: r_Cdmin_top, r_ke0_top            ! vertical physics: top/bottom drag coef.
   !
   USE lib_mpp, ONLY: ctl_stop, ctl_nam
   USE fldread, ONLY: FLD, FLD_N
   USE in_out_manager ! I/O manager
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   isf_stp, isf_init, isf_nam  ! routine called in sbcmod and divhor

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: isfstp.F90 11876 2019-11-08 11:26:42Z mathiot $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
 
   SUBROUTINE isf_stp( kt, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isf_stp  ***
      !!
      !! ** Purpose : compute total heat flux and total fwf due to ice shelf melt
      !!
      !! ** Method  : For each case (parametrisation or explicity cavity) :
      !!              - define the before fields
      !!              - compute top boundary layer properties
      !!                (in case of parametrisation, this is the 
      !!                 depth range model array used to compute mean far fields properties)
      !!              - compute fluxes
      !!              - write restart variables
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time step
      INTEGER, INTENT(in) ::   Kmm   ! ocean time level index
      !
      INTEGER :: jk                              ! loop index
#if defined key_qco
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ze3t   ! 3D workspace
#endif
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('isf')
      !
      !=======================================================================
      ! 1.: compute melt and associated heat fluxes in the ice shelf cavities
      !=======================================================================
      !
      IF ( ln_isfcav_mlt ) THEN
         !
         ! 1.1: before time step 
         IF ( kt /= nit000 ) THEN 
            risf_cav_tsc_b (:,:,:) = risf_cav_tsc (:,:,:)
            fwfisf_cav_b(:,:)      = fwfisf_cav(:,:)
         END IF
         !
         ! 1.2: compute misfkb, rhisf_tbl, rfrac (deepest level, thickness, fraction of deepest cell affected by tbl)
         rhisf_tbl_cav(:,:) = rn_htbl * mskisf_cav(:,:)
#if defined key_qco
         DO jk = 1, jpk
            ze3t(:,:,jk) = e3t(:,:,jk,Kmm)
         END DO 
         CALL isf_tbl_lvl( ht(:,:), ze3t           , misfkt_cav, misfkb_cav, rhisf_tbl_cav, rfrac_tbl_cav )
#else
         CALL isf_tbl_lvl( ht(:,:),  e3t(:,:,:,Kmm), misfkt_cav, misfkb_cav, rhisf_tbl_cav, rfrac_tbl_cav )
#endif
         !
         ! 1.3: compute ice shelf melt
         CALL isf_cav( kt, Kmm, risf_cav_tsc, fwfisf_cav )
         !
      END IF
      ! 
      !=================================================================================
      ! 2.: compute melt and associated heat fluxes for not resolved ice shelf cavities
      !=================================================================================
      !
      IF ( ln_isfpar_mlt ) THEN
         !
         ! 2.1: before time step 
         IF ( kt /= nit000 ) THEN 
            risf_par_tsc_b(:,:,:) = risf_par_tsc(:,:,:)
            fwfisf_par_b  (:,:)   = fwfisf_par  (:,:)
         END IF
         !
         ! 2.2: compute misfkb, rhisf_tbl, rfrac (deepest level, thickness, fraction of deepest cell affected by tbl)
         ! by simplicity, we assume the top level where param applied do not change with time (done in init part)
         rhisf_tbl_par(:,:) = rhisf0_tbl_par(:,:)
#if defined key_qco
         DO jk = 1, jpk
            ze3t(:,:,jk) = e3t(:,:,jk,Kmm)
         END DO
         CALL isf_tbl_lvl( ht(:,:), ze3t           , misfkt_par, misfkb_par, rhisf_tbl_par, rfrac_tbl_par )
#else
         CALL isf_tbl_lvl( ht(:,:),  e3t(:,:,:,Kmm), misfkt_par, misfkb_par, rhisf_tbl_par, rfrac_tbl_par )
#endif
         !
         ! 2.3: compute ice shelf melt
         CALL isf_par( kt, Kmm, risf_par_tsc, fwfisf_par )
         !
      END IF
      !
      !==================================================================================
      ! 3.: output specific restart variable in case of coupling with an ice sheet model
      !==================================================================================
      !
      IF ( ln_isfcpl .AND. lrst_oce ) CALL isfcpl_rst_write(kt, Kmm)
      !
      IF( ln_timing )   CALL timing_stop('isf')
      !
   END SUBROUTINE isf_stp

   
   SUBROUTINE isf_init( Kbb, Kmm, Kaa )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isfstp_init  ***
      !!
      !! ** Purpose :   Initialisation of the ice shelf public variables
      !!
      !! ** Method  :   Read the namisf namelist, check option compatibility and set derived parameters
      !!
      !! ** Action  : - read namisf parameters
      !!              - allocate memory
      !!              - output print
      !!              - ckeck option compatibility
      !!              - call cav/param/isfcpl init routine
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Kmm, Kaa   ! ocean time level indices
      !!----------------------------------------------------------------------
      !
      ! constrain: l_isfoasis need to be known
      !
      CALL isf_nam()                                              ! Read namelist
      !
      CALL isf_alloc()                                            ! Allocate public array
      !
      CALL isf_ctl()                                              ! check option compatibility
      !
      IF( ln_isfcav ) CALL isf_load( Kmm, risfload )              ! compute ice shelf load
      !
      ! terminate routine now if no ice shelf melt formulation specify
      IF( ln_isf ) THEN
         !
         IF( ln_isfcav_mlt )   CALL isf_cav_init()                ! initialisation melt in the cavity
         !
         IF( ln_isfpar_mlt )   CALL isf_par_init()                ! initialisation parametrised melt
         !
         IF( ln_isfcpl     )   CALL isfcpl_init( Kbb, Kmm, Kaa )  ! initialisation ice sheet coupling
         !
      END IF
         
  END SUBROUTINE isf_init

  
  SUBROUTINE isf_ctl()
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isf_ctl  ***
      !!
      !! ** Purpose :   output print and option compatibility check
      !!
      !!----------------------------------------------------------------------
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'isf_init : ice shelf initialisation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namisf :'
         !
         WRITE(numout,*) '   ice shelf cavity (open or parametrised)  ln_isf = ', ln_isf
         WRITE(numout,*)
         !
         IF ( ln_isf ) THEN
            WRITE(numout,*) '      Add debug print in isf module           ln_isfdebug     = ', ln_isfdebug
            WRITE(numout,*)
            WRITE(numout,*) '      melt inside the cavity                  ln_isfcav_mlt   = ', ln_isfcav_mlt
            IF ( ln_isfcav_mlt) THEN
               WRITE(numout,*) '         melt formulation                         cn_isfcav_mlt= ', TRIM(cn_isfcav_mlt)
               WRITE(numout,*) '         thickness of the top boundary layer      rn_htbl      = ', rn_htbl
               WRITE(numout,*) '         gamma formulation                        cn_gammablk  = ', TRIM(cn_gammablk) 
               IF ( TRIM(cn_gammablk) .NE. 'spe' ) THEN 
                  WRITE(numout,*) '         gammat coefficient                       rn_gammat0   = ', rn_gammat0  
                  WRITE(numout,*) '         gammas coefficient                       rn_gammas0   = ', rn_gammas0  
                  WRITE(numout,*) '         top background ke used (from namdrg_top) rn_vtide**2  = ', rn_vtide**2
                  WRITE(numout,*) '         top drag coef.    used (from namdrg_top) rn_Cd0       = ', r_Cdmin_top
               END IF
            END IF
            WRITE(numout,*) ''
            !
            WRITE(numout,*) '      ice shelf melt parametrisation          ln_isfpar_mlt   = ', ln_isfpar_mlt
            IF ( ln_isfpar_mlt ) THEN
               WRITE(numout,*) '         isf parametrisation formulation         cn_isfpar_mlt   = ', TRIM(cn_isfpar_mlt)
            END IF
            WRITE(numout,*) ''
            !
            WRITE(numout,*) '      Coupling to an ice sheet model          ln_isfcpl       = ', ln_isfcpl
            IF ( ln_isfcpl ) THEN
               WRITE(numout,*) '         conservation activated ln_isfcpl_cons     = ', ln_isfcpl_cons
               WRITE(numout,*) '         number of call of the extrapolation loop  = ', nn_drown
            ENDIF
            WRITE(numout,*) ''
            !
         ELSE
            !
            IF ( ln_isfcav ) THEN
               WRITE(numout,*) ''
               WRITE(numout,*) '   W A R N I N G: ice shelf cavities are open BUT no melt will be computed or read from file !'
               WRITE(numout,*) ''
            END IF
            !
         END IF

         IF (ln_isfcav) THEN
            WRITE(numout,*) '      Ice shelf load method                   cn_isfload        = ', TRIM(cn_isfload)
            WRITE(numout,*) '         Temperature used to compute the ice shelf load            = ', rn_isfload_T
            WRITE(numout,*) '         Salinity    used to compute the ice shelf load            = ', rn_isfload_S
         END IF
         WRITE(numout,*) ''
         FLUSH(numout)

      END IF
      !

      !---------------------------------------------------------------------------------------------------------------------
      ! sanity check  ! issue ln_isfcav not yet known as well as l_isfoasis  => move this call in isf_stp ?
      ! melt in the cavity without cavity
      IF ( ln_isfcav_mlt .AND. (.NOT. ln_isfcav) ) &
          &   CALL ctl_stop('ice shelf melt in the cavity activated (ln_isfcav_mlt) but no cavity detected in domcfg (ln_isfcav), STOP' )
      !
      ! ice sheet coupling without cavity
      IF ( ln_isfcpl .AND. (.NOT. ln_isfcav) ) &
         &   CALL ctl_stop('coupling with an ice sheet model detected (ln_isfcpl) but no cavity detected in domcfg (ln_isfcav), STOP' )
      !
      IF ( ln_isfcpl .AND. ln_isfcpl_cons .AND. ln_linssh ) &
         &   CALL ctl_stop( 'The coupling between NEMO and an ice sheet model with the conservation option does not work with the linssh option' )
      !
      IF ( l_isfoasis .AND. .NOT. ln_isf ) CALL ctl_stop( ' OASIS send ice shelf fluxes to NEMO but NEMO does not have the isf module activated' )
      !
      IF ( l_isfoasis .AND. ln_isf ) THEN
         !
         CALL ctl_stop( 'namelist combination ln_cpl and ln_isf not tested' )
         !
         ! NEMO coupled to ATMO model with isf cavity need oasis method for melt computation 
         IF ( ln_isfcav_mlt .AND. TRIM(cn_isfcav_mlt) /= 'oasis' ) CALL ctl_stop( 'cn_isfcav_mlt = oasis is the only option availble if fwf send by oasis' )
         IF ( ln_isfpar_mlt .AND. TRIM(cn_isfpar_mlt) /= 'oasis' ) CALL ctl_stop( 'cn_isfpar_mlt = oasis is the only option availble if fwf send by oasis' )
         !
         ! oasis melt computation not tested (coded but not tested)
         IF ( ln_isfcav_mlt .OR. ln_isfpar_mlt ) THEN
            IF ( TRIM(cn_isfcav_mlt) == 'oasis' ) CALL ctl_stop( 'cn_isfcav_mlt = oasis not tested' )
            IF ( TRIM(cn_isfpar_mlt) == 'oasis' ) CALL ctl_stop( 'cn_isfpar_mlt = oasis not tested' )
         END IF
         !
         ! oasis melt computation with cavity open and cavity parametrised (not coded)
         IF ( ln_isfcav_mlt .AND. ln_isfpar_mlt ) THEN
            IF ( TRIM(cn_isfpar_mlt) == 'oasis' .AND. TRIM(cn_isfcav_mlt) == 'oasis' ) CALL ctl_stop( 'cn_isfpar_mlt = oasis and cn_isfcav_mlt = oasis not coded' )
         END IF
         !
         ! compatibility ice shelf and vvl
         IF( .NOT. ln_vvl_zstar .AND. ln_isf ) CALL ctl_stop( 'Only vvl_zstar has been tested with ice shelf cavity' )
         !
      END IF
   END SUBROUTINE isf_ctl

   
   SUBROUTINE isf_nam
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isf_nam  ***
      !!
      !! ** Purpose :   Read ice shelf namelist cfg and ref
      !!
      !!----------------------------------------------------------------------
      INTEGER               :: ios                  ! Local integer output status for namelist read
      !!----------------------------------------------------------------------
      NAMELIST/namisf/ ln_isf        ,                                                           & 
         &             cn_gammablk   , rn_gammat0    , rn_gammas0    , rn_htbl, sn_isfcav_fwf,   &
         &             ln_isfcav_mlt , cn_isfcav_mlt , sn_isfcav_fwf ,                           &
         &             ln_isfpar_mlt , cn_isfpar_mlt , sn_isfpar_fwf ,                           &
         &             sn_isfpar_zmin, sn_isfpar_zmax, sn_isfpar_Leff,                           &
         &             ln_isfcpl     , nn_drown      , ln_isfcpl_cons, ln_isfdebug, rn_vtide,    &
         &             cn_isfload    , rn_isfload_T  , rn_isfload_S  , cn_isfdir  ,              &
         &             rn_isfpar_bg03_gt0
      !!----------------------------------------------------------------------
      !
      READ  ( numnam_ref, namisf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namisf in reference namelist' )
      !
      READ  ( numnam_cfg, namisf, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namisf in configuration namelist' )
      IF(lwm) WRITE ( numond, namisf )

   END SUBROUTINE isf_nam
   !!
   !!======================================================================
END MODULE isfstp
