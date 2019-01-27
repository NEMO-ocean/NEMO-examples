MODULE trazdf
   !!==============================================================================
   !!                 ***  MODULE  trazdf  ***
   !! Ocean active tracers:  vertical component of the tracer mixing trend
   !!==============================================================================
   !! History :  1.0  !  2005-11  (G. Madec)  Original code
   !!            3.0  !  2008-01  (C. Ethe, G. Madec)  merge TRC-TRA
   !!            4.0  !  2017-06  (G. Madec)  remove explict time-stepping option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_zdf       : Update the tracer trend with the vertical diffusion
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables 
   USE domvvl         ! variable volume
   USE phycst         ! physical constant
   USE zdf_oce        ! ocean vertical physics variables
   USE sbc_oce        ! surface boundary condition: ocean
   USE ldftra         ! lateral diffusion: eddy diffusivity
   USE ldfslp         ! lateral diffusion: iso-neutral slope 
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends: tracer trend manager
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_zdf       ! called by step.F90
   PUBLIC   tra_zdf_imp   ! called by trczdf.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: trazdf.F90 10572 2019-01-24 15:37:13Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_zdf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf  ***
      !!
      !! ** Purpose :   compute the vertical ocean tracer physics.
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   jk   ! Dummy loop indices
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztrdt, ztrds   ! 3D workspace
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('tra_zdf')
      !
      IF( kt == nit000 )  THEN
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) 'tra_zdf : implicit vertical mixing on T & S'
         IF(lwp)WRITE(numout,*) '~~~~~~~ '
      ENDIF
      !
      IF( neuler == 0 .AND. kt == nit000 ) THEN   ;   r2dt =      rdt   ! at nit000, =   rdt (restarting with Euler time stepping)
      ELSEIF( kt <= nit000 + 1           ) THEN   ;   r2dt = 2. * rdt   ! otherwise, = 2 rdt (leapfrog)
      ENDIF
      !
      IF( l_trdtra )   THEN                  !* Save ta and sa trends
         ALLOCATE( ztrdt(jpi,jpj,jpk) , ztrds(jpi,jpj,jpk) )
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal)
      ENDIF
      !
      !                                      !* compute lateral mixing trend and add it to the general trend
      CALL tra_zdf_imp( kt, nit000, 'TRA', r2dt, tsb, tsa, jpts ) 

!!gm WHY here !   and I don't like that !
      ! DRAKKAR SSS control {
      ! JMM avoid negative salinities near river outlet ! Ugly fix
      ! JMM : restore negative salinities to small salinities:
!!$      WHERE( tsa(:,:,:,jp_sal) < 0._wp )   tsa(:,:,:,jp_sal) = 0.1_wp
!!gm

      IF( l_trdtra )   THEN                      ! save the vertical diffusive trends for further diagnostics
         DO jk = 1, jpkm1
            ztrdt(:,:,jk) = ( ( tsa(:,:,jk,jp_tem)*e3t_a(:,:,jk) - tsb(:,:,jk,jp_tem)*e3t_b(:,:,jk) ) &
               &          / (e3t_n(:,:,jk)*r2dt) ) - ztrdt(:,:,jk)
            ztrds(:,:,jk) = ( ( tsa(:,:,jk,jp_sal)*e3t_a(:,:,jk) - tsb(:,:,jk,jp_sal)*e3t_b(:,:,jk) ) &
              &           / (e3t_n(:,:,jk)*r2dt) ) - ztrds(:,:,jk)
         END DO
!!gm this should be moved in trdtra.F90 and done on all trends
         CALL lbc_lnk_multi( 'trazdf', ztrdt, 'T', 1. , ztrds, 'T', 1. )
!!gm
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_zdf, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_zdf, ztrds )
         DEALLOCATE( ztrdt , ztrds )
      ENDIF
      !                                          ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' zdf  - Ta: ', mask1=tmask,               &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( ln_timing )   CALL timing_stop('tra_zdf')
      !
   END SUBROUTINE tra_zdf

 
   SUBROUTINE tra_zdf_imp( kt, kit000, cdtype, p2dt, ptb, pta, kjpt ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf_imp  ***
      !!
      !! ** Purpose :   Compute the after tracer through a implicit computation
      !!     of the vertical tracer diffusion (including the vertical component 
      !!     of lateral mixing (only for 2nd order operator, for fourth order 
      !!     it is already computed and add to the general trend in traldf) 
      !!
      !! ** Method  :  The vertical diffusion of a tracer ,t , is given by:
      !!          difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(t) )
      !!      It is computed using a backward time scheme (t=after field)
      !!      which provide directly the after tracer field.
      !!      If ln_zdfddm=T, use avs for salinity or for passive tracers
      !!      Surface and bottom boundary conditions: no diffusive flux on
      !!      both tracers (bottom, applied through the masked field avt).
      !!      If iso-neutral mixing, add to avt the contribution due to lateral mixing.
      !!
      !! ** Action  : - pta  becomes the after tracer
      !!---------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000   ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype   ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt     ! number of tracers
      REAL(wp)                             , INTENT(in   ) ::   p2dt     ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb      ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta      ! in: tracer trend ; out: after tracer field
      !
      INTEGER  ::  ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::  zrhs, zzwi, zzws ! local scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zwi, zwt, zwd, zws
      !!---------------------------------------------------------------------
      !
      !                                               ! ============= !
      DO jn = 1, kjpt                                 !  tracer loop  !
         !                                            ! ============= !
         !  Matrix construction
         ! --------------------
         ! Build matrix if temperature or salinity (only in double diffusion case) or first passive tracer
         !
         IF(  ( cdtype == 'TRA' .AND. ( jn == jp_tem .OR. ( jn == jp_sal .AND. ln_zdfddm ) ) ) .OR.   &
            & ( cdtype == 'TRC' .AND. jn == 1 )  )  THEN
            !
            ! vertical mixing coef.: avt for temperature, avs for salinity and passive tracers
            IF( cdtype == 'TRA' .AND. jn == jp_tem ) THEN   ;   zwt(:,:,2:jpk) = avt(:,:,2:jpk)
            ELSE                                            ;   zwt(:,:,2:jpk) = avs(:,:,2:jpk)
            ENDIF
            zwt(:,:,1) = 0._wp
            !
            IF( l_ldfslp ) THEN            ! isoneutral diffusion: add the contribution 
               IF( ln_traldf_msc  ) THEN     ! MSC iso-neutral operator 
                  DO jk = 2, jpkm1
                     DO jj = 2, jpjm1
                        DO ji = fs_2, fs_jpim1   ! vector opt.
                           zwt(ji,jj,jk) = zwt(ji,jj,jk) + akz(ji,jj,jk)  
                        END DO
                     END DO
                  END DO
               ELSE                          ! standard or triad iso-neutral operator
                  DO jk = 2, jpkm1
                     DO jj = 2, jpjm1
                        DO ji = fs_2, fs_jpim1   ! vector opt.
                           zwt(ji,jj,jk) = zwt(ji,jj,jk) + ah_wslp2(ji,jj,jk)
                        END DO
                     END DO
                  END DO
               ENDIF
            ENDIF
            !
            ! Diagonal, lower (i), upper (s)  (including the bottom boundary condition since avt is masked)
            IF( ln_zad_Aimp ) THEN         ! Adaptive implicit vertical advection
               DO jk = 1, jpkm1
                  DO jj = 2, jpjm1
                     DO ji = fs_2, fs_jpim1   ! vector opt. (ensure same order of calculation as below if wi=0.)
                        zzwi = - p2dt * zwt(ji,jj,jk  ) / e3w_n(ji,jj,jk  )
                        zzws = - p2dt * zwt(ji,jj,jk+1) / e3w_n(ji,jj,jk+1)
                        zwd(ji,jj,jk) = e3t_a(ji,jj,jk) - zzwi - zzws   &
                           &                 + p2dt * ( MAX( wi(ji,jj,jk  ) , 0._wp ) - MIN( wi(ji,jj,jk+1) , 0._wp ) )
                        zwi(ji,jj,jk) = zzwi + p2dt *   MIN( wi(ji,jj,jk  ) , 0._wp )
                        zws(ji,jj,jk) = zzws - p2dt *   MAX( wi(ji,jj,jk+1) , 0._wp )
                    END DO
                  END DO
               END DO
            ELSE
               DO jk = 1, jpkm1
                  DO jj = 2, jpjm1
                     DO ji = fs_2, fs_jpim1   ! vector opt.
                        zwi(ji,jj,jk) = - p2dt * zwt(ji,jj,jk  ) / e3w_n(ji,jj,jk)
                        zws(ji,jj,jk) = - p2dt * zwt(ji,jj,jk+1) / e3w_n(ji,jj,jk+1)
                        zwd(ji,jj,jk) = e3t_a(ji,jj,jk) - zwi(ji,jj,jk) - zws(ji,jj,jk)
                    END DO
                  END DO
               END DO
            ENDIF
            !
            !! Matrix inversion from the first level
            !!----------------------------------------------------------------------
            !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
            !
            !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
            !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
            !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
            !        (        ...               )( ...  ) ( ...  )
            !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
            !
            !   m is decomposed in the product of an upper and lower triangular matrix.
            !   The 3 diagonal terms are in 3d arrays: zwd, zws, zwi.
            !   Suffices i,s and d indicate "inferior" (below diagonal), diagonal
            !   and "superior" (above diagonal) components of the tridiagonal system.
            !   The solution will be in the 4d array pta.
            !   The 3d array zwt is used as a work space array.
            !   En route to the solution pta is used a to evaluate the rhs and then 
            !   used as a work space array: its value is modified.
            !
            DO jj = 2, jpjm1        !* 1st recurrence:   Tk = Dk - Ik Sk-1 / Tk-1   (increasing k)
               DO ji = fs_2, fs_jpim1            ! done one for all passive tracers (so included in the IF instruction)
                  zwt(ji,jj,1) = zwd(ji,jj,1)
               END DO
            END DO
            DO jk = 2, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwt(ji,jj,jk-1)
                  END DO
               END DO
            END DO
            !
         ENDIF 
         !         
         DO jj = 2, jpjm1           !* 2nd recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
            DO ji = fs_2, fs_jpim1
               pta(ji,jj,1,jn) = e3t_b(ji,jj,1) * ptb(ji,jj,1,jn) + p2dt * e3t_n(ji,jj,1) * pta(ji,jj,1,jn)
            END DO
         END DO
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  zrhs = e3t_b(ji,jj,jk) * ptb(ji,jj,jk,jn) + p2dt * e3t_n(ji,jj,jk) * pta(ji,jj,jk,jn)   ! zrhs=right hand side
                  pta(ji,jj,jk,jn) = zrhs - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) * pta(ji,jj,jk-1,jn)
               END DO
            END DO
         END DO
         !
         DO jj = 2, jpjm1           !* 3d recurrence:    Xk = (Zk - Sk Xk+1 ) / Tk   (result is the after tracer)
            DO ji = fs_2, fs_jpim1
               pta(ji,jj,jpkm1,jn) = pta(ji,jj,jpkm1,jn) / zwt(ji,jj,jpkm1) * tmask(ji,jj,jpkm1)
            END DO
         END DO
         DO jk = jpk-2, 1, -1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  pta(ji,jj,jk,jn) = ( pta(ji,jj,jk,jn) - zws(ji,jj,jk) * pta(ji,jj,jk+1,jn) )   &
                     &             / zwt(ji,jj,jk) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !                                            ! ================= !
      END DO                                          !  end tracer loop  !
      !                                               ! ================= !
   END SUBROUTINE tra_zdf_imp

   !!==============================================================================
END MODULE trazdf
