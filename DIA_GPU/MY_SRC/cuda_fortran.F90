#ifdef key_gpu
MODULE cuda_fortran
    USE cudafor
    USE dom_oce
    USE lib_mpp
    CONTAINS
      ATTRIBUTES(global) &
      SUBROUTINE dia_hsb_kernel(surf, e3t, surf_ini, e3t_ini, ts, hc_loc_ini, sc_loc_ini,     &
                            & tmask, tmask_ini, zwrkv, zwrkh, zwrks, zwrk, jpi, jpj, jpk, jpt, Kmm)
            !
            INTEGER, VALUE :: jpi, jpj, jpk, jpt                                                                         ! MPI process sub-domain and stream tile
            INTEGER, VALUE :: Kmm                                                                                   !ocean time level indice
            REAL(8)        :: surf(jpi, jpj), e3t(jpi, jpj, jpk, jpt), surf_ini(jpi, jpj), e3t_ini(jpi, jpj, jpk), ts(jpi, jpj, jpk, 2, jpt)
            REAL(8)        :: hc_loc_ini(jpi, jpj, jpk), sc_loc_ini(jpi, jpj, jpk), tmask(jpi, jpj, jpk), tmask_ini(jpi, jpj, jpk)
            REAL(8)        :: zwrkv(jpi, jpj, jpk), zwrkh(jpi, jpj, jpk), zwrks(jpi, jpj, jpk), zwrk(jpi, jpj, jpk)
            !

            INTEGER        :: i, j, k                       ! dummy indexes

            !
            i = blockDim%x * (blockIdx%x -1) + threadIdx%x
            j = blockDim%y * (blockIdx%y -1) + threadIdx%y
            k = blockDim%z * (blockIdx%z -1) + threadIdx%z
            !
            IF ( (i .le. jpi ) .AND. (j .le. jpj) .AND. (k .le. jpk-1) ) THEN
                zwrkv(i, j, k) =   surf    (i,j) * e3t    (i,j,k,Kmm) * tmask    (i,j,k)      &
                &                - surf_ini(i,j) * e3t_ini(i,j,k)     * tmask_ini(i,j,k)
                
                zwrkh(i, j, k) = ( surf(i, j) * e3t(i, j, k, Kmm) * ts(i, j, k, 1, Kmm) - surf_ini(i, j) * hc_loc_ini(i, j, k) )
                
                zwrks(i, j, k) = ( surf(i, j) * e3t(i, j, k, Kmm) * ts(i, j, k, 2, Kmm) - surf_ini(i, j) * sc_loc_ini(i, j, k) )
                
                zwrk(i, j, k)  =   surf(i, j) * e3t(i, j, k, Kmm) * tmask(i, j, k)
            END IF
      END SUBROUTINE dia_hsb_kernel

!     ATTRIBUTES(global) &
!     SUBROUTINE dia_hsb_kernel1d(surf, e3t_n, surf_ini, e3t_ini, &
!        &  tsn, hc_loc_ini, sc_loc_ini, tmask, zwrkv, zwrkh, zwrks, zwrk, jpi, jpj, jpk)
!        IMPLICIT NONE
!        REAL(kind=8) :: surf(:, :), e3t_n(:, :, :), surf_ini(:, :), e3t_ini(:, :, :), tsn(:, :, :, :), &
!                            & hc_loc_ini(:, :, :), sc_loc_ini(:, :, :), tmask(:, :, :), zwrkv(:, :, :),      &
!                            & zwrkh(:, :, :), zwrks(:, :, :), zwrk(:, :, :)
!        REAL(kind=8), SHARED       :: sdata(*)
!        INTEGER, VALUE             :: jpi, jpj, jpk
!        INTEGER                    :: i, si, ti, globsize
!
!        globsize = jpi*jpj*jpk
!
!        i  = blockDim%x * (blockIdx%x-1) + threadIdx%x
!        ti = threadIdx%x
!        si = MOD(i, jpk)+1 !Shared memory indexing
!
!        !Compute volume (zwrk), volume variation (zwrkv), heat deviation (zwrkh)
!        !and salinity deviation (zwrks)
!        IF ( i .le. globsize ) THEN
!            sdata(ti) = surf(si) * e3t_n(i)
!            zwrkv(i) = ( sdata(ti) - surf_ini(si) * e3t_ini(i) ) * tmask(i) * surf(si)
!            zwrkh(i) = ( sdata(ti) * tsn(i) - surf_ini(si) * hc_loc_ini(i) ) * tmask(i) * surf(si)
!            zwrks(i) = ( sdata(ti) * tsn(globsize + i) - surf_ini(si) * sc_loc_ini(i) ) &
!                    * tmask(i) * surf(si)
!            zwrk(i) = sdata(ti) * tmask(i) * surf(si)
!        END IF
!     END SUBROUTINE dia_hsb_kernel1d


       ATTRIBUTES(global) &
       SUBROUTINE filter_cuda(ptab, mask, jpi, jpj, jpk)
           !
           REAL(kind=8)         :: ptab(jpi,jpj,jpk) , mask(jpi,jpj)
           !
           INTEGER, VALUE  :: jpi, jpj, jpk                 ! MPI process sub-domain
           INTEGER         :: i, j, k                 ! dummy indexes
           !
           i = blockDim%x * (blockIdx%x -1) + threadIdx%x
           j = blockDim%y * (blockIdx%y -1) + threadIdx%y
           k = blockDim%z * (blockIdx%z -1) + threadIdx%z
           !
           tile = 1
           IF ( (i .le. jpi ) .AND. (j .le. jpj) .AND. (k .le. jpk-1) ) THEN
               ptab(i, j, k)  = ptab(i, j, k) * mask(i, j)
           END IF
       END SUBROUTINE filter_cuda

       ATTRIBUTES(global) SUBROUTINE array3dto1d(d_inp,d_out,ipi, ipj, ipk)
           !
           REAL(kind=8)                 :: d_inp(:,:,:)
           REAL(kind=8), intent(out)    :: d_out(:)
           !
           INTEGER, VALUE  :: ipi, ipj, ipk                 ! MPI process sub-domain
           INTEGER         :: i, j, k                       ! dummy indexes
           !
           i = blockDim%x * (blockIdx%x -1) + threadIdx%x
           j = blockDim%y * (blockIdx%y -1) + threadIdx%y
           k = blockDim%z * (blockIdx%z -1) + threadIdx%z
           !
           IF ( (i .le. ipi ) .AND. (j .le. ipj) .AND. (k .le. ipk) ) THEN
               d_out(i + ipi*(j-1) + ipi*ipj*(k-1) )  = d_inp(i, j, k)
           END IF
       END SUBROUTINE array3dto1d

       !Knuth's trick
       ATTRIBUTES(device) &
       SUBROUTINE DDPDD_d(ydda, yddb)
           IMPLICIT NONE
           COMPLEX(kind=8), INTENT(in   ) :: ydda              !Scalar to add
           COMPLEX(kind=8), INTENT(inout) :: yddb              !Total sum
           REAL(kind=8)                   :: zerr, zt1, zt2
       
           zt1  = REAL(ydda) + REAL(yddb)
           zerr = zt1 - REAL(ydda)
           zt2  = ( (REAL(yddb) - zerr) + (REAL(ydda) - (zt1 - zerr)) ) &
                 & + AIMAG(ydda) + AIMAG(yddb)
           !The result is t1 + t2 after normalization
           yddb = CMPLX( zt1 + zt2, zt2 - ((zt1 + zt2) - zt1), 8 )
       END SUBROUTINE DDPDD_d       
END MODULE
#endif
