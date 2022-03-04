#ifdef key_gpu
MODULE gpu_manager
    USE lib_mpp        ! distributed memory computing library
    USE cudafor
    USE par_oce
    !!======================================================================
    !!                       ***  MODULE  gpu_manager  ***
    !! GPU CUDA utilities : Defines run parameters for GPGPU routines
    !!=====================================================================
    !! History :   1.0  !  2019-10  (M. Faria)   original code
    !!----------------------------------------------------------------------

    !!----------------------------------------------------------------------
    IMPLICIT NONE
    PUBLIC
    !!----------------------------------------------------------------------
    !!                   gpu configuration parameters
    !!----------------------------------------------------------------------
    INTEGER          :: cudaistat, usecudev
    INTEGER          :: nodesize
    INTEGER          :: cudadev
    INTEGER          :: cusocket
CONTAINS
    SUBROUTINE setdevice()
        nodesize  = 1                                                    !Node number of MPI processes, CTE-Power default 40
        cusocket  = 1                                                    !Number of NVlink or PCI-E nodes, CTE-Power default 2
        cudadev   = 1                                                    !Number of GPUs on a node, CTE-Power default 4
        usecudev  = MOD( mpprank /  ( nodesize / cusocket )  , cudadev ) !Bind MPI RANK to GPU, consider balanced NVlink and constant number of devices per node 
        cudaistat = cudaSetDevice( usecudev  )                           !Set which GPU use. Ex for 4 GPUs. cpus_per_socket / sockets_per_node for "bind-to core"
        print *, 'CUDA istat ', cudaistat, 'rank', mpprank, 'cuda device', cudadev 
    END SUBROUTINE setdevice

END MODULE gpu_manager
#endif
