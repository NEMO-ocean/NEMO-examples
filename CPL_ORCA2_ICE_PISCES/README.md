#  CPL_ORCA2_ICE_PISCES

This configuration is similar to ORCA2_ICE_PISCES but can be coupled, via OASIS3-MCT (v3 or v4) to a toy atmosphere. Sources of the toy atmosphere model are located at tools/TOYATM. Both NEMO test configuration and toy sources must be linked with OASIS library. In addition to ORCA2_ICE_PISCES input files, files from tools/TOYATM/EXP directory (OASIS auxiliary and parameter files) must be copied to the working directory. Coupled system must be launched following MPI MPMD syntax (e.g. mpirun -n 20 nemo.exe : -n 20 toyatm.exe).

