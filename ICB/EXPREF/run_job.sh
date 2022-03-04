#!/bin/bash 
#OAR -n nemo_sette
#OAR -l core=2,walltime=00:05:00
#OAR -O nemo_sette.o%jobid%
#OAR -E nemo_sette.e%jobid%
#OAR --project tipaccs

. ~/bin/load_intelmodule.bash

  export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
  export OMP_NUM_THREADS=1
  cd $PBS_O_WORKDIR
  export XIO_HOME=$xios_path
#
  echo " ";
  export OMP_NUM_THREADS=1
  O_PER_NODE=2
  X_PER_NODE=2
  OCORES=2
  XCORES=0
  if [ $OCORES -le 32 ] ; then O_PER_NODE=$OCORES; fi
  export SETTE_DIR=/bettik/mathiotp/NEMO/NEMO_dev/branches/2020/tickets_icb_1900/sette

  export EXE_DIR=/bettik/mathiotp/NEMO/NEMO_dev/branches/2020/tickets_icb_1900/tests/ICB/EXP00/
  ulimit -c unlimited
  ulimit -s unlimited
#
# end of set up
###############################################################
#
# change to the working directory 
#
  cd $EXE_DIR
  echo Directory is `pwd`
  
  if [ $XCORES -gt 0 ]; then
#
#  Run MPMD case
#
     #XIOS will run on a separate node so will run in parallel queue
     if [ ! -f ./xios_server.exe ] && [ -f ${XIO_HOME}/bin/xios_server.exe ]; then
        cp ${XIO_HOME}/bin/xios_server.exe .
     fi
     if [ ! -f ./xios_server.exe ]; then
        echo "./xios_server.exe not found"
        echo "run aborted"
        exit
     fi
       echo time aprun -b -n $XCORES -N $X_PER_NODE ./xios_server.exe : -n $OCORES -N $O_PER_NODE ./nemo
            time aprun -b -n $XCORES -N $X_PER_NODE ./xios_server.exe : -n $OCORES -N $O_PER_NODE ./nemo
#
  else
#
# Run SPMD case
#
    echo time mpirun -n $OCORES ./nemo
         time mpirun -n $OCORES ./nemo
  fi
#
