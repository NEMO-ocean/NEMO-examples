#!/bin/bash
#set -vx
# ncmax $var_nm $fl_nm : What is maximum of variable?
function ncmax { ncap2 -O -C -v -s "foo=${1}.max();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
# ncmin $var_nm $fl_nm : What is minimum of variable?
function ncmin { ncap2 -O -C -v -s "foo=${1}.min();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
# ncmdn $var_nm $fl_nm : What is median of variable?
function ncmdn { ncap2 -O -C -v -s "foo=gsl_stats_median_from_sorted_data(${1}.sort());print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }

##
## simple report generator for the test case
##

##
## Variables which may need to be adapted to your experiment:
##
#  RUNDIR = directory where the test case is executed: contains all outputs
RUNDIR=/gpfswork/rech/omr/romr001/OUT/CPLTESTCASE/2020-03-31120816
#  NB_NEMO_IT = expected total number of NEMO iterations
NB_NEMO_IT=160
#  NB_OASIS_OUTFILES = number of debug.root.0* OASIS output files
NB_OASIS_OUTFILES=2
##
## END of variables to be checked - Nothing need to be changed below
##
# check if directory is here
  if [ ! -d $RUNDIR ]; then
    printf "%-27s %s %s\n" $RUNDIR  "directory does not exist. Check RUNDIR variable in script. Stop" 
    return
  fi

cd $RUNDIR

echo " "
echo "Check results of test case in directory: " `pwd`
echo " "
##
## Check if OASIS execution has been successful
##
echo "      OASIS successful (true if OASIS outputs in debug.root.0? includes SUCCESSFUL RUN)  : "
count=0
for file in debug.root.0*
do
  echo $file ; grep "SUCCESSFUL RUN" $file
  count=`expr $count + 1`
done
echo "OASIS success checked on $count files"
[ $count = $NB_OASIS_OUTFILES ] && echo true || echo false
##
## Check if NEMO execution has been sucessful
##
echo " "
echo "      NEMO execution is successful if the run.stat file contains one line for each of NB_NEMO_IT iterations, indicating they have indeed been computed"
 if [ ! -f ./run.stat ]; then
   echo " the run.stat file does not exist: NEMO did not end its first time step"
   echo " NEMO UNSUCESSFUL. Stop"
   return
 fi
echo "From run.stat NEMO output file, NEMO has executed the 160 time steps:"
nemo_iterations=`wc -l ./run.stat | awk {'print $1'} `; [ $nemo_iterations = $NB_NEMO_IT ] && echo true || echo false

##
## Check mean value of sst field seen by toyatm
##
 if [ ! -f ./ATSSTSST_toyatm_01.nc ]; then
   echo " the ATSSTSST_toyatm_01.nc file does not exist: the test is not successful"
   echo " Test case UNSUCESSFUL. Stop"
   return
 fi
echo " "
echo "Examining ATSSTSST variable sea surface temperature as seen by toyatm, unit is degree Kelvin (min. should be around 271., max. around 302., median around 280.)" 
ASSTmin=`ncmin ATSSTSST ATSSTSST_toyatm_01.nc`
ASSTmax=`ncmax ATSSTSST ATSSTSST_toyatm_01.nc` 
ASSTmed=`ncmdn  ATSSTSST ATSSTSST_toyatm_01.nc`
echo "Minimum value of ATSSTSST variable in ATSSTSST_toyatm_01.nc file = " $ASSTmin
echo "Maximum value of ATSSTSST variable in ATSSTSST_toyatm_01.nc file = " $ASSTmax
echo "Median value of ATSSTSST variable in ATSSTSST_toyatm_01.nc file = "  $ASSTmed
MINMAX=0
if [ ${ASSTmin%%.*} -lt 270 -o ${ASSTmax%%.*} -gt 310 ]; then
echo " Min. or max. values of ATSSTSST do not look reasonable. Check the test again "
MINMAX=1
fi
##
## Summary
##
echo " "  
if [ $count = $NB_OASIS_OUTFILES ] && [ $nemo_iterations = $NB_NEMO_IT ] && [ $MINMAX = 0 ] 
then
  echo " The run looks very succesful!"
  echo " Have a look to the ASTSSTSST.nc file (using ncview for example): sea surface temperatures as seen by the toyatm and compare it to the reference file (CPL/ref_ATSSTSST_last_time_step.jpg) "
  echo " Units are in degrees Kelvin and it will confirm the test is successful"
  echo " "
else
  echo "The test case is unsuccessful. Check all inputs and outputs"
fi
