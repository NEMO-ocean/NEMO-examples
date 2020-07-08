#!/bin/bash

################################################################
#
# Script to launch a set of STATION_ASF simulations
#
# L. Brodeau, 2020
#
################################################################

# What directory inside "tests" actually contains the compiled "nemo.exe" for STATION_ASF ?
TC_DIR="STATION_ASF2"

expdir=`basename ${PWD}`; # we expect "EXPREF" or "EXP00" normally...

# NEMOGCM root directory:
NEMO_ROOT_DIR="/home/laurent/DEV/NEMO/NEMOGCM_trunk"

# NEMOGCM root directory where to fetch compiled STATION_ASF nemo.exe:
SASF_WRK_DIR="${NEMO_ROOT_DIR}/tests/${TC_DIR}"

# Directory where to run the simulation:
PROD_DIR="${HOME}/tmp/STATION_ASF_gh"


####### End of normal user configurable section #######

#================================================================================

SASF_REF_DIR=`pwd | sed -e "s|/${expdir}||g"` ; # STATION_ASF directory where to work from, normally in here!
if [ ! -d ${SASF_REF_DIR} ]; then echo " Mhhh, no EXPREF directory ${SASF_REF_DIR} !"; exit; fi

# NEMO executable to use is:
NEMO_EXE="${SASF_WRK_DIR}/BLD/bin/nemo.exe"
if [ ! -f ${NEMO_EXE} ]; then echo " Mhhh, no compiled 'nemo.exe' found into `dirname ${NEMO_EXE}` !"; exit; fi

DATA_IN_DIR="${SASF_REF_DIR}/input_data" ; # Directory containing sea-surface + atmospheric input data
if [ ! -d ${DATA_IN_DIR} ]; then echo "PROBLEM!!! => did not find directory 'input_data' with input forcing..."; exit; fi

SASF_EXPREF=${SASF_REF_DIR}/${expdir}  ; # STATION_ASF EXPREF directory from which to use namelists and XIOS xml files...
if [ ! -d ${SASF_EXPREF} ]; then echo " Mhhh, no ${expdir} directory ${SASF_EXPREF} !"; exit; fi


echo "###########################################################"
echo "#        S T A T I O N   A i r  -  S e a   F l u x        #"
echo "###########################################################"
echo
echo "  * NEMO reference root directory is: ${NEMO_ROOT_DIR}"
echo "  * STATION_ASF work directory is: ${SASF_WRK_DIR}"
echo "       ==> NEMO EXE to use: ${NEMO_EXE}"
echo
echo "  * Input forcing data into: ${DATA_IN_DIR}"
echo "  * Production will be done into: ${PROD_DIR}"
echo "  * Directory in which namelists and xml files are fetched:"
echo "       ==> ${SASF_EXPREF}"
echo

mkdir -p ${PROD_DIR}

rsync -avP ${NEMO_EXE}          ${PROD_DIR}/

for ff in "context_nemo.xml" "file_def_nemo-oce.xml" "iodef.xml"; do
    if [ ! -f ${SASF_EXPREF}/${ff} ]; then echo " Mhhh, ${ff} not found into ${SASF_EXPREF} !"; exit; fi
    rsync -avPL ${SASF_EXPREF}/${ff} ${PROD_DIR}/
done

# Getting reference/defaults files from reference NEMO distro:
rdir="${NEMO_ROOT_DIR}/cfgs/SHARED"
for ff in "domain_def_nemo.xml" "grid_def_nemo.xml" "field_def_nemo-oce.xml" "namelist_ref"; do
    if [ ! -f ${rdir}/${ff} ]; then echo " Mhhh, ${ff} not found into ${rdir} !"; exit; fi
    ln -sf ${rdir}/${ff} ${PROD_DIR}/.
done


# Copy forcing to work directory:
rsync -avP ${DATA_IN_DIR}/Station_PAPA_50N-145W*.nc ${PROD_DIR}/

for CASE in "ECMWF" "COARE3p6" "NCAR" "ECMWF-noskin" "COARE3p6-noskin"; do

    echo ; echo
    echo "============================="
    echo " Going for ${CASE} experiment"
    echo "============================="
    echo

    scase=`echo "${CASE}" | tr '[:upper:]' '[:lower:]'`

    rm -f ${PROD_DIR}/namelist_cfg
    rsync -avPL ${SASF_EXPREF}/namelist_${scase}_cfg ${PROD_DIR}/namelist_cfg

    cd ${PROD_DIR}/
    echo
    echo "Launching NEMO !"
    ./nemo.exe 1>out_nemo.out 2>err_nemo.err
    echo "Done!"
    echo

    # Moving output files:
    mkdir -p output
    mv -f STATION_ASF-${CASE}_*_grid*.nc output/
    
    # Saving logs:
    mkdir -p ${CASE}_log
    mv -f *.out *.err ocean.output output.namelist.dyn ${CASE}_log/

done
