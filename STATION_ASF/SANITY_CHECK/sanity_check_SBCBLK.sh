#!/bin/bash

################################################################
#
# Script to launch a set of STATION_ASF simulations
#
#  Ocean-only version....
#
# L. Brodeau, 2020
#
#
# Note: the current setup only performs an ocean-only test (no sea-ice),
#       even though the current script already recongnises if the nemo.exe
#       has been compiled with "key_si3". For now it will just test the
#       fluxes over open water.
#
################################################################

CONFIG="STATION_ASF" ; # general name of the configuration
CONFIG_BLD="STATION_ASF2" ; # name of config as build in NEMO... (directory inside "tests" actually contains the compiled test-case?)

# Atmo + SSX forcing to use and sea-ice support?
FORCING="IDEALIZED" ; i_sea_ice=0 ; SFORC="input_output_VALIDATION_IDEALIZED.nc"


# Root directory NEMOGCM reference depository where to fetch compiled STATION_ASF nemo.exe + default namelists:
NEMO_REF_DIR="`dirname ${PWD} | sed -e 's|/tests/STATION_ASF||g'`" ; # that should normally do the trick!

# NEMOGCM root directory where to fetch compiled STATION_ASF nemo.exe:
SASF_WRK_DIR="${NEMO_REF_DIR}/tests/${CONFIG_BLD}"

# DATA_IN_DIR => Directory containing sea-surface + atmospheric forcings:
DATA_IN_DIR="${NEMO_REF_DIR}/tests/${CONFIG}/input_data"

# Directory where to run the simulation:
PROD_DIR="${HOME}/tmp/${CONFIG}_test"

# MPI launch command for your system
# Even though we are running on 1 proc here, NEMO has likely been compiled with MPI,
# if it's not the case then just leave "MPI_LAUNCH" void...
MPI_LAUNCH="mpirun -n 1"

####### End of normal user configurable section #######
#================================================================================



# Should the `analyze_output.py` script provide more output ???
cmore=""
if [ "$1" != "" ]; then
    if [ "$1" = "more" ]; then
        cmore="m"
    else
        echo "Only 'more' is accepted as argument for script `basename $0` !"
        exit
    fi
fi

rm -f SBCBLK.success SBCBLK.fail

HERE=`pwd`

echo

# Is the nemo.exe of STATION_ASF compiled with sea-ice support (SI3) ?
i_si3=0
FCPP="${SASF_WRK_DIR}/cpp_${CONFIG_BLD}.fcm"
if [ ! -f ${FCPP} ]; then echo " Mhhh, we could not find 'cpp_STATION_ASF.fcm' into `dirname ${FCPP}` !"; exit; fi
ca=`cat ${FCPP} | grep 'key_si3'`

if [ "${ca}" = "" ]; then
    echo " *** NEMO was NOT compiled with sea-ice support (SI3) !" ; echo
    if [ ${i_sea_ice} -eq 1 ]; then
        echo " ===> so you cannot request ice-air flux computation !"
        echo "     ===> please set i_sea_ice=0 or compile STATION_ASF with CPP key 'key_si3' !"
        echo ; exit
    fi
else
    echo " *** NEMO was apparently compiled with sea-ice support (SI3) !" ; echo
    i_si3=1
fi


DIR_NL=${FORCING}/oce ; # directory where to find the namelists...

# NEMO executable to use is:
NEMO_EXE="${SASF_WRK_DIR}/BLD/bin/nemo.exe"
if [ ! -f ${NEMO_EXE} ]; then echo " Mhhh, no compiled 'nemo.exe' found into `dirname ${NEMO_EXE}` !"; exit; fi

SASF_EXPREF="${NEMO_REF_DIR}/tests/STATION_ASF/EXPREF" ; # STATION_ASF EXPREF directory from which to fetch and XIOS xml files...

CFGS_SHARED="${NEMO_REF_DIR}/cfgs/SHARED"
if [ ! -d ${CFGS_SHARED} ]; then echo "PROBLEM!!! => could not find directory ${CFGS_SHARED} !"; exit; fi

if [ ! -d ${DATA_IN_DIR} ]; then echo "PROBLEM!!! => could not find directory 'input_data' with input forcing..."; exit; fi

cdt_cmpl="`\ls -l ${NEMO_EXE} | cut -d' ' -f 6,7,8`"

echo "###########################################################"
echo "#        S T A T I O N   A i r  -  S e a   F l u x        #"
echo "###########################################################"
echo
echo "  * NEMO reference root directory is: ${NEMO_REF_DIR}"
echo "  * STATION_ASF work directory is: ${SASF_WRK_DIR}"
echo "       ==> NEMO EXE to use: ${NEMO_EXE}"
echo "           ==> compiled: ${cdt_cmpl} !"
echo
echo "  * Input forcing data into: ${DATA_IN_DIR}"
echo "  * Production will be done into: ${PROD_DIR}"
echo "  * Directory in which namelists and xml files are fetched:"
echo "       ==> ${SASF_EXPREF}"
echo
sleep 2

mkdir -p ${PROD_DIR}

rsync -avP ${NEMO_EXE}          ${PROD_DIR}/
ln -sf     ${NEMO_EXE}          ${PROD_DIR}/nemo.exe.lnk


# XIOS xml file
################

list_xml_ref="field_def_nemo-oce.xml domain_def_nemo.xml grid_def_nemo.xml"
list_xml_cfg="iodef.xml file_def_nemo-oce.xml"
fcntxt="context_nemo_OCE.xml"
if [ ${i_sea_ice} -eq 1 ]; then
    list_xml_ref+=" field_def_nemo-ice.xml"
    list_xml_cfg+=" file_def_nemo-ice.xml"
    fcntxt="context_nemo_OCE+ICE.xml"
fi

# The "context_nemo.xml" file:
if [ ! -f ${SASF_EXPREF}/${fcntxt} ]; then echo " Mhhh, ${fcntxt} not found into ${SASF_EXPREF} !"; exit; fi
rsync -avPL ${SASF_EXPREF}/${fcntxt} ${PROD_DIR}/context_nemo.xml

# All remaining "*.xml" files:
for ff in ${list_xml_cfg} ; do
    if [ ! -f ${SASF_EXPREF}/${ff} ]; then echo " Mhhh, ${ff} not found into ${SASF_EXPREF} !"; exit; fi
    rsync -avPL ${SASF_EXPREF}/${ff} ${PROD_DIR}/
done
for ff in ${list_xml_ref} ; do
    if [ ! -f ${CFGS_SHARED}/${ff} ]; then echo " Mhhh, ${ff} not found into ${CFGS_SHARED} !"; exit; fi
    rsync -avPL ${CFGS_SHARED}/${ff} ${PROD_DIR}/
done


# Copy forcing to work directory:
echo ; echo "Forcing files to use:"
\ls -l ${DATA_IN_DIR}/${SFORC}
echo
rsync -avP ${DATA_IN_DIR}/${SFORC} ${PROD_DIR}/
echo; echo

for CASE in "NCAR" "ECMWF" "COARE3p0" "COARE3p6" "ANDREAS"; do

        echo ; echo
        echo "======================================================================="
        echo " Going for experiment: ${CASE} bulk param. with ${FORCING} forcing "
        echo "======================================================================="
        echo

        scase=`echo "${CASE}" | tr '[:upper:]' '[:lower:]'`

        rm -f ${PROD_DIR}/namelist_cfg
        fnml="${HERE}/${DIR_NL}/namelist_${scase}_cfg"
        if [ ! -f ${fnml} ]; then echo " Mhhh, test namelist ${fnml} not found !"; exit; fi

        echo "   ===> namelist to use is: ${DIR_NL}/namelist_${scase}_cfg !"; echo
        
        # The namelists:
        rsync -avPL ${HERE}/${DIR_NL}/namelist_${scase}_cfg        ${PROD_DIR}/namelist_cfg
        rsync -avPL ${CFGS_SHARED}/namelist_ref                    ${PROD_DIR}/namelist_ref
    
        cd ${PROD_DIR}/
        echo
        echo "Launching NEMO ! (==> ${MPI_LAUNCH} ./nemo.exe)"
        ${MPI_LAUNCH} ./nemo.exe 1>out_nemo.out 2>err_nemo.err
        echo "Done!"
        echo

        # Moving output files:
        mkdir -p output
        mv -f ${CONFIG}-${CASE}_${FORCING}_*_grid*.nc output/

        # Saving logs:
        mkdir -p ${CASE}_${FORCING}_log
        mv -f *.out *.err ocean.output output.namelist.dyn ${CASE}_${FORCING}_log/

done




# Now we can compare the results with sanity range!

echo ; echo ; echo
echo " *** Now time for sanity-check of heat flux and wind stress components ***"
echo

cd ${HERE}
python3 ./analyze_output.py ${DATA_IN_DIR}/${SFORC} ${PROD_DIR}/output ${cmore}


