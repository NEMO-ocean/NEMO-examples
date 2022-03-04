#!/bin/bash

################################################################
#
# Script to launch a set of STATION_ASF simulations
#
#  Ocean-only version....
#
# L. Brodeau, 2020
#
################################################################

CONFIG="STATION_ASF" ; # general name of the configuration
CONFIG_BLD="STATION_ASF2" ; # name of config as build in NEMO... (directory inside "tests" actually contains the compiled test-case?)

# Atmo + SSX forcing to use and sea-ice support?
#     => FORCING: name of forcing
#     => i_sea_ice: whether to compute fluxes over sea-ice as well
#     => SFORC: string sufficient to copy relevant files as in "*${SFORC}*.nc"
FORCING="PAPA"                 ; i_sea_ice=0 ; SFORC="Station_PAPA_50N-145W"
#FORCING="ERA5_NorthGreenland" ; i_sea_ice=1 ; SFORC="ERA5_NorthGreenland_surface_84N_-36E_1h" ; # WITH ice/air flux computation
#FORCING="ERA5_NorthGreenland" ; i_sea_ice=0 ; SFORC="ERA5_NorthGreenland_surface_84N_-36E_1h" ; # "ERA5_arctic" WITHOUT ice/air flux computation


# Root directory NEMOGCM reference depository where to fetch compiled STATION_ASF nemo.exe + default namelists:
NEMO_REF_DIR="`dirname ${PWD} | sed -e 's|/tests/STATION_ASF||g'`" ; # that should normally do the trick!

# NEMOGCM root directory where to fetch compiled STATION_ASF nemo.exe:
SASF_WRK_DIR="${NEMO_REF_DIR}/tests/${CONFIG_BLD}"

# DATA_IN_DIR => Directory containing sea-surface + atmospheric forcings:
DATA_IN_DIR="${NEMO_REF_DIR}/tests/${CONFIG}/input_data"

# Directory where to run the simulation:
PROD_DIR="${HOME}/tmp/${CONFIG}"

# MPI launch command for your system
# Even though we are running on 1 proc here, NEMO has likely been compiled with MPI,
# if it's not the case then just leave "MPI_LAUNCH" void...
MPI_LAUNCH="mpirun -n 1"

####### End of normal user configurable section #######
#================================================================================

iconv_1d=1
if [ "`which ncks`" = "" ]; then
    echo
    echo "WARNING: you do not seem to have NCO installed here... You should!"
    echo "      => anyway! will do without, but output netCDF files will remaine 3x3 in space :( "
    echo
    iconv_1d=0
    sleep 3
fi

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

# Ok things are a bit different whether we compute fluxes over open ocean only or open ocean + sea-ice

DIR_F=`echo ${FORCING} | cut -d'_' -f1` ; # "ERA5_blabla" => "ERA5" !

if [ ${i_sea_ice} -eq 1 ]; then
    # OPEN-OCEAN/AIR + SEA-ICE/AIR flux computation
    #   - simpler if we only use one Open-ocean/air algo (ECMWF)
    LIST_OA_ALGOS="ECMWF" ;          # list of air-sea algorithms to test
    LIST_IA_ALGOS="AN05 LG15 LU12 CSTC" ; # list of air-ice algorithms to test
    DIR_NL=${DIR_F}/oce+ice ; # where to fetch the namelists from...
else
    # Only OPEN-OCEAN/AIR flux computation
    LIST_OA_ALGOS="NCAR ECMWF COARE3p6 ANDREAS"; # list of air-sea algorithms to test
    LIST_IA_ALGOS=""
    DIR_NL=${DIR_F}/oce ; # where to fetch the namelists from...
fi
if [ ! -d ${DIR_NL} ]; then echo " Mhhh, ${DIR_NL} not found !"; exit; fi



# NEMO executable to use is:
NEMO_EXE="${SASF_WRK_DIR}/BLD/bin/nemo.exe"
if [ ! -f ${NEMO_EXE} ]; then echo " Mhhh, no compiled 'nemo.exe' found into `dirname ${NEMO_EXE}` !"; exit; fi

SASF_EXPREF=`pwd`  ; # STATION_ASF EXPREF directory from which to use namelists and XIOS xml files...

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
ln -sf      ${NEMO_EXE}          ${PROD_DIR}/nemo.exe.lnk


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
\ls -l ${DATA_IN_DIR}/*${SFORC}*.nc
echo
rsync -avPL ${DATA_IN_DIR}/*${SFORC}*.nc ${PROD_DIR}/
echo; echo


if [ "${LIST_IA_ALGOS}" = "" ]; then LIST_IA_ALGOS="none"; fi

echo

for OA_A in ${LIST_OA_ALGOS}; do
    for IA_A in ${LIST_IA_ALGOS}; do

        echo        
        CASE="${OA_A}"
        if [ ${i_sea_ice} -eq 1 ] && [ "${IA_A}" != "none" ]; then CASE="${OA_A}-${IA_A}"; fi

        echo ; echo
        echo "======================================================================="
        echo " Going for experiment: ${CASE} bulk param. with ${FORCING} forcing "
        echo "======================================================================="
        echo

        scase=`echo "${CASE}" | tr '[:upper:]' '[:lower:]'`

        rm -f ${PROD_DIR}/namelist_cfg
        fnml="${SASF_EXPREF}/${DIR_NL}/namelist_${scase}_cfg"
        if [ ! -f ${fnml} ]; then echo " Mhhh, test namelist ${fnml} not found !"; exit; fi

        echo "   ===> namelist to use is: ${DIR_NL}/namelist_${scase}_cfg !"; echo
        
        # The namelists:
        #rsync -avPL ${fnml} ${PROD_DIR}/namelist_cfg
        sed -e "s|<FORCING>|${FORCING}|g" -e "s|<SFORC>|${SFORC}|g" ${fnml} > ${PROD_DIR}/namelist_cfg

        rsync -avPL ${CFGS_SHARED}/namelist_ref                    ${PROD_DIR}/namelist_ref
        if [ ${i_sea_ice} -eq 1 ]; then
            rsync -avPL ${SASF_EXPREF}/namelist_ice_cfg ${PROD_DIR}/namelist_ice_cfg
            rsync -avPL ${CFGS_SHARED}/namelist_ice_ref ${PROD_DIR}/namelist_ice_ref
        fi
    
        cd ${PROD_DIR}/
        echo
        echo "Launching NEMO ! (==> ${MPI_LAUNCH} ./nemo.exe)"
        ${MPI_LAUNCH} ./nemo.exe 1>out_nemo.out 2>err_nemo.err
        echo "Done!"
        echo

        # Moving output files:
        mkdir -p output
        mv -f ${CONFIG}-${CASE}_${FORCING}_*_grid*.nc output/
        if [ ${i_sea_ice} -eq 1 ]; then mv -f ${CONFIG}-${CASE}_${FORCING}_*_icemod.nc output/; fi

        # Saving logs:
        mkdir -p ${CASE}_${FORCING}_log
        mv -f *.out *.err ocean.output output.namelist.dyn ${CASE}_${FORCING}_log/

        if [ ${iconv_1d} -eq 1 ]; then
        # Making 3x3 to 1 !
        cd output/
        list=`\ls ${CONFIG}*${CASE}_${FORCING}_*.nc | grep -v '_restart_'| grep -v '_1p.nc'`
        for ff in ${list}; do
            echo
            fn=`echo ${ff} | sed -e "s|.nc|_1p.nc|g"`
            CMD="ncks -O -d x,1 -d y,1 ${ff} -o ${fn}"
            echo " *** ${CMD}"; ${CMD}
            echo
        done
        fi

    done

done
