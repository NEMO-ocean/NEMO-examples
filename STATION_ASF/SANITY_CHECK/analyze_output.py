#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
#######################################################################################
# This script analyzes the output of STATION_ASF test-case with IDEALIZED forcing
# in order to test the validity of computed fluxes and bulk transfer coefficient.
# Beside an explicit standard output message, a result file is spawned:
#    * SBCBLK.success => the test passed !
#    * SBCBLK.fail    => the test failed !
#
#   Brodeau, 2020
#
########################################################################################

import sys
from os import path as path
import math
import numpy as nmp
from netCDF4 import Dataset

l_t_shift = False  ; # because time interp. is set to FALSE into "&namsbc_blk" of NEMO...
#                    #  ==> so time array is shifted by 30 minutes but fluxes are the same (persitence of input fields ???)

l_alg = [ 'ECMWF' , 'NCAR' , 'COARE3p0', 'COARE3p6', 'ANDREAS' ]
nb_alg     =    len(l_alg)


# Variables to check:
l_var_rf = [  'Qsen'   ,   'Qlat'   ,  'Qlw'   ,   'Tau',   'Cd'  ,   'Ce'   ] ; # In forcing file "IDEALIZED/input_output_VALIDATION_IDEALIZED.nc"
l_var_ot = [ 'qsb_oce' , 'qla_oce'  , 'qlw_oce',  'taum', 'Cd_oce', 'Ce_oce' ] ; # names in the NEMO/STATION_ASF output file (check file_def_oce.xml)
nb_var  = len(l_var_rf)

dir_figs='.'
size_fig=(13,8)
#fig_ext='png'
fig_ext='svg'


rDPI=100.


class fclrs:
    OKGR = '\033[92m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'



# Getting arguments:
narg = len(sys.argv)
if not narg in [3,4]:
    print('Usage: '+sys.argv[0]+' <forcing_+_validation_file> <NEMO-STATION_ASF_output_directory> (<m> for more/debug)'); sys.exit(0)
cf_rf = sys.argv[1]
cdir_out = sys.argv[2]
l_more = False
if narg==4:
    l_more = ( sys.argv[3] in ['m','M'] )
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt



    
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Populating and checking existence of files to be read
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def chck4f(cf):
    cmesg = 'ERROR: File '+cf+' does not exist !!!'
    if not path.exists(cf): print(cmesg) ; sys.exit(0)

print('\n')

# Input forcing/valid file:
chck4f(cf_rf)
id_rf = Dataset(cf_rf)
vt = id_rf.variables['time'][:]
cunit_t = id_rf.variables['time'].units
id_rf.close()
Nt_rf = len(vt)
vtime_rf = nmp.zeros(Nt_rf); vtime_rf[:]  = vt[:]
del vt
print(' *** in forcing/valid file, "time" is in "'+cunit_t+'", Nt = '+str(Nt_rf)+'\n')


# STATION_ASF output files (1 file per algorithm):
cf_nemo = []
for ja in range(nb_alg):
    cfi = cdir_out+'/STATION_ASF-'+l_alg[ja]+'_IDEALIZED_1h_20200101_20200105_gridT.nc'
    chck4f(cfi)
    cf_nemo.append(cfi)
print('\n *** NEMO/STATION_ASF output files we are goin to check:')
for ja in range(nb_alg): print(cf_nemo[ja])
print('\n')
#-----------------------------------------------------------------

# Getting time array from the first file:
id_nm = Dataset(cf_nemo[0])
vt = id_nm.variables['time_counter'][:]
cunit_t = id_nm.variables['time_counter'].units ; print(' "time_counter" is in "'+cunit_t+'"')
id_nm.close()
Nt = len(vt)
vtime_nm = nmp.zeros(Nt); vtime_nm[:]  = vt[:]
del vt

if Nt != Nt_rf-1: print('ERROR: the two files do not agree in terms of record lengrth: '+Nt_rf-1+' vs '+Nt)

print('\n *** Excellent! We are going to look at surface fluxes under '+str(Nt)+' different scenarios of air-sea stability/wind-speed conditions...\n')

# 30 minute shift, just like NEMO
vtime = nmp.zeros(Nt)
if l_t_shift:
    vtime[:] = 0.5*(vtime_rf[:-1] + vtime_rf[1:])
else:
     vtime[:] =     vtime_rf[:-1]
# Debug:
#for jt in range(3):
#    print(' Ref.  , Nemo "', vtime_rf[jt], vtime_nm[jt], vtime[jt])
#sys.exit(0)




##
IREPORT = nmp.zeros((nb_alg,nb_var), dtype=int)




# Loop on the fields to control...
###################################

jv=0
for cv in l_var_rf:
    cv_rf_m = cv+'_mean'
    cv_rf_t = cv+'_tol'
    cv_nemo = l_var_ot[jv]
    print('\n\n ==== Checking variable '+cv_nemo+' against '+cv_rf_m+' !')

    F_rf_m = nmp.zeros( Nt )
    F_rf_t = nmp.zeros( Nt )
    F_nemo = nmp.zeros((Nt,nb_alg))

    wnd_rf = nmp.zeros( Nt ) ; #DEBUG
    tzt_rf = nmp.zeros( Nt ) ; #DEBUG
    qzt_rf = nmp.zeros( Nt ) ; #DEBUG
    sst_rf = nmp.zeros( Nt ) ; #DEBUG
    
    id_rf = Dataset(cf_rf)
    trf_m   = id_rf.variables[cv_rf_m][:] ; # Nt+1
    trf_t   = id_rf.variables[cv_rf_t][:] ; # Nt+1
    trf_wnd = id_rf.variables['wndspd'][:,1,1] ; # DEBUG
    trf_tzt = id_rf.variables['t_air'][:,1,1] ; # DEBUG
    trf_qzt = id_rf.variables['rh_air'][:,1,1] ; # DEBUG
    trf_sst = id_rf.variables['sst'][:,1,1] ; # DEBUG
    id_rf.close()
    
    if l_t_shift:
        # 30 minute shift, just like NEMO
        F_rf_m[:] = 0.5 * (trf_m[:-1] + trf_m[1:])
        F_rf_t[:] = 0.5 * (trf_t[:-1] + trf_t[1:])
    else:
        F_rf_m[:] =        trf_m[:-1]
        F_rf_t[:] =        trf_t[:-1]

        wnd_rf[:] =        trf_wnd[:-1]        
        tzt_rf[:] =        trf_tzt[:-1]
        qzt_rf[:] =        trf_qzt[:-1]
        sst_rf[:] =        trf_sst[:-1]        
    
    
    for ja in range(nb_alg):
        
        calgo = l_alg[ja]

        print(' *** '+calgo+' => '+cf_nemo[ja])

        id_nemo = Dataset(cf_nemo[ja])
        F_nemo[:,ja] = id_nemo.variables[cv_nemo][:,1,1] ; # it's 3x3 spatial domain, taking middle point !
        id_nemo.close()


        if l_more:
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # How does all this look on on a figure?
            cfig = l_var_rf[jv]+'_'+calgo+'.'+fig_ext
            print(' *** will plot '+cfig)
            fig = plt.figure(num=1, figsize=size_fig, facecolor='w', edgecolor='k')
            ax1 = plt.axes([0.08, 0.25, 0.9, 0.7])
            #
            plt.plot(vtime, F_nemo[:,ja], label='NEMO['+calgo+']', zorder=1)
            plt.plot(vtime, F_rf_m[:], color='k', label='MEAN REF!', zorder=10)
            # +- rtol enveloppe:
            plt.fill_between(vtime, F_rf_m[:]-F_rf_t[:], F_rf_m[:]+F_rf_t[:], alpha=0.2)
            #
            ax1.grid(color='k', linestyle='-', linewidth=0.3)
            plt.legend(loc='best', ncol=1, shadow=True, fancybox=True)
            plt.savefig(cfig, dpi=int(rDPI), transparent=False)
            plt.close(1)
            print('')
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        # Does the field look okay with respect to reference +- tolerance ?

        l_overshoot_a = nmp.any( F_nemo[:,ja] > F_rf_m[:]+F_rf_t[:] )
        l_overshoot_b = nmp.any( F_nemo[:,ja] < F_rf_m[:]-F_rf_t[:] )

        if l_overshoot_a: print(fclrs.FAIL+'\n ***** BAD overshoot + for '+calgo+' for variable '+cv+' !!!\n'+fclrs.ENDC )
            
        if l_overshoot_b: print(fclrs.FAIL+'\n ***** BAD overshoot - for '+calgo+' for variable '+cv+' !!!'+fclrs.ENDC )

        if l_overshoot_a or l_overshoot_b:
            print(fclrs.FAIL+'\n ***** TEST NOT PASSED FOR '+calgo+' for variable '+cv+' !!!\n'+fclrs.ENDC )
        else:
            # We're all good !
            IREPORT[ja,jv] = 1
            if l_more: print(fclrs.OKGR+'\n ***** TEST PASSED FOR '+calgo+' for variable '+cv+' :D !!!\n'+fclrs.ENDC )
                    
    jv=jv+1



l_ok =  nmp.sum(IREPORT[:,:]) == nb_var*nb_alg



if l_ok:
    ctxt      = 'PASSED'
    ccol      = fclrs.OKGR
    cf_report = 'SBCBLK.success'
else:
    ctxt      = 'FAILED'
    ccol      = fclrs.FAIL
    cf_report = 'SBCBLK.fail'
    

print(ccol+'\n\n ############ FINAL REPORT ############\n'+fclrs.ENDC)

f = open(cf_report, 'w')

f.write("### Sanity-check report for SBCBLK generated via 'STATION_ASF/EXP00/sbcblk_sanity_check.sh'\n\n")

for ja in range(nb_alg):
    calgo = l_alg[ja]
    
    if nmp.sum(IREPORT[ja,:]) == nb_var:
        # Success for this algo
        cbla = ' ***** Algorithm "'+calgo+'" PASSED sanity check !!!\n\n'
        print(fclrs.OKGR+cbla+fclrs.ENDC ) ; f.write(cbla)
    else:
        # Algo FAILS!
        cbla = ' ***** Algorithm "'+calgo+'" FAILED sanity check !!!\n'
        print(fclrs.FAIL+cbla+fclrs.ENDC ) ; f.write(cbla)
        (idx_fail,) = nmp.where(IREPORT[ja,:]==0)
        for jv in idx_fail:
            cbla = '   ==> on variable '+l_var_rf[jv]+' !\n\n'
            print(fclrs.FAIL+cbla+fclrs.ENDC ) ; f.write(cbla)

# Conclusion:
clist=''
for cc in l_var_ot: clist=clist+cc+', '
cbla = ' Test performed on the following NEMO prognostic variables:\n   ==> '+clist[:-2]+'\n'
print(ccol+cbla+fclrs.ENDC ) ; f.write(cbla+'\n')

cbla = '    ####################################\n'  +\
       '    ###   TEST '+ctxt+' FOR SBCBLK !   ###\n'+\
       '    ####################################\n\n'

print(ccol+cbla+fclrs.ENDC ) ; f.write(cbla)

f.close()
