#!/usr/bin/env python
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

# Post-diagnostic of STATION_ASF /  L. Brodeau, 2019

import sys
from os import path as path
import math
import numpy as nmp
from netCDF4 import Dataset
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


reload(sys)
sys.setdefaultencoding('utf8')

cstation = 'PAPA'

cy1     = '2018' ; # First year
cy2     = '2018' ; # Last year

jt0 = 0

dir_figs='.'
size_fig=(13,8)
size_fig0=(12,10)
fig_ext='svg'

clr_red = '#AD0000'
clr_sat = '#ffed00'
clr_mod = '#008ab8'

rDPI=100.

L_ALGOS = [ 'COARE3p6' , 'ECMWF'   , 'NCAR' ]
l_color = [  '#ffed00' , '#008ab8' , '0.4'  ] ; # colors to differentiate algos on the plot
l_width = [     3      ,    2      ,  1     ] ; # line-width to differentiate algos on the plot
l_style = [    '-'     ,   '-'     , '--'   ] ; # line-style


# Variables to compare for A GIVEN algorithm
###############################################
#L_VNEM0 = [ 



# Variables to compare between algorithms
############################################
L_VNEM  = [   'Cd_oce'  ,   'Ce_oce'  ,   'qla'     ,     'qsb'     ,     'qt'     ,   'qlw'     ,  'taum'     ,    'dt_skin'         ]
L_VARO  = [     'Cd'    ,     'Ce'    ,   'Qlat'    ,    'Qsen'     ,     'Qnet'   ,   'Qlw'     ,  'Tau'      ,    'dT_skin'         ] ; # name of variable on figure
L_VARL  = [ r'$C_{D}$'  , r'$C_{E}$'  , r'$Q_{lat}$', r'$Q_{sens}$' , r'$Q_{net}$' , r'$Q_{lw}$' , r'$|\tau|$' , r'$\Delta T_{skin}$' ] ; # name of variable in latex mode
L_VUNT  = [     ''      ,     ''      , r'$W/m^2$'  , r'$W/m^2$'    , r'$W/m^2$'   , r'$W/m^2$'  , r'$N/m^2$'  ,      'K'             ]
L_VMAX  = [    0.0075   ,    0.005    ,     75.     ,     75.       ,    800.      ,     25.     ,    1.2      ,       0.7            ]
L_VMIN  = [    0.0005   ,    0.0005   ,   -250.     ,   -125.       ,   -400.      ,   -150.     ,    0.       ,      -0.7            ]
L_ANOM  = [   False     ,   False     ,   True      ,    True       ,    True      ,    True     ,   True      ,      False           ]


nb_algos = len(L_ALGOS) ; print(nb_algos)

# Getting arguments:
narg = len(sys.argv)
if narg != 2:
    print 'Usage: '+sys.argv[0]+' <DIR_OUT_SASF>'; sys.exit(0)
cdir_data = sys.argv[1]



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Populating and checking existence of files to be read
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def chck4f(cf):
    cmesg = 'ERROR: File '+cf+' does not exist !!!'
    if not path.exists(cf): print cmesg ; sys.exit(0)

###cf_in = nmp.empty((), dtype="S10")
cf_in = []
for ja in range(nb_algos):
    cfi = cdir_data+'/output/'+'STATION_ASF-'+L_ALGOS[ja]+'_1h_'+cy1+'0101_'+cy2+'1231_gridT.nc'
    chck4f(cfi)
    cf_in.append(cfi)
print('Files we are goin to use:')
for ja in range(nb_algos): print(cf_in[ja])
#-----------------------------------------------------------------


# Getting time array from the first file:
id_in = Dataset(cf_in[0])
vt = id_in.variables['time_counter'][jt0:]
cunit_t = id_in.variables['time_counter'].units ; print(' "time_counter" is in "'+cunit_t+'"')
id_in.close()
nbr = len(vt)




vtime = nmp.zeros(nbr)

vt = vt + 1036800. + 30.*60. # BUG!??? don't get why false in epoch to date conversion, and yet ncview gets it right!
for jt in range(nbr): vtime[jt] = mdates.epoch2num(vt[jt])

ii=nbr/300
ib=max(ii-ii%10,1)
xticks_d=int(30*ib)


rat = 100./float(rDPI)
params = { 'font.family':'Open Sans',
           'font.size':       int(15.*rat),
           'legend.fontsize': int(15.*rat),
           'xtick.labelsize': int(15.*rat),
           'ytick.labelsize': int(15.*rat),
           'axes.labelsize':  int(16.*rat)
}
mpl.rcParams.update(params)
font_inf = { 'fontname':'Open Sans', 'fontweight':'normal', 'fontsize':18.*rat }
font_x   = { 'fontname':'Open Sans', 'fontweight':'normal', 'fontsize':15.*rat }





# First for each algorithm we compare some input vs out put variables:                                                                       

# t_skin
vtemp_in = [ 'sst'            ,               'theta_zu'                 ,               't_skin'                   ] ; #,               'theta_zt'                
vtemp_lb = [ 'SST (bulk SST)' , r'$\theta_{zu}$ (pot. air temp. at 10m)' ,           'Skin temperature'             ] ; #, r'$\theta_{zt}$ (pot. air temp. at 2m)' 
vtemp_cl = [  'k'             ,                clr_mod                   ,                clr_red                   ] ; #,                clr_mod                  
vtemp_lw = [   3              ,                   1.3                    ,                   0.7                    ] ; #,                 2                       
vtemp_ls = [     '-'          ,                   '-'                    ,                   '-'                    ] ; #,                 '-'                     

ntemp = len(vtemp_in)

xxx = nmp.zeros((nbr,ntemp))

for ja in range(nb_algos):
    #
    # Temperatures...
    id_in = Dataset(cf_in[ja])
    for jv in range(ntemp):
        xxx[:,jv] = id_in.variables[vtemp_in[jv]][jt0:,1,1] # only the center point of the 3x3 spatial domain!
    id_in.close()

    fig = plt.figure(num = 1, figsize=size_fig0, facecolor='w', edgecolor='k')
    ax1 = plt.axes([0.07, 0.2, 0.9, 0.75])
    ax1.set_xticks(vtime[::xticks_d])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    plt.xticks(rotation='60', **font_x)

    for jv in range(ntemp):
        plt.plot(vtime, xxx[:,jv], '-', color=vtemp_cl[jv], linestyle=vtemp_ls[jv], linewidth=vtemp_lw[jv], label=vtemp_lb[jv], zorder=10)

    ax1.set_ylim(0., 17.) ; ax1.set_xlim(vtime[0],vtime[nbr-1])
    plt.ylabel(r'Temperature [$^{\circ}$C]')

    ax1.grid(color='k', linestyle='-', linewidth=0.3)
    plt.legend(loc='best', ncol=1, shadow=True, fancybox=True)
    ax1.annotate('Algo: '+L_ALGOS[ja]+', station: '+cstation, xy=(0.4, 1.), xycoords='axes fraction',  bbox={'facecolor':'w', 'alpha':1., 'pad':10}, zorder=50, **font_inf)
    plt.savefig('01_temperatures_'+L_ALGOS[ja]+'.'+fig_ext, dpi=int(rDPI), transparent=False)
    plt.close(1)


del xxx



# Now we compare output variables from bulk algorithms between them:

nb_var = len(L_VNEM)

xF  = nmp.zeros((nbr,nb_algos))
xFa = nmp.zeros((nbr,nb_algos))


for jv in range(nb_var):
    print('\n *** Treating variable: '+L_VARO[jv]+' !')

    for ja in range(nb_algos):
        #
        id_in = Dataset(cf_in[ja])
        xF[:,ja] = id_in.variables[L_VNEM[jv]][jt0:,1,1] # only the center point of the 3x3 spatial domain!
        if ja == 0: cvar_lnm = id_in.variables[L_VNEM[jv]].long_name
        id_in.close()

    fig = plt.figure(num = jv, figsize=size_fig, facecolor='w', edgecolor='k')
    ax1 = plt.axes([0.08, 0.25, 0.9, 0.7])
    ax1.set_xticks(vtime[::xticks_d])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    plt.xticks(rotation='60', **font_x)

    for ja in range(nb_algos):
        plt.plot(vtime, xF[:,ja], '-', color=l_color[ja], linestyle=l_style[ja], linewidth=l_width[ja], label=L_ALGOS[ja], zorder=10+ja)

    ax1.set_ylim(L_VMIN[jv], L_VMAX[jv]) ; ax1.set_xlim(vtime[0],vtime[nbr-1])
    plt.ylabel(L_VARL[jv]+' ['+L_VUNT[jv]+']')

    ax1.grid(color='k', linestyle='-', linewidth=0.3)
    plt.legend(loc='best', ncol=1, shadow=True, fancybox=True)
    ax1.annotate(cvar_lnm+', station: '+cstation, xy=(0.3, 1.), xycoords='axes fraction',  bbox={'facecolor':'w', 'alpha':1., 'pad':10}, zorder=50, **font_inf)
    plt.savefig(L_VARO[jv]+'.'+fig_ext, dpi=int(rDPI), transparent=False)
    plt.close(jv)

    if L_ANOM[jv]:

        for ja in range(nb_algos): xFa[:,ja] = xF[:,ja] - nmp.mean(xF,axis=1)

        if nmp.sum(xFa[:,:]) == 0.0:
            print('     Well! Seems that for variable '+L_VARO[jv]+', choice of algo has no impact a all!')
            print('          ==> skipping anomaly plot...')

        else:

            # Want a symetric y-range that makes sense for the anomaly we're looking at:
            rmax = nmp.max(xFa) ; rmin = nmp.min(xFa)
            rmax = max( abs(rmax) , abs(rmin) )
            romagn = math.floor(math.log(rmax, 10)) ; # order of magnitude of the anomaly  we're dealing with
            rmlt = 10.**(int(romagn)) / 2.
            yrng = math.copysign( math.ceil(abs(rmax)/rmlt)*rmlt , rmax)
            #print 'yrng = ', yrng ;  #sys.exit(0)

            fig = plt.figure(num = 10+jv, figsize=size_fig, facecolor='w', edgecolor='k')
            ax1 = plt.axes([0.08, 0.25, 0.9, 0.7])
            ax1.set_xticks(vtime[::xticks_d])
            ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
            plt.xticks(rotation='60', **font_x)

            for ja in range(nb_algos):
                plt.plot(vtime, xFa[:,ja], '-', color=l_color[ja], linewidth=l_width[ja], label=L_ALGOS[ja], zorder=10+ja)

            ax1.set_ylim(-yrng,yrng) ; ax1.set_xlim(vtime[0],vtime[nbr-1])
            plt.ylabel(L_VARL[jv]+' ['+L_VUNT[jv]+']')
            ax1.grid(color='k', linestyle='-', linewidth=0.3)
            plt.legend(bbox_to_anchor=(0.45, 0.2), ncol=1, shadow=True, fancybox=True)
            ax1.annotate('Anomaly of '+cvar_lnm, xy=(0.3, 0.97), xycoords='axes fraction',  bbox={'facecolor':'w', 'alpha':1., 'pad':10}, zorder=50, **font_inf)
            plt.savefig(L_VARO[jv]+'_anomaly.'+fig_ext, dpi=int(rDPI), transparent=False)
            plt.close(10+jv)




