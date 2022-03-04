#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##########################################################################################
# Post-diagnostic of STATION_ASF for air-sea fluxes (over open ocean)
#
#  L. Brodeau, 2020
##########################################################################################

import sys
from os import path, listdir
import argparse as ap
from math import floor, ceil, copysign, log
import numpy as nmp
from netCDF4 import Dataset,num2date
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

dir_figs='.'
size_fig=(13,8.5)
size_fig0=(12,10)

clr_red = '#AD0000'
clr_mod = '#008ab8'

rDPI=100.

#l_color = [ '0.85'    ,  '#ffed00' , '#008ab8' , '0.4'  ] ; # colors to differentiate algos on the plot
#l_width = [  4        ,     3      ,    2      ,  1     ] ; # line-width to differentiate algos on the plot
#l_style = [ '-'       ,    '-'     ,   '-'     , '--'   ] ; # line-style

#ffed00: yellow ON
#E8A727: ornage

l_color = [ '0.3' , '#E8A727', '0.1'  , '#008ab8'  ] ; # colors to differentiate algos on the plot
l_width = [   2   ,    2     ,  1.5   ,    2       ] ; # line-width to differentiate algos on the plot
l_style = [  '-'  ,    '-'   , '--'   ,   '-'      ] ; # line-style


# Variables to compare between algorithms
############################################
crealm = 'open-ocean'
L_VNEM = [   'Cd_oce'  ,   'Ce_oce'  ,   'qla_oce' , 'qsb_oce'     ,     'qt_oce' ,   'qlw_oce' ,  'taum'     ,    'dt_skin'         ]
L_VARO = [     'Cd'    ,     'Ce'    ,   'Qlat'    ,    'Qsen'     ,     'Qnet'   ,   'Qlw'     ,  'Tau'      ,    'dT_skin'         ] ; # name of variable on figure
L_VARL = [ r'$C_{D}$'  , r'$C_{E}$'  , r'$Q_{lat}$', r'$Q_{sens}$' , r'$Q_{net}$' , r'$Q_{lw}$' , r'$|\tau|$' , r'$\Delta T_{skin}$' ] ; # name of variable in latex mode
L_VUNT = [     ''      ,     ''      , r'$W/m^2$'  , r'$W/m^2$'    , r'$W/m^2$'   , r'$W/m^2$'  , r'$N/m^2$'  ,      'K'             ]
L_BASE = [    0.0005   ,    0.0005   ,      5.     ,      5.       ,       5      ,      5.     ,    0.05      ,       0.05           ]
L_PREC = [      3      ,      3      ,      0      ,      0        ,       0      ,      0      ,     2       ,        3             ]
L_ANOM = [   False     ,   False     ,   True      ,    True       ,    True      ,    True     ,   True      ,      False           ]
L_MAXT = [  10000.   ,     10000.,  10000.   ,  10000.      ,  10000.    ,  10000.      , 10000.     ,      1.5    ]
L_MINT = [    0.001  ,     0.001 , -10000.   , -10000.      , -10000.    , -10000.      ,-10000.     ,   -10000.   ]

# About STATION_ASF output files to read:
cpref = 'STATION_ASF-'          ; np = len(cpref)
csuff = '_gridT.nc'             ; ns = len(csuff)
cclnd = '_1h_YYYY0101_YYYY1231' ; nc = len(cclnd)


################## ARGUMENT PARSING / USAGE ################################################################################################
parser = ap.ArgumentParser(description='Generate pixel maps of a given scalar.')
#
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-d', '--dirout' , required=True,                 help='Path to (production) directory where STATION_ASF was run')
requiredNamed.add_argument('-f', '--forcing', required=True, default="PAPA", help='Name of forcing (ex: PAPA, ERA5_arctic')
#
parser.add_argument('-C', '--conf',   default="STATION_ASF",  help='specify NEMO config (ex: STATION_ASF)')
parser.add_argument('-s', '--ystart', default="2018",         help='specify first year of experiment (ex: 2018)')
parser.add_argument('-e', '--yend',   default="2018",         help='specify last  year of experiment (ex: 2018)')
#
parser.add_argument('-t', '--itype',   default="png",         help='specify the type of image you want to create (ex: png, svg, etc.)')
#parser.add_argument('-l', '--lev' , type=int, default=0,    help='specify the level to use if 3D field (default: 0 => 2D)')
#parser.add_argument('-I', '--ice' , action='store_true',    help='draw sea-ice concentration layer onto the field')
#
args = parser.parse_args()
#
cdir_data = args.dirout
cforcing  = args.forcing
#
CONF      = args.conf
cy1       = args.ystart
cy2       = args.yend
#
fig_ext   = args.itype
#jk    = args.lev
#lshow_ice = args.ice
#
#print(''); print(' *** cdir_data = ', cdir_data); print(' *** cforcing  = ', cforcing)
#print(' *** CONF = ', CONF); print(' *** cy1 = ', cy1); print(' *** cy2 = ', cy2)
###############################################################################################################################################



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Populating and checking existence of files to be read
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

dir_out = cdir_data+'/output'
ldir = listdir(dir_out)

cf_in    = []
list_exp = []
list_frc = []
for fn in ldir:
    fpn = dir_out+'/'+fn
    if path.isfile(dir_out+'/'+fn):
        if fn[:np]==cpref and fn[-ns:]==csuff and cforcing in fn:
            print('\n file: '+fn)
            clab = fn[np:-nc-ns]
            [ cexp, cfrc ] = str.split(clab, '_', 1)
            print('  ===> Experiment = '+cexp+', Forcing = '+cfrc)
            list_exp.append(cexp)
            list_frc.append(cfrc)
            cf_in.append(fpn)
nbf = len( set(list_frc) )
if not nbf == 1:
    print('PROBLEM: we found files for more that one forcing: ', set(list_frc))
    sys.exit(0)

nb_exp = len(list_exp)


print('\n\nThere are '+str(nb_exp)+' experiments to compare:')
for ja in range(nb_exp): print('  * '+list_exp[ja]+'\n'+'     ==> '+cf_in[ja]+'\n')

if nb_exp > len(l_color):
    print('PROBLEM: the max number of experiments for comparison is '+str(len(l_color))+' for now...')
    sys.exit(0)


#-----------------------------------------------------------------


def round_bounds( x1, x2,  base=5, prec=3 ):
    rmin =  base * round( floor(float(x1)/base), prec )
    rmax =  base * round(  ceil(float(x2)/base), prec )
    return rmin, rmax


# Getting time array from the first file:
id_in = Dataset(cf_in[0])
vt = id_in.variables['time_counter'][:]
cunit_t = id_in.variables['time_counter'].units ; print(' "time_counter" is in "'+cunit_t+'"')
id_in.close()
Nt = len(vt)

vtime = num2date(vt, units=cunit_t) ; # something human!
vtime = vtime.astype(dtype='datetime64[D]')

ii=Nt/300
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
vtemp_in = [ 'sst'           ,               'theta_zu'               ,               'theta_zt'               ,     't_skin'       ]
vtemp_lb = [ 'SST (bulk SST)', r'$\theta_{zu}$ (pot. air temp. at zu)', r'$\theta_{zt}$ (pot. air temp. at zt)', 'Skin temperature' ]
vtemp_cl = [  'k'            ,                clr_mod                 ,                'purple'                ,      clr_red       ]
vtemp_lw = [   3             ,                   1.3                  ,                   1.3                  ,         0.7        ]
vtemp_ls = [     '-'         ,                   '-'                  ,                   '--'                 ,         '-'        ]
ntemp = len(vtemp_in)

xxx = nmp.zeros((Nt,ntemp))

for ja in range(nb_exp):
    #
    # Temperatures...
    id_in = Dataset(cf_in[ja])
    for jv in range(ntemp):
        xxx[:,jv] = id_in.variables[vtemp_in[jv]][:,1,1] # only the center point of the 3x3 spatial domain!
    id_in.close()

    fig = plt.figure(num = 1, figsize=size_fig0, facecolor='w', edgecolor='k')
    ax1 = plt.axes([0.07, 0.2, 0.9, 0.75])
    ax1.set_xticks(vtime[::xticks_d])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    plt.xticks(rotation='60', **font_x)

    for jv in range(ntemp):
        plt.plot(vtime, xxx[:,jv], '-', color=vtemp_cl[jv], linestyle=vtemp_ls[jv], linewidth=vtemp_lw[jv], label=vtemp_lb[jv], zorder=10)

    idx_okay = nmp.where( nmp.abs(xxx) < 1.e+10 )
    fmin, fmax = round_bounds( nmp.min(xxx[idx_okay]) , nmp.max(xxx[idx_okay]), base=5, prec=0 )
    ax1.set_ylim(fmin, fmax) ; ax1.set_xlim(vtime[0],vtime[Nt-1])
    plt.ylabel(r'Temperature [$^{\circ}$C]')

    ax1.grid(color='k', linestyle='-', linewidth=0.3)
    plt.legend(loc='best', ncol=1, shadow=True, fancybox=True)
    ax1.annotate('Algo: '+list_exp[ja]+', station: '+cforcing, xy=(0.5, 1.), xycoords='axes fraction', ha='center',  bbox={'facecolor':'w', 'alpha':1., 'pad':10}, zorder=50, **font_inf)
    plt.savefig('01_temperatures_'+list_exp[ja]+'_'+cforcing+'.'+fig_ext, dpi=int(rDPI), transparent=False)
    plt.close(1)


del xxx



# Now we compare output variables from bulk algorithms between them:

nb_var = len(L_VNEM)

xF  = nmp.zeros((Nt,nb_exp))
xFa = nmp.zeros((Nt,nb_exp))


for jv in range(nb_var):
    print('\n *** Treating variable: '+L_VARO[jv]+' !')

    for ja in range(nb_exp):
        #
        id_in = Dataset(cf_in[ja])
        xF[:,ja] = id_in.variables[L_VNEM[jv]][:,1,1]  ; # only the center point of the 3x3 spatial domain!
        if ja == 0: cvar_lnm = id_in.variables[L_VNEM[jv]].long_name
        id_in.close()
        #
        id_toolarge, = nmp.where( xF[:,ja] > L_MAXT[jv] ) # 
        xF[id_toolarge,ja] = L_MAXT[jv]
        id_toosmall, = nmp.where( xF[:,ja] < L_MINT[jv] ) ; #print("id_toosmall =", id_toosmall)
        xF[id_toosmall,ja] = L_MINT[jv]

    idx_okay = nmp.where( nmp.abs(xF) < 1.e+10 )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fig = plt.figure(num = jv, figsize=size_fig, facecolor='w', edgecolor='k')
    ax1 = plt.axes([0.083, 0.23, 0.9, 0.7])
    ax1.set_xticks(vtime[::xticks_d])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    plt.xticks(rotation='60', **font_x)

    for ja in range(nb_exp):
        fplot = nmp.ma.masked_where( xF[:,ja]==0., xF[:,ja] )
        plt.plot(vtime, fplot, '-', color=l_color[ja], \
                 linestyle=l_style[ja], linewidth=l_width[ja], label=list_exp[ja], alpha=0.6 ) #zorder=10+ja)

    fmin, fmax = round_bounds( nmp.min(xF[idx_okay]) , nmp.max(xF[idx_okay]), base=L_BASE[jv], prec=L_PREC[jv])
    ax1.set_ylim(fmin, fmax) ; ax1.set_xlim(vtime[0],vtime[Nt-1])
    plt.ylabel(L_VARL[jv]+' ['+L_VUNT[jv]+']')

    ax1.grid(color='k', linestyle='-', linewidth=0.3)
    plt.legend(loc='best', ncol=1, shadow=True, fancybox=True)
    ax1.annotate(cvar_lnm+', station: '+cforcing, xy=(0.5, 1.04), xycoords='axes fraction', \
                 ha='center', bbox={'facecolor':'w', 'alpha':1., 'pad':10}, \
                 zorder=50, **font_inf)
    plt.savefig(L_VARO[jv]+'_'+cforcing+'_'+crealm+'.'+fig_ext, dpi=int(rDPI), transparent=False)
    plt.close(jv)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




    def symetric_range( pmin, pmax ):
        # Returns a symetric f-range that makes sense for the anomaly of "f" we're looking at...
        from math import floor, copysign, log, ceil
        zmax = max( abs(pmax) , abs(pmin) )
        romagn = floor(log(zmax, 10)) ; # order of magnitude of the anomaly  we're dealing with
        rmlt = 10.**(int(romagn)) / 2.
        frng = copysign( ceil(abs(zmax)/rmlt)*rmlt , zmax)
        return frng


    if L_ANOM[jv]:

        for ja in range(nb_exp): xFa[:,ja] = xF[:,ja] - nmp.mean(xF,axis=1)

        if nmp.sum(nmp.abs(xFa[:,:])) == 0.0:
            print('     Well! Seems that for variable '+L_VARO[jv]+', choice of algo has no impact a all!')
            print('          ==> skipping anomaly plot...')

        else:

            yrng = symetric_range( nmp.min(xFa) , nmp.max(xFa) )

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            fig = plt.figure(num = 10+jv, figsize=size_fig, facecolor='w', edgecolor='k')
            ax1 = plt.axes([0.09, 0.23, 0.9, 0.7])
            ax1.set_xticks(vtime[::xticks_d])
            ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
            plt.xticks(rotation='60', **font_x)

            for ja in range(nb_exp):
                fplot = nmp.ma.masked_where( xFa[:,ja]==0., xFa[:,ja] )
                plt.plot(vtime, fplot, '-', color=l_color[ja], \
                         linewidth=l_width[ja], label=list_exp[ja], alpha=0.6) #, zorder=10+ja)

            ax1.set_ylim(-yrng,yrng) ; ax1.set_xlim(vtime[0],vtime[Nt-1])
            plt.ylabel(L_VARL[jv]+' ['+L_VUNT[jv]+']')
            ax1.grid(color='k', linestyle='-', linewidth=0.3)
            plt.legend(bbox_to_anchor=(0.45, 0.2), ncol=1, shadow=True, fancybox=True)
            ax1.annotate('Anomaly of '+cvar_lnm+', station: '+cforcing, xy=(0.5, 1.04), xycoords='axes fraction', \
                         ha='center', bbox={'facecolor':'w', 'alpha':1., 'pad':10}, \
                         zorder=50, **font_inf)
            plt.savefig(L_VARO[jv]+'_'+cforcing+'_anomaly_'+crealm+'.'+fig_ext, dpi=int(rDPI), transparent=False)
            plt.close(10+jv)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
