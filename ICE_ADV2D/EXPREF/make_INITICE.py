#!/usr/bin/python

import os,sys
from netCDF4 import Dataset as netcdf
import numpy as np
import matplotlib.pyplot as plt
from math import exp
from math import ceil

resname=''

# input file
fcoord='mesh_mask.nc'

# output file
fflx='initice.nc'

print '   creating init ice file  ' +fflx

# Reading coordinates file
nccoord=netcdf(fcoord,'r')
nav_lon=nccoord.variables['nav_lon']
nav_lat=nccoord.variables['nav_lat']
time_counter=1
LON1= nav_lon.shape[1]
LAT1= nav_lon.shape[0]
print 'nav_lon.shape[1]' ,nav_lon.shape[1]
print 'LON1 ', LON1
print 'LAT1 ', LAT1

# Creating INITICE netcdf file
nc=netcdf(fflx,'w')
nc.createDimension('y',LAT1)
nc.createDimension('x',LON1)
nc.createDimension('time_counter',None)    # Setting dimension size to 0 or None makes it unlimited. 

cdflon=nc.createVariable('nav_lon','f',('y','x'))
cdflat=nc.createVariable('nav_lat','f',('y','x'))
cdftimecounter=nc.createVariable('time_counter','f',('time_counter'))

# ati : Fraction of open waters in sea ice - units %
# hti : Sea ice thickness - units m
# hts : Snow thickness - units m
# smi : Sea ice salinity:
# tmi : Sea ice internal temperature - units K
# tsu : Sea ice surface temperature - units K
#
# Take constant values from namelist &namiceini of NEMO
rn_hti_ini=2.0
rn_hts_ini=0.2            #  initial real snow thickness (m)
rn_ati_ini=0.9            #  initial ice concentration   (-)
rn_smi_ini=6.3            #  initial ice salinity     (g/kg)
rn_tmi_ini=270.           #  initial ice/snw temperature (K)
rn_tsu_ini=270.           #  initial sea ice temperature (K)
#
cdfati=nc.createVariable('ati','f',('time_counter','y','x'))
cdfati.units='Percentage'
cdfati.long_name='Sea ice concentration'
cdfhti=nc.createVariable('hti','f',('time_counter','y','x'))
cdfhti.long_name='Sea ice thickness'
cdfhti.units='m'
cdfhts=nc.createVariable('hts','f',('time_counter','y','x'))
cdfhts.long_name='Snow thickness'
cdfhts.units='m'
cdfsmi=nc.createVariable('smi','f',('time_counter','y','x'))
cdfsmi.long_name='Sea ice salinity'
cdfsmi.units='pss'
cdftmi=nc.createVariable('tmi','f',('time_counter','y','x'))
cdftmi.long_name='Sea ice internal temperature'
cdftmi.units='Kelvin'
cdftsu=nc.createVariable('tsu','f',('time_counter','y','x'))
cdftsu.long_name='Sea ice surface temperature'
cdftsu.units='Kelvin'

cdflon[:,:]=nav_lon[:,:]
cdflat[:,:]=nav_lat[:,:]
cdftimecounter[0]=1

# Fill fields
#print 'cdfati[:,1]', cdfati[:,1]  -> 32 values

# Add a gaussian for sea ice thickness here
cdfhti[:,:,:]=0.
cdfhts[:,:,:]=0.
cdfati[:,:,:]=0.
cdfsmi[:,:,:]=0.
cdftmi[:,:,:]=rn_tmi_ini
cdftsu[:,:,:]=rn_tsu_ini

# --------------------------------------
# for basin=99x99km with dx=1km ; dy=3km
#sigx=-0.04
#sigy=-0.04*9.
#xshift=50.-1.
#yshift=17.-1.
#dlat=7
#dlon=21

# --------------------------------------
# for basin=99x99km with dx=3km ; dy=3km
sigx=-0.04
sigy=-0.04
xshift=50.-1.
yshift=50.-1.
dlat=21
dlon=21

# --- gaussian and square experiment ---
for y in np.arange(dlat,LAT1-dlat,1) :
    for x in np.arange(dlon,LON1-dlon,1) :
        cdfhti[:,y,x] = rn_hti_ini*exp(sigx*(x-xshift)**2)*exp(sigy*(y-yshift)**2)
        cdfhts[:,y,x] = rn_hts_ini*exp(sigx*(x-xshift)**2)*exp(sigy*(y-yshift)**2)
        cdfati[:,y,x] = rn_ati_ini
        cdfsmi[:,y,x] = rn_smi_ini

nc.close()
nccoord.close()

#sys.exit()
