#!/usr/bin/python

import os,sys
from netCDF4 import Dataset as netcdf
import numpy as np

resname=''

# input file
fcoord='mesh_mask.nc'

# output file
fflx='initice.nc'

print '   creating init ice file  ' +fflx

# Reading coordinates file
nccoord=netcdf(fcoord,'r')
nav_lon=nccoord.variables['glamt'][0,:,:]
nav_lat=nccoord.variables['gphit'][0,:,:]
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
rn_hts_ini=0.1            #  initial real snow thickness (m)
rn_ati_ini=0.99           #  initial ice concentration   (-)
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

cdfati[:,:,:]=rn_ati_ini
cdfhts[:,:,:]=rn_hts_ini
cdfsmi[:,:,:]=rn_smi_ini
cdftmi[:,:,:]=rn_tmi_ini
cdftsu[:,:,:]=rn_tsu_ini

# --- add noise in initial thickness to help nucleation of deformation ---
#for y in np.arange(1,1000,1) :
#    for x in np.arange(1,1000,1) :
#        cdfhti[:,y,x] = rn_hti_ini+0.02*rn_hti_ini*np.random.uniform(-1,1)
cdfhti[:,:,:]=rn_hti_ini

nc.close()
nccoord.close()

#sys.exit()
