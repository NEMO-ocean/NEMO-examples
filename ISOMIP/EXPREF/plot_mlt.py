from netCDF4 import Dataset
import numpy as np
from numpy import ma
import argparse
import matplotlib.pyplot as plt
import matplotlib

parser = argparse.ArgumentParser()
parser.add_argument("-f" , metavar='file_name'   , help="names of input files" , type=str  , nargs="+", required=True )
parser.add_argument("-v" , metavar='var_name'    , help="variable list"        , type=str  , nargs=1  , required=True )
args = parser.parse_args()

# read mesh_mask
ncid   = Dataset('mesh_mask.nc')
lat2d  = ncid.variables['gphit'    ][  :,:].squeeze()
lon2d  = ncid.variables['glamt'    ][  :,:].squeeze()
msk    = ncid.variables['tmaskutil'][0,:,:].squeeze()
ncid.close()

plt.figure(figsize=np.array([210,210]) / 25.4)

# read psi.nc
ncid   = Dataset(args.f[0])
var2d  = ncid.variables[args.v[0]][-1,:,:].squeeze() 
var2dm = ma.masked_where(msk==0.0,var2d)
# convert in m/y
var2dm = var2dm * 86400 * 365 / 1e3
ncid.close()
 
# define colorbar 
vlevel=np.arange(-1.6,1.8,0.2)
pcol = plt.contourf(lon2d,lat2d,var2dm,levels=vlevel,extend='both')
vlevel=np.arange(-1.6,1.8,0.4)
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
plt.contour(lon2d,lat2d,var2dm,levels=vlevel,colors='k')
plt.grid()
plt.title('melt rate ISOMIP (m/y)')
plt.ylabel('Latitude',fontsize=14)
plt.xlabel('Longitude',fontsize=14)
cbar = plt.colorbar(pcol, ticks=vlevel)
cbar.ax.tick_params(labelsize=14)

# save figure
plt.savefig('mlt.png', format='png', dpi=300)

plt.show()

