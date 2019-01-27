from netCDF4 import Dataset
import numpy as np
from numpy import ma
import argparse
import matplotlib.pyplot as plt

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
# convert in Sv
var2dm = var2dm / 1e6
ncid.close()
 
# define colorbar 
vlevel=np.arange(0.00,0.36,0.02)
pcol = plt.contourf(lon2d,lat2d,var2dm,levels=vlevel)
plt.clf()

# plot contour
ax = plt.subplot(1, 1, 1)
ax.contour(lon2d,lat2d,var2dm,levels=vlevel)
ax.grid()
ax.set_title('PSI ISOMIP (Sv)')
ax.set_ylabel('Latitude',fontsize=14)
ax.set_xlabel('Longitude',fontsize=14)

# plot colorbar
plt.subplots_adjust(left=0.1,right=0.89, bottom=0.1, top=0.89, wspace=0.1, hspace=0.1)
cax = plt.axes([0.91, 0.1, 0.02, 0.79])
cbar= plt.colorbar(pcol, ticks=vlevel, cax=cax)
cbar.ax.tick_params(labelsize=14)

# save figure
plt.savefig('psi.png', format='png', dpi=300)

plt.show()

