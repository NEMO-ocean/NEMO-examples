from netCDF4 import Dataset
import numpy as np
import argparse
import matplotlib.pyplot as plt

# read argument
parser = argparse.ArgumentParser()
parser.add_argument("-f" , metavar='file_name'   , help="names of input files" , type=str  , nargs=1  , required=True )
parser.add_argument("-v" , metavar='var_name'    , help="variable list"        , type=str  , nargs=1  , required=True )
args = parser.parse_args()

# read mesh mask
ncid  = Dataset('mesh_mask.nc')
vx2d  = ncid.variables['gphit'  ][0,:,0].squeeze()
vx2dv = ncid.variables['gphiv'  ][0,:,0].squeeze()
y2d   = ncid.variables['gdepw_0'][0,:,:,1].squeeze()*-1
y2dt  = ncid.variables['gdept_0'][0,:,:,1].squeeze()*-1
msk   = ncid.variables['tmask'  ][0,:,:,1].squeeze()
ncid.close()

# build x 2d array
x2d=y2d*0.0
x2dv=y2d*0.0
for jk in range(0,y2d.shape[0]):
   x2d[jk,:]=vx2d[:]
   x2dv[jk,:]=vx2d[:]

plt.figure(figsize=np.array([210,210]) / 25.4)

# read data and mask it
ncid   = Dataset(args.f[0])
var2d  = ncid.variables[args.v[0]][-1,:,:,:].squeeze() 
var2dm = var2d[:,:]
var2dm[msk==0] = -1
ncid.close()

# define colorbar
vlevel=np.arange(0,0.13,0.01)  
pcol = plt.contourf(x2d,y2d,var2dm,levels=vlevel)
plt.clf()

# plot contour
ax = plt.subplot(1, 1, 1)
ax.contour(x2dv,y2dt,var2dm,levels=vlevel)
ax.grid()
ax.set_title('MOC ISOMIP (Sv)')
ax.set_ylabel('Depth (m)',fontsize=14)
ax.set_xlabel('Latitude',fontsize=14)

# plot colorbar
plt.subplots_adjust(left=0.1,right=0.89, bottom=0.1, top=0.89, wspace=0.1, hspace=0.1)
cax = plt.axes([0.91, 0.1, 0.02, 0.79])
cbar= plt.colorbar(pcol, ticks=vlevel, cax=cax)
cbar.ax.tick_params(labelsize=14)

# save figure
plt.savefig('moc.png', format='png', dpi=300)

plt.show()

