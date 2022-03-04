#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Draw a plot (or and animation)
    of passive tracer at the bottom
    in the DOME experiment
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#
# Parent grid data:
ncid = Dataset('DOME_grid_T.nc')
lon0 = ncid.variables['nav_lon_grid_T'][:, :]
lat0 = ncid.variables['nav_lat_grid_T'][:, :]
work = ncid.variables['btra'][:, :, :]
zos = ncid.variables['zos'][:, :, :]
ncid.close()
#
(jpt0, jpj0, jpi0) = np.shape(work)
mask = np.where((zos == 0.), True, False)
tra0 = np.ma.array(work, mask=mask, hard_mask=True)

#
# Child grid data:
ncid = Dataset('1_DOME_grid_T.nc')
lon1 = ncid.variables['nav_lon_grid_T'][:, :]
lat1 = ncid.variables['nav_lat_grid_T'][:, :]
work = ncid.variables['btra'][:, :, :]
zos = ncid.variables['zos'][:, :, :]
sp = ncid.variables['Agrif_sponge'][:, :, :]
ncid.close()
#
(jpt1, jpj1, jpi1) = np.shape(work)
mask = np.where((zos == 0.), True, False)
tra1 = np.ma.array(work, mask=mask, hard_mask=True)

# Get resolution in km:
res0 = np.abs(lon0[0, 1] - lon0[1, 0])
res1 = np.abs(lon1[0, 1] - lon1[0, 0])

# Shift lon, lat in order to have a pixel centred around the right location:
lon0 = lon0 - res0/2.
lat0 = lat0 - res0/2.
lon1 = lon1 - res1/2.
lat1 = lat1 - res1/2.

# Indexes to skip ghost zone:
nghost = 3 # + int(res0/res1)*2 - 1
imin = nghost + 1
jmin = nghost + 1
imax = jpi1 - nghost - 1
jmax = jpj1 - nghost - 1
#
fig, ax = plt.subplots()
fig.set_tight_layout(True)
mycmap = plt.cm.nipy_spectral
plt.set_cmap(mycmap)
i = 0 
label = ' Bottom tracer concentration day: {00}'.format(i+1)
pcol = plt.pcolor(lon0, lat0, np.squeeze(tra0[i, :, :]),
                  vmin=0.01, vmax=1., edgecolor='0.9', cmap=mycmap)
pcol = plt.pcolor(lon1[jmin:jmax, imin:imax], lat1[jmin:jmax, imin:imax], np.squeeze(tra1[i, jmin:jmax, imin:imax]),
                  vmin=0.01, vmax=1., edgecolor='0.4', cmap=mycmap)
# plt.contour(lon1[jmin:jmax, imin:imax], lat1[jmin:jmax, imin:imax], np.squeeze(sp[i, jmin:jmax, imin:imax]))

plt.axis('scaled')
# plt.xlim((-800., 250.))
# plt.ylim((-500.,  50.))
plt.xlim((-1400., 200.))
plt.ylim((-300.,  50.))
plt.gcf().set_size_inches(6, 3)

plt.ylabel('Y (km)', fontsize=14)
plt.xlabel('X (km)', fontsize=14)
cbar = plt.colorbar(pcol, shrink=0.5)
cbar.ax.tick_params(labelsize=14)

plt.title(label)


def update(t):
    txt = ' Bottom tracer concentration day: {00}'.format(t+1)
    print('Process day', t+1)
    pcol = plt.pcolor(lon0, lat0, np.squeeze(tra0[t, :, :]),
                      vmin=0.01, vmax=1., edgecolor='0.9', cmap=mycmap)
    pcol = plt.pcolor(lon1[jmin:jmax, imin:imax], lat1[jmin:jmax, imin:imax], np.squeeze(tra1[t, jmin:jmax, imin:imax]),
                      vmin=0.01, vmax=1., edgecolor='0.4', cmap=mycmap)
    plt.title(txt)
    return pcol


anim = FuncAnimation(fig, update, frames=np.arange(1, 39), interval=150)
anim.save('DOME_anim_btra.gif', dpi=80, writer='imagemagick')
# plt.savefig('DOME_bottom_tracer.png')
# plt.show()
