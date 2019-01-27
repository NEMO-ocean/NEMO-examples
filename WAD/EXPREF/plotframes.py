#!/usr/bin/env python
################################################################################## 
# Utility to produce animation frames directly from the WAD_TEST_CASE output.
# It can be used without arguments, i.e.:
#     python plotframes.py
# in which case a frame is created every tenth time level from SSH data extracted from
# WAD_1ts_00010101_00010101_grid_T.nc along the centre of the basin (j=17). A closed 
# basin is assumed and frames are named wadfr0000.png etc. Bathymetry information is 
# extracted from the mesh_mask.nc file. The frames are annotated with a timestamp that
# assumes an 18s baroclinic timestep.
# 
# All these settings can be overridden with command-line arguments. See:
#     python plotframes.py -h 
# for details. For example:
#     python plotframes.py  -nt 300 -stride 30 -froot mywad
# 
# Two major variations are also supported for specific test cases:
#     python plotframes.py -obc
# plots the right-hand side of the basin as an open boundary (test case 7) and:
#     python plotframes.py -use_sal
# colours each gridcell according to its salinity value (test case 6)
################################################################################## 
import os, sys
from argparse import ArgumentParser
import numpy as np
import netCDF4

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from netCDF4 import Dataset

#
# Turn off the unhelpful warning about open figures
#
matplotlib.rcParams['figure.max_open_warning'] = 0

if __name__ == '__main__':
    parser = ArgumentParser(description=
    """
    produce frames for the animation of results from the WAD_TEST_CASES.
    Mostly this can be run without arguments but command line arguments may be
    used to override defaults. These are necessary in cases with open boundaries
    (e.g. nn_wad_test=7) and cases where it is desired to show variations in salinity
    (e.g. nn_wad_test=6).
     e.g. plotframes.py -tfile <T-grid file> -bfile  <bathymetry file>
    -froot <root name for frames>
    -nt <maximum number of time frames to process>
    -stride <stride through time frames> 
    -rdt <length of baroclinic timestep (s)> 
    -obc
    -use_sal
    """)
    parser.add_argument('-tfile',dest='tfile',help='T-grid file if not WAD_1ts_00010101_00010101_grid_T.nc', default='WAD_1ts_00010101_00010101_grid_T.nc')
    parser.add_argument('-bfile',dest='bfile',help='Bathymetry file if not mesh_mask.nc', default='mesh_mask.nc')
    parser.add_argument('-j',dest='jrow',help='jrow; j-row to extract and plot (default: 17 (fortran index))', type=int, default=16)
    parser.add_argument('-froot',dest='froot',help='froot; root name for frames (default: wadfr)', default='wadfr')
    parser.add_argument('-nt',dest='nfmax',help='nfmax; maximum number of frames to produce', type=int, default=None)
    parser.add_argument('-stride',dest='tinc',help='tinc; stride through time frames (default: 10)', type=int, default=10)
    parser.add_argument('-rdt',dest='rdt',help='rdt; length of baroclinic timestep (s) (default: 18.0)', type=float, default=18.0)
    parser.add_argument('-obc',help='Right-hand side boundary is open', action="store_true")
    parser.add_argument('-use_sal',help='colour polygons according to salinity variations', action="store_true")
    args   = parser.parse_args()
    tfile  = args.tfile
    bfile  = args.bfile
    jrow   = args.jrow
    froot  = args.froot
    nfmax  = args.nfmax
    stride = args.tinc
    rdt    = args.rdt
    obc    = args.obc
    use_sal= args.use_sal


fw = Dataset(tfile)
ssh = fw.variables['sossheig'][:,jrow,:]
vot = fw.variables['sosaline'][:,jrow,:]
if use_sal:
 sal = fw.variables['vosaline'][:,:,jrow,:]
 nz = sal.shape[1]
fw.close()

fw = Dataset(bfile)
bat = fw.variables['ht_wd'][0,jrow,:]
if use_sal:
 mbat = fw.variables['mbathy'][0,jrow,:]
fw.close()

#print "ssh"
#print ssh.shape
#print "bat"
#print bat.shape 
#print "vot"
#print vot.shape 
nt = ssh.shape[0]
nx = ssh.shape[1]
#print nx,nt

bat = -1.*bat
batmin = np.amin(bat)
batmax = np.amax(bat)
brange = batmax - batmin
tol = 0.1*brange
#print batmin,batmax,' ho'
if obc:
 batrhs = batmin
else:
 batrhs = batmax

if nfmax is None:
  nfmax = nt

nf = 0
ntmax = np.minimum(nt,nfmax*stride)

if not use_sal:
#
# plot solid single colour polygons just showing ssh variation
#
  for t in range(0,ntmax,stride):
   wadfr = froot+"{:0>4d}.png".format(nf)
   nf = nf + 1
  
   tfac = rdt/3600.0
   t24  = np.int(np.mod(t*tfac,24))
   dy   = np.int(t*tfac/24.0)
   mn   = np.int(np.rint((np.mod(t*tfac,24) - t24 )*60))
   hour  = "t={:0>2d}:{:0>2d}:{:0>2d} ".format(dy,t24,mn)
   hour2  = "  (days:hrs:mins)"
   batpts = np.zeros((nx+4,2))
   sshpts = np.zeros((2*nx,2))
   votpts = np.zeros((nx,2))
  
   for pt in range(nx):
      batpts[pt+2,0] = pt
      batpts[pt+2,1] = bat[pt]
      sshpts[pt,0] = pt
      sshpts[pt,1] = ssh[t,pt]
      votpts[pt,0] = pt
      votpts[pt,1] = np.minimum(36.,vot[t,pt])
      votpts[pt,1] = np.maximum(30.0,votpts[pt,1])
      votpts[pt,1] = batmin +0.2*brange + (votpts[pt,1]-30.)*brange/6.0
  
   batpts[nx+1,1] = batrhs
   batpts[nx+2,0] = nx-1
   batpts[nx+2,1] = batrhs
   batpts[nx+3,0] = nx-1
   batpts[nx+3,1] = batmin
   batpts[0,0] = 0.0
   batpts[0,1] = batmin
   batpts[1,0] = 0.0
   batpts[1,1] = batmax
   batpts[2,1] = batmax
   sshpts[nx-1,0] = nx-1
   sshpts[nx-1,1] = sshpts[nx-2,1]
   sshpts[0,0] = 0.0
   votpts[nx-2,0] = nx-1
   votpts[nx-1,0] = nx-1
   votpts[0,0] = 0.0
   votpts[nx-1,1] = batmax + tol
   votpts[0,1] = batmax + tol
   sshpts[0,1] = sshpts[1,1]
   for pt in range(nx):
      sshpts[pt+nx,0]=batpts[nx+1-pt,0]
      sshpts[pt+nx,1]=batpts[nx+1-pt,1]
   
  
   xs, ys = zip(*votpts)
   fig, ax = plt.subplots()
  
   patches = []
   polygon = Polygon(batpts, True)
   patches.append(polygon)
   p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=1.0)
   p.set_facecolors(['#f1a9a9'])
  
   patches = []
   polygon2 = Polygon(sshpts, True)
   patches.append(polygon2)
   p2 = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=1.0)
   p2.set_facecolors(['#44a1ff'])
   
# Maximum depth set here to -10m
   ax.set_ylim([-10., 6.0])
   ax.set_xlim([0., 51.0])
   ax.add_collection(p2)
   ax.add_collection(p)
   ax.plot(xs,ys, '--', color='black', ms=10)
  
   plt.annotate(hour,xy=(2,batmin+0.1*brange))
   plt.annotate(hour2,xy=(2,batmin+0.05*brange))
   plt.savefig(wadfr)

else:
#
# plot each gridcell coloured according to its salinity value
#
  for t in range(0,ntmax,stride):
   wadfr = froot+"{:0>4d}.png".format(nf)
   nf = nf + 1

   tfac = rdt/3600.0
   t24  = np.int(np.mod(t*tfac,24))
   dy   = np.int(t*tfac/24.0)
   mn   = np.int(np.rint((np.mod(t*tfac,24) - t24 )*60))
   hour  = "t={:0>2d}:{:0>2d}:{:0>2d} ".format(dy,t24,mn)
   hour2  = "  (days:hrs:mins)"
   batpts = np.zeros((nx+4,2))
   votpts = np.zeros((nx,2))
   salpts = np.zeros((nx*nz,6,2))
   salmin = 28.
   salmax = 37.
   salrange = salmax - salmin
   faccol = np.zeros((nx*nz))
   cl = 0
   for pt in range(nx):
      batpts[pt+2,0] = pt
      batpts[pt+2,1] = bat[pt]
      votpts[pt,0] = pt
      votpts[pt,1] = np.minimum(35.,vot[t,pt])
      votpts[pt,1] = np.maximum(30.0,votpts[pt,1])
      votpts[pt,1] = batmin +0.2*brange + (votpts[pt,1]-30.)*brange/6.0
   batpts[nx+1,1] = batmax
   batpts[nx+2,0] = nx-1
   batpts[nx+2,1] = batmax
   batpts[nx+3,0] = nx-1
   batpts[nx+3,1] = batmin
   batpts[0,0] = 0.0
   batpts[0,1] = batmin
   batpts[1,0] = 0.0
   batpts[1,1] = batmax
   batpts[2,1] = batmax
   votpts[nx-2,0] = nx-1
   votpts[nx-1,0] = nx-1
   votpts[0,0] = 0.0
   votpts[nx-1,1] = batmax + tol
   votpts[0,1] = batmax + tol
  
   cl = 0
   for pt in range(nx):
      mz = np.maximum(1,mbat[pt])
      im1 = np.maximum(pt-1,1)
      ip1 = np.minimum(pt+1,nx-1)
      dz  = (ssh[t,pt] - batpts[pt+2,1] )/mz
      dz1 = 0.5*(ssh[t,pt] + ssh[t,im1] - batpts[pt+2,1] - batpts[im1+2,1] )/mz
      dz2 = 0.5*(ssh[t,pt] + ssh[t,ip1] - batpts[pt+2,1] - batpts[ip1+2,1] )/mz
      dz  = np.maximum(dz ,0.0)
      dz1 = np.maximum(dz1,0.0)
      dz2 = np.maximum(dz2,0.0)
      bat1 = 0.5*( batpts[pt+2,1] + batpts[im1+2,1] )
      bat2 = 0.5*( batpts[pt+2,1] + batpts[ip1+2,1] )
      ptm = np.maximum(pt-0.5,0.0)
      ptx = np.minimum(pt+0.5,nx-1)
      for z in range(mz):
        if ( sal[t,mz-1-z,pt] > 0.0 ):
          salpts[cl,0,0] = pt
          salpts[cl,1,0] = ptm
          salpts[cl,2,0] = ptm
          salpts[cl,3,0] = pt
          salpts[cl,4,0] = ptx
          salpts[cl,5,0] = ptx
          salpts[cl,0,1] = batpts[pt+2,1] +dz*z
          salpts[cl,1,1] = bat1 +dz1*z
          salpts[cl,2,1] = bat1 +dz1*(z+1)
          salpts[cl,3,1] = batpts[pt+2,1] +dz*(z+1)
          salpts[cl,4,1] = bat2 +dz2*(z+1)
          salpts[cl,5,1] = bat2 +dz2*z
          faccol[cl] = 100*(sal[t,mz-1-z,pt] - salmin) / salrange
          faccol[cl] = np.maximum(faccol[cl],0.0)
          faccol[cl] = np.minimum(faccol[cl],100.0)
          cl = cl + 1
   
  
   votpts2 = votpts[2:nx-4,:]
   xs, ys = zip(*votpts2)
   fig, ax = plt.subplots()
   patches = []
  
   for pt in range(cl-1):
      polygon = Polygon(salpts[pt,:,:], True)
      patches.append(polygon)
  
   p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=1.0)
   
   p.set_array(faccol)
   p.set_edgecolor('face')
  
   patches = []
   polygon = Polygon(batpts, True)
   patches.append(polygon)
   p2 = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=1.0)
   p2.set_facecolors(['#f1a9a9'])
  
# Maximum depth set here to -8m (suitable for test case 6 only)
   ax.set_ylim([-8., 6.0])
   ax.set_xlim([0., 51.0])
   ax.add_collection(p)
   ax.add_collection(p2)
   ax.plot(xs,ys, '--', color='black', ms=10)
  
   plt.annotate(hour,xy=(2,batmin+0.1*brange))
   plt.annotate(hour2,xy=(2,batmin+0.05*brange))
   plt.savefig(wadfr)
