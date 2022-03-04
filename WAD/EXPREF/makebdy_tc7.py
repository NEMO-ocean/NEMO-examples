from netCDF4 import Dataset
import numpy as np

pathout = "bdyssh_tc7"

nx = 23 
ny = 1
nt = 24
ndays=4

#-------------------------------------------------------
# Create bdyssh_tc7_m01d01.nc, bdyssh_tc7_m01d02.nc etc.
#-------------------------------------------------------

pathstart="{}_m12d30.nc".format(pathout)
for nd in range(ndays):
 print pathstart
 ssh = np.zeros((nt,ny,nx))
 for nnt in range(nd*nt,(nd+1)*nt):
   tx = 2.5*np.cos((3.141592654/6.0)*(nnt))
   print nnt, tx
   for nnx in range(nx):
     for nny in range(ny):
       ssh[nnt-nd*nt,nny,nnx] = tx


 fo = Dataset(pathstart, 'w', format='NETCDF4')
 nxo = fo.createDimension('x', nx)
 nyo = fo.createDimension('y', ny)
 nto = fo.createDimension('t', None)
 ssho = fo.createVariable('sshbdy', 'f4',('t','y','x'))

 ssho[:,:,:] = ssh[:,:,:]
 ssho.long_name = 'bdy ssh boundary condition'
 ssho.standard_name = 'sshbdy'
 ssho.units = 'm'
#
 fo.close()
 pathstart="{}_m01d{:0>2d}.nc".format(pathout,nd+1)

#-------------------------------------------------------
# Create bdyuv_tc7_m01d01.nc, bdyuv_tc7_m01d02.nc etc.
# u is -(1/H)*d(ssh)/dt; v =0.0
#-------------------------------------------------------

pathout = "bdyuv_tc7"


pathstart="{}_m12d30.nc".format(pathout)
for nd in range(ndays):
 print pathstart
 u = np.zeros((nt,ny,nx))
 for nnt in range(nd*nt,(nd+1)*nt):
  tx = 2.5*(3.141592654/6.0)*np.sin((3.141592654/6.0)*(nnt+1.0))/10.0
  print nnt, tx
  for nnx in range(nx):
    for nny in range(ny):
       u[nnt-nd*nt,nny,nnx] = tx

 v = np.zeros((nt,ny,nx))

 fo = Dataset(pathstart, 'w', format='NETCDF4')
 nxo = fo.createDimension('x', nx)
 nyo = fo.createDimension('y', ny)
 nto = fo.createDimension('t', None)
 uo = fo.createVariable('ubdy', 'f4',('t','y','x'))
 vo = fo.createVariable('vbdy', 'f4',('t','y','x'))

 uo[:,:,:] = u[:,:,:]
 uo.long_name = 'bdy u boundary condition'
 uo.standard_name = 'ubdy'
 uo.units = 'm/s'
#
 vo[:,:,:] = v[:,:,:]
 vo.long_name = 'bdy v boundary condition'
 vo.standard_name = 'vbdy'
 vo.units = 'm/s'
#
 fo.close()
 pathstart="{}_m01d{:0>2d}.nc".format(pathout,nd+1)
