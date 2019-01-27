ncks -O -d y,17,17,1 WAD_1ts_00010101_00010101_grid_T.nc sshtime.nc
ncks -4 -A -v gdepw_0,ht_wd -d y,17,17,1 -d z,10,10,1 -C mesh_mask.nc sshtime.nc

python2.7 matpoly2.py sshtime.nc
animate wadfr*.png
