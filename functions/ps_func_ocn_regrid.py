#! /usr/bin/env python
#coding=utf-8
#--------------------------------------------------------------------
# for micom regriding,               2014 0815 by M. Shen
#      modify, creating netcdf file, 2015 0227 by M. Shen
# load module
#
import time
start_time = time.time()
import sys
import os
import numpy as np
import scipy.io.netcdf as nio
import scipy.interpolate as interp
from netCDF4 import Dataset as NCDataset
# from Scientific.IO.NetCDF import NetCDFFile as NCDataset  # No Scientific on Hexagon
print("--- %s s is used for loading libs! ---" % (time.time() - start_time))
#time.sleep(1.1)
#--------------------------------------------------------------------
# for testing 
#
# from the projected coordinate system (EPSG 3857) to a geographic coordinate system (EPSG 4326).
#
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#import matplotlib.mlab as mlab
#np.set_printoptions(threshold=np.nan)
#--------------------------------------------------------------------
# define parameters
#
if len(sys.argv) < 2:
  print "The file name and variable name are necessary!!!"
  sys.exit()

DataFilename = sys.argv[1]
vname = sys.argv[2]
#SourceBackupNC = sys.argv[3]

#vname = "var2"
#DataFilename = "/Users/earnestshen/temp/pSUMO_mpiom_22000101_22001231.nc"
#SourceBackupNC = "IGT15_var3Dlev.nc"
if sys.platform == "linux2" :
  #GridPath = "/work/earnest/nn9039k/noresm/inputdata/ocn/micom/tnx2v1/20130206"
  GridPath = "/work/shared/noresm/inputdata/ocn/micom/tnx2v1/20130206"

if sys.platform == "darwin" :
  GridPath = "/Users/earnestshen/Research/Data/OBData/T31/GridInfo"

FGrid = GridPath + "/grid.nc"
GData = nio.netcdf_file(FGrid)
grids_lon = GData.variables['plon'].data
grids_lat = GData.variables['plat'].data
temp = np.reshape(grids_lon, grids_lon.size)
grids = np.zeros((grids_lon.size,0))
grids = np.append( grids, np.atleast_2d(temp).transpose(), 1)
temp = np.reshape(grids_lat, grids_lat.size)
grids = np.append( grids, np.atleast_2d(temp).transpose(), 1)
# local NC file
tempNCfn = "tempIGT_" + vname + ".nc"
#os.system("cp " + SourceBackupNC + " " + tempNCfn)

#  grids = readnetgrid(FGrid,'oces')
## griddata can be used 
## WOA grid
#RLW = -179.5
#RLE = 179.5
#RLS = -89.5
#RLN = 89.5
# EN4 grid
RLW = -180
RLE = 180
RLS = -90
RLN = 90

dlon = 1
dlat = 1
Long = np.arange(RLW,RLE+dlon,dlon)
Lati = np.arange(RLS,RLN+dlat,dlat)[::-1]
mLo, mLa = np.meshgrid(Long, Lati)
interpMethod = 'linear' # linear, cubic, nearest
#--------------------------------------------------------------------
# Read target file
#
VData = nio.netcdf_file(DataFilename)
VarC = VData.variables[vname]
Vtime = VData.variables["time"].data
missV = VarC._attributes['_FillValue']
data = np.float32(VarC.data.copy())
data[data==missV] = np.nan
try : 
  print "Re-range the data ..."
  data = ( data * VarC.scale_factor ) + VarC.add_offset
except : 
  print "No need to re-range the data ..."


print "ID: ", ("--- %s s is used for parameters defining! ---" % (time.time() - start_time))

fData = NCDataset(DataFilename, 'r')
varData = fData.variables[vname]
varTData = fData.variables['time']

#tData.variables['var1'][:] = np.float32(regriddata)

#--------------------------------------------------------------------
# regriding
#
# regriddata = zeros(size(mlon',1),size(mlat',2),size(Lev,1),12);
# for checking: grids [38,4] is (179.6613,-0.1405)
if len(data.shape) == 4 : # 3D case
  regriddata = np.float32(np.zeros((data.shape[0], data.shape[1], mLo.shape[0], mLo.shape[1])))
  varDepth = fData.variables['depth']
if len(data.shape) == 3 :
  regriddata = np.float32(np.zeros((data.shape[0], mLo.shape[0], mLo.shape[1])))

#for im in range(data.shape[0]) : # for month
for im in range(1) : # for month
  if len(data.shape) == 4 : 
    for iD in range(data.shape[1]) : # for diff depth
      temp = data[im,iD,:,:].copy().reshape(grids_lon.size)
      #temp = data[im,data.shape[1]-(iD+1),::-1,:].copy().reshape(grids_lon.size)
      temp_T = np.atleast_2d(temp).transpose() 
      temp_griddata = interp.griddata(grids, temp_T, (mLo, mLa), method=interpMethod)
      temp_put = temp_griddata.reshape(mLo.shape[0], mLo.shape[1]).copy()
      temp_put[np.isnan(temp_put)] = missV
      regriddata[im,iD,:,:] = np.float32(temp_put)
  if len(data.shape) == 3 :
    temp = data[im,::-1,:].copy().reshape(grids_lon.size)
    temp_T = np.atleast_2d(temp).transpose()
    temp_griddata = interp.griddata(grids, temp_T, (mLo, mLa), method=interpMethod)
    temp_put = temp_griddata.reshape(mLo.shape[0], mLo.shape[1]).copy()
    temp_put[np.isnan(temp_put)] = missV
    regriddata[im,:,:] = np.float32(temp_put)
  print "ID: ", ("--- %s s is used! ---" % (time.time() - start_time))
#  print data[im,0,60,4], 'v.s.', regriddata[im,39,64,120]
#  print 'C1: ', data.shape, ' C2: ', temp_put.shape, 'Check 2: ', temp.shape
#--------------------------------------------------------------------
# There is no Scientific on Hexagon
#
#print 'copy target -> ', tempNCfn
tData = NCDataset(tempNCfn, 'w')
#time = tData.variables['time'] 
#lats = tData.variables['lat']
#lons = tData.variables['lon']
tData.createDimension('time',None)
tData.createDimension('latitude',Lati.size)
tData.createDimension('longitude',Long.size)
time = tData.createVariable('time','float32',('time',))
lats = tData.createVariable('latitude','float32',('latitude',))
lons = tData.createVariable('longitude','float32',('longitude',))
lats.units = 'degrees_north'
lons.units = 'degrees_east'
for attr in varTData.ncattrs() : 
  exec "time." + attr + " = varTData." + attr 

time[:] = varTData 
lats[:] = Lati
lons[:] = Long
#temp = tData.variables['var1']
#for attr in varData.ncattrs() :
#  if attr != u'_FillValue' and attr != u'add_offset' and attr != u'add_offset': 
#    exec "temp." + attr + " = varData." + attr

if len(regriddata.shape) == 4 :
  tData.createDimension('depth',varDepth.shape[0])
  depth = tData.createVariable('depth','float32',('depth',))
  depth.units = 'm'
  for attr in varDepth.ncattrs() :
    exec "depth." + attr + " = varDepth." + attr
  depth[:] = varDepth
  temp = tData.createVariable(vname,'float32',('time','depth','latitude','longitude'), fill_value=missV)
if len(data.shape) == 3 :
		temp = tData.createVariable(vname,'float32',('time','latitude','longitude'), fill_value=missV)
#regriddata[np.isnan(regriddata)] = np.float32(missV)
#temparr = regriddata[:,::-1,:]
#temp[:] = temparr
temp[:] = regriddata

fData.close()
tData.sync()
tData.close()
GData.close()
#--------------------------------------------------------------------
# regriding test
#
#ax = plt.gca()
#def func(x, y):
#  return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2+x+y


#grid_x, grid_y = np.mgrid[0:1:10j, 0:2:20j]
#points = np.zeros((100,0))
#point_temp = np.random.rand(100,1)
#points = point_temp
#point_temp = np.random.rand(100,1)*2.0
#points = np.append( points, point_temp, axis=1)
#points = np.random.rand(100, 2)
#points[:,[1]] = points[:,[1]].__imul__(2)
#values = func(points[:,0], points[:,1])
#grid_z2 = interp.griddata(points, values, (grid_x, grid_y), method='cubic')
# figure simple verify
#ax = plt.gca()
#clevs = np.arange(np.nanmin(grid_z2),np.nanmax(grid_z2),0.05)
#pvar = regriddata[0,39,:,:].copy()
#clevs = np.arange(np.floor(np.nanmin(pvar)),np.floor(np.nanmax(pvar)),1)
#CS = plt.contour(mLo, mLa, pvar, linewidths=2, level=clevs, cmap=cm.jet)

#plt.show()

os.system("mv " + tempNCfn + " " + DataFilename[:-8:])
os.system("rm " + DataFilename)
print " ", "Done!!"
