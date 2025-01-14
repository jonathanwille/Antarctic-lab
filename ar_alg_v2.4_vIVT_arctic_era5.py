
#from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import shapefile
import netCDF4 as nc
import cartopy

from pylab import *
#from mpl_toolkits.basemap import Basemap, maskoceans, addcyclic, shiftgrid
from scipy.interpolate import griddata
import io
from matplotlib import animation
import sys
import operator
import glob, os
import pickle
import csv
np.set_printoptions(threshold=sys.maxsize)
#import pdb; pdb.set_trace()



#Read in the 12 files containing the monthly values from 1979_2017

month = ['01','02','03','04','05','06','07','08','09','10','11','12']



file_test = "/net/thermo/atmosdyn/era5/cdf/2023/12/Z20231231_21"
fh = Dataset(file_test, mode='r',format="NETCDF4")
lat = fh.variables['lat'][:]  

lat_plot = fh.variables['lat'][:]
lon_plot = fh.variables['lon'][:]
lon = fh.variables['lon'][:]  
#lon = lon - 180
lonvar = fh.variables['lon']
latvar = fh.variables['lat']
lat_scan = lat[254:351]  #The latitude range of consideration is from 37 to 85
nlon = lon.size
nlat = lat.size
lon2,lat2 = np.meshgrid(lon,lat_plot)
res_lon = lon[1] - lon[0]
res_lat = lat[1] - lat[0]
origin_lat,origin_lon = lat2[0,0],lon2[0,0]
lat_stepsize = lat2[1,0] - lat2[0,0]
lon_stepsize = lon2[0,1] - lon2[0,0]

timevar = fh.variables['time']


b = 0 #counter for the month
for mon in month: #loop through the files; the MAR melting data and the AR indicies 
	print(mon)

	file1 = "/home/jwille/xenon_data/IGE_data/vIVT_IWV/ERA5/coarse/80_23/vIVT/lat_rst_arctic_vIVT_ERA5_2DIVT_"+mon+"_80_23.nc"
	print(file1)
	#Read in data

	fh1 = Dataset(file1, mode='r',format="NETCDF4")
	lat_test = fh1.variables['lat'][:]  
	lon_test = fh1.variables['lon'][:]  
	ivt_full = fh1.variables['vIVT'][:,:,:]
	ivt = ivt_full[:,14:111,:]   #The latitude range of consideration is from 37 to 85
	ivt[ivt < 0] = 0


	fh1.close()

	#Run the algorithm

	lat_index_list = []
	lon_index_list = []
	landfall_idx = []
	ivt_per  = np.percentile(ivt, 98, axis=0) #Calc the 98th percentile of vIVT at all points

	indices = np.where(ivt > ivt_per)
	time = indices[0]
	y = lat_scan[indices[1]]
	x = lon[indices[2]]

	river_idx = []
	landfall_idx = []
	timesteps = np.arange(0,len(ivt[:,0,0]))

	for timestep in timesteps :
		timestep_idx = np.where(timestep == time)[0]
		print(timestep)
		y_idx = y[timestep_idx]
		x_idx = x[timestep_idx]

		y_splitted_temp = np.split(y_idx, np.where(np.diff(y_idx) > res_lat)[0] +1)
		x_splitted_temp = np.split(x_idx, np.where(np.diff(y_idx) > res_lat)[0] +1)



		y_longest = max(y_splitted_temp, key=len)
		x_longest = max(x_splitted_temp, key=len)


		x_reverse = []
		y_reverse = []



		try:
			if y_longest.max() - y_longest.min() > 20 :
				reverse_grid = np.arange(min(x_longest),max(x_longest)+0.5,res_lon)
				for i in reverse_grid:
					x_index_reverse = np.where(x_longest == i)
					x_reverse = np.concatenate((x_reverse,x_longest[x_index_reverse]))
					y_reverse = np.concatenate((y_reverse,y_longest[x_index_reverse]))

		except ValueError:
			pass
		x_splitted = np.split(x_reverse, np.where(np.diff(x_reverse) > 20)[0] +1)
		y_splitted = np.split(y_reverse, np.where(np.diff(x_reverse) > 20)[0] +1)



		try:
			if x_splitted[0][0]+360 - x_splitted[-1][-1] < 20:
				x_splitted[-1] = np.concatenate((x_splitted[-1],x_splitted[0]))
				x_splitted = np.delete(x_splitted,0,0)

				y_splitted[-1] = np.concatenate((y_splitted[-1],y_splitted[0]))
				y_splitted = np.delete(y_splitted,0,0)

	
		except IndexError:
			pass		




		x_shape = []
		y_shape = []

		x_shape_landfall = []
		y_shape_landfall = []

		for i in range(0,len(y_splitted)):
			x_reverse2 = []
			y_reverse2 = []

			x_final = []
			y_final = []
	
			n = 0

			try:	
				reverse_grid2 = np.arange(min(y_splitted[i]),max(y_splitted[i])+0.5,res_lat)
				for j in reverse_grid2:
					y_index_reverse2 = np.where(y_splitted[i] == j)
					x_reverse2 = np.concatenate((x_reverse2,x_splitted[i][y_index_reverse2]))
					y_reverse2 = np.concatenate((y_reverse2,y_splitted[i][y_index_reverse2]))


			except ValueError:
				pass

			try:
				y_splitted_final = np.split(y_reverse2, np.where(np.diff(y_reverse2) > res_lat)[0] +1)
				x_splitted_final = np.split(x_reverse2, np.where(np.diff(y_reverse2) > res_lat)[0] +1)

				for z in range(0,len(y_splitted_final)):
					if y_splitted_final[z].max() - y_splitted_final[z].min() > 20:
						y_final = np.concatenate((y_final, y_splitted_final[z]))
						x_final = np.concatenate((x_final, x_splitted_final[z]))
				
				x_shape = np.concatenate((x_shape,x_final))
				y_shape = np.concatenate((y_shape,y_final))



			except ValueError:
				pass

		lat_index = (y_shape - origin_lat) / lat_stepsize 
		lon_index = (x_shape - origin_lon) / lon_stepsize 
		#	print(y_shape)
			#print(lat_index)
		lat_index = lat_index.astype(int)
		lon_index = lon_index.astype(int)
			#print(y_shape)
			#print(lat_index)
		if(len(lat_index) > 0):
			lat_index_list.append(lat_index.data)
		else:
			lat_index_list.append(lat_index)
	

		if(len(lon_index) > 0):
			lon_index_list.append(lon_index.data)
		else:
			lon_index_list.append(lon_index)


	ivt_full = None
	ivt = None
	np.save('/home/jwille/xenon_data/IGE_data/ar_catalogues/era5/coord_idx/'+mon+'_80_23_arc_vivt_lat_idx_v2.4',lat_index_list)
	np.save('/home/jwille/xenon_data/IGE_data/ar_catalogues/era5/coord_idx/'+mon+'_80_23_arc_vivt_lon_idx_v2.4',lon_index_list)



