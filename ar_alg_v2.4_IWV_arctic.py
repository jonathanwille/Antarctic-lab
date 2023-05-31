
#from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

import netCDF4 as nc


from pylab import *
#from mpl_toolkits.basemap import Basemap, maskoceans, addcyclic, shiftgrid

import io
from matplotlib import animation
import sys
import operator
import glob, os
import pickle
import csv
np.set_printoptions(threshold=sys.maxsize)
#import pdb; pdb.set_trace()



#Read in the 12 files containing the monthly precipitable water values from 1980_2018



#os.chdir("/media/willej/storage1/data/ERAI/pwat/month/79_17")
datevar_list = []
month = ['01','02','03','04','05','06','07','08','09','10','11','12']
#month = ['04','05','06','07','08','09','10','11','12']
file3 = "/home/jwille/xenon_data/IGE_data/land_sea_masks/merra2_ocn_land_fraction.nc"
fh3 = Dataset(file3, mode='r',format="NETCDF4")

lat = fh3.variables['lat'][:] 
lon = fh3.variables['lon'][:] 
lon = lon + 180 
lonvar = fh3.variables['lon']
latvar = fh3.variables['lat']
#lat_scan = fh3.variables['lat'][::-1] 
#lat_scan = lat[255:341]  #The latitude range of consideration is from 37.5 to 80
#lat_scan = lat[265:351]  #The latitude range of consideration is from 42.5 to 85
lat_scan = lat[255:351]  #The latitude range of consideration is from 37.5 to 85


#lon_scan = lon + 180
nlon = lon.size
nlat = lat.size


lon2,lat2 = np.meshgrid(lon,lat)
res_lon = lon[1] - lon[0]
res_lat = lat[1] - lat[0]



#lat_lon_land = np.column_stack((lat_land, lon_land))

origin_lat,origin_lon = lat2[0,0],lon2[0,0]
lat_stepsize = lat2[1,0] - lat2[0,0]
lon_stepsize = lon2[0,1] - lon2[0,0]


#The minimum latitude that equals 0(ocean) for each longitude value 


b = 0 #counter for the month

for mon in month: #loop through the files; the MAR melting data and the AR indicies 
	print(mon)

	#file1 = '/media/willej/8A309E89309E7C3F/data/vIVT_IWV/MERRA2/80_20/ivt/lat_rst_arctic_ARTMIP_MERRA_vIVT_'+mon+'_80_20.nc'
	file1 = '/home/jwille/xenon_data/IGE_data/vIVT_IWV/MERRA2/80_22/IWV/lat_rst_arctic_ARTMIP_MERRA_IWV_'+mon+'_80_22.nc'


#for file in sorted(glob.glob("*")):
#	print(file)


#Read in data

	fh1 = Dataset(file1, mode='r',format="NETCDF4")
	#fh2 = Dataset(file2, mode='r',format="NETCDF4")

	lat_test = fh1.variables['lat'][:]  
	lon_test = fh1.variables['lon'][:]  
	iwv_full = fh1.variables['IWV'][:,:,:]


	iwv = iwv_full[:,15:111,:]   #The latitude range of consideration is from 37.5 to 85


	fh1.close()

	lat_index_list = []
	lon_index_list = []
#Run the algorithm



	iwv_per  = np.percentile(iwv, 98, axis=0) #Calc the 98th percentile of pwat at all points
	


	indices = np.where(iwv > iwv_per)
	time = indices[0]
	y = lat_scan[indices[1]]
	x = lon[indices[2]]


	#timesteps = unique(time) # unique timesteps, avoids if(timestep != prev_timestep) loop
	timesteps = np.arange(0,len(iwv[:,0,0]))
	#timesteps = np.arange(0,323)
	river_idx = []

	for timestep in timesteps :
		timestep_idx = np.where(timestep == time)[0]
		print(timestep)
		y_idx = y[timestep_idx]
		x_idx = x[timestep_idx]
		#print(y_idx)

		y_splitted_temp = np.split(y_idx, np.where(np.diff(y_idx) > res_lat)[0] +1)
		x_splitted_temp = np.split(x_idx, np.where(np.diff(y_idx) > res_lat)[0] +1)

	
		#print(y_splitted_temp)

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

	np.save('/home/jwille/xenon_data/IGE_data/ar_catalogues/merra2/coord_idx/iwv/'+mon+'_80_22_arc_iwv_lat_idx_v2.4',lat_index_list)
	np.save('/home/jwille/xenon_data/IGE_data/ar_catalogues/merra2/coord_idx/iwv/'+mon+'_80_22_arc_iwv_lon_idx_v2.4',lon_index_list)






