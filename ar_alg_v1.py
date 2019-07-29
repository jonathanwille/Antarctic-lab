
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

#Define the boundaries of the quadrants and halves in relation to ERA-Interim data. 
#There are three different lines where the quadrant/half names are inputed
#Also the AR time index output files have to be changed based on the quadrant/half used

EAIS = slice(0,120)
DDU = slice(120,240)
WAIS = slice(240,360)
wed = slice(360,480)
EAIS_DDU = slice(0,240)
WAIS_wed = slice(240,480)
WAIS_wed2 = slice(210,480)
amery = slice(65,185)
pen = slice(320,440)
all = slice(0,480)
law = slice(90,210)
WAIS_DDU = slice(120,360)
DDU2 = slice(100,280) 


#Read in the 12 files containing the monthly precipitable water values from 1979_2017


os.chdir("/media/willej/storage1/data/ERAI/pwat/month/79_17")
#os.chdir("/media/willej/8A309E89309E7C3F/data/ERAI/pwat_full")
datevar_list = []
mon = ['01','02','03','04','05','06','07','08','09','10','11','12']


file2 = "/media/willej/storage1/data/ERAI/land_sea.nc"
fh2 = Dataset(file2, mode='r',format="NETCDF4")
lsm = fh2.variables['lsm'][0,::-1,WAIS_wed]
lat_mask = fh2.variables['latitude'][::-1]  


#The minimum latitude that equals 0(ocean) for each longitude value 
coast_idx = np.argmin(lsm,axis=0)


b = 0 #counter for the month
for file in sorted(glob.glob("*")):
	print(file)


#Read in data

	fh1 = Dataset(file, mode='r',format="NETCDF4")

#	lat = fh1.variables['lat'][20:61]  #The latitude range of consideration is from -45 to -75
	lat = fh1.variables['lat'][20:67]  #The latitude range of consideration is from -45 to -79.5
	lon = fh1.variables['lon'][WAIS_wed]  #The longitude range of consideration is from 0 to 89.25 for the DM-PE quadrent /EAIS
	pwat = fh1.variables['pwat'][:,20:67,WAIS_wed]
	time2=fh1.variables['time'][:]
	units=fh1.variables['time'].units
	



	try :
		cal = fh1.variables['time'].calendar
	except AttributeError : 
	# Attribute doesnt exist
		cal = 'standard' 
	datevar=(nc.num2date(time2,units=units,calendar=cal))
	fh1.close()
#Run the algorithm

	pwat_per  = np.percentile(pwat, 98, axis=0) #Calc the 98th percentile of pwat at all points
	
	indices = np.where(pwat > pwat_per)
	time = indices[0]
	y = lat[indices[1]]
	x = lon[indices[2]]
	x2 = indices[2]


	timesteps = unique(time) # unique timesteps, avoids if(timestep != prev_timestep) loop

	river_idx = []
	lat_min_idx = []

	for timestep in timesteps :
		timestep_idx = np.where(timestep == time)[0]
		y_idx = y[timestep_idx]
		x_idx = x[timestep_idx]
		x2_idx = x2[timestep_idx]

		if not (y_idx[:-1] - y_idx[1:] > 3).any(): # checks the difference between successor/predecessor latitudes at the array level (no loop)
			if y_idx.max() - y_idx.min() > 20 :
				river_idx.append(timestep) 
				i = np.argmin(y_idx)
				t = np.argmax(y_idx)
					#print(timestep)
				#	print(i)
					#print(x2[i])
					#print(coast_idx[x2[i]])
					#print(lat_mask[coast_idx[x2[i]]])
					#print(y_idx[i])
					#print(y_idx[t])
					#print(y_idx[i])
#					print(x_idx[i])

				if lat_mask[coast_idx[x2[i]]] >= y_idx[i]:
					lat_min_idx.append(timestep)
					y_idx = np.asarray(y_idx)
					x_idx = np.asarray(x_idx)
					#	with open('/home/willej/data/ERAI/river_coord/'+str(datevar[timestep])+'_lat_WAIS_wed_pwat_min_idx_old.dat', 'w') as file:
					#		file.write('\n'.join(str(idx) for idx in y_idx))
					#	with open('/home/willej/data/ERAI/river_coord/'+str(datevar[timestep])+'_lon_WAIS_wed_pwat_min_idx_old.dat', 'w') as file:
					#		file.write('\n'.join(str(idx) for idx in x_idx))
					#	np.save('/home/willej/data/ERAI/river_coord/'+str(datevar[timestep])+'_lat_WAIS_wed_river_min_idx.dat',y_idx)
					#	np.save('/home/willej/data/ERAI/river_coord/'+str(datevar[timestep])+'_lon_WAIS_wed_river_min_idx.dat',x_idx)

			
#	print(river_idx)
#	print(lat_min_idx)
#	print(datevar[lat_min_idx])
	datevar_list.append(datevar[lat_min_idx])


######################################## Creation fichiers textes
#Create two text files: One for the AR time indices for al rivers and one for ARs that make landfall
#Also I haven't been able to figure out how to combine these two indicies into one text file
	with open('/home/willej/Documents/data/ERAI/'+mon[b]+'_79_17_WAIS_wed_river_idx_v1_3.dat', 'w') as file:
		file.write('\n'.join(str(idx) for idx in river_idx))
	with open('/home/willej/Documents/data/ERAI/'+mon[b]+'_79_17_WAIS_wed_river_min_idx_v1_3.dat', 'w') as file:
		file.write('\n'.join(str(idx) for idx in lat_min_idx))

	print(mon[b]) #output the month to check if the file matches the correct month
	b = b + 1
#v2 98 and 1 degrees
#v3 98 and 3 degrees
#v_old 98, 5 cont, and 2 long

	#delete all the stored variables in case of a future problem
	del river_idx
	del fh1
	del lat
	del lon
	del pwat
	del time2
	del units
	del pwat_per
	del y
	del time
	del indices
	del lat_min_idx


datevar_list = np.asarray(datevar_list)
print(len(datevar_list))
#with open('/home/willej/data/ERAI/ar_dates/wed_ar_dates.dat', 'w') as file:
#	file.write('\n'.join(str(idx) for idx in datevar_list))

print(datevar_list)

























