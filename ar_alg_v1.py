#Atmospheric river detection algorithm using IWV. Version 1.0.0

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import shapefile
import netCDF4 as nc
import cartopy
from pylab import *
from scipy.interpolate import griddata
import io
from matplotlib import animation
import sys
import operator
import glob, os
import pickle
import csv


#Define the boundaries of the quadrants and halves in relation to reanalysis dataset
#There are three different lines where the quadrant/half names are inputed
#Line 43
#Line 57
#Line 58

#Also the AR time index output files have to be changed based on the quadrant/half used

EAIS = slice(0,120) #The DML-PEL quadrant 
DDU = slice(120,240) #The Wilkes Land quadrant 
WAIS = slice(240,360) #The WAIS quadrant
wed = slice(360,480) # The AP-Weddell quadrant
EAIS_DDU = slice(0,240) # East Antarctica half
WAIS_wed = slice(240,480) # West Antarctica half

#Read in the 12 files containing the monthly precipitable water values from 1979-2017
os.chdir("/media/willej/storage1/data/ERAI/pwat/month/79_17")

datevar_list = [] #empty list later used to store the AR dates
mon = ['01','02','03','04','05','06','07','08','09','10','11','12']

file2 = "/media/willej/storage1/data/ERAI/land_sea.nc"  #Land/sea mask file
fh2 = Dataset(file2, mode='r',format="NETCDF4")
lsm = fh2.variables['lsm'][0,::-1,WAIS_wed] #Reverse the land mask looking from the 90 S to north
lat_mask = fh2.variables['latitude'][::-1]  


#The minimum latitude that equals 0(ocean) for each longitude value 
coast_idx = np.argmin(lsm,axis=0)

b = 0 #counter for the month
for file in sorted(glob.glob("*")):
	print(file)

#Read in data
	fh1 = Dataset(file, mode='r',format="NETCDF4")
	lat = fh1.variables['lat'][20:67]  #The latitude range of consideration is from -45 to -79.5
	lon = fh1.variables['lon'][WAIS_wed] 
	pwat = fh1.variables['pwat'][:,20:67,WAIS_wed]
	time2=fh1.variables['time'][:]
	units=fh1.variables['time'].units
	
#Extract the date information
	try :
		cal = fh1.variables['time'].calendar
	except AttributeError : 
	# Attribute doesnt exist
		cal = 'standard' 
	datevar=(nc.num2date(time2,units=units,calendar=cal))
	fh1.close()

#Begin the algorithm

	pwat_per  = np.percentile(pwat, 98, axis=0) #Calc the 98th percentile of pwat at all points
	
	#Create indicies where the pwat exceeds the 98th percentile
	indices = np.where(pwat > pwat_per) #Create indicies where the pwat exceeds the 98th percentile
	time = indices[0] #The time indices
	y = lat[indices[1]] #The latitude numbers
	x = lon[indices[2]] #The longitude numbers
	x2 = indices[2] #The longitude indices

	timesteps = unique(time) # unique timesteps, avoids if(timestep != prev_timestep) loop

	river_idx = [] #This list will store all AR detections
	landfall_idx = [] #This list will store AR landfalls only 

	for timestep in timesteps : #Cycle through the timesteps in each monthly file
		timestep_idx = np.where(timestep == time)[0] #organize the indices by time
		y_idx = y[timestep_idx]
		x_idx = x[timestep_idx]
		x2_idx = x2[timestep_idx]

		if not (y_idx[:-1] - y_idx[1:] > 1).any(): # checks the difference between successor/predecessor latitudes at the array level (no loop)
			if y_idx.max() - y_idx.min() > 20 : #The difference between the minimum and maximum latitude
				river_idx.append(timestep) #Attach the timestep to the all AR list
				i = np.argmin(y_idx)
				t = np.argmax(y_idx)
	
				if lat_mask[coast_idx[x2[i]]] >= y_idx[i]: #If the minimum latitude is at or below the latitude of the coastline
					landfall_idx.append(timestep) #Then the AR is considered a landfall
					
	datevar_list.append(datevar[landfall_idx])


######################################## Creation fichiers textes
#Create two text files: One for the AR time indices for all rivers and one for ARs that make landfall

	with open('/home/willej/Documents/data/ERAI/'+mon[b]+'_79_17_WAIS_wed_all_idx_v1.dat', 'w') as file:
		file.write('\n'.join(str(idx) for idx in river_idx))
	with open('/home/willej/Documents/data/ERAI/'+mon[b]+'_79_17_WAIS_wed_landfall_idx_v1.dat', 'w') as file:
		file.write('\n'.join(str(idx) for idx in landfall_idx))

	print(mon[b]) #output the month to check if the file matches the correct month
	b = b + 1

#Output the list of AR landfall dates

datevar_list = np.asarray(datevar_list)
with open('/home/willej/data/ERAI/ar_dates/wed_ar_dates.dat', 'w') as file:
	file.write('\n'.join(str(idx) for idx in datevar_list))

print(datevar_list)






