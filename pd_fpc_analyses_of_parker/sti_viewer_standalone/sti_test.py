"""
scatt_light.py
May 2016
P. Fagrelius


Script for running sbigcam_sti.py for measuring the scattered light from the Mayall telescope

"""

import numpy as np
#import matplotlib.pyplot as plt
import sys, os, time
from sbigcam_sti import SBIGCam

from astropy.time import Time
from astropy.coordinates import (EarthLocation,SkyCoord, get_sun, AltAz, Angle)
import astropy.units as u
from astropy.utils.data import download_file
from astropy.coordinates.angle_utilities import angular_separation
import ephem
from datetime import datetime


#####INPUTS######

data_dir = '/home/msdos/SBIG/'
exptime = 1000 #ms
rest = 5 #sec 
DARK = False
num = 1 #'cont'= continuous
name = 'sti_test_' #filename identified, will include time
shutt = True  #True shutter closes, False shutter stays open 
window_mode=False
################


def get_moon(time):
	location = EarthLocation(lat=(31+57/60+50.28/3600),lon=(-111+35./60+59.41/3600.),height=2098) #kitt_peak
	time = Time(time)

	moon = ephem.Moon()
	obs = ephem.Observer()
	obs.lat = location.latitude.to_string(u.deg, sep=':')
	obs.lon = location.longitude.to_string(u.deg, sep=':')
	obs.elevation = location.height.to(u.m).value
	obs.date = time.datetime

	moon.compute(obs)
	m_RA = float(moon.ra)
	m_DEC = float(moon.dec)
	return(m_RA, m_DEC)


##### On-Mountain Script ########
camera = SBIGCam()
camera.verbose = True

print('The images will be saved in ',dir)

if not camera.open_camera():
	print("Can't establish conneciton to camera")
	sys.exit()

#print('The exposure time is for each exposure is ',exptime,' ms')
camera.set_exposure_time(exptime)

if DARK:
	print("This will be a dark image. If that's not what you intended, change teh value in scatt_light.py")
	camera.set_dark(True)
else:
	camera.set_dark(False)

if num=='lin_test':
	EXP = [100,200,300,400,500,1000,1500]
	for i in EXP:
		camera.set_exposure_time(i)
		take_image = camera.start_exposure(shutter=shutt)
		image, start = take_image
		filename = (name+"%04d"+'_'+datetime.now().strftime("%H_%M_%S.%f")+'.fits') % (i,) 
		camera.write_fits(image,data_dir+'/'+filename,i,start,get_moon(start))
	camera.close_camera()
		
elif num!='cont':
	N = int(num)
	while N>0:
		print("%02d" % (N,))
		if window_mode:
			try:
				camera.set_window_mode(top=100,left=100) 
				print('Window mode set')
			except:
				print('Window mode error') 
		take_image = camera.start_exposure(shutter=shutt)
		image, start = take_image
		#print(start)
		filename = (name+"%02d"+'_'+datetime.now().strftime("%H_%M_%S.%f")+'.fits') % (N,)
		camera.write_fits(image,data_dir+'/'+filename,exptime,start,get_moon(start))
		N-=1
	camera.close_camera()
elif num == 'cont':
	print('To stop continous exposures, press Ctrl-C')
	try:
		while True:
			print(datetime.now().strftime("%H_%M_%S.%f"))
			take_image = camera.start_exposure(shutter=shutt)
			image, start = take_image
			filename = name+datetime.now().strftime("%H_%M_%S.%f")+'.fits'
			camera.write_fits(image,data_dir+'/'+filename,exptime,start,get_moon(start))
			time.sleep(rest)
	except KeyboardInterrupt:
		camera.close_camera()
            
        
        




