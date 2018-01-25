#/bin/python

import os
import sys
import numpy as np
import time
import Pyro4
import csv

from DOSlib.application import Application
from DOSlib.discovery import discoverable
from DOSlib.util import dos_parser
from DOSlib.advertise import Seeker

exptime = float(sys.argv[1])
num_exp = 10

levels = {0.5:{'dc':95,'current':55},1.0:{'dc':45,'current':35},2.0:{'dc':25,'current':15},5.0:{'dc':10,'current':2},10.:{'dc':9,'current':1}}
led_current = float(levels[exptime]['current'])
fid_dc = float(levels[exptime]['dc'])

#Setup some new folders
base_dir = '/data/fvc/'
csv_name = str(exptime)+'_'+time.strftime('%d%m%y_%H:%M')
full_dir = base_dir + csv_name
csv_dir = base_dir + 'csv_files/'+csv_name
print("Files will be saved on desifvc at /data/fvc/. The list of files created during this test will be saved on this machine in this file: %s.csv"%(csv_dir))

role1 = 'ILLUMINATOR'
role2 = 'FVC'
role3 = 'PC0'
led = None
fvc = None
pc0 = None

#Identify Applications
s = Seeker('-dos-','DOStest')
while led == None:
    s.seek()
    if role1 in s.devices:
        led = Pyro4.Proxy(s.devices[role1]['pyro_uri'])
        print('ILLUMINATOR connected')
    else:
        print('Not connecting to the ILLUMINATOR application')
    time.sleep(1)
    
while fvc == None:
    s.seek()
    if role2 in s.devices:
        fvc = Pyro4.Proxy(s.devices[role2]['pyro_uri'])
        print('FVC connected')
    else:
        print('Not connecting to the FVC application')
    time.sleep(1)

ss = Seeker('-dos-','PetalControl')
while pc0 == None:
    ss.seek()
    if role3 in ss.devices:
        pc0 = Pyro4.Proxy(ss.devices[role3]['pyro_uri'])
        print('PC0 connected')
    else:
        print('Not connecting to the PC0 application')
    time.sleep(1)



fvc.set(exptime=exptime)
fvc.set(image_dir=full_dir)

#Prepare for measurements
led.set(channel=1)
led.set(led='off')
pc0.set_fiducial(20000,0)
print('Everything is turned off')
error = fvc.calibrate_bias(0)
if error == 'SUCCESS':
    print('Background image taken successfully')
else:
    print('Bias image was not taken')

filenames = []
 
led.set(channel=1)
led.set(current=led_current) 
led.set(led='on')
pc0.set_fiducial(20000,fid_dc)
pc0.set_fiducial(1114,fid_dc*2.)
print('Fiducials and fibers turned back on') 
error = fvc.calibrate_image()
if error == 'SUCCESS':
    print('Autotune successful')
    for i in range(num_exp):
        measure = fvc.take_exposure()
        if type(measure) == 'SUCCESS':
            print("Measure  was successful")
            fname = fvc.get('image_name')
            filenames.append(fname)
        else:
            print(measure)
else:
   print("Autotune was not successful")
f=open(csv_dir+'.csv','w')
w=csv.writer(f)
for f in filenames:
    w.writerow([f])
    n = os.path.splitext(f)[0]
    w.writerow([n+'.red'])
    w.writerow([n+'.pos'])
    w.writerow([n+'.cat'])
