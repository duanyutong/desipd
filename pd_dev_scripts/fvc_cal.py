import numpy as np
import time
import Pyro4
import sys,os

from DOSlib.application import Application
from DOSlib.discovery import discoverable
from DOSlib.util import dos_parser
from DOSlib.advertise import Seeker

"""
Script to turn off lights, take calibration image, turn lights back on, and then take calibrate_image and measure to make sure everything is working

"""
#Fiducials that need to be turned off
bad_fids = [] #[1111,1113,1103]
bad_zero = [] #[0,0,0]
#Fiducials that need to be brighter than usualZZ
dim_fid = 1114
etime =2.0
led_current = 20.
fid_dc = 100. 
data_dir = '/data/fvc/'
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

fvc.set(exptime=etime)
#Calibration Sequence
#led.set(channel=2)
#led.set(led='off')
#led.set(channel=3)
#led.set(led='off')
#led.set(channel=4)
#led.set(led='off')
led.set(channel=1)
led.set(led='off')
pc0.set_fiducial(20000,0)
print('Everything is turned off')
error = fvc.calibrate_bias(0)
if error == 'SUCCESS':
    print('Background image taken successfully')
else:
    print('Bias image was not taken')


led.set(channel=1)
led.set(current=20)
led.set(led='on')
pc0.set_fiducial(20000,50)
pc0.set_fiducial(dim_fid,100)
pc0.set_fiducials(bad_fids,bad_zero)
print('Fiducials and fibers turned back on')
#set_target = fvc.set(target_file=t_file)
#print("Setting target file: ", set_target)
error = fvc.calibrate_image()
if error == 'SUCCESS':
    print('Autotune successful')
    measure = fvc.measure()
    if type(measure) == dict:
        print("Calibration was successful")
    else:
        print(measure)
        print('Measure was not successful')
else:
    print('Autotune not successful')



