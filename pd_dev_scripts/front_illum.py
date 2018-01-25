import numpy as np
import time
import Pyro4
import sys,os

from DOSlib.application import Application
from DOSlib.discovery import discoverable
from DOSlib.util import dos_parser
from DOSlib.advertise import Seeker


"""Script to turn on all front illuminating LEDs

Input arguments: led_mode fvc_image
led_mode = on, off
fvc_image = take_image, None
"""
#Script instructions
etime = 1.0
led_current = 20
data_dir = '/data/fvc/'
fvc_machine = 'msdos@desifvc.kpno.noao.edu'
fvc_name='frontillum.fits'
fvc_address = fvc_machine+':'+data_dir+fvc_name


#Read in Arguments
led_mode = sys.argv[1]
fvc_image = sys.argv[2]

# initialization done when starting petal application
petal_id = 0
role1 = 'ILLUMINATOR'
role2 = 'FVC'
role3 = 'PC0'
led = None
fvc = None
pc0 = None

# Find Illuminator and FVC applications
s = Seeker('-dos-','DOStest')
while led == None:
    s.seek()
    if role1 in s.devices:
        led = Pyro4.Proxy(s.devices[role1]['pyro_uri'])
        print('ILLUMINATOR connected')
    time.sleep(1)
    
while fvc == None:
    s.seek()
    if role2 in s.devices:
        fvc = Pyro4.Proxy(s.devices[role2]['pyro_uri'])
        print('FVC connected')
    else:
        print('Not connecting to anything using Pyro')
    time.sleep(1)

ss = Seeker('-dos-','PetalControl')
while pc0 == None:
    ss.seek()
    if role3 in ss.devices:
        pc0 = Pyro4.Proxy(ss.devices[role3]['pyro_uri'])
        print('PC0 connected')
    else:
        print('PCO not connecting')
    
#Turn on/off LEDs
if led_mode.lower() == 'on':
    print("Turning on LEDs with %d current" % led_current)
    led.set(channel=2)
    led.set(current=led_current)
    led.set(led='on')
    led.set(channel=3)
    led.set(current=led_current)
    led.set(led='on')
    led.set(channel=4)
    led.set(current=led_current)
    led.set(led='on')
elif led_mode.lower() == 'off':
    print("Turning off Front Illuminating LEDs")
    led.set(channel=2)
    led.set(led='off')
    led.set(channel=3)
    led.set(led='off')
    led.set(channel=4)
    led.set(led='off')


# Take image with FVC
if fvc_image.lower() == 'take_image':
    fvc.set(exptime=float(etime))
    fvc.set(image_name=fvc_name)
    print("Taking %d second FVC exposure" % etime)
    fvc.take_exposure()
else: 
	pass

    
def sh(script):
    os.system("bash -c '%s'" % script)

sh("scp %s ." % fvc_address)
sh("ds9 frontillum.fits")
sh("rm frontillum.fits")




