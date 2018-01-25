# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 18:20:47 2016

@author: Duan Yutong (dyt@physics.bu.edu, dyt@lbl.gov)

Script to take FPC images via DOS

"""

import time
import Pyro4
import os
#from DOSlib.application import Application
#from DOSlib.discovery import discoverable
#from DOSlib.util import dos_parser
from DOSlib.advertise import Seeker

id_dataset = '20160923'
image_dir = os.path.join('/data/images/fpc', id_dataset)
expid = 0
expno = 100
exptime = 0.001
local_copy = True
role1 = 'FPC'
fpc = None

#Identify Applications
s = Seeker('-dos-','DOStest')
#while led == None:
#    print('Connecting to ILLUMINATOR')
#    s.seek()
#    if role1 in s.devices:
#        led = Pyro4.Proxy(s.devices[role1]['pyro_uri'])
#        print('ILLUMINATOR connected')
#    time.sleep(1)

while fpc == None:
    print('Connecting to FPC...')
    s.seek()
    if role1 in s.devices:
        fpc = Pyro4.Proxy(s.devices[role1]['pyro_uri'])
        print('FPC connected')
    time.sleep(1)

fpc.set(image_dir = image_dir)

print('Dataset to be saved in: {}'.format(image_dir))
for i in range(expno):

#    expid = expid + 1
#    print('Exposing bias frame {} of {}'.format(expid, expno*3))
#    fpc.expose(expid, exptype="bias",exptime=0,local_copy=local_copy)
#    expid = expid + 1
#    print('Exposing dark frame {} of {}'.format(expid, expno*3))
#    fpc.expose(expid, exptype="dark",exptime=exptime,local_copy=local_copy)
    expid = expid + 1
    print('Exposing light frame {} of {}'.format(expid, expno*3))
    fpc.expose(expid, exptype="object",exptime=exptime,local_copy=local_copy)



