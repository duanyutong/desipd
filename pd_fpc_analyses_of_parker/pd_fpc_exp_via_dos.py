# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 18:20:47 2016
Script to take FPC images via DOS
@author: givoltage
"""

# import numpy as np
import time
import Pyro4
# import sys
import os
#import matplotlib.pyplot as plt
# from skimage import measure
# import cv2
# from astropy.io import fits
# import astropy.units as units
# from astropy.stats import mad_std
# from astropy.stats import sigma_clipped_stats
#from astropy.visualization import LogStretch, SqrtStretch
#from astropy.visualization.mpl_normalize import ImageNormalize
# from skimage import measure
#from photutils.background import Background
#from photutils import daofind
# from photutils import CircularAperture
#from photutils import detect_sources
#from photutils.utils import random_cmap
#from photutils import source_properties, properties_table
#from photutils import EllipticalAperture
# from scipy import stats
#import sys
# import glob
# from photutils.morphology import centroid_com, centroid_1dg, centroid_2dg


#from DOSlib.application import Application
#from DOSlib.discovery import discoverable
#from DOSlib.util import dos_parser
from DOSlib.advertise import Seeker

#path_dataset = '/data/images/fpc/linear/expreq_0.3'
#exposureEnabled = True
#reductionEnabled = True
#dirLocal = r'/data/images/'
# flatEnabled = False
#datasetDirLocal = os.path.join(dirLocal, datasetID)
#masterSaveDir = os.path.join(datasetDirLocal, '')

id_dataset = '20160827'
image_dir = os.path.join('/data/images/fpc', id_dataset)
expid = 0
expno = 1000
exptime = 0.001
role1 = 'FPC2'
local_copy = True
#role1 = 'ILLUMINATOR'
# role2 = 'FPC'
#led = None
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
    print('Connecting to FPC')
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



