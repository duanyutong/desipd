# -*- coding: utf-8 -*-
"""
Taken from FPC code using astropy.photutils

Created on Mon May 16 14:41:22 2016

@author: Duan Yutong (dyt@lbl.gov)
"""

# TODO: pipeline for fits files to be sent back to us for analysis
# 	- to do at Kitt Peak, should be easy
# TODO: calibration for conversion from instrumental flux and magnitude ->
#	apparent magnitude
#	- need known reference
# TODO: add support for dithering pattern of telescope or positioners in order
#	to find position which maximises intensity
#
# 	- fit with a 2d elliptical gaussian kernal?
# 	- in what way, if any, do we want to visualise this dataset and results?
# TODO: Write position information to fits header
# TODO: check if FVC code is useful for FPC
#	- pending repository access

# TODO: [x] image reduction of object, bias, dark, flat
#	- these raw images are taken by camera controller, part of ICS
#	new SBIG camera arrived, accessories to be ordered
# TODO: [x] camera controller and interface definition. ICS/TCS integration?
#	- done and is part of ICS; FPC analysis is separate from ICS
#	for now just read in and analyse files from specified directory
# TODO: [x] SBIG camera crash problem when saturated?
#	- not a concern for the new camera
# TODO: [x] address negative intensity values after background subtraction
#	- not a concern since background subtraction brings median/mean to zero
#	and will always have a spread around 0 which necessarily goes negative
#	either manually set to zero or do nothing (built-in to astropy/photutils)

#%% import modules

import os
# import numpy as np
#import matplotlib.pyplot as plt
# from skimage import measure
# import cv2
from astropy.io import fits
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
import glob
# from photutils.morphology import centroid_com, centroid_1dg, centroid_2dg

#%% read in fits file

# test files
# directory
# parentdir = r'K:\Google Drive\BUPHY\AS203 Principles of Astronomy II\night lab\zeta gem\as203-2016_spring-images'
parentDir = '/Volumes/data/images/fpc/'
# use os.sep after drive letter
datasetID = 'flats'
datasetDir = os.path.join(parentDir, datasetID)

# filename match patterns
path_pattern_bias   = os.path.join(datasetDir, '*bias*.fit*')
path_pattern_dark   = os.path.join(datasetDir, '*dark*.fit*')
path_pattern_flat   = os.path.join(datasetDir, '*flat*.fit*')
path_pattern_object = os.path.join(datasetDir, '*object*.fit*')

# perform search and create full pathnames
paths_bias   = glob.glob(path_pattern_bias)
paths_dark   = glob.glob(path_pattern_dark)
paths_flat   = glob.glob(path_pattern_flat)
paths_object = glob.glob(path_pattern_object)

# dtype = hdulist[0].data.dtype 	# get data type # always float16
# get numbers of files for each type
for i in range(len(paths_bias)):
    hdulist = fits.open(paths_bias[i], mode='update')
    hdulist[0].header['OBSTYPE'] = 'zero'
    hdulist.close()
for i in range(len(paths_dark)):
    hdulist = fits.open(paths_dark[i], mode='update')
    hdulist[0].header['OBSTYPE'] = 'dark'
    hdulist.close()
for i in range(len(paths_flat)):
    hdulist = fits.open(paths_flat[i], mode='update')
    hdulist[0].header['OBSTYPE'] = 'flat'
    hdulist.close()
for i in range(len(paths_object)):
    hdulist = fits.open(paths_object[i], mode='update')
    hdulist[0].header['OBSTYPE'] = 'object'
    hdulist.close()
