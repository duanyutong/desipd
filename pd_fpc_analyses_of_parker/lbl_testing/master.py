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
import numpy as np
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
from scipy import stats
#import sys
import glob
# from photutils.morphology import centroid_com, centroid_1dg, centroid_2dg

#%% read in fits file

# test files
# directory
# parentdir = r'K:\Google Drive\BUPHY\AS203 Principles of Astronomy II\night lab\zeta gem\as203-2016_spring-images'
parentDir = os.path.join(os.path.expanduser("~"), r'Google Drive\BUPHY\AS203 Principles of Astronomy II\night lab\zeta gem\as203-2016_spring-images')
# use os.sep after drive letter
datasetDirs = glob.glob(os.path.join(parentDir, '2016*'))
for datasetDir in datasetDirs:
    # os.chdir(dir)
    # os.listdir(dir)
    
    # filename match patterns
    
    print('Now processing dataset: ' + datasetDir)
    path_pattern_bias   = os.path.join(datasetDir, '*bias*.fit*')
    path_pattern_dark   = os.path.join(datasetDir, '*dark*.fit*')
    path_pattern_flat   = os.path.join(datasetDir, '*flat*.fit*')
    path_pattern_object = os.path.join(datasetDir, '*object*.fit*')
    path_pattern_light  = os.path.join(datasetDir, '*light*.fit*')
    
    # perform search and create full pathnames
    paths_bias   = glob.glob(path_pattern_bias)
    paths_dark   = glob.glob(path_pattern_dark)
    paths_flat   = glob.glob(path_pattern_flat)
    paths_object = glob.glob(path_pattern_object)
    if len(paths_object) == 0:
        paths_object  = glob.glob(path_pattern_light)
    
    # in order to preallocate arrays
    hdulist = fits.open(paths_bias[0])
    shape = hdulist[0].data.shape 	# get image dimension
    datatype = np.float32
    # dtype = hdulist[0].data.dtype 	# get data type # always float16
    # get numbers of files for each type
    num_bias   = len(paths_bias)
    num_dark   = len(paths_dark)
    num_flat   = len(paths_flat)
    num_object = len(paths_object)
    
    # preallocate arrays
    bias3D   = np.empty([shape[0], shape[1], num_bias],   dtype = datatype)
    dark3D   = np.empty([shape[0], shape[1], num_dark],   dtype = datatype)
    flat3D   = np.empty([shape[0], shape[1], num_flat],   dtype = datatype)
    object3D = np.empty([shape[0], shape[1], num_object], dtype = datatype)
    # fill arrays with fits data
    for i in range(len(paths_bias)):
    	bias3D[:,:,i]   = fits.open(paths_bias[i])  [0].data
    for i in range(len(paths_dark)):
    	dark3D[:,:,i]   = fits.open(paths_dark[i])  [0].data
    for i in range(len(paths_flat)):
    	flat3D[:,:,i]   = fits.open(paths_flat[i])  [0].data
    for i in range(len(paths_object)):
    	object3D[:,:,i] = fits.open(paths_object[i])[0].data
    
    #%% image reduction (from object, bias, dark, flat) - for FPC only
    
    # master bias from  median combining
    m_bias = np.median(bias3D, axis=2)
    
    # master dark from  median combining, -m_bias
    m_dark = np.median(dark3D, axis=2) - m_bias
    
    # master flat
    for i in range(len(paths_flat)):    # rescale each flat frame by the mode
    	flat3D[:,:,i] = flat3D[:,:,i] - m_bias - m_dark   # remove m_bias, m_dark
    	[mode, _] = stats.mode(flat3D[:,:,i], axis=None)
    	flat3D[:,:,i] = flat3D[:,:,i]/mode
    m_flat = np.median(flat3D, axis=2)   # median combine all flat frames
    m_flat_normalised = m_flat / np.mean(m_flat)   # normalised by the mean
    
    # master object from median combining, -m_bias, -m_dark
    m_object = np.median(object3D, axis=2) - m_bias - m_dark
    
    # divide out normalised flat to obtain final master image
    master = m_object / m_flat_normalised
    
    #debugging
    #maxindex = np.argmax(m_bias)
    #maxindex = np.unravel_index(maxindex, (1020,1530))
    #print(bias3D[maxindex[0],maxindex[1],:])
    #print(dark3D[maxindex[0],maxindex[1],:])
    #print(flat3D[maxindex[0],maxindex[1],:])
    #print(object3D[maxindex[0],maxindex[1],:])
    
    # save master images
    saveDir = os.path.join(datasetDir, 'master')
    # save master images
    if not os.path.exists(saveDir):	# check and create master folder
        os.makedirs(saveDir)
    # enter master folder
    # os.chdir(savedir)
    # create hdu objects
    hdu_m_bias            = fits.PrimaryHDU(m_bias)
    hdu_m_dark            = fits.PrimaryHDU(m_dark)
    hdu_m_flat            = fits.PrimaryHDU(m_flat)
    hdu_m_flat_normalised = fits.PrimaryHDU(m_flat_normalised)
    hdu_m_object          = fits.PrimaryHDU(m_object)
    hdu_master            = fits.PrimaryHDU(master)
    # write hdu to fits
    hdu_m_bias.writeto(
        os.path.join(saveDir, 'master_bias.fits'),            clobber = True)
    hdu_m_dark.writeto(
        os.path.join(saveDir, 'master_dark.fits'),            clobber = True)
    hdu_m_flat.writeto(
        os.path.join(saveDir, 'master_flat.fits'),            clobber = True)
    hdu_m_flat_normalised.writeto(
        os.path.join(saveDir, 'master_flat_normalised.fits'), clobber = True)
    hdu_m_object.writeto(
        os.path.join(saveDir, 'master_object.fits'),          clobber = True)
    hdu_master.writeto(
        os.path.join(saveDir, 'master.fits'),                 clobber = True)
        
    print('Finished processing dataset: ' + datasetDir)
