# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 19:29:25 2016

@author: givoltage

functions for FPC/positioner dithering script

"""

import numpy as np
import time
import Pyro4
# import sys
import os
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


#from DOSlib.application import Application
#from DOSlib.discovery import discoverable
#from DOSlib.util import dos_parser
from DOSlib.advertise import Seeker

''' sample values'''
datasetID = '20160825_01'
dirLocal = r'/data/images/fpc/'
expno = 9
exptime = 0.1

def expose_fpc_dataset(datasetID, dirLocal, expno, exptime):
    
    '''
    it takes a dataset to create each master image. the dataset is taken when
    nothing is moving.
     
    each dataset is saved in dirLocal/datasetID/. datasetID should increase by
    1 for each pos movement.
    
    expno is the number of exposures in each sequence
    
    exptime is used for dark and object frames
    '''

    datasetDirLocal = os.path.join(dirLocal, datasetID)
    expidini = 1 # initial value for exposure id
    
    #role1 = 'ILLUMINATOR'
    role2 = 'FPC'
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
        if role2 in s.devices:
            fpc = Pyro4.Proxy(s.devices[role2]['pyro_uri'])
            print('FPC connected')
        time.sleep(1)
        
    fpc.set(image_dir = datasetDirLocal)
    
    print('Taking dataset to be saved in: {}'.format(datasetDirLocal))
    expids = np.arange(expidini, expidini + expno, 1, dtype=int)
    for expid in expids:
        print('Exposing bias frame {}'.format(expid))
        fpc.expose(expid, exptype="bias",exptime=0,local_copy=False)
        write_fits(fpc.get_data(expid), expid)
    expids = expids + expno
    print(expids)
    for expid in expids:
        print('Exposing dark frame {}'.format(expid))
        fpc.expose(expid, exptype="dark",exptime=exptime,local_copy=False)
        write_fits(fpc.get_data(expid), expid)
    expids = expids + expno
    for expid in expids:
        print('Exposing object frame {}'.format(expid))
        fpc.expose(expid, exptype="object",exptime=exptime,local_copy=False)
        write_fits(fpc.get_data(expid), expid)
    
def master_from_header(datasetID, dirLocal):
    
    '''
    data reduction function gives a single reduced master image from a dataset
    '''
    flatEnabled = False
    masterSaveDir = os.path.join(dirLocal, datasetID)
    
    #%% read in fits file
    
    # get all fits files
    path_pattern = os.path.join(masterSaveDir, '*.fit*')
    paths = glob.glob(path_pattern)
    obstypes = []
    for path in paths:
        hdulist = fits.open(path)
        obstype = hdulist[0].header['OBSTYPE']
        obstypes.append(obstype)
    frames_dict = {'path': paths}
    frames_dict['filename'] = [os.path.basename(path) for path in paths]
    frames_dict['obstype'] = obstypes
    
    # perform search and create full pathnames based on header keyword
    paths_bias   = [frames_dict['path'][i] for i in range(len(paths))
                        if frames_dict['obstype'][i] == (
                                                         'ZERO'
                                                         or 'zero'
                                                         or 'BIAS'
                                                         or 'bias'
                                                         )
                    ]
    paths_dark   = [frames_dict['path'][i] for i in range(len(paths)) 
                        if frames_dict['obstype'][i] == (
                                                         'DARK'
                                                         or 'dark'
                                                         )
                    ]
    paths_flat   = [frames_dict['path'][i] for i in range(len(paths)) 
                        if frames_dict['obstype'][i] == (
                                                         'FLAT' 
                                                         or 'flat'        
                                                         or 'DOME FLAT'
                                                         or 'dome flat'
                                                         or 'SKY FLAT'
                                                         or 'sky flat'
                                                         )
                     ]
    paths_object = [frames_dict['path'][i] for i in range(len(paths)) 
                        if frames_dict['obstype'][i] == (
                                                         'OBJECT'
                                                         or 'object'
                                                         or 'LIGHT'
                                                         or 'light'
                                                         )
                    ]
    
    # in order to preallocate arrays
    hdulist = fits.open(paths_bias[0])
    shape = hdulist[0].data.shape 	# get image dimension
    # dtype = hdulist[0].data.dtype 	# get data type # always float16
    datatype = np.float32   # let's be more precise
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
    
    # master bias from median combining
    m_bias = np.median(bias3D, axis=2)
    
    # master dark from median combining, -m_bias
    m_dark = np.median(dark3D, axis=2) - m_bias
    
    # master object from median combining, -m_bias, -m_dark
    m_object = np.median(object3D, axis=2) - m_bias - m_dark
    
    # master flat
    if flatEnabled and num_flat >= 3: 
        for i in range(len(paths_flat)):    # rescale each flat frame by the mode
            flat3D[:,:,i] = fits.open(paths_flat[i]) [0].data - m_bias - m_dark
            [mode, _] = stats.mode(flat3D[:,:,i], axis=None)
            flat3D[:,:,i] = flat3D[:,:,i]/mode
        m_flat = np.median(flat3D, axis=2)   # median combine all flat frames
        m_flat_normalised = m_flat / np.mean(m_flat)   # normalised by the mean
        # divide out normalised flat to obtain final master image
        master = m_object / m_flat_normalised    
    else:
        master = m_object
        
    #debugging
    #maxindex = np.argmax(m_bias)
    #maxindex = np.unravel_index(maxindex, (1020,1530))
    #print(bias3D[maxindex[0],maxindex[1],:])
    #print(dark3D[maxindex[0],maxindex[1],:])
    #print(flat3D[maxindex[0],maxindex[1],:])
    #print(object3D[maxindex[0],maxindex[1],:])
    
    # save master images
    if not os.path.exists(masterSaveDir):	# check and create master folder
        os.makedirs(masterSaveDir)
    
    # os.chdir(savedir)    # enter master folder
    # create hdu objects
    hdu_m_bias   = fits.PrimaryHDU(m_bias)
    hdu_m_dark   = fits.PrimaryHDU(m_dark)
    hdu_m_object = fits.PrimaryHDU(m_object)
    hdu_master   = fits.PrimaryHDU(master)
    # write exposure time for flux/magnitude calculation
    expreq_bias = fits.open(paths_bias[0])[0].header['EXPREQ']
    expreq      = fits.open(paths_object[0])[0].header['EXPREQ']
    print('EXPREQ for bias is: {}'.format(expreq_bias))
    hdu_m_bias.header  ['EXPREQ'] = expreq_bias
    hdu_m_dark.header  ['EXPREQ'] = expreq
    hdu_m_object.header['EXPREQ'] = expreq
    hdu_master.header  ['EXPREQ'] = expreq
    if 'EXPTIME' in fits.open(paths_object[0])[0].header:
        exptime_bias = fits.open(paths_bias[0])[0].header['EXPTIME']
        exptime      = fits.open(paths_object[0])[0].header['EXPTIME']
        hdu_m_bias.header  ['EXPTIME'] = exptime_bias
        hdu_m_dark.header  ['EXPTIME'] = exptime
        hdu_m_object.header['EXPTIME'] = exptime
        hdu_master.header  ['EXPTIME'] = exptime
    
    # write hdu to fits
    hdu_m_bias.writeto(
        os.path.join(masterSaveDir, 'master_bias.fits'),   clobber = True)
    hdu_m_dark.writeto(
        os.path.join(masterSaveDir, 'master_dark.fits'),   clobber = True)
    hdu_m_object.writeto(
        os.path.join(masterSaveDir, 'master_object.fits'), clobber = True)
    hdu_master.writeto(
        os.path.join(masterSaveDir, 'master.fits'),        clobber = True)
    
    if flatEnabled and num_flat >= 3: 
        hdu_m_flat            = fits.PrimaryHDU(m_flat)
        hdu_m_flat_normalised = fits.PrimaryHDU(m_flat_normalised)
        hdu_m_flat.writeto(
            os.path.join(masterSaveDir, 'master_flat.fits'),            clobber = True)
        hdu_m_flat_normalised.writeto(
            os.path.join(masterSaveDir, 'master_flat_normalised.fits'), clobber = True)
