# -*- coding: utf-8 -*-
""" FPC code using astropy.photutils
    Created on Mon May 16 14:41:22 2016
    @author: Duan Yutong (dyt@lbl.gov)
    
    ### bias
    
    Artemis Capture, the software that comes with the SBIG camera, does
    not support taking bias images. Other image capture software, e.g. MaximDL,
    may support it. However, bias is not an aboslute must for this script to run.
    
    Program will automatically detect filenames containing "bias".
    
    ### Dark
    
    When taking one or a series of exposures, Artemis Capture provides the
    option to take a "dark" frame and subtract it from the object frame. This
    subtraction takes care of both dark current and bias offset. So, if using
    Artemis Capture, it is recommended that this subtraction method always be
    used for most accurate correction and convenience, since dark is very
    temperature-senstiive.
    
    Filenames should contain "subtracted" if they are already dark-subtracted.
    Script will detect this keyword in the filename.
    
    ### Flat
    
    Flat is not required but highly recommended to correct for pixel gain
    variations.
    
    Script will detect existence of flats.
    
    ### LEDoff
    
    If using other software, we could take all 4 frame types -
    bias, dark, flat, and object - follow the standard CCD reduction
    procedure, and calculate the background using sigma-clipping.
    
    Alternatively, we could also turn off the back-illuminated LED to get a
    background frame which is a quick cheat that only works for ProtoDESI and
    is otherwise unfeasible for actual stars. With this method, the only thing
    needed afterwards is flat-field correction (auto-detect, optional).
    
    ### Decision logic // all flags start with "ccd"
    
    *   if ledoff images exist, check preference of data reduction method
    *   if ledoff is preferred, subtracted ledoff background, and optionally
        divide by flats if they exist
    *   if ledoff DNE, or is not preferred, check existence of bias, flat, and
        "subtracted" images.
        
        if all images are dark-subtracted, no dark is used
        if either object or flat frames are not dark-subtracted, script will
        take all dark frames to make a master dark, assuming temperature is
        constant for all types of frames, and apply the master dark
    
    Any type of frame has to have at least 3 images to be considered valid.
    
    ### Before running script
    
    Specify the file naming scheme in the customisation section, and path
    structure.

"""

# TODO: event display for a single location (this script treats single location)
# TODO: calibration for conversion from instrumental flux and magnitude ->
#   apparent magnitude
#   - need known reference
# TODO: add support for dithering pattern of telescope or positioners in order
#   to find position which maximises intensity
#   - fit with a 2d elliptical gaussian kernal?
#   - in what way, if any, do we want to visualise this dataset and results?
# TODO: Write position information/coordinates to fits header
# TODO: check if FVC code is useful for FPC
# TODO: expand the section where standard reduction is preferred
# TODO: pipeline for fits files to be sent back to us for analysis
#   - to do at Kitt Peak, should be easy

# TODO: [x] image reduction of object, bias, dark, flat
#   - these raw images are taken by camera controller, part of ICS
#   new SBIG camera arrived, accessories to be ordered
# TODO: [x] camera controller and interface definition. ICS/TCS integration?
#   - done and is part of ICS; FPC analysis is separate from ICS
#   for now just read in and analyse files from specified directory
# TODO: [x] SBIG camera crash problem when saturated?
#   - not a concern for the new camera
# TODO: [x] address negative intensity values after background subtraction
#   - not a concern since background subtraction brings median/mean to zero
#   and will always have a spread around 0 which necessarily goes negative
#   either manually set to zero or do nothing (built-in to astropy/photutils)

#%% import modules

import os
import numpy as np
import matplotlib.pyplot as plt
import datetime
# from skimage import measure
# import cv2
from astropy.io import fits
# import astropy.units as units
# from astropy.stats import mad_std
# from astropy.stats import sigma_clipped_stats
from astropy.visualization import LogStretch, SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
# from skimage import measure
from photutils.background import Background
#from photutils import daofind
# from photutils import CircularAperture
from photutils import detect_sources
from photutils.utils import random_cmap
from photutils import source_properties, properties_table
from photutils import EllipticalAperture
from scipy import stats
from sklearn.metrics import mean_squared_error
import sys
import glob
# from photutils.morphology import centroid_com, centroid_1dg, centroid_2dg

#%% Customise settings before running

# x position
xpos = '0502'
# flags
ccdPreferLEDoff = True
ccdUseFlat      = True
datatype = np.float32
# filename patterns
name_pattern_bias              = '*bias*.fit*'
name_pattern_dark              = '*dark*.fit*'
name_pattern_flat_subtracted   = '*flat_subtracted*.fit*'
name_pattern_flat              = '*flat*.fit*'
name_pattern_object_subtracted = '*object_subtracted*.fit*'
name_pattern_object            = '*object*.fit*'
name_pattern_ledoff_subtracted = '*ledoff_subtracted*.fit*'
name_pattern_ledoff            = '*ledoff*.fit*'

# directory
parentDir = r'K:\Google Drive\DESI\protoDESI\images\fpc_data'   # directory
# parentDir = os.path.join(os.path.expanduser("~"), r'Google Drive\BUPHY\AS203 Principles of Astronomy II\night lab\zeta gem\as203-2016_spring-images')
# use os.sep after drive letter
datasetDir = os.path.join(parentDir)

#%% read in fits file and set flags

os.chdir(datasetDir)
os.listdir(datasetDir)

# filename match patterns
path_pattern_bias              = os.path.join(datasetDir, name_pattern_bias)
path_pattern_dark              = os.path.join(datasetDir, name_pattern_dark)
path_pattern_flat_subtracted   = os.path.join(datasetDir, name_pattern_flat_subtracted)
path_pattern_flat              = os.path.join(datasetDir, name_pattern_flat)
path_pattern_object_subtracted = os.path.join(datasetDir, name_pattern_object_subtracted)
path_pattern_object            = os.path.join(datasetDir, name_pattern_object)
path_pattern_ledoff_subtracted = os.path.join(datasetDir, name_pattern_ledoff_subtracted)
path_pattern_ledoff            = os.path.join(datasetDir, name_pattern_ledoff)

# perform search and create full pathnames
paths_bias              = glob.glob(path_pattern_bias)
paths_dark              = glob.glob(path_pattern_dark)
paths_flat_subtracted   = glob.glob(path_pattern_flat_subtracted)
paths_flat              = glob.glob(path_pattern_flat)
paths_object_subtracted = glob.glob(path_pattern_object_subtracted)
paths_object            = glob.glob(path_pattern_object)
paths_ledoff_subtracted = glob.glob(path_pattern_ledoff_subtracted)
paths_ledoff            = glob.glob(path_pattern_ledoff)

# set flags, all initial boolean value false
ccdBiasExist             = False
ccdDarkExist             = False
ccdFlatSubtractedExist   = False
ccdFlatExist             = False
ccdObjectSubtractedExist = False
ccdObjectExist           = False
ccdLEDoffSubtractedExist = False
ccdLEDoffExist           = False

# modify flags by filename keyword detection
if len(paths_bias) >=3:     # bias
    ccdBiasExist = True
if len(paths_dark) >= 3:    # dark
    ccdDarkExist = True
if len(paths_flat_subtracted) >=3:  # flat-subtracted
    ccdFlatSubtractedExist = True
elif len(paths_flat) >=3:
    ccdFlatExist = True
if len(paths_object_subtracted) >=3:    # flat-subtracted
    ccdObjectSubtractedExist = True
elif len(paths_object) >=3:
    ccdObjectExist = True
if len(paths_ledoff_subtracted) >=3:    # flat-subtracted
    ccdLEDoffSubtractedExist = True
elif len(paths_ledoff) >=3:
    ccdLEDoffExist = True

#%% CCD image reduction

#  if LEDoff background is preferred and LEDoff images do exist

if ccdPreferLEDoff and (ccdLEDoffSubtractedExist or ccdLEDoffExist):
    
    # if image already comes dark-subtracted
    if ccdLEDoffSubtractedExist:
        
        # in order to preallocate arrays, get some info
        hdulist = fits.open(paths_ledoff_subtracted[0])
        shape = hdulist[0].data.shape   # get image dimension
        # get numbers of files for each type
        num_ledoff_subtracted = len(paths_ledoff_subtracted)
        # preallocate arrays
        ledoff_subtracted3D = np.empty([shape[0], shape[1], num_ledoff_subtracted], dtype = datatype)
        # fill arrays with fits data
        for i in range(len(paths_ledoff_subtracted)):
            ledoff_subtracted3D[:,:,i] = fits.open(paths_ledoff_subtracted[i]) [0].data
        # master ledoff_subtracted from median combining
        master_ledoff_subtracted = np.median(ledoff_subtracted3D, axis=2)
    
    # if ledoff images have seperate dark
    elif ccdLEDoffExist and ccdDarkExist: # LEDoff and dark both exist
    
        # in order to preallocate arrays, get some info
        hdulist = fits.open(paths_ledoff[0])
        shape = hdulist[0].data.shape   # get image dimension
        # get numbers of files for each type
        num_dark  = len(paths_dark)
        num_ledoff = len(paths_ledoff)
        # preallocate arrays
        dark3D = np.empty([shape[0], shape[1], num_dark],  dtype = datatype)
        ledoff3D = np.empty([shape[0], shape[1], num_ledoff], dtype = datatype)
        # fill arrays with fits data
        for i in range(len(paths_dark)):
            dark3D[:,:,i]  = fits.open(paths_dark[i]) [0].data
        for i in range(len(paths_ledoff)):
            ledoff3D[:,:,i]  = fits.open(paths_ledoff[i]) [0].data
        master_dark = np.median(dark3D, axis=2)
        master_ledoff = np.median(ledoff3D, axis=2)
        master_ledoff_subtracted = master_ledoff - master_dark
        
    # now master_ledoff_subtracted has been created in any case
    # if object comes dark-subtracted
    if ccdObjectSubtractedExist:
        
        # in order to preallocate arrays, get some info
        hdulist = fits.open(paths_object_subtracted[0])
        shape = hdulist[0].data.shape   # get image dimension
        # get numbers of files for each type
        num_object_subtracted = len(paths_object_subtracted)
        # preallocate arrays
        object_subtracted3D = np.empty([shape[0], shape[1], num_object_subtracted], dtype = datatype)
        # fill arrays with fits data
        for i in range(len(paths_object_subtracted)):
            object_subtracted3D[:,:,i] = fits.open(paths_object_subtracted[i]) [0].data
        # master ledoff_subtracted from median combining
        master_object_subtracted = np.median(object_subtracted3D, axis=2)
        
    # if object has seperate dark
    elif ccdObjectExist and ccdDarkExist:
        
        # in order to preallocate arrays, get some info
        hdulist = fits.open(paths_object[0])
        shape = hdulist[0].data.shape   # get image dimension
        # get numbers of files for each type
        num_dark  = len(paths_dark)
        num_object = len(paths_object)
        # preallocate arrays
        dark3D = np.empty([shape[0], shape[1], num_dark],  dtype = datatype)
        object3D = np.empty([shape[0], shape[1], num_object], dtype = datatype)
        # fill arrays with fits data
        for i in range(len(paths_dark)):
            dark3D[:,:,i]  = fits.open(paths_dark[i]) [0].data
        for i in range(len(paths_object)):
            object3D[:,:,i]  = fits.open(paths_object[i]) [0].data
        master_dark = np.median(dark3D, axis=2)
        master_object = np.median(object3D, axis=2)
        master_object_subtracted = master_object - master_dark

    # now master_object_subtracted has been created in any case
    # finally master using LEDoff
    # subtracted means free of LEDoff background
    master_subtracted = master_object_subtracted - master_ledoff_subtracted
    
    #%% still LEDoff preferred - optional: flat-field correction
    if ccdUseFlat:
        
        # if flat comes dark-subtracted
        if ccdFlatSubtractedExist:
            
            # in order to preallocate arrays, get some info
            hdulist = fits.open(paths_flat_subtracted[0])
            shape = hdulist[0].data.shape   # get image dimension
            num_flat_subtracted = len(paths_flat_subtracted)
            # preallocate arrays
            flat_subtracted3D  = np.empty([shape[0], shape[1], num_flat_subtracted], dtype = datatype)
            for i in range(len(paths_flat_subtracted)):
                print('now processing i =', i)
                flat_subtracted3D[:,:,i]  = fits.open(paths_flat_subtracted[i]) [0].data
                # rescale each flat frame by the mode
                [mode, _] = stats.mode(flat_subtracted3D[:,:,i], axis=None)
                flat_subtracted3D[:,:,i] = flat_subtracted3D[:,:,i]/mode
            # master flat, median combine all flats
            master_flat_subtracted = np.median(flat_subtracted3D, axis=2)  
            # normalised by the mean
            master_flat_subtracted_normalised = master_flat_subtracted / np.mean(master_flat_subtracted)
        
        # if flat has seperate dark
        elif ccdFlatExist and ccdDarkExist:
            # in order to preallocate arrays, get some info
            hdulist = fits.open(paths_flat_subtracted[0])
            shape = hdulist[0].data.shape   # get image dimension
            # get numbers of files for each type
            num_dark  = len(paths_dark)
            num_flat = len(paths_flat)
            # preallocate arrays
            dark3D = np.empty([shape[0], shape[1], num_dark], dtype = datatype)
            flat_subtracted3D = np.empty([shape[0], shape[1], num_flat], dtype = datatype)
            for i in range(len(paths_dark)):
                dark3D[:,:,i]  = fits.open(paths_dark[i]) [0].data
            master_dark = np.median(dark3D, axis=2)
            for i in range(len(paths_flat)):
                flat_subtracted3D[:,:,i]  = fits.open(paths_flat[i]) [0].data - master_dark
                # rescale each flat frame by the mode
                [mode, _] = stats.mode(flat_subtracted3D[:,:,i], axis=None)
                flat_subtracted3D[:,:,i] = flat_subtracted3D[:,:,i]/mode
            # master flat, median combine all flats
            master_flat_subtracted = np.median(flat_subtracted3D, axis=2)  
            # normalised by the mean
            master_flat_subtracted_normalised = master_flat_subtracted / np.mean(master_flat_subtracted)
            
        # now master_flat_subtracted_normalised has been created in any case
        # create normalised master_subtracted using flat-field
        master_subtracted_normalised = master_subtracted / master_flat_subtracted_normalised

#%% No LEDoff is available, then standard reduction is preferred
# this section needs to be expanded to accomodate other frame type combinations

elif ccdBiasExist and ccdDarkExist and ccdObjectExist:
    
    # if flat does not exist, reduce images without flat
    if not ccdFlatExist:
        
        # in order to preallocate arrays, get some info
        hdulist = fits.open(paths_bias[0])
        shape = hdulist[0].data.shape   # get image dimension
        # dtype = hdulist[0].data.dtype     # get data type # always float16
        # get numbers of files for each type
        num_bias   = len(paths_bias)
        num_dark   = len(paths_dark)
        num_object = len(paths_object)
        
        # preallocate arrays
        bias3D  = np.empty([shape[0], shape[1], num_bias],  dtype = datatype)
        dark3D  = np.empty([shape[0], shape[1], num_dark],  dtype = datatype)
        object3D = np.empty([shape[0], shape[1], num_object], dtype = datatype)
        
        # fill arrays with fits data
        for i in range(len(paths_bias)):
            bias3D[:,:,i]  = fits.open(paths_bias[i]) [0].data
        for i in range(len(paths_dark)):
            dark3D[:,:,i]  = fits.open(paths_dark[i]) [0].data
        for i in range(len(paths_object)):
            object3D[:,:,i] = fits.open(paths_object[i])[0].data
        
        # master bias from  median combining
        master_bias = np.median(bias3D, axis=2)
        
        # master dark from  median combining, -master_bias
        master_dark = np.median(dark3D, axis=2) - master_bias

        # master object from median combining, -master_bias, -master_dark
        master_object = np.median(object3D, axis=2)
        
        # master
        master = master_object - master_bias - master_dark

        #debugging
        #maxindex = np.argmax(master_bias)
        #maxindex = np.unravel_index(maxindex, (1020,1530))
        #print(bias3D[maxindex[0],maxindex[1],:])
        #print(dark3D[maxindex[0],maxindex[1],:])
        #print(flat3D[maxindex[0],maxindex[1],:])
        #print(object3D[maxindex[0],maxindex[1],:])
        
    # otherwise, use flat-field correction to normalise
    elif ccdUseFlat:
        num_flat = len(paths_flat)
        flat3D = np.empty([shape[0], shape[1], num_flat],  dtype = datatype)
        for i in range(len(paths_flat)):
            # remove master_bias, master_dark
            flat3D[:,:,i]  = fits.open(paths_flat[i]) [0].data -master_bias -master_dark
            [mode, _ ] = stats.mode(flat3D[:,:,i], axis = None)
            flat3D[:,:,i] = flat3D[:,:,i] / mode
        master_flat = np.median(flat3D, axis=2)   # median combine all flat frames
        master_flat_normalised = master_flat / np.mean(master_flat)   # normalised by the mean
        # divide out normalised flat to obtain final master image
        master_normalised = master_object / master_flat_normalised
        
#%% save all master images
    
# check and create master folder
os.chdir(datasetDir)
if not os.path.exists('master'):
    os.makedirs('master')
os.chdir(os.path.join(datasetDir, 'master'))

### the folloiwng must exist and be saved

# save master_object
hdu_master_object = fits.PrimaryHDU(master_object) # create hdu objects
hdu_master_object.writeto('master_object.fits', clobber = True) # write hdu to fits

### the following may or may not exist

# save _master_object_subtracted
if 'master_object_subtracted' in locals():
    hdu_master_object_subtracted = fits.PrimaryHDU(master_object_subtracted) # create hdu objects
    hdu_master_object_subtracted.writeto('master_object_subtracted.fits', clobber = True) # write hdu to fits

# save master_ledoff
if 'master_ledoff' in locals():
    hdu_master_ledoff = fits.PrimaryHDU(master_ledoff)
    hdu_master_ledoff.writeto('master_ledoff.fits', clobber = True)

# save master_ledoff_subtracted
if 'master_ledoff_subtracted' in locals():
    hdu_master_ledoff_subtracted = fits.PrimaryHDU(master_ledoff_subtracted)
    hdu_master_ledoff_subtracted.writeto('master_ledoff_subtracted.fits', clobber = True)

# save master_bias
if 'master_bias' in locals():
    hdu_master_dark = fits.PrimaryHDU(master_bias)
    hdu_master_dark.writeto('master_bias.fits', clobber = True)
    
# save master_dark
if 'master_dark' in locals():
    hdu_master_dark = fits.PrimaryHDU(master_dark)
    hdu_master_dark.writeto('master_dark.fits', clobber = True)
    
# save master_flat_subtracted_normalised
if 'master_flat' in locals():
    hdu_master_flat = fits.PrimaryHDU(master_flat_normalised)
    hdu_master_flat.writeto('master_flat_normalised.fits', clobber = True)
    hdu_master_flat_normalised = fits.PrimaryHDU(master_flat_normalised)
    hdu_master_flat_normalised.writeto('master_flat_normalised.fits', clobber = True)
    
# save master_flat_subtracted_normalised
if 'master_flat_subtracted' in locals():
    hdu_master_flat_subtracted = fits.PrimaryHDU(master_flat_subtracted)
    hdu_master_flat_subtracted.writeto('master_flat_subtracted.fits', clobber = True)
    hdu_master_flat_subtracted_normalised = fits.PrimaryHDU(master_flat_subtracted_normalised)
    hdu_master_flat_subtracted_normalised.writeto('master_flat_subtracted_normalised.fits', clobber = True) 

# save master_subtracted
if 'master_subtracted' in locals():
    hdu_master_subtracted = fits.PrimaryHDU(master_subtracted)
    hdu_master_subtracted.writeto('master_subtracted.fits', clobber = True)

# save master_subtracted_normalised
if 'master_subtracted_normalised' in locals():
    hdu_master_subtracted_normalised  = fits.PrimaryHDU(master_subtracted_normalised)
    hdu_master_subtracted_normalised.writeto('master_subtracted_normalised.fits', clobber = True)
    
# save master_normalised
if 'master' in locals():
    hdu_master = fits.PrimaryHDU(master)
    hdu_master.writeto('master.fits', clobber = True)
    
# save master_normalised
if 'master_normalised' in locals():
    hdu_master_normalised = fits.PrimaryHDU(master_normalised)
    hdu_master_normalised.writeto('master_normalised.fits', clobber = True)


#%% Save master.log containing parameters and date and time

masterLogText = \
    'UTC time:'                                                        +'\n'+\
    '    ' + str(datetime.datetime.utcnow())                           +'\n'+\
    ''                                                                 +'\n'+\
    'Local time:'                                                      +'\n'+\
    '    ' + str(datetime.datetime.now())                              +'\n'+\
    ''                                                                 +'\n'+\
    'Master images were created under the following conditions:'       +'\n'+\
    '    ccdPreferLEDoff          = ' + repr(ccdPreferLEDoff)          +'\n'+\
    '    ccdUseFlat               = ' + repr(ccdUseFlat)               +'\n'+\
    '    ccdBiasExist             = ' + repr(ccdBiasExist)             +'\n'+\
    '    ccdDarkExist             = ' + repr(ccdDarkExist)             +'\n'+\
    '    ccdFlatSubtractedExist   = ' + repr(ccdFlatSubtractedExist)   +'\n'+\
    '    ccdFlatExist             = ' + repr(ccdFlatExist)             +'\n'+\
    '    ccdObjectSubtractedExist = ' + repr(ccdObjectSubtractedExist) +'\n'+\
    '    ccdObjectExist           = ' + repr(ccdObjectExist)           +'\n'+\
    '    ccdLEDoffSubtractedExist = ' + repr(ccdLEDoffSubtractedExist) +'\n'+\
    '    ccdLEDoffExist           = ' + repr(ccdLEDoffExist)
    
print(masterLogText + '\n')
print('Saving master.log to', datasetDir)
masterLog = open('master.log', 'w')
masterLog.write(masterLogText)
masterLog.close()

sys.exit() # end after image reduction

#%% aperture photometry from source segmentation

# determine threshold for background detection

# if LEDoff was used, get threshold from LEDoff/background
if 'master_ledoff_subtracted' in locals():
    ledoff_pred = np.mean(master_ledoff_subtracted) * np.ones(master_ledoff_subtracted.shape)
    mse = mean_squared_error(master_ledoff_subtracted, ledoff_pred)  
    rmse = np.sqrt(mse)
    threshold = 7.0 * rmse
    
# if no LEDoff was used, background subtraction is needed
else:
    # create preliminary mask    
    """ make_source_mask not yet available in photutils v0.2.1
        wait for v0.3 release
    """
    #from photutils import make_source_mask
    #masterMask = make_source_mask(master, snr=2, npixels=5, dilate_size=11)
    
    # background subtraction
    """ create 2D image of background and background rms and 
        apply sigma-clipping to each region in the low-res  
        background map to get mean, median, and std/rms. 
        sigma-clipping is the most widely used method though not as 
        good as using mask; still superior to robust standard 
        deviation using median absolute deviation (MAD-STD)
    """
    
    # [mean, median, std] = sigma_clipped_stats(master, sigma=3.0, iters=5)
    bkg = Background(master_normalised, (100, 100), filter_shape=(3, 3), method='median')
    # bkg = Background(master, (50, 50), filter_size=(3, 3), method='median')
    # plt.imshow(bkg.background, norm=normalisation, origin='lower', cmap=plt.cm.gray)
    plt.imshow(bkg.background, origin='lower', cmap=plt.cm.gray)
    [fig, ax] = plt.subplots(figsize=(8, 8))
    # make background-substracted image
    master_subtracted = master_normalised - bkg.background
    hdu_subtracted = fits.PrimaryHDU(master_subtracted)
    # save background subtracted image
    hdu_subtracted.writeto('master_subtracted.fits', clobber = True)
    # plot
    plt.imshow(master_subtracted, origin='lower', cmap=plt.cm.gray)
    
    # segmentation at a given sigma level, for regional properties
    threshold = 5.0 * bkg.background_rms # since data is background-subtracted

# perform segmentation whether flat was available or not
if 'master_subtracted_normalised' in locals():
    segm = detect_sources(master_subtracted_normalised, threshold, npixels=5)
elif 'master_subtracted' in locals():
    segm = detect_sources(master_subtracted, threshold, npixels=5)
print(segm.labels)
cmapRand = random_cmap(segm.max+1, random_state=12345)
plt.imshow(segm, origin='lower', cmap=cmapRand)

# measure regional source properties from segmentation
# the centroid is from image moments, already intensity-weighted
if 'bkg' in locals():
    props = source_properties(master_subtracted, segm,
        error = bkg.background_rms, background = bkg.background)
else:
    props = source_properties(master_subtracted_normalised, segm,
        error = master_ledoff_subtracted - np.mean(master_ledoff_subtracted),
        background = master_ledoff_subtracted)
        
# instrumental magnitude = -2.5 * log10(flux)
for i in range(len(props)):
    props[i].mag_instr = -2.5 * np.log10(props[i].source_sum)
propsTableColumns = ['id', 'xcentroid', 'ycentroid', 'area', 'max_value',
    'source_sum', 'mag_instr']
# there are other properties available, see list of SourceProperties
# http://photutils.readthedocs.io/en/latest/api/photutils.segmentation.SourceProperties.html#photutils.segmentation.SourceProperties
    
propsTable = properties_table(props, columns = propsTableColumns)
print(propsTable)

#%% plots for visualisation

apertures = []
for prop in props:
    position = (prop.xcentroid.value, prop.ycentroid.value)
    a = prop.semimajor_axis_sigma.value * 3.0
    b = prop.semiminor_axis_sigma.value * 3.0
    theta = prop.orientation.value
    apertures.append(EllipticalAperture(position, a, b, theta=theta))
norm = ImageNormalize(stretch=SqrtStretch())
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 18))

if 'bkg' in locals():
    ax1.imshow(master_subtracted, origin='lower', cmap='Greys_r', norm=norm)
else:
    ax1.imshow(master_subtracted_normalised, origin='lower', cmap='Greys_r', norm=norm)
ax2.imshow(segm, origin='lower', cmap=cmapRand)
for aperture in apertures:
    aperture.plot(color='blue', lw=1.5, alpha=0.5, ax=ax1)
    aperture.plot(color='white', lw=1.5, alpha=1.0, ax=ax2)
    