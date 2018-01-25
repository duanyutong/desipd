# -*- coding: utf-8 -*-
"""
Created on Mon May 16 14:41:22 2016

@author: Duan Yutong (dyt@lbl.gov)
"""

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


#%% function test values

path_dataset = r'C:\Users\givoltage\Downloads\20160615'
path_flat = None
save_dir_rel = ''
clobber = True

#%% set image directory

# test files
# directory

# os.chdir(dir)
# os.listdir(dir)
# parentdir = r'K:\Google Drive\BUPHY\zeta gem\as203-2016_spring-images'
# parentDir = os.path.join(os.path.expanduser("~"), r'Downloads')
# use os.sep after drive letter
# datasetID = '20160615'

#%% function for data reduction

def reduce_dataset(path_dataset,
                   path_flat = None,
                   save_dir_rel = '',
                   clobber = False):

    '''
    To reduce a dataset with corresponding flat frames, specify flats in kw.
    To reduce a dataset without flat, leave path_flat as none.
    To reduce a dataset of flat images only, treat it as a dataset w/o flat.

    '''

    #%% this is a naive version without flats

    def reduce_dataset_woflat(path_dataset,
                              save_dir_rel = '',
                              clobber = True,
                              dataset_is_flat = False,
                              suppress_return = True):

        # check and create master folder
        dir_save = os.path.join(path_dataset, save_dir_rel)
        if not os.path.exists(dir_save):
            os.makedirs(dir_save)

        # read in fits file
        path_pattern = os.path.join(path_dataset, '*.fit*')
        paths = glob.glob(path_pattern)
        # exclude existing master files
        paths = [path for path in paths if 'master' not in path]

        # read obstype from all headers
        obstypes = []
        for path in paths:
            hdulist = fits.open(path)
            obstype = hdulist[0].header['OBSTYPE']
            obstypes.append(obstype)
        frames_dict = {'path': paths}
        # frames_dict['filename'] = [os.path.basename(path) for path in paths]
        frames_dict['obstype'] = obstypes

        # perform search based on obstype
        paths_bias = [frames_dict['path'][i] for i in range(len(paths))
            if frames_dict['obstype'][i] == (
                'ZERO' or 'zero' or 'BIAS' or 'bias')]
        paths_dark = [frames_dict['path'][i] for i in range(len(paths))
            if frames_dict['obstype'][i] == (
                'DARK' or 'dark')]
        if dataset_is_flat:
            paths_object = [frames_dict['path'][i] for i in range(len(paths))
                if frames_dict['obstype'][i] == (
                'FLAT' or 'flat' or 'DOME FLAT' or 'dome flat' or 'SKY FLAT'
                or 'sky flat')]
        else:
            paths_object = [frames_dict['path'][i] for i in range(len(paths))
                if frames_dict['obstype'][i] == (
                    'OBJECT' or 'object' or 'LIGHT' or 'light')]

        # get info in order to preallocate arrays
        hdulist = fits.open(paths_bias[0])
        shape = hdulist[0].data.shape 	# get image dimension
        # dtype = hdulist[0].data.dtype 	# get data type # always float16
        datatype = np.float32   # let's just have the utmost precision
        # get numbers of files for each type
        num_bias = len(paths_bias)
        num_dark = len(paths_dark)
        num_object = len(paths_object)

        # preallocate arrays
        bias3D = np.empty([shape[0], shape[1], num_bias], dtype = datatype)
        dark3D = np.empty([shape[0], shape[1], num_dark], dtype = datatype)
        object3D = np.empty([shape[0], shape[1], num_object], dtype = datatype)

        # fill arrays with fits data
        for i in range(num_bias):
            bias3D[:, :, i] = fits.open(paths_bias[i]) [0].data
        for i in range(num_dark):
            dark3D[:, :, i] = fits.open(paths_dark[i]) [0].data

        # master bias from median combining
        m_bias = np.median(bias3D, axis=2)
        # master dark from median combining, -m_bias
        m_dark = np.median(dark3D, axis=2) - m_bias
        # create hdu objects
        hdu_m_bias = fits.PrimaryHDU(m_bias)
        hdu_m_dark = fits.PrimaryHDU(m_dark)

        # # write exposure time in header for flux/magnitude calculation
        expreq_bias = fits.open(paths_bias[0]) [0].header['EXPREQ']
        expreq_dark = fits.open(paths_dark[0]) [0].header['EXPREQ']
        hdu_m_bias.header ['EXPREQ'] = expreq_bias
        hdu_m_dark.header ['EXPREQ'] = expreq_dark
        # if exptime is present, add it, too
        if 'EXPTIME' in fits.open(paths_dark[0]) [0].header:
            exptime_bias = fits.open(paths_bias[0]) [0].header['EXPTIME']
            exptime = fits.open(paths_dark[0]) [0].header['EXPTIME']
            hdu_m_bias.header ['EXPTIME'] = exptime_bias
            hdu_m_dark.header ['EXPTIME'] = exptime

        # write hdu to fits
        hdu_m_bias.writeto(
            os.path.join(dir_save, 'master_bias.fits'), clobber = clobber)
        hdu_m_dark.writeto(
            os.path.join(dir_save, 'master_dark.fits'), clobber = clobber)

        # master light which can be object or flat
        if dataset_is_flat:

            paths_flat = paths_object
            num_flat = len(paths_flat)
            flat3D = np.empty([shape[0], shape[1], num_flat], dtype = datatype)

            for i in range(len(paths_flat)):
                # subtract bias and dark
                flat3D[:,:,i] = fits.open(paths_flat[i]) [0].data \
                    - m_bias - m_dark
                [mode, _ ] = stats.mode(flat3D[:, :, i], axis=None)
                # each flat is divide by mode
                flat3D[:,:,i] = object3D[:,:,i]/mode
            # median combine all flat frames
            m_flat = np.median(flat3D, axis=2)
            # normalised by the mean of master flat
            m_flat_normalised = m_flat / np.mean(m_flat)
            # create hdu objects
            hdu_m_flat = fits.PrimaryHDU(m_flat)
            hdu_m_flat_normalised = fits.PrimaryHDU(m_flat_normalised)

            # write exposure time in header for flux/magnitude calculation
            expreq_flat = fits.open(paths_flat[0]) [0].header['EXPREQ']
            hdu_m_flat.header ['EXPREQ'] = expreq_flat
            hdu_m_flat_normalised.header ['EXPREQ'] = expreq_flat
            # if exptime is present, add it, too
            if 'EXPTIME' in fits.open(paths_flat[0]) [0].header:
                exptime = fits.open(paths_object[0]) [0].header['EXPTIME']
                hdu_m_flat.header ['EXPTIME'] = exptime
                hdu_m_flat_normalised.header ['EXPTIME'] = exptime

            # write hdu to fits
            hdu_m_flat.writeto(
                os.path.join(dir_save, 'master_flat.fits'),
                clobber = clobber)
            hdu_m_flat_normalised.writeto(
                os.path.join(dir_save, 'master_flat_normalised.fits'),
                clobber = clobber)

            if not suppress_return:
                return m_flat_normalised

        else: # dataset is object, not flat

            for i in range(len(paths_object)):
                object3D[:, :, i] = fits.open(paths_object[i]) [0].data
            # master object from median combining, -m_bias, -m_dark
            m_object = np.median(object3D, axis=2) - m_bias - m_dark
            # create hdu objects
            hdu_m_object = fits.PrimaryHDU(m_object)

            # write exposure time in header for flux/magnitude calculation
            expreq_object = fits.open(paths_object[0]) [0].header['EXPREQ']
            hdu_m_object.header ['EXPREQ'] = expreq_object
            # if exptime is present, add it, too
            if 'EXPTIME' in fits.open(paths_object[0]) [0].header:
                exptime = fits.open(paths_object[0]) [0].header['EXPTIME']
                hdu_m_object.header ['EXPTIME'] = exptime
            # write hdu to fits
            hdu_m_object.writeto(os.path.join(dir_save, 'master_object.fits'),
                                 clobber = clobber)

            if not suppress_return:
                return m_object

    #%% data reduction actually starts here

    # determine whether there is flat available
    if path_flat is None:
        # no flat avilable, just reduce the given dataset
        reduce_dataset_woflat(path_dataset,
                              save_dir_rel = save_dir_rel,
                              clobber = clobber,
                              dataset_is_flat = False,
                              suppress_return = True)
    else:
        # then we gotta reduce the dataset of flats
        # and then do flat field correction
        m_object = reduce_dataset_woflat(path_dataset,
                                         save_dir_rel = save_dir_rel,
                                         clobber = clobber,
                                         dataset_is_flat = False,
                                         suppress_return = False)
        m_flat_normalised = reduce_dataset_woflat(path_flat,
                                                  save_dir_rel = save_dir_rel,
                                                  clobber = clobber,
                                                  dataset_is_flat = True,
                                                  suppress_return = False)

        # flat field correction
        master = m_object / m_flat_normalised
        hdu_master = fits.PrimaryHDU(master)
        hdu_master.writeto(
            os.path.join(path_dataset, save_dir_rel, 'master.fits'),
            clobber = clobber)
