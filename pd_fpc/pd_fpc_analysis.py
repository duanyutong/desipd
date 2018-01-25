# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 00:48:30 2016

@author: Duan Yutong (dyt@physics.bu.edu, dyt@lbl.gov)

"""

import os
import logging
import pickle
import numpy as np
from scipy import interpolate
from skimage.measure import marching_cubes
from dateutil.parser import parse

import matplotlib
matplotlib.use('Agg') # in case no backend was automatically chosen
import matplotlib.pyplot as plt
plt.ioff() # turn off interative mode to suppress plot output in console
from matplotlib.patches import Circle
# from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import art3d, Axes3D
# from matplotlib.legend_handler import HandlerLine2D
from matplotlib.dates import date2num, DateFormatter
from matplotlib.backends.backend_pdf import PdfPages

import astropy
from astropy.io import fits
from photutils import EllipticalAperture

os.chdir(os.path.dirname(os.path.abspath(__file__)))
from pd_fpc_photometry import segmentation_photometry, aper_photometry

#%% FPC analysis class for data reduction and photometry

class FPCAnalysis:

    def __init__(self, analysis_name):

        """
        Use appropriate ANALYSIS_NAME including: field name/ID, date and time.

        """

        # set constants
        if os.name == 'nt':
            self.DIR_FITS = r'D:\fits'
            self.DIR_ANALYSIS = r'D:\fpc_analysis'
        elif os.name == 'posix':
            self.DIR_FITS = '/global/project/projectdirs/' \
                            +'desi/protodesi/images/fits/'
            self.DIR_ANALYSIS = '/global/project/projectdirs/' \
                                +'desi/protodesi/images/fpc_analysis/'
            #self.DIR_ANALYSIS = '/global/u2/d/dyt/fpc_analysis/'

        self.DUMP_PICKLE = False
        self.CLOBBER = False
        self.ANALYSIS_NAME = analysis_name
        self.create_dir(self.ANALYSIS_NAME)
        self.DIR_SAVE = os.path.join(self.DIR_ANALYSIS, self.ANALYSIS_NAME)
        self.DATATYPE = np.float32
        self.SHAPE = None

        # fibre parameters
        self.FIBRE_IDS = {'3001', '3002', '3003'}  # fibre names in a set
        self.FIBRE_BEST = '3002' # best fibre which has the highest througput
        self.FIBRE_COORDS = {'3001': (1196, 1058), # when telescope at zenith
                             '3002': (1855, 1059),
                             '3003': (1200, 1374)}
        self.FIBRE_COORDS_FLAT =  {'3001': (1202, 1089), # when at white spot
                                   '3002': (1860, 1090),
                                   '3003': (1205, 1404)}
        self.FIBRE_N = len(self.FIBRE_IDS)
        self.FIBRE_DIST_THRESHOLD = 40 # should be 10 px if agree, 38 if not
        self.FIBRE_RADIUS = 53 # in microns
        # define approximate isophotal ellipses for each object
        self.ISOPHOTAL_EXT = 5 # approximate isophotal extent
        self.FIBRE_ELLIPTICAL_APERTURES = {
            '3001':
                {'a': 13*self.ISOPHOTAL_EXT,
                 'b': 10*self.ISOPHOTAL_EXT,
                 'theta':  0.055},
            '3002':
                {'a': 12*self.ISOPHOTAL_EXT,
                 'b':  9*self.ISOPHOTAL_EXT,
                 'theta': -0.095},
            '3003':
                {'a': 17*self.ISOPHOTAL_EXT,
                 'b': 12*self.ISOPHOTAL_EXT,
                 'theta': -0.125}}
        self.FIBRE_ELLIPTICAL_APERTURES_FLAT = {
            '3001':
                {'a': 12.4090*self.ISOPHOTAL_EXT,
                 'b':  9.1608*self.ISOPHOTAL_EXT,
                 'theta': 0.063960888},
            '3002':
                {'a': 11.3185*self.ISOPHOTAL_EXT,
                 'b':  8.7111*self.ISOPHOTAL_EXT,
                 'theta': -0.059023251},
            '3003':
                {'a': 16.6673*self.ISOPHOTAL_EXT,
                 'b': 11.8841*self.ISOPHOTAL_EXT,
                 'theta': -0.135632276}}
        # plate scale from in unit of microns/arcsec
        self.PLATE_SCALE = 57.3
        self.Z_CIRCLE = -1500 # vertical position of bottom plane in 3d plot
        # dither sequence grid indices
        self.GRID_INDEX_TEL = np.array([
            [11, 10,  9, 24, 23],
            [12,  2,  1,  8, 22],
            [13,  3,  0,  7, 21],
            [14,  4,  5,  6, 20],
            [15, 16, 17, 18, 19]])

        self.GRID_INDEX_POS = np.array([
            [19, 18, 17, 16, 15],
            [20,  6,  5,  4, 14],
            [21,  7,  0,  3, 13],
            [22,  8,  1,  2, 12],
            [23, 24,  9, 10, 11]])

        # Photometry parameters
        self.BKG_SIGMA = 1.5
        self.SOURCE_SNR = 1.01
        self.FWHM_KERNEL = 25
        self.X_SIZE_KERNEL = 100    # default is 8*sigma
        self.Y_SIZE_KERNEL = 80     # default is 8*sigma

        # initialise dictionary storing dataset metadata
        self.datasets={}

        # setup logger
        self.logger = logging.getLogger()
        logging.basicConfig()
        try:
            self.logger.removeHandler(self.logger.handlers[0])
        except:
            # this is the first logger instance with no handle attached
            pass
        logpath = os.path.join(self.DIR_SAVE,
                               'fpc_analysis-' + self.ANALYSIS_NAME + '.log')
        logging.basicConfig(filename = logpath,
                            level = logging.INFO,
                            format = '%(asctime)s %(name)s '
                                        + '[%(levelname)s]: %(message)s',
                            datefmt= '%Y/%m/%d %H:%M:%S')
        print('Log is being saved to:', logpath)

    def create_dir(self, dir_rel):

        """
        given relative path, create folder for this analysis under DIR_ANALYSIS

        """

        dir_abs = os.path.join(self.DIR_ANALYSIS, dir_rel)
        if not os.path.exists(dir_abs):
            os.makedirs(dir_abs)
            print('Created new directory:', dir_abs)
        else:
            print('Directory already exists:', dir_abs)

    def expid_to_filename(self, expid_input):

        if isinstance(expid_input, (int, np.integer)):
            filename = 'PROTODESI_' + str(expid_input).zfill(8) + '.fits'
            return filename

        else:
            filenames = [None]*len(expid_input)
            for (i, expid) in enumerate(expid_input):
                filenames[i] = self.expid_to_filename(expid)
            return filenames

    def expid_to_filename_reduced(self, expid_input):

        if isinstance(expid_input, (int, np.integer)):
            filename = 'PROTODESI_' \
                        + str(expid_input).zfill(8) \
                        + '-reduced.fits'
            return filename

        else:
            filenames = [None]*len(expid_input)
            for (i, expid) in enumerate(expid_input):
                filenames[i] = self.expid_to_filename_reduced(expid)
            return filenames

    def expid_to_path(self, expid_input):

        if isinstance(expid_input, (int, np.integer)):
            path = os.path.join(self.DIR_FITS,
                                self.expid_to_filename(expid_input))
            return path
        else:
            paths_abs = [None]*len(expid_input)
            for (i, expid) in enumerate(expid_input):
                paths_abs[i] = self.expid_to_path(expid)
            return paths_abs

    def expid_to_path_reduced(self, expid_input):

        if isinstance(expid_input, (int, np.integer)):
            path = os.path.join(self.DIR_SAVE,
                                self.expid_to_filename_reduced(expid_input))
            return path
        else:
            paths_abs = [None]*len(expid_input)
            for (i, expid) in enumerate(expid_input):
                paths_abs[i] = self.expid_to_path_reduced(expid)
            return paths_abs

    def filter_expids(self, exptypes_input, exptype_output,
                      expid_initial, expid_end):

        """
        if exposures were taken in alternating exptypes, this gives you a
        filtered list of a single type

        """

        expids = np.arange(expid_initial, expid_end+1, 1)
        expids_filtered = [expids[i] for i in range(len(expids))
            if i % len(exptypes_input) == exptypes_input.index(exptype_output)]

        return expids_filtered

    def specify_dataset(self, dataset_type , expid_initial, expid_final,
                        dataset_name, bias_name = None, dark_name = None,
                        tile_id = None,
                        focus = None,
                        dither_mode = None,
                        dither_grid = None,
                        dither_step = None):

        """
        accepted dataset_type:
            'bias stack': all frames are bias to be combined
            'dark stack': all frames are dark to be combined
            'image stack': repeated light frames to be combined and subtracted
            'image sequence': light frames that need bias dark subtraction only

        for image stack and image sequence, specify the corresponding bias
        and dark names

        """

        expids = list(np.arange(expid_initial, expid_final+1, 1))
        self.datasets[dataset_name]={
            'dataset_type': dataset_type,
            'expids': expids,
            'dataset_name': dataset_name,
            'bias_name': bias_name,
            'dark_name': dark_name,
            'tile_id': tile_id,
            'focus': focus,
            'dither_mode': dither_mode,
            'dither_grid': dither_grid,
            'dither_step': dither_step}

    def check_consistency(self, expids):

        """
        check if all expids are consistent

        all expids in a dataset should have same exptime and targets
        (except for bias)

        more checks to be added, incomplete

        """

        if len(expids) <= 0:
            self.logger.info('Empty expids list: {}'.format(expids))
            return True

        paths = self.expid_to_path(expids)
        exptimes = []

        for path in paths:
            try:
                header_fpc = fits.open(path)['FPC'].header
            except:
                header_fpc = fits.open(path)[0].header

            if header_fpc['OBSTYPE'].lower() in ['bias', 'zero']:
                exptime = 0
            else:
                exptime = header_fpc['EXPREQ']
                exptimes.append(exptime)

        if all(exptime == exptimes[0] for exptime in exptimes):
            self.logger.info('Consistency check passed for {}.'.format(expids))
            return True
        else:
            self.logger.error('Exptime is inconsistent for {}'.format(expids))
            self.logger.error('Exptimes are: {}'.format(exptimes))
            return False

    def read_dataset_header(self, dataset_name):

        """
        check consistency and read header
        more checks to be added, now it's just exptime

        """

        # check if all FPC image and dark frames have the same exptime
        dataset_type = self.datasets[dataset_name]['dataset_type']

        if dataset_type == 'bias stack':
            self.logger.info('Checking consistency for bias stack: {}...'
                .format(dataset_name))
            expids_check = self.datasets[dataset_name]['expids']
            if not self.check_consistency(expids_check):
                self.logger.info('Consistency check failed for bias stack {}'
                                 .format(dataset_name))

        elif dataset_type == 'dark stack':
            self.logger.info('Checking consistency for dark stack: {}...'
                .format(dataset_name))
            expids_check = self.datasets[dataset_name]['expids']

            if not self.check_consistency(expids_check):
                self.logger.info('Consistency check failed for '
                                 + 'dark stack {}'.format(dataset_name))

        elif dataset_type == 'image stack':
            self.logger.info('Checking consistency for image stack: {}...'
                .format(dataset_name))
            dark_name = self.datasets[dataset_name]['dark_name']
            expids_check = self.datasets[dataset_name]['expids'] \
                            + self.datasets[dark_name]['expids']

            if not self.check_consistency(expids_check):
                self.logger.info('Consistency check failed for '
                                 + 'image stack {}'.format(dataset_name))

        elif dataset_type == 'image sequence':
            self.logger.info('Checking consistency for image sequence: {}...'
                .format(dataset_name))
            dark_name = self.datasets[dataset_name]['dark_name']
            expids_check = self.datasets[dataset_name]['expids'] \
                            + self.datasets[dark_name]['expids']
            if not self.check_consistency(expids_check):
                self.logger.info('Consistency check failed for '
                                 + 'image sequence {}'.format(dataset_name))
        else:
            self.logger.error('Invalid exptype for dataset {}'
                .format(dataset_name))
            return False

        # if everything checks out and no False has been returned, read header
        expids = self.datasets[dataset_name]['expids']
        for expid in expids:
            self.datasets[dataset_name][expid] = {}
            try:
                hdrpri = fits.open(self.expid_to_path(expid))['PRIMARY'].header
                hdrfpc = fits.open(self.expid_to_path(expid))['FPC'].header
            except:
                hdrpri = fits.open(self.expid_to_path(expid))[0].header
                hdrfpc = fits.open(self.expid_to_path(expid))[0].header
            self.datasets[dataset_name][expid]['header_primary'] = hdrpri
            self.datasets[dataset_name][expid]['header_fpc'] = hdrfpc

        # add expreq from FPC to exptime
        if hdrfpc['OBSTYPE'].lower() in ['bias', 'zero']:
            self.datasets[dataset_name]['exptime'] = 0
            self.logger.info('Exptime set to 0 for bias stack {}'
                                 .format(dataset_name))
        else:
            exptime = \
                self.datasets[dataset_name][expids[0]]['header_fpc']['EXPREQ']
            self.datasets[dataset_name]['exptime'] = exptime
            self.logger.info('Exptime set to {}s for dataset {}'
                             .format(exptime, dataset_name))

    def read_datasets_header(self):

        # determine data shape
        self.set_shape()

        # read dataset header in order: bias, dark, image
        for dataset_name in self.datasets.keys():
            if self.datasets[dataset_name]['dataset_type'] == 'bias stack':
                self.read_dataset_header(dataset_name)
        for dataset_name in self.datasets.keys():
            if self.datasets[dataset_name]['dataset_type'] == 'dark stack':
                self.read_dataset_header(dataset_name)
        for dataset_name in self.datasets.keys():
            if self.datasets[dataset_name]['dataset_type'] not in \
                ['bias stack', 'dark stack']:
                self.read_dataset_header(dataset_name)


    def set_shape(self):

        dataset_name = list(self.datasets.keys())[0]
        path = self.expid_to_path(self.datasets[dataset_name]['expids'][0])
        try:
            self.SHAPE = fits.open(path)['FPC'].data.shape
        except:
            self.SHAPE = fits.open(path)[0].data.shape
        self.logger.info('Data shape set as: {}'.format(self.SHAPE))
        return True

    def reduce_stack(self, dataset_name):

        """
        Several scenarios:
            * bias stack to be combined
            * dark stack to be combined, using existing master bias
            * image stack to be combined, using existing master bias and dark
        """
        dataset_type = self.datasets[dataset_name]['dataset_type']

        if dataset_type == 'bias stack':
            return self.reduce_bias_stack(dataset_name)
        elif dataset_type == 'dark stack':
            return self.reduce_dark_stack(dataset_name)
        elif dataset_type == 'image stack':
            return self.reduce_image_stack(dataset_name)

    def reduce_bias_stack(self, dataset_name):

        self.logger.info('Reducing bias stack: {}...'.format(dataset_name))
        filename_master = dataset_name + '-master_bias.fits'
        filename_error = dataset_name + '-master_error_bias.fits'
        path_save_prefix = os.path.join(self.DIR_SAVE, dataset_name)
        self.datasets[dataset_name]['filename_master'] = filename_master
        self.datasets[dataset_name]['filename_error'] = filename_error
        # check if master bias already exists        
        if filename_master in os.listdir(self.DIR_SAVE) \
            and filename_error in os.listdir(self.DIR_SAVE) \
            and not self.CLOBBER:
            self.logger.info('Master bias files already exist. '
                             + 'Overwriting not allowed.')
            return True

        # get bias paths
        expids = self.datasets[dataset_name]['expids']
        paths_bias = self.expid_to_path(expids)
        # get numbers of files for each type
        num_bias = len(paths_bias)
        self.logger.info('Reading {} bias frames...'.format(num_bias))
        # preallocate array
        bias_stack = np.empty([self.SHAPE[0], self.SHAPE[1], num_bias],
                              dtype = self.DATATYPE)
        # fill arrays with fits data
        for i in range(num_bias):
            bias_stack[:, :, i] = fits.open(paths_bias[i]) ['FPC'].data
        # master bias from median combining
        m_bias = np.median(bias_stack, axis=2)
        m_bias_error = np.std(bias_stack, axis=2)

        # save
        hdu_m_bias = fits.PrimaryHDU(m_bias)
        hdu_m_bias_error = fits.PrimaryHDU(m_bias_error)
        hdu_m_bias.writeto(path_save_prefix+'-master_bias.fits', 
                            clobber = True)
        hdu_m_bias_error.writeto(path_save_prefix+'-master_error_bias.fits', 
                            clobber = True)
        self.logger.info('Master bias saved to {}'.format(self.DIR_SAVE))

    def reduce_dark_stack(self, dataset_name):

        self.logger.info('Reducing dark stack: {}...'.format(dataset_name))
        filename_master = dataset_name + '-master_dark.fits'
        filename_error = dataset_name + '-master_error_dark.fits'
        path_save_prefix = os.path.join(self.DIR_SAVE, dataset_name)
        self.datasets[dataset_name]['filename_master'] = filename_master
        self.datasets[dataset_name]['filename_error'] = filename_error
        # check if master dark already exists
        if filename_master in os.listdir(self.DIR_SAVE) \
            and filename_error in os.listdir(self.DIR_SAVE) \
            and not self.CLOBBER:
            self.logger.info('Master dark files already exist. '
                             + 'Overwriting not allowed.')
            return True

        # read master bias
        bias_name = self.datasets[dataset_name]['bias_name']
        bias_filename = self.datasets[bias_name]['filename_master']
        self.logger.info('Reading master bias {}...'.format(bias_filename))
        m_bias = fits.open(os.path.join(self.DIR_SAVE, bias_filename)) [0].data

        # get dark paths
        expids = self.datasets[dataset_name]['expids']
        paths_dark = self.expid_to_path(expids)
        # get numbers of files for each type
        num_dark = len(paths_dark)
        self.logger.info('Reading {} dark frames...'.format(num_dark))
        # preallocate array
        dark_stack = np.empty([self.SHAPE[0], self.SHAPE[1], num_dark],
                              dtype = self.DATATYPE)
        # fill arrays with fits data
        for i in range(num_dark):
            dark_stack[:, :, i] = fits.open(paths_dark[i]) ['FPC'].data \
                                    - m_bias
        # master dark from median combining -m_bias
        m_dark = np.median(dark_stack, axis=2)
        m_dark_error = np.std(dark_stack, axis=2)

        # save
        hdu_m_dark = fits.PrimaryHDU(m_dark)
        hdu_m_dark_error = fits.PrimaryHDU(m_dark_error)
        hdu_m_dark.writeto(path_save_prefix+'-master_dark.fits', 
                            clobber = True)
        hdu_m_dark_error.writeto(path_save_prefix+'-master_error_dark.fits', 
                            clobber = True)
        self.logger.info('Master bias saved to {}'.format(self.DIR_SAVE))

    def reduce_image_stack(self, dataset_name):

        self.logger.info('Reducing image stack: {}...'.format(dataset_name))
        if self.datasets[dataset_name]['dataset_type'] != 'image stack':
            self.logger.error('Wrong dataset type.')
            return False       
        
        filename_master = dataset_name + '-master_image.fits'
        filename_error = dataset_name + '-master_error_image.fits'
        path_save_prefix = os.path.join(self.DIR_SAVE, dataset_name)
        self.datasets[dataset_name]['filename_master'] = filename_master
        self.datasets[dataset_name]['filename_error'] = filename_error
        # check if master image already exists
        if filename_master in os.listdir(self.DIR_SAVE) \
            and filename_error in os.listdir(self.DIR_SAVE) \
            and not self.CLOBBER:
            self.logger.info('Master image files already exist. '
                             + 'Overwriting not allowed.')
            return True

        # read master bias
        bias_name = self.datasets[dataset_name]['bias_name']
        bias_filename = self.datasets[bias_name]['filename_master']
        self.logger.info('Reading master bias {}...'.format(bias_filename))
        m_bias = fits.open(os.path.join(self.DIR_SAVE, bias_filename)) [0].data

        # read master dark
        dark_name = self.datasets[dataset_name]['dark_name']
        dark_filename = self.datasets[dark_name]['filename_master']
        self.logger.info('Reading master dark {}...'.format(bias_filename))
        m_dark = fits.open(os.path.join(self.DIR_SAVE, dark_filename)) [0].data

        expids = self.datasets[dataset_name]['expids']
        paths_image = self.expid_to_path(expids)
        # get numbers of files for each type
        num_image = len(paths_image)
        self.logger.info('Reading {} image frames...'.format(num_image))
        # preallocate array
        image_stack = np.empty([self.SHAPE[0], self.SHAPE[1], num_image],
                              dtype = self.DATATYPE)
        for i in range(len(paths_image)):
            try:
                image_stack[:, :, i] = fits.open(paths_image[i])['FPC'].data \
                                        - m_bias - m_dark
            except:
                image_stack[:, :, i] = fits.open(paths_image[i])[0].data \
                                        - m_bias - m_dark
        # master image from median combining, -m_bias, -m_dark
        m_image = np.median(image_stack, axis=2)
        m_error = np.std(image_stack, axis=2)
        # create hdu objects
        hdu_m_image = fits.PrimaryHDU(m_image)
        hdu_m_error = fits.PrimaryHDU(m_error)

        # write exptime to header
        hdu_m_image.header['EXPREQ'] = self.datasets[dataset_name]['exptime']
        hdu_m_image.header['EXPTIME'] = self.datasets[dataset_name]['exptime']

        hdu_m_image.writeto(path_save_prefix+'-master_image.fits', 
                            clobber = True)
        hdu_m_error.writeto(path_save_prefix+'-master_error_image.fits', 
                            clobber = True)
        self.logger.info('Master images saved to {}'.format(self.DIR_SAVE))

    def reduce_image_sequence(self, dataset_name):

        """
        input dataset is either:

        1. sequence of 5x5 continuous image frames from a dither pattern, or
        2. sequence of stability test repeated image exposures

        both of which cases need bias and dark subtraction

        """

        self.logger.info('Reducing image sequence: {}...'.format(dataset_name))

        # read master bias
        bias_name = self.datasets[dataset_name]['bias_name']
        bias_filename = self.datasets[bias_name]['filename_master']
        self.logger.info('Reading master bias {}...'.format(bias_filename))
        m_bias = fits.open(os.path.join(self.DIR_SAVE, bias_filename)) [0].data

        # read master dark
        dark_name = self.datasets[dataset_name]['dark_name']
        dark_filename = self.datasets[dark_name]['filename_master']
        self.logger.info('Reading master dark {}...'.format(dark_filename))
        m_dark = fits.open(os.path.join(self.DIR_SAVE, dark_filename)) [0].data

        expids = self.datasets[dataset_name]['expids']


        # paths_image = self.expid_to_path(expids)
        for i, expid in enumerate(expids):
        #for i, path_image in enumerate(paths_image):
            filename = self.expid_to_filename(expid)
            filename_save = self.expid_to_filename_reduced(expid)

            # check if reduced image already exists
            if filename_save in os.listdir(self.DIR_SAVE) and not self.CLOBBER:
                self.logger.info('{} already exists. Overwriting not allowed.'
                    .format(filename_save))
            else:
                self.logger.info('Reducing {} of {} in image sequence: {}...'
                    .format(i+1, len(expids), filename))
                # read and reduce image data
                hdulist = fits.open(self.expid_to_path(expid))
                hdulist['FPC'].data = hdulist['FPC'].data - m_bias - m_dark
                hdulist['FPC'].header['REDUCED'] = True
                hdulist['FPC'].header['COMMENT'] = 'This image has been ' \
                    + 'reduced by bias and dark subtraction'
                # save
                path_save = os.path.join(self.DIR_SAVE, filename_save)
                hdulist.writeto(path_save, clobber = self.CLOBBER)
                self.logger.info('Reduced {} of {} saved to {}'
                    .format(i+1, len(expids), path_save))

        self.logger.info('All {} reduced images in sequence saved to {}'
            .format(len(expids), self.DIR_SAVE))
        
        # prepare error image by summing error bias and error dark
        self.logger.info('Creating master error image for reduced frames...')
        bias_name = self.datasets[dataset_name]['bias_name']
        dark_name = self.datasets[dataset_name]['dark_name']
        path_error_bias = os.path.join(
                               self.DIR_SAVE, 
                               self.datasets[bias_name]['filename_error'])
        path_error_dark = os.path.join(
                               self.DIR_SAVE,
                               self.datasets[dark_name]['filename_error'])
        self.datasets[dataset_name]['filename_error'] = \
                dataset_name+'-master_error_image.fits'
        if dataset_name+'-master_error_image.fits' \
            in os.listdir(self.DIR_SAVE) \
            and not self.CLOBBER:
            self.logger.info('{} already exisits. Overwriting not allowed.'
                             .format(dataset_name+'-master_error_image.fits'))
            return True
        else:
            m_error = np.power( 
                               np.square(fits.open(path_error_bias)[0].data) 
                               + np.square(fits.open(path_error_dark)[0].data), 
                               1/2)
                    
            path_save_error = os.path.join(self.DIR_SAVE,
                                    dataset_name+'-master_error_image.fits')
            hdu_m_error = fits.PrimaryHDU(m_error)
            hdu_m_error.writeto(path_save_error,  clobber = self.CLOBBER)
            self.logger.info('Master error image saved to {}'
                             .format(path_save_error))

    def photometry(self, dataset_name):

        '''
        Performs photometry, compiles source properties, cleans up 
        self.datasets, and saves dict to pickle
        '''

        dataset_type = self.datasets[dataset_name]['dataset_type']

        if dataset_type == 'image sequence':

            expids = self.datasets[dataset_name]['expids']
            paths_reduced = self.expid_to_path_reduced(expids)

            # check if all reduced image files exist
            checklist = [os.path.isfile(path_reduced)
                            for path_reduced in paths_reduced]
            if all(checklist):
                self.logger.info('All reduced images exist and '
                                 + 'ready for photometry.')
            else:
                self.logger.error('Not all reduced image files exist: {}'
                    .format(checklist))
                return False
                
            # read in master error image for all image frames in the sequence
            filename_error = self.datasets[dataset_name]['filename_error']
            path_error = os.path.join(self.DIR_SAVE, filename_error)

            for i, expid in enumerate(expids):

                self.logger.info(('Performing photometry {} of {}, ' \
                    + 'exposure {}...')
                    .format(i+1, len(expids), expid))
                [segm, props_list] = segmentation_photometry(
                                        paths_reduced[i],
                                        path_error_abs = path_error,
                                        logger = self.logger,
                                        bkg_sigma = self.BKG_SIGMA,
                                        source_snr = self.SOURCE_SNR,
                                        fwhm_kernel = self.FWHM_KERNEL,
                                        x_size_kernel = self.X_SIZE_KERNEL,
                                        y_size_kernel = self.Y_SIZE_KERNEL,
                                        dump_pickle = self.DUMP_PICKLE,
                                        clobber = self.CLOBBER)
                # do not use segm.nlabels as segm may not exist
                # it is None if no source is detected
                self.logger.info('Detected {} sources from exposure {}.'
                    .format(len(props_list), expid)) 
                
                # analyse source properties for each expid
                self.analyse_photometry_results(
                    dataset_name, expid, props_list)
                # set variables to 0 to save memory
                del segm, props_list
         
            # make fibre photometry data into lists
            self.compile_sequence(dataset_name)

        elif dataset_type == 'image stack':

            filename_master = self.datasets[dataset_name]['filename_master']
            filename_error = self.datasets[dataset_name]['filename_error']
            path_master = os.path.join(self.DIR_SAVE, filename_master)
            path_error = os.path.join(self.DIR_SAVE, filename_error)
            [segm, props_list] = segmentation_photometry(
                        path_master,
                        path_error_abs = path_error,
                        logger = self.logger,
                        bkg_sigma = self.BKG_SIGMA,
                        source_snr = self.SOURCE_SNR,
                        fwhm_kernel = self.FWHM_KERNEL,
                        x_size_kernel = self.X_SIZE_KERNEL,
                        y_size_kernel = self.Y_SIZE_KERNEL,
                        dump_pickle = self.DUMP_PICKLE,
                        clobber = self.CLOBBER)
            # master image photometry results always saved under expid = 0
            self.analyse_photometry_results(dataset_name, 0, props_list)
            
            # del segm, props_list

    def analyse_photometry_results(self, dataset_name, expid, phot_input):
        
        dataset_type = self.datasets[dataset_name]['dataset_type']
        
        if dataset_type == 'image sequence':
            # organise results in SourceProperties into dictionary
            self.logger.info('Analysing photometry results of '
                + 'image sequence {}, exposure {}...'
                .format(dataset_name, expid))
                
        elif dataset_type == 'image stack':
            self.logger.info('Analysing photometry results of '
                + 'matser image: {}...'
                .format(self.datasets[dataset_name]['filename_master']))
                
        if type(phot_input) is list:
        # input photometry results are a list of source properties
            props_list = phot_input
            
            # fibre identification
            fibres_identified = set({})
            # initialise dict for the given expid
            if expid in self.datasets[dataset_name].keys():
                if type(self.datasets[dataset_name][expid]) is not dict:
                    self.datasets[dataset_name][expid] = {}
                else:
                    # this dict already exists. keep it intact
                    pass
            else:
                # expid not yet in keys, initialise it as a dict
                self.datasets[dataset_name][expid] = {}
    
            for fibre_id in self.FIBRE_IDS:
    
                self.logger.info('Looking for fibre {} in source properties...'
                    .format(fibre_id))
                # initialise dict for the given fibre_id under expid
                self.datasets[dataset_name][expid][fibre_id] = {}
    
                for props in props_list:
                    # clean up props object, remove large, unuseful attributes
                    props = self.clean_up_props(props)
                    # check distance with known fibre position
                    # prop.centroid is in (y, x) format
                    dist = np.linalg.norm(
                                np.fliplr([props.centroid.value])[0]
                                - self.FIBRE_COORDS[fibre_id])
                                
                    if dist < self.FIBRE_DIST_THRESHOLD:
                        
                        self.logger.info(('Fibre {} identified in exposure '
                                + '{}, distance to known centroid: {} pixels')
                                .format(fibre_id, expid, dist))
                        fibres_identified.add(fibre_id)
                        # read props and save it in dictionary
                        self.datasets[dataset_name][expid][fibre_id]['props'] \
                            = props
    
    #                    # remove unwanted props to save memory and disk space
    #                    if self.EXCLUDE_ATTRIBUTES:
    #                        for attribute in self.ATTRIBUTES_EXCLUDED:
    #                            setattr(self.datasets[dataset_name][expid] \
    #                                    [fibre_id]['props'],
    #                                    attribute,
    #                                    None)
    
                        # remove the props of identified source from props list
                        props_list.remove(props)
                        # stop searching for current fibre
                        break
                    else:
                        # detected source is not the fibre we are looking for
                        self.logger.info('Cannot match {} to known fibre {},' \
                                    ' distance to known centroid: {} pixels'
                                    .format(np.fliplr(
                                                [props.centroid.value])[0],
                                            fibre_id,
                                            dist))
            
            del props_list, phot_input
            
            self.datasets[dataset_name][expid]['fibres_identified'] = \
                fibres_identified
            self.logger.info('A total of {} fibres identified: {}'
                .format(len(fibres_identified), fibres_identified))
            
            # measure missing fibres
            if fibres_identified == self.FIBRE_IDS:
                self.logger.info(('All fibres have been identified in ' 
                                  + 'exposure {}')
                                  .format(expid))
            else:
                self.logger.warning('Fibres still missing: {}'
                                .format(self.FIBRE_IDS - fibres_identified))
                self.measure_missing_sources(dataset_name, expid)
            
        else:
        # input photometry results are a table from aperture photometry
            phot_table = phot_input
            fibres_identified = \
                self.datasets[dataset_name][expid]['fibres_identified']
            
            for fibre_id in self.FIBRE_IDS:
                
                self.logger.info(('Looking for missing fibre {} in '
                    + 'aperture photometry table...')
                    .format(fibre_id))
                    
                for row in range(len(phot_table)):
                    
                    centre = np.array((phot_table[row]['xcenter'], 
                              phot_table[row]['ycenter']))
                    dist = np.linalg.norm(centre - self.FIBRE_COORDS[fibre_id])
                    
                    if dist < self.FIBRE_DIST_THRESHOLD:
                        
                        self.logger.info(('Fibre {} identified in exposure '
                                + '{}, distance to known centroid: {} pixels')
                                .format(fibre_id, expid, dist))
                        fibres_identified.add(fibre_id)
                        # read phot_table and save it in dictionary
                        self.datasets[dataset_name][expid][fibre_id] \
                            ['phot_table'] = phot_table[row]
                        break
                    else:
                        # this row is not the fibre we are looking for
                        self.logger.info('Cannot match {} to known fibre {},' \
                                    ' distance to known centroid: {} pixels'
                                    .format(centre, fibre_id, dist))
            del phot_input, phot_table
            self.datasets[dataset_name][expid]['fibres_identified'] = \
                fibres_identified
                
    def measure_missing_sources(self, dataset_name, expid):
        
        # if not all three fibres detected, perform aperture photometry
        # using pre-defined apertures
        fibres_identified = \
            self.datasets[dataset_name][expid]['fibres_identified']
        if len(fibres_identified) < 3:
            fibres_missing = self.FIBRE_IDS - fibres_identified
            # prepare apertures for missing fibres
            apertures = []
            for fibre_id in fibres_missing:
                position = self.FIBRE_COORDS[fibre_id]
                a = self.FIBRE_ELLIPTICAL_APERTURES[fibre_id]['a']
                b = self.FIBRE_ELLIPTICAL_APERTURES[fibre_id]['b']
                theta = \
                    self.FIBRE_ELLIPTICAL_APERTURES[fibre_id]['theta']
                apertures.append(
                    EllipticalAperture(position, a, b, theta=theta))
            # perform aperture photometry
            if self.datasets[dataset_name]['dataset_type'] == 'image sequence':
                path = self.expid_to_path_reduced(expid)
            elif self.datasets[dataset_name]['dataset_type'] == 'image stack':
                path = os.path.join(self.DIR_SAVE, 
                                self.datasets[dataset_name]['filename_master'])
            else:
                self.logger.error('Invalid dataset type.')
            phot_table = aper_photometry(
                            path,
                            apertures,
                            logger = self.logger,
                            bkg_sigma = self.BKG_SIGMA,
                            dump_pickle = self.DUMP_PICKLE,
                            clobber = self.CLOBBER)
            # analyse aperture photometry table
            self.analyse_photometry_results(
                dataset_name, expid, phot_table)
    
        # check if all fibres are detected or measured with
        # pre-defined apertures for this expid
        fibres_identified = \
            self.datasets[dataset_name][expid]['fibres_identified']
        if fibres_identified == self.FIBRE_IDS:
            self.logger.info(('All fibres have been measured in ' 
                              + 'exposure {}')
                              .format(expid))
                    
        else:
            fibres_missing = self.FIBRE_IDS - fibres_identified
            self.logger.error('Not all fibres have been measured. \
                Fibres missing: {}'
                .format(fibres_missing))
            return False

    def compile_sequence(self, dataset_name):

        '''
        put source properties from all expids together into lists,
        and store under self.datasets[dataset_name][fibre_id]

        '''

        self.logger.info('Compiling source properties for image sequence {}...'
            .format(dataset_name))

        dataset_type = self.datasets[dataset_name]['dataset_type']

        if dataset_type != 'image sequence':
            self.logger.error('Wrong method to compile photometry, '
                + 'not a sequence.')
            return False

#        # find the first exposure that has at least a fibre identified
#        for expid in self.datasets[dataset_name]['expids']:
#            if len(self.datasets[dataset_name][expid]['fibres_identified']) > 0:
#                fibre_id = self.datasets[dataset_name][expid] \
#                            ['fibres_identified'][0]
#                props = self.datasets[dataset_name][expid][fibre_id]['props']
#        attributes = 
        
        # get list of attributes available
        attributes = []
        for expid in self.datasets[dataset_name]['expids']:
            for fibre_id in self.FIBRE_IDS:
                try:
                    props = self.datasets[dataset_name] \
                                [expid][fibre_id]['props']
                    attributes = attributes + dir(
                        self.datasets[dataset_name][expid][fibre_id]['props'])
                    self.logger.info(('Attributes acquired from '
                                     + 'Exposure {} '
                                     + 'Fibre {} '
                                     + 'source propertiess: {}')
                                     .format(expid, fibre_id, attributes))
                    break
                except:
                    pass
                
        # make sure all table column names are included in attributes as well
        for expid in self.datasets[dataset_name]['expids']:
            for fibre_id in self.FIBRE_IDS:
                try:
                    phot_table = self.datasets[dataset_name] \
                                    [expid][fibre_id]['phot_table']
                    attributes = attributes + phot_table.colnames
                    self.logger.info(('Attributes acquired from '
                                     + 'Exposure {} '
                                     + 'Fibre {} '
                                     + 'aperture photometry table: {}')
                                     .format(expid, fibre_id, attributes))
                    break
                except:
                    pass
                
        # start compiling photometry
        for fibre_id in self.FIBRE_IDS:
            self.datasets[dataset_name][fibre_id] = {}

            for attribute in attributes:
                values = []
                for expid in self.datasets[dataset_name]['expids']:
                    try:
                        # read props for sources identified
                        props = self.datasets[dataset_name][expid] \
                                [fibre_id]['props']
                        prop_type = type(getattr(props, attribute))
                        if prop_type == astropy.units.quantity.Quantity:
                            # read value only, strip units
                            value = getattr(props, attribute).value
                        else:
                            value = getattr(props, attribute)
                    except:
                        # if this fibre is missing, read its photometry table
                        # assuming the attribute is present
                        try:
                            value = self.datasets[dataset_name][expid] \
                                        [fibre_id]['phot_table'][attribute]
                        except:
                            # if attribute is not present in phot table, skip
                            value = np.nan
                    
                    # detailed info on what keys and values are copied
#                    self.logger.info(('Fibre: {}, '
#                        + 'attribute: {}, '
#                        + 'expid: {}, '
#                        + 'value: {}')
#                        .format(fibre_id, attribute, expid, value))
                        
                    values.append(value)
                    
                # add list of values to dictionary
                self.datasets[dataset_name][fibre_id][attribute] = values
        
        # also get date and time of observation from header
        self.datasets[dataset_name]['dates'] = \
            [parse(self.datasets[dataset_name]
             [expid]['header_primary']['DATE-OBS'])
             for expid in self.datasets[dataset_name]['expids']]

#        try:
#            deltara = hdrpri['DELTARA']
#            deltadec = hdrpri['DELTADEC']
#        except:
#            deltara = None
#            deltadec = None
#        self.datasets[dataset_name]['exptime'] = hdrfpc['EXPREQ']
#
#        self.datasets[dataset_name]['DATE-OBS'] = hdrfpc['DATE-OBS']
#        self.datasets[dataset_name]['TIME-OBS'] = hdrfpc['TIME-OBS']
#        self.datasets[dataset_name]['MJD-OBS'] = hdrfpc['MJD-OBS']
#        self.datasets[dataset_name]['DELTARA'] = deltara
#        self.datasets[dataset_name]['DELTADEC'] = deltadec

        # check if all expids are included
        for key in self.datasets[dataset_name][fibre_id].keys():
            if len(self.datasets[dataset_name][fibre_id][key]) \
                != len(self.datasets[dataset_name]['expids']):
                self.logger.error('Keys not compiled completely for'
                    + 'fibre {} propertiess attribute: {}'
                    .format(fibre_id, key))
                return False

    def compile_stacks(self, dataset_names):
        
        '''
        compile photometry results for batch of image stacks
        
        '''
        
        for dataset_name in dataset_names:
            
            self.logger.info('Compiling photometry tables '
                             + 'for image stack {}...'
                             .format(dataset_name))
            
            dataset_type = self.datasets[dataset_name]['dataset_type']
            if dataset_type != 'image stack':
                self.logger.error('Wrong method to compile photometry, '
                    + 'not an image stack.')
                return False
            
            # get list of attributes available
            # attributes = ['flux']
            for expid in self.datasets[dataset_name]['expids']:
                for fibre_id in self.FIBRE_IDS:
                    try:
                        props = self.datasets[dataset_name] \
                                    [expid][fibre_id]['props']
                        attributes = dir(self.datasets[dataset_name] 
                                         [expid][fibre_id]['props'])
                        break
                    except:
                        pass
                    
            # make sure all table column names are included in attributes too
            for expid in self.datasets[dataset_name]['expids']:
                for fibre_id in self.FIBRE_IDS:
                    try:
                        phot_table = self.datasets[dataset_name] \
                                        [expid][fibre_id]['phot_table']
                        attributes.append(phot_table.colnames)
                        break
                    except:
                        pass
                    
            # start compiling photometry
            for fibre_id in self.FIBRE_IDS:
                self.datasets[dataset_name][fibre_id] = {}
    
                for attribute in attributes:
                    # values = []
                    expid = 0
                    # for expid in self.datasets[dataset_name]['expids']:
                    try:
                        # read props for sources identified
                        props = self.datasets[dataset_name][expid] \
                                [fibre_id]['props']
                        prop_type = type(getattr(props, attribute))
                        if prop_type == astropy.units.quantity.Quantity:
                            # read value only, strip units
                            value = getattr(props, attribute).value
                        else:
                            value = getattr(props, attribute)
                    except:
                        # if this fibre is missing, read its photometry table
                        # assuming the attribute is present
                        try:
                            value = self.datasets[dataset_name][expid] \
                                        [fibre_id]['phot_table'][attribute]
                        except:
                            # if attribute is not present in phot table, skip
                            value = np.nan
                    
                    # detailed info on what keys and values are copied
    #                    self.logger.info(('Fibre: {}, '
    #                        + 'attribute: {}, '
    #                        + 'expid: {}, '
    #                        + 'value: {}')
    #                        .format(fibre_id, attribute, expid, value))
                        
                    # values.append(value)
                    # add list of values to dictionary
                    self.datasets[dataset_name][fibre_id][attribute] = value

    def clean_up_props(self, props):

        '''
        remove unwated attributes from SourceProperties object,
        beginnging with '_', which are input data

        dir(props) for a complete list

        '''

        attributes_to_remove = [attribute for attribute in dir(props)
            if attribute[0] == '_' and attribute[0:2] != '__']
        for attribute in attributes_to_remove:
            setattr(props, attribute, None)

        return props

    def dump_dict(self):
        '''
        save dictionary to pickle for later analysis

        '''
        path = os.path.join(self.DIR_SAVE, self.ANALYSIS_NAME + '-pickle.obj')
        self.logger.info('Dumping dictionary to pickle: {}...'.format(path))
        file = open(path, 'wb') # 'wb' means write and binay, i.e. overwriting
        pickle.dump(self.datasets, file)
        self.logger.info('Dictionary dumped to: {}'.format(path))

    def load_dict(self):

        path = os.path.join(self.DIR_SAVE, self.ANALYSIS_NAME + '-pickle.obj')
        self.logger.info('Loading dictionary from pickle: {}...'.format(path))
        file = open(path, 'rb') # read binary mode
        self.datasets = pickle.load(file)
        self.logger.info('Dictionary loaded from: {}'.format(path))

        # return True

    def find_array_centroid(self, array):

        '''
        simple function that calculates weighted centroid of an array
        origin is assumed to be the first element at top left corner

        '''

        if type(array) == np.ndarray:
            if array.ndim == 2:

                array = np.float64(array)
                x = range(0, array.shape[0])
                y = range(0, array.shape[1])
                (xarray, yarray) = np.meshgrid(x, y)

                xcentroid = (xarray*array).sum() / array.sum()
                ycentroid = (yarray*array).sum() / array.sum()

                return (xcentroid, ycentroid)
        else:
            print('Input must be a numpy array.')
            return False
            
    def plot_flats(self, dataset_names):
        
        '''
        plot flux vs exptime for dome flats
        input is a list of dataset names
        
        '''
        
        exptime =  [self.datasets[dataset_name]['exptime']
                    for dataset_name in dataset_names]
        source_sum = {}
        source_sum_err = {}
        flux = {}
        flux_err = {}

        for fibre_id in self.FIBRE_IDS:
            source_sum[fibre_id] = []
            source_sum_err[fibre_id] = []
            flux[fibre_id] = []
            flux_err[fibre_id] = []
            for dataset_name in dataset_names:
                try:
                    s = self.datasets[dataset_name] \
                            [0][fibre_id]['props']['source_sum']
                    se = self.datasets[dataset_name] \
                            [0][fibre_id]['props']['source_sum_err']
                    f = self.datasets[dataset_name] \
                            [0][fibre_id]['props']['flux']
                    fe = self.datasets[dataset_name] \
                            [0][fibre_id]['props']['flux_err']
                except:
                    self.logger.error(('Source proerties from '
                        + 'segmentation photometry unavailable for '
                        + 'flat dataset {}')
                        .format(dataset_name))
                    # we just get empty lists
                    
                source_sum[fibre_id].append(s)
                source_sum_err[fibre_id].append(se)
                flux[fibre_id].append(f)
                flux_err[fibre_id].append(fe)
        
        # calculate relative fibre throughput
        thrpt_rel = {}
        thrpt_rel_err = {}
        for fibre_id in self.FIBRE_IDS:
            thrpt_rel[fibre_id] = np.divide(flux[fibre_id], 
                                            flux[self.FIBRE_BEST])
            thrpt_rel_err[fibre_id] = np.std(thrpt_rel[fibre_id])
            thrpt_rel[fibre_id] = np.mean(thrpt_rel[fibre_id])

        self.logger.info('Creating flat linearity plots for datasets {}...'
                         .format(dataset_names))
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12, 6))
        # specify properties for three fibres
        colours = ['r', 'g', 'b']
        markers = ['-o', '-x', '-*']
#        legendmap1 = {}
#        legendmap2 = {}
    
        for i, fibre_id in enumerate(sorted(self.FIBRE_IDS)):
#            [line1] = ax1.plot(exptime, source_sum[fibre_id], 
#                c=colours[i], marker=markers[i], label='Fibre #'+fibre_id)
            ax1.errorbar(exptime, source_sum[fibre_id], 
                         yerr=source_sum_err[fibre_id], 
                         fmt=markers[i], markersize=4,
                         c=colours[i], ecolor=colours[i],
                         label='Fibre #'+fibre_id)
#            [line2] = ax2.plot(exptime, flux[fibre_id],
#                c=colours[i], marker=markers[i], label='Fibre #'+fibre_id)
            ax2.errorbar(exptime, flux[fibre_id],
                         yerr=flux_err[fibre_id], 
                         fmt=markers[i], markersize=4,
                         c=colours[i], ecolor=colours[i],
                         label='Fibre #'+fibre_id)
#            legendmap1[line1] = HandlerLine2D(numpoints=3)
#            legendmap2[line2] = HandlerLine2D(numpoints=3)
        ax1.grid()
        ax1.set_xlabel(r'Exposure Time $t_\mathrm{exp}/\mathrm{s}$')
        ax1.set_xlim(0, np.amax(exptime)+1)
        ax1.set_ylabel(r'Fibre Intensity Sum $\sum(I-B)/\mathrm{ADU}$')
        yformatter = matplotlib.ticker.ScalarFormatter()
        yformatter.set_powerlimits((1, 4))
        ax1.set_ylim(0, 5e7)
        ax1.yaxis.set_major_formatter(yformatter)
        ax1.legend(loc=2)
        ax1.set_title('Fibre Throughput and FPC Linearity')
        #ax1.legend(handler_map=legendmap1, loc=2)
        
        
        ax2.set_title('Fibre Throughput and FPC Linearity')
        ax2.set_xlabel(r'Exposure Time $t_\mathrm{exp}/\mathrm{s}$')
        ax2.set_xlim(0, np.amax(exptime)+0.5)
        ax2.set_ylabel(r'Exptime Averaged Flux '
                       + r'$f/\mathrm{ADU} \cdot \mathrm{s}^{-1}$')
        ax2.set_xlim(0, np.amax(exptime)+1)
        ax2.set_ylim(0, 5e6)
        ax2.yaxis.set_major_formatter(yformatter)
        ax2.legend(loc=2)
        #ax2.legend(handler_map=legendmap2, loc=2)
        ax2.grid()
        
        # add text for fibre throughput
        textstr = 'Relative Throughput'
        for fibre_id in sorted(self.FIBRE_IDS):
            textstr = textstr + '\n#{0}: {1:.3f} $\pm$ {2:.3f}'.format(
                        fibre_id, 
                        thrpt_rel[fibre_id], 
                        thrpt_rel_err[fibre_id])
        textbbox = {'boxstyle': 'square',
                    'facecolor': 'white',
                    'alpha': 0.5}
        ax1.text(0.96, 0.04, textstr, fontsize=12, linespacing=1.5,
                 transform=ax1.transAxes, bbox=textbbox,
                 ha='right', va='bottom')
        ax2.text(0.96, 0.04, textstr, fontsize=12, linespacing=1.5,
                 transform=ax2.transAxes, bbox=textbbox,
                 ha='right', va='bottom')
        
        # padding to avoid overlap between axes
        plt.tight_layout()
        # save
        fig.savefig(os.path.join(
                        self.DIR_SAVE, 
                        self.ANALYSIS_NAME+'-throughput_linearity.png'), 
                    dpi=600)
        pp = PdfPages(os.path.join(
                        self.DIR_SAVE, 
                        self.ANALYSIS_NAME+'-throughput_linearity.pdf'))
        pp.savefig(fig, dpi=600)
        pp.close()
        
        plt.close('all')
        self.logger.info('Flat linearity plots saved to {}'
                 .format(self.DIR_SAVE))
        
    def plot_stability(self, dataset_name):
        
        tile = self.datasets[dataset_name]['tile_id']
        dates = self.datasets[dataset_name]['dates']
        x = date2num(dates)
        flux = {}
        flux_err = {}
        flux_mean = {}
        flux_mean_err = {}
        for fibre_id in self.FIBRE_IDS:
            flux[fibre_id] = self.datasets[dataset_name][fibre_id]['flux']
            flux_err[fibre_id] = self.datasets[dataset_name] \
                                    [fibre_id]['flux_err']
            flux_mean[fibre_id] = np.mean(flux[fibre_id])
            flux_mean_err[fibre_id] = np.std(flux[fibre_id])
        flux_max = np.amax([np.amax(flux[key]) for key in flux.keys()])
        self.logger.info('Creating stability plots for dataset {}...'
              .format(dataset_name))
        (fig, ax) = plt.subplots(figsize = (8, 6))
        colours = ['r', 'g', 'b']
        markers = ['-o', '-x', '-*'] #markers = [None, None, None] # 
        for i, fibre_id in enumerate(sorted(self.FIBRE_IDS)):
#            [line] = ax.plot(x, flux[fibre_id], 
#                c=colours[i], marker=markers[i], label='Fibre #'+fibre_id)
            ax.errorbar(x, flux[fibre_id], yerr=flux_err[fibre_id],
                        fmt=markers[i], markersize=4,
                        c=colours[i], ecolor=colours[i], 
                        label='Fibre #'+fibre_id)
        ax.grid()
        ax.set_title('Target {} Stability Test'.format(tile))
        ax.set_xlabel(r'Datetime $t/\mathrm{s}$')
        ax.xaxis.set_major_formatter(DateFormatter('%m/%d %H:%M:%S'))
        yformatter = matplotlib.ticker.ScalarFormatter()
        yformatter.set_powerlimits((1, 4))
        ax.yaxis.set_major_formatter(yformatter)
        ax.set_ylabel(r'Exptime Averaged Flux '
                       + r'$f/\mathrm{ADU} \cdot \mathrm{s}^{-1}$')
        ax.set_ylim(0, flux_max/0.7)
        ax.legend(loc=2)
        fig.autofmt_xdate()
        
        # add text for fibre throughput
        textstr = r'Mean Flux $\bar{f}/10^6 \mathrm{ADU} \cdot s^{-1}$'
        for fibre_id in sorted(self.FIBRE_IDS):
            textstr = textstr + '\n#{0}: {1:.3f} $\pm$ {2:.3f}'.format(
                        fibre_id, 
                        flux_mean[fibre_id]/1e6, 
                        flux_mean_err[fibre_id]/1e6)
        textbbox = {'boxstyle': 'square',
                    'facecolor': 'white',
                    'alpha': 0.5}
        ax.text(0.96, 0.04, textstr, fontsize=12, linespacing=1.5,
                 transform=ax.transAxes, bbox =textbbox,
                 ha='right', va='bottom')

        # padding to avoid overlap between axes
        plt.tight_layout()
        # save
        fig.savefig(os.path.join(self.DIR_SAVE, 
                                 dataset_name+'-stability_plot.png'), 
                    dpi=600)
        pp = PdfPages(os.path.join(self.DIR_SAVE, 
                                   dataset_name+'-stability_plot.pdf'))
        pp.savefig(fig, dpi=600)
        pp.close()
        
        plt.close('all')
        self.logger.info('Stability plots saved to {}'
                         .format(self.DIR_SAVE))

    def plot_dither_2d(self, dataset_name):
        
        self.logger.info('Creating 2D dither plots for dataset {}...'
              .format(dataset_name))

        tile = self.datasets[dataset_name]['tile_id']
        grid = self.datasets[dataset_name]['dither_grid']
        step = self.datasets[dataset_name]['dither_step']
        dither_mode = self.datasets[dataset_name]['dither_mode']

        if grid == (5, 5):
            if dither_mode == 'telescope':
                grid_index = self.GRID_INDEX_TEL
            if dither_mode == 'positioner':
                grid_index = self.GRID_INDEX_POS
                
        # determine x, y max for the coordinates
        xmax = (grid[0]-1)/2*step
        ymax = (grid[1]-1)/2*step
        # construct x and y labels
        x = np.linspace(-xmax, xmax, grid[0], endpoint=True)
        y = np.linspace(-ymax, ymax, grid[0], endpoint=True)
        xlabels = [str(value) for value in x]
        ylabels = list(reversed([str(value) for value in y])) # origin on top
        # x and y label positioners on matrix plot, must be from -2 to 2
        x = np.linspace(-(grid[0]-1)/2, (grid[0]-1)/2, grid[0], endpoint=True)
        y = np.linspace(-(grid[0]-1)/2, (grid[0]-1)/2, grid[0], endpoint=True)
                                
        for fibre_id in self.FIBRE_IDS:
            path_save_prefix = os.path.join(self.DIR_SAVE,
                                dataset_name
                                + '-' +dither_mode + '_dither'
                                + '-fibre_#' + fibre_id)
            # read in z value list
            z_list = np.int64(self.datasets[dataset_name][fibre_id]['flux'])
            expids_list = np.int64(self.datasets[dataset_name]['expids'])
            # fill in z array and put z values in their right place
            z = np.int64(np.empty(grid))
            expids = np.int64(np.empty(grid))
            for i in range(grid[0]):
                for j in range(grid[1]):
                    index = grid_index[i,j]
                    z[i,j] = z_list[index]
                    expids[i,j] = expids_list[index]

            # python matrices are row major 
            # but plot images are effectively column major
            zlabels = np.transpose(z)
            expids = np.transpose(expids)
            
            # calculates intensity-weighted centroid
            (xcentroid, ycentroid) = self.find_array_centroid(z)
            self.logger.info('Weighted centroid in array calculated as: {}'
                             .format((xcentroid, ycentroid)))
            # plot grid
            fig1, ax1 = plt.subplots(figsize=(8, 6))
            mat = ax1.matshow(z, cmap=plt.cm.coolwarm)
            # plot centroid
            ax1.scatter(xcentroid, ycentroid, 
                     marker='+', c='green', s=600, alpha=1)
            for i in range(grid[0]):
                for j in range(grid[1]):
                    ax1.text(x[i]+(grid[0]-1)/2, y[j]+(grid[1]-1)/2, 
                             '#' + str(expids[i,j]) + 
                             '\n' + str(zlabels[i,j]),
                             size=8, ha='center',va='center')
            ax1.xaxis.set_ticks_position('bottom')
            #ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax1.set_xticks(x + (grid[0]-1)/2)
            ax1.set_xticklabels(xlabels)
            ax1.set_xlabel(r'Relative RA Offset $\alpha/$"')
            ax1.set_yticks(y + (grid[0]-1)/2)
            ax1.set_yticklabels(ylabels)
            ax1.set_ylabel(r'Relative DEC offset $\delta/$"')
            ax1.set_title(r'Target {} Fibre #{} {} Dither Test'
                          .format(tile, fibre_id, dither_mode.title()))
            colorbar = fig1.colorbar(mat, aspect=13)
            colorbar.set_label(r'Exptime Averaged Flux '
                               + r'$f/\mathrm{ADU}\cdot\mathrm{s}^{-1}$')
            colorbar.formatter.set_powerlimits((0, 4))
            fig1.tight_layout()
            # ax1.xaxis.set_ticks_position('bottom')
            
            # save fig1
            fig1.savefig(path_save_prefix
                         + '-plot2d_matrix.png',
                         dpi=800, pad_inches=0)
            pp1 = PdfPages(path_save_prefix
                           + '-plot2d_matrix.pdf')
            pp1.savefig(fig1, dpi=800)
            pp1.close()
    
            # plot 2d interpolated heat map
            for interp in ['bilinear', 'spline36', 'sinc', 'lanczos']:
                fig2, ax2 = plt.subplots(figsize=(8, 6))
                im = ax2.imshow(z, cmap=plt.cm.coolwarm, 
                                interpolation=interp)
                ax2.scatter(xcentroid, ycentroid, 
                         marker='+', c='green', s=600, alpha=1)
                for i in range(grid[0]):
                    for j in range(grid[1]):
                        ax2.text(x[i]+(grid[0]-1)/2, y[j]+(grid[1]-1)/2, 
                                 '#' + str(expids[i,j]) + 
                                 '\n' + str(zlabels[i,j]),
                                 size=8, ha='center',va='center')
                ax2.xaxis.set_ticks_position('bottom')
                #ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                ax2.set_xticks(x + (grid[0]-1)/2)
                ax2.set_xticklabels(xlabels)
                ax2.set_xlabel(r'Relative RA Offset $\alpha/$"')
                ax2.set_yticks(y + (grid[0]-1)/2)
                ax2.set_yticklabels(ylabels)
                ax2.set_ylabel(r'Relative DEC offset $\delta/$"')
                ax2.set_title(r'Target {} Fibre #{} {} Dither Test'
                              .format(tile, fibre_id, dither_mode.title()))
                colorbar = fig2.colorbar(im, aspect=13)
                colorbar.set_label(r'Exptime Averaged Flux '
                            + r'$f/\mathrm{ADU}\cdot\mathrm{s}^{-1}$')
                colorbar.formatter.set_powerlimits((0, 4))
                fig2.tight_layout()
                
                # save fig2
                fig2.savefig(path_save_prefix 
                             + '-plot2d_' + interp + '.png',
                             dpi=800, pad_inches=0)
                pp2 = PdfPages(path_save_prefix 
                               + '-plot2d_' + interp + '.pdf',)
                pp2.savefig(fig2, dpi=800)
                pp2.close()
                
                plt.close('all')
                
        self.logger.info('2D dither plots saved to {}...'
                         .format(self.DIR_SAVE))

    def plot_dither_3d(self, dataset_name):
        
        self.logger.info('Creating 3D dither plots for dataset {}...'
                         .format(dataset_name))
        tile = self.datasets[dataset_name]['tile_id']
        grid = self.datasets[dataset_name]['dither_grid']
        step = self.datasets[dataset_name]['dither_step']
        dither_mode = self.datasets[dataset_name]['dither_mode']
        # radius of fibre circles on plot
        circle_radius = self.FIBRE_RADIUS/self.PLATE_SCALE

        if grid == (5, 5):
            if dither_mode == 'telescope':
                grid_index = self.GRID_INDEX_TEL
            if dither_mode == 'positioner':
                grid_index = self.GRID_INDEX_POS
        else:
            self.logger.error('Grid indices undefined for: {}'.format(grid))

        # vertically flip grid index array becaus python index starts from top
        # contrary to natural orientation of coordinate system
        grid_index = np.flipud(grid_index)

        # determine x, y max for the coordinates
        xmax = (grid[0]-1)/2*step
        ymax = (grid[1]-1)/2*step

        # construct x and y coordinates
        x = np.linspace(-xmax, xmax, grid[0], endpoint=True)
        y = np.linspace(-ymax, ymax, grid[0], endpoint=True)
        (x, y) = np.meshgrid(x, y)
        z_circle = self.Z_CIRCLE
        expids = np.int64(np.empty(grid))
        # x and y coordinates for finer interpolated plot
        xnew = np.linspace(-xmax, xmax, grid[0]*10, endpoint=True)
        ynew = np.linspace(-ymax, ymax, grid[0]*10, endpoint=True)
        (xnew, ynew) = np.meshgrid(xnew, ynew)

        for fibre_id in self.FIBRE_IDS:

            # read in z value list
            z_list = np.int64(self.datasets[dataset_name][fibre_id]['flux'])
            expids_list = np.int64(self.datasets[dataset_name]['expids'])

            # fill in z array and put z values in their right place
            z = np.int64(np.empty(grid))
            for i in range(grid[0]):
                for j in range(grid[1]):
                    index = grid_index[i,j]
                    z[i,j] = z_list[index]
                    expids[i,j] = expids_list[index]

            # calculates intensity-weighted centroid
            (xcentroid, ycentroid) = self.find_array_centroid(np.flipud(z))
            # shift the origin to centre of matrix, and then rescale w/ step
            xcentroid = (xcentroid - (grid[0]-1)/2) * step
            ycentroid = ((grid[0]-1)/2 - ycentroid) * step
            centroid_l1 = [[xcentroid-step/2, xcentroid+step/2],
                           [ycentroid, ycentroid,],
                           [z_circle, z_circle]]
            centroid_l2 = [[xcentroid, xcentroid],
                           [ycentroid-step/2, ycentroid+step/2],
                           [z_circle, z_circle]]

            # define ticks (original data) to be used for interpolation
            tck = interpolate.bisplrep(x, y, z, s=0)
            # create interpolated z data
            # take transpose because interpolate routine assumes 0, 0 at
            # upper left corner, not botttom left
            znew = np.transpose(interpolate.bisplev(xnew[0,:], ynew[:,0], tck))
            
            for elevation in [20, 40, 60, 80]:
                
                # === plot z data ===
                fig1 = plt.figure(figsize=(8, 6))
                ax1 = fig1.gca(projection='3d')
                # add circles representing fibres and texts inside circles
                circles = []
                for i in range(grid[0]):
                    for j in range(grid[1]):
                        circle_centre = (x[i,j], y[i,j])
                        circles.append(Circle(circle_centre, circle_radius,
                            linewidth=0.5, edgecolor='lime', facecolor='none',
                            alpha=0.7, antialiased=True, zorder=-1))
                        ax1.text(x[i,j], y[i,j], z_circle,
                                 '#' + str(expids[i,j]) + '\n' + str(z[i,j]),
                                  zdir='x', size=4, zorder=0,
                                  ha='center',va='center')
                for circle in circles:
                    ax1.add_patch(circle)
                    art3d.pathpatch_2d_to_3d(circle, z=z_circle, zdir='z')
                for cent_line in [centroid_l1, centroid_l2]:
                    ax1.plot(cent_line[0], cent_line[1], cent_line[2],
                         color='red', linewidth=1, alpha=0.6)
                # add surface plot on top of circles
                surf = ax1.plot_surface(x, y, z, zorder = 2,
                                rstride=1, cstride=1, cmap=plt.cm.coolwarm,
                                linewidth=0.1, antialiased=True, alpha=0.6)
                ax1.scatter(x, y, z, zorder=3, alpha=0.6,
                            marker='o', s=7, c='red', edgecolors='face')
                ax1.set_xlim(-xmax - circle_radius, xmax + circle_radius)
                ax1.set_ylim(-ymax - circle_radius, ymax + circle_radius)
                ax1.set_zlim(z_circle, z.max())
                zformatter = matplotlib.ticker.ScalarFormatter()
                zformatter.set_powerlimits((1, 4))
                ax1.zaxis.set_major_formatter(zformatter)
                # ax2.scatter(xcentroid, ycentroid, z_circle, zorder=1,
                #             marker='+', s=50, c='red', alpha=1, zdir='z')
                ax1.text(-xmax-step/2, ymax-step/2, z.max()*0.8,
                         'Weighted Centroid \n ({0:.2f}, {1:.2f})'
                         .format(xcentroid, ycentroid),
                         zorder=0)
                ax1.set_title('Target {} Fibre #{} {} Dither Test'
                              .format(tile, fibre_id, dither_mode.title()))
                ax1.set_xlabel(r'Relative RA Offset $\alpha/$"')
                ax1.set_ylabel(r'Relative DEC offset $\delta/$"')
                ax1.set_zlabel(r'Exptime Averaged Flux '
                               + r'$f/\mathrm{ADU}\cdot\mathrm{s}^{-1}$',
                               labelpad=18)
                ax1.view_init(elev=elevation, azim=280)
                colorbar = fig1.colorbar(surf, shrink=0.6, aspect=13)
                colorbar.set_label(r'Exptime Averaged Flux '
                                   + r'$f/\mathrm{ADU}\cdot\mathrm{s}^{-1}$')
                colorbar.formatter.set_powerlimits((0, 4))
                
                fig1.tight_layout()
                fig1.show()
    
                # === plot interplated data ===
                fig2 = plt.figure(figsize=(8, 6))
                ax2 = fig2.gca(projection='3d')
                # add circles representing fibres and texts inside circles
                circles = []
                for i in range(grid[0]):
                    for j in range(grid[1]):
                        circle_centre = (x[i,j], y[i,j])
                        circles.append(Circle(circle_centre, circle_radius,
                            linewidth=0.5, edgecolor='lime', facecolor='none',
                            alpha=0.7, antialiased=True, zorder=-1))
                        ax2.text(x[i,j], y[i,j], z_circle,
                                 '#'+ str(expids[i,j]) + '\n' + str(z[i,j]),
                                  zdir='x', size=4, zorder=0,
                                  ha='center',va='center')
                for circle in circles:
                    ax2.add_patch(circle)
                    art3d.pathpatch_2d_to_3d(circle, z=z_circle, zdir='z')
                for cent_line in [centroid_l1, centroid_l2]:
                    ax2.plot(cent_line[0], cent_line[1], cent_line[2],
                         color='red', linewidth=1, alpha=0.6)
                # add surface plot on top of circles
                surf = ax2.plot_surface(xnew, ynew, znew, zorder=2,
                                rstride=1, cstride=1, cmap=plt.cm.coolwarm,
                                linewidth=0.1, antialiased=True, alpha=0.6)
                ax2.scatter(x, y, z, zorder=3,
                            marker='o', s=7, c='red', edgecolors='face', 
                            alpha=0.6)
                # ax2.scatter(xcentroid, ycentroid, z_circle, zorder=1,
                #             marker='+', s=50, c='red', alpha=1, zdir='z')
                ax2.text(-xmax-step/2, ymax-step/2, z.max()*0.8,
                         'Weighted Centroid \n ({0:.2f}, {1:.2f})'
                         .format(xcentroid, ycentroid),
                         zorder=0)
                ax2.set_xlim(-xmax - circle_radius, xmax + circle_radius)
                ax2.set_ylim(-ymax - circle_radius, ymax + circle_radius)
                ax2.set_zlim(z_circle, z.max())
                ax2.zaxis.set_major_formatter(zformatter)
                ax2.set_title('Target {} Fibre #{} {} Dither Test'
                    .format(tile, fibre_id, dither_mode.title()))
                ax2.set_xlabel(r'Relative RA Offset $\alpha/$"')
                ax2.set_ylabel(r'Relative DEC offset $\delta/$"')
                ax2.set_zlabel(r'Exptime Averaged Flux '
                            + r'$f/\mathrm{ADU} \cdot \mathrm{s}^{-1}$',
                            labelpad=18)
                ax2.view_init(elev=elevation, azim=280)
                fig2.tight_layout()
                #ax.zaxis.set_major_locator(LogLocator(10))
                colorbar = fig2.colorbar(surf, shrink=0.6, aspect=13)
                colorbar.set_label(r'Exptime Averaged Flux '
                                   + r'$f/\mathrm{ADU}\cdot\mathrm{s}^{-1}$')
                colorbar.formatter.set_powerlimits((0, 4))
                fig2.show()

                # === save plots ===
                path_save_prefix = os.path.join(self.DIR_SAVE,
                                                dataset_name
                                                + '-' +dither_mode + '_dither'
                                                + '-fibre_#' + fibre_id)
    
                fig1.savefig(path_save_prefix 
                             + '-plot3d_scatter_elev' 
                             + str(elevation) +'.png',
                             dpi=800, pad_inches=0)
                fig2.savefig(path_save_prefix 
                             + '-plot3d_surface_interp_elev' 
                             + str(elevation) +'.png',
                             dpi=800, pad_inches=0)
    
                pp1 = PdfPages(path_save_prefix 
                               + '-plot3d_scatter_elev' 
                               + str(elevation) +'.pdf')
                pp1.savefig(fig1, dpi=800)
                pp1.close()
    
                pp2 = PdfPages(path_save_prefix 
                               + '-plot3d_surface_interp_elev' 
                               + str(elevation) +'.pdf')
                pp2.savefig(fig2, dpi=800)
                pp2.close()
                
                plt.close('all')
        
        self.logger.info('3D dither plots saved to {}'
                         .format(self.DIR_SAVE))
    
    def plot_dither_focus(self, datasets):
        
        '''
        3D density distribution of a few dither sequences that differ by focus
        
        '''
        
        self.logger.info('Creating 3D dither focus plots for datasets {}...'
                         .format(repr(datasets)))
        # assuming the following parameters used are the same for all datasets
        dataset_name = datasets[0]
        tile = self.datasets[dataset_name]['tile_id']
        grid = self.datasets[dataset_name]['dither_grid']
        step = self.datasets[dataset_name]['dither_step']
        dither_mode = self.datasets[dataset_name]['dither_mode']

        if grid == (5, 5):
            if dither_mode == 'telescope':
                grid_index = self.GRID_INDEX_TEL
            if dither_mode == 'positioner':
                grid_index = self.GRID_INDEX_POS
        else:
            self.logger.error('Grid indices undefined for: {}'.format(grid))
            
        # determine x, y max for the coordinates
        xmax = (grid[0]-1)/2*step
        ymax = (grid[1]-1)/2*step
        
        # vertically flip grid index array becaus python index starts from top
        # contrary to natural orientation of coordinate system
        grid_index = np.flipud(grid_index)
        x = np.linspace(-xmax, xmax, grid[0], endpoint=True)
        y = np.linspace(-ymax, ymax, grid[0], endpoint=True)
        z = [self.datasets[dataset_name]['focus'] for dataset_name in datasets]
        flux = {}
        (x, y, z) = np.meshgrid(x, y, z)
        for fibre_id in self.FIBRE_IDS:
            # for each fibre, flux is a 3d array, with focus as 3rd dimension
            flux[fibre_id] = np.empty([grid[0], grid[1], len(datasets)])
            for k, dataset_name in enumerate(datasets):
                # read in flux as a list from 25 expids
                flux_list = np.int64(self.datasets[dataset_name]
                                     [fibre_id]['flux'])
                # for each slice, i.e. each focus k, fill 2d entries
                for i in range(grid[0]):
                    for j in range(grid[1]):
                        index = grid_index[i,j]
                        flux[fibre_id][i, j, k] = flux_list[index]
            # start plotting

    def plot_dither_focus_marching_cubes(self, datasets):
        
        '''
        3D density distribution of a few dither sequences that differ by focus
        failed attempt -
        marching_cubes doesn't interpolate to create a smooth surface
        
        '''
        
        self.logger.info('Creating 3D dither focus plots for datasets {}...'
                         .format(repr(datasets)))
        # assuming the following parameters used are the same for all datasets
        dataset_name = datasets[0]
        tile = self.datasets[dataset_name]['tile_id']
        grid = self.datasets[dataset_name]['dither_grid']
        step = self.datasets[dataset_name]['dither_step']
        dither_mode = self.datasets[dataset_name]['dither_mode']

        if grid == (5, 5):
            if dither_mode == 'telescope':
                grid_index = self.GRID_INDEX_TEL
            if dither_mode == 'positioner':
                grid_index = self.GRID_INDEX_POS
        else:
            self.logger.error('Grid indices undefined for: {}'.format(grid))
            
        # determine x, y max for the coordinates
        xmax = (grid[0]-1)/2*step
        ymax = (grid[1]-1)/2*step
        
        # vertically flip grid index array becaus python index starts from top
        # contrary to natural orientation of coordinate system
        grid_index = np.flipud(grid_index)
        x = np.linspace(-xmax, xmax, grid[0], endpoint=True)
        y = np.linspace(-ymax, ymax, grid[0], endpoint=True)
        z = [self.datasets[dataset_name]['focus'] for dataset_name in datasets]
        flux = {}
        vertices = {}
        faces = {}
        for fibre_id in self.FIBRE_IDS:
            
            # for each fibre, flux is a 3d array, with focus as 3rd dimension
            flux[fibre_id] = np.empty([grid[0], grid[1], len(datasets)])
            for k, dataset_name in enumerate(datasets):
                # read in flux as a list from 25 expids
                flux_list = np.int64(self.datasets[dataset_name]
                                     [fibre_id]['flux'])
                # for each slice, i.e. each focus k, fill 2d entries
                for i in range(grid[0]):
                    for j in range(grid[1]):
                        index = grid_index[i,j]
                        flux[fibre_id][i, j, k] = flux_list[index]
            # calculate isosurface 
            (vertices[fibre_id], faces[fibre_id]) = marching_cubes(
                flux[fibre_id],
                level = np.mean(flux[fibre_id]),
                spacing = (0.1, 0.1, 0.1))
            # plot
            fig = plt.figure(figsize=(8, 6))
            ax = fig.gca(projection='3d')
            ax.plot_trisurf(vertices[fibre_id][:, 0], vertices[fibre_id][:,1], 
                            faces[fibre_id], vertices[fibre_id][:, 2],
                            cmap=plt.cm.spectral, lw=0.1)
            ax.set_title('Target {} Fibre #{} {} Dither+Focus Test'
                         .format(tile, fibre_id, dither_mode.title()))
            ax.set_xlabel(r'Relative RA Offset $\alpha/$"')
            ax.set_ylabel(r'Relative DEC offset $\delta/$"')
            ax.set_zlabel(r'Telescope Focus')
            fig.savefig(os.path.join(self.DIR_SAVE, 
                                     'dither_focus-'+fibre_id+'-3d_plot.png'), 
                        dpi=1000)
            pp = PdfPages(os.path.join(self.DIR_SAVE, 
                                     'dither_focus-'+fibre_id+'-3d_plot.pdf'))
            pp.savefig(fig, dpi=1000)
            pp.close()
            plt.show
            
        #plt.close('all')
            
#%% if script is run as main
if __name__ == '__main__':
    
    import os
    import glob
    import multiprocessing
    
    def runfile(file):
        exec(open(file).read())
    
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    filelist = glob.glob('pd_analysis-*.py')
    print('Multiprocessing scripts {}'.format(repr(filelist)))
    with multiprocessing.Pool(8) as pool:
        pool.map(runfile, filelist)
