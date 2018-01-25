#%% -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 22:49:36 2016

@author: Duan Yutong (dyt@physics.bu.edu)

Reduce: enter relative path of parent folder

"""

import os
import sys
import functools
import glob
import stat
import numpy as np
import scipy
from astropy.io import fits
import paramiko
import pyqtgraph as pg
from PyQt4 import QtCore, QtGui, uic
#from sklearn.metrics import mean_squared_error

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

from photutils.background import Background
from photutils import detect_sources
from photutils.utils import random_cmap
from photutils import source_properties, properties_table

#%% Settings

autosync = False
printEnabled = True
writeEnabled = False

DEVICE_LIST = ['sti', 'fvc', 'fpc']
# hosts
DOMAIN_STI = 'desisti.kpno.noao.edu'
DOMAIN_FVC = 'desifvc.kpno.noao.edu'
DOMAIN_FPC = 'desirpi1.kpno.noao.edu'
USERNAME = 'msdos'
PASSWORD = 'MS-d0s'

# paths
DIR_REMOTE_STI = r'/home/msdos/SBIG/'
DIR_REMOTE_FVC = r'/data/images/fvc/'
DIR_REMOTE_FPC = r'/data/images/'
DIR_LOCAL_STI = r'/Volumes/data/images/sti/'
DIR_LOCAL_FVC = r'/Volumes/data/images/fvc/'
DIR_LOCAL_FPC = r'/Volumes/data/images/fpc/'
#dirLocal = os.path.join(os.path.expanduser("~"), 'Downloads', 'data', 'fpc')
#id_test = '20160825_0001'
#id_dataset = '1'
TEXT_BROWSER_FONT = 'Courier'
TEXT_BROWSER_FONT_SIZE = 13

#%% Qt setup

pg.mkQApp()
uiPath = os.path.dirname(os.path.abspath(__file__))
uiFile = os.path.join(uiPath, 'fpc_gui.ui')
[WindowTemplate, TemplateBaseClass] = uic.loadUiType(uiFile)

#%% function isdir

def isdir(sftp, path):
    try:
        return stat.S_ISDIR(sftp.stat(path).st_mode)
    except IOError:
        #Path does not exist, so by definition not a directory
        return False

#%% function open_sftp

def open_sftp(domain, user, pw):

    '''
    Given ssh info, return an SFTP session

    '''
    print('SSHing to {} ...'.format(domain))
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(domain, username=user, password=pw)

    print('Establishing SFTP...')
    sftp = ssh.open_sftp()

    return sftp

#%% function sftp_get_all

def sftp_get_all(sftp, path_remote, path_local):

    '''
    Transfer all files recursively via SFTP

    '''

    if isdir(sftp, path_remote):
        # if remote path is a directory
        # check that remote path ends with / to ensure os seperator
        # inserted is consistent across platform
        if path_remote[-1] != '/':
            path_remote = path_remote+'/'
        # create lists of names and full paths
        names = sftp.listdir(path_remote)
        paths = []
        for name in names:
            paths.append(os.path.join(path_remote, name))
        # copy each remote path in paths
        for i in range(len(paths)):
            name = names[i]
            path = paths[i]

            # if system is windows, : in filename is not allowed
            # replace with underscore
            if os.path.sep == '\\':
                name_local = name.replace(':', '_')
            elif os.path.sep == '/':
                name_local = name

            destination = os.path.join(path_local, name_local)
            print('Syncing {}'.format(path))
            if isdir(sftp, path):
                #print('{} is a directory'.format(name))
                # create local directory
                if not os.path.exists(destination):	# check and create
                    os.makedirs(destination)
            sftp_get_all(sftp, path, destination)
    else:
        # remote path is a file
        if os.path.isfile(path_local):
            # file with same name already exists
            pass
        else:
            print('Transferring file {} ...'.format(path_remote))
            sftp.get(path_remote, path_local,
                callback = sftp_progress
                )

#%% function sftp_new_data_checker

def sftp_new_data_checker(sftp, path_remote, path_local):

    '''
    check for new data on controllers recursively and compare to data server
    '''

    if isdir(sftp, path_remote):
        # if remote path is a directory
        # check that remote path ends with / to ensure compatibility
        # with linux file server
        if path_remote[-1] != '/':
            path_remote = path_remote+'/'
        # create lists of names and full paths
        names = sftp.listdir(path_remote)
        paths = []
        for name in names:
            paths.append(os.path.join(path_remote, name))
        # copy each remote path in paths
        for i in range(len(paths)):
            name = names[i]
            path = paths[i]

            # if system is windows, : in filename is not allowed
            # replace with underscore
            if os.path.sep == '\\':
                name_local = name.replace(':', '_')
            elif os.path.sep == '/':
                name_local = name

            destination = os.path.join(path_local, name_local)
            print('Checking {}'.format(path))
            if isdir(sftp, path):
                #print('{} is a directory'.format(name))
                if not os.path.exists(destination):
                    # folder does not exist, new data available
                    return True
                else:
                    return sftp_new_data_checker(sftp, path, destination)
    else:
        # remote path is a file
        if os.path.isfile(path_local):
            # file with same name already exists
            return False
        else:
            # file does not already exist
            return True

#%% function sftp_progress

def sftp_progress(transferred, total):

    '''
    Show SFTP transfer progress

    '''
    percentage = transferred/total*100
    total_mb = total/(1024**2)
    print('{0:.2f}%, Total {1:.2f} MB \r'.format(percentage, total_mb))
#        self.ui.statusbar.showMessage('{0:.2f}%, Total {1:.2f} MB \r'.format(percentage, total_mb))

#%% Listening Host for data sync

class ListeningHost(QtCore.QThread):

    '''
    Listening host takes an SFTP session object and monitors new data
    All it does is changing the newDataAvailable flag and emits a notification

    '''
    emitter = QtCore.pyqtSignal(object)

    def __init__(self, sftp, path_remote, path_local):

        self.newDataAvailable = False
        self.sftp = sftp
        self.pathRemote = path_remote
        self.pathLocal = path_local
        QtCore.QThread.__init__(self)

    def run(self):

        while self.newDataAvailable is False:

#           print('Listening host is actively checking for new data...')
            self.newDataAvailable = sftp_new_data_checker(
                            self.sftp, self.pathRemote, self.pathLocal)
            self.emitter.emit(None)

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

    #%% worker function that gets wrapped

    def reduce_dataset_worker(path_dataset,
                              save_dir_rel = '',
                              clobber = False,
                              dataset_is_flat = False,
                              suppress_return = True):

        # if any master file already exists, and clobber = False, then break
        filenames_combined = '\t'.join(os.listdir(path_dataset))
        if 'master' in filenames_combined.lower() and clobber == False :
            print('Master files already exist. Overwriting not allowed.')
            return

        print('Reducing directory {}...'.format(path_dataset))
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
        print('Read {} bias frames'.format(num_bias))
        print('Read {} dark frames'.format(num_dark))
        print('Read {} object/flat frames'.format(num_object))

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
        print('Master bias saved to {}'.format(
            os.path.join(dir_save, 'master_bias.fits')))
        hdu_m_dark.writeto(
            os.path.join(dir_save, 'master_dark.fits'), clobber = clobber)
        print('Master dark saved to {}'.format(
            os.path.join(dir_save, 'master_dark.fits')))

        # master light which can be object or flat
        if dataset_is_flat:

            paths_flat = paths_object
            num_flat = len(paths_flat)
            flat3D = np.empty([shape[0], shape[1], num_flat], dtype = datatype)

            for i in range(len(paths_flat)):
                # subtract bias and dark
                flat3D[:,:,i] = fits.open(paths_flat[i]) [0].data \
                    - m_bias - m_dark
                [mode, _ ] = scipy.stats.mode(flat3D[:, :, i], axis=None)
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
            print('Master flat saved to {}'.format(
                os.path.join(dir_save, 'master_flat.fits')))

            hdu_m_flat_normalised.writeto(
                os.path.join(dir_save, 'master_flat_normalised.fits'),
                clobber = clobber)
            print('Master flat normalised saved to {}'.format(
                os.path.join(dir_save, 'master_flat_normalised.fits')))

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

            print('Master object saved to {}'.format(
                os.path.join(dir_save, 'master_object.fits')))
            if not suppress_return:
                return m_object

    #%% data reduction actually starts here

    # determine whether there is flat available
    if path_flat is None:

        print('No flat provided, reducing given dataset only...')
        reduce_dataset_worker(path_dataset,
                              save_dir_rel = save_dir_rel,
                              clobber = clobber,
                              dataset_is_flat = False,
                              suppress_return = True)
    else:
        # then we gotta reduce the dataset of flats
        # and then do flat field correction
        print('Flat provided, using flat at {}'.format(path_flat))
        print('Reducing object dataset')
        m_object = reduce_dataset_worker(path_dataset,
                                         save_dir_rel = save_dir_rel,
                                         clobber = clobber,
                                         dataset_is_flat = False,
                                         suppress_return = False)
        print('Reducing flat dataset')
        m_flat_normalised = reduce_dataset_worker(path_flat,
                                                  save_dir_rel = save_dir_rel,
                                                  clobber = clobber,
                                                  dataset_is_flat = True,
                                                  suppress_return = False)

        print('Flat field correction...')
        master = m_object / m_flat_normalised
        hdu_master = fits.PrimaryHDU(master)
        hdu_master.writeto(
            os.path.join(path_dataset, save_dir_rel, 'master.fits'),
            clobber = clobber)
        print('Master image with flat field correction saved to {}'.format(
            os.path.join(path_dataset, save_dir_rel, 'master.fits')))

#%% Qt MainWindow Class

class MainWindow(TemplateBaseClass):

    def __init__(self):

        TemplateBaseClass.__init__(self)
        # super(MainWindow, self).__init__()

        # create the main window
        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        main.ui.path_flat = None
        #self.ui.plotBtn.clicked.connect(self.plot)
        self.show()
        self.setWindowTitle('ProtoDESI Image Viewer and FPC Analysis')

        # establish ssh and sftp, used by listner and sync functions
        self.establish_sftp('all')
        # update file list
        self.update_file_list('all')
        # create and start listeners
        self.create_listener('all')

        # initialise interface components

#        figure_photometry = Figure()
#        self.addmpl(figure_photometry)
#
#        # show first image initially by default
#        if len(self.ui.listRawSTI) > 0:
#            path_file_abs = self.ui.listRaw.item(0).text()
#            self.update_raw(path_file_abs, 'sti')
#
#        if len(self.ui.listRawFVC) > 0:
#            path_file_abs = self.ui.listRaw.item(0).text()
#            self.update_raw(path_file_abs, 'fvc')
#
#        if len(self.ui.listRawFPC) > 0:
#            path_file_abs = self.ui.listRaw.item(0).text()
#            self.update_raw(path_file_abs, 'fpc')

    def msg(self, text):

        self.ui.statusbar.showMessage(text)
        if printEnabled:
            print(text)

    def create_listener(self, device):

        if device == 'all':

            for each_device in DEVICE_LIST:
                self.create_listener(each_device)

        else:

            self.msg('Listening host initialising for {}...'.format(device))

            if device == 'sti':

                self.listeningHostSTI = ListeningHost(
                    self.sftpSTI, DIR_REMOTE_STI, DIR_LOCAL_STI)
                self.listeningHostSTI.emitter.connect(self.on_receiving)
                self.listeningHostSTI.start()

            elif device == 'fvc':

                self.listeningHostFVC = ListeningHost(
                    self.sftpFVC, DIR_REMOTE_FVC, DIR_LOCAL_FVC)
                self.listeningHostFVC.emitter.connect(self.on_receiving)
                self.listeningHostFVC.start()

            elif device == 'fpc':

                self.listeningHostFPC = ListeningHost(
                    self.sftpFPC, DIR_REMOTE_FPC, DIR_LOCAL_FPC)
                self.listeningHostFPC.emitter.connect(self.on_receiving)
                self.listeningHostFPC.start()

        self.msg('Listening host active for {}.'.format(device))

    def restart_listener(self, device):

        if device == 'all':

            for each_device in DEVICE_LIST:
                self.restart_listener(each_device)

        elif device == 'sti':

            self.listeningHostSTI.newDataAvailable = False
            self.listeningHostSTI.start()

        elif device == 'fvc':

            self.listeningHostFVC.newDataAvailable = False
            self.listeningHostFVC.start()

        elif device == 'fpc':

            self.listeningHostFPC.newDataAvailable = False
            self.listeningHostFPC.start()

    def addmpl(self, fig):

        self.ui.canvas = FigureCanvas(fig)
        self.ui.segViewLayout.addWidget(self.ui.canvas)
        self.ui.canvas.draw()
        self.ui.toolbar = NavigationToolbar(self.ui.canvas,
                                        self.ui.segView, coordinates=True)
        self.ui.segViewLayout.addWidget(self.ui.toolbar)

    def rmmpl(self):
        self.ui.segViewLayout.removeWidget(self.ui.canvas)
        self.ui.canvas.close()
        self.ui.segViewLayout.removeWidget(self.ui.toolbar)
        self.ui.toolbar.close()

    def establish_sftp(self, device):

        if device == 'all':

            for each_device in DEVICE_LIST:
                self.establish_sftp(each_device)

        elif device == 'sti':
            self.sftpSTI = open_sftp(DOMAIN_STI, USERNAME, PASSWORD)

        elif device == 'fvc':
            self.sftpFVC = open_sftp(DOMAIN_FVC, USERNAME, PASSWORD)

        elif device == 'fpc':
            self.sftpFPC = open_sftp(DOMAIN_FPC, USERNAME, PASSWORD)

        self.msg('SFTP to {} established.'.format(device))

    def sync(self, device):

        '''
        responds to user click on sync button, and perform actual sync

        '''

        # determine device
        if device == 'all':

            for each_device in DEVICE_LIST:
                self.sync(each_device)

        elif device == 'sti':

            listener = self.listeningHostSTI
            sftp = self.sftpSTI
            domain = DOMAIN_STI
            path_remote = DIR_REMOTE_STI
            path_local = DIR_LOCAL_STI

        elif device == 'fvc':

            listener = self.listeningHostFVC
            sftp = self.sftpFVC
            domain = DOMAIN_FVC
            path_remote = DIR_REMOTE_FVC
            path_local = DIR_LOCAL_FVC

        elif device == 'fpc':

            listener = self.listeningHostFPC
            sftp = self.sftpFPC
            domain = DOMAIN_FPC
            path_remote = DIR_REMOTE_FPC
            path_local = DIR_LOCAL_FPC

        # perform sync
        self.msg('Syncing new data from {}' + '{}'
                    .format(domain, path_remote))

        if listener.newDataAvailable:

            if writeEnabled:
                sftp_get_all(sftp, path_remote, path_local)
                self.msg('Sync complete for {}.'.format(device))
            else:
                self.msg('Write to disk disabled.')
            # update file list
            self.update_file_list(device)
            # restart listener
            self.restart_listener(device)

        else:
            self.msg('No new data available from {}.'.format(device))

    def on_receiving(self):

        # receiver for listener emitter
        if self.listeningHostSTI.newDataAvailable:
            self.msg('New ST-i data available. Sync Now.')
            if autosync:
                self.msg('Autosync is on. Auto-syncing ST-i...')
                self.sync('sti')
        else:
            self.msg('No new ST-i data available.')

        if self.listeningHostFVC.newDataAvailable:
            self.msg('New FVC data available. Sync Now.')
            if autosync:
                self.msg('Autosync is on. Auto-syncing FVC...')
                self.sync('fvc')
        else:
            self.msg('No new FVC data available.')

        if self.listeningHostFPC.newDataAvailable:
            self.msg('New FPC data available. Sync Now.')
            if autosync:
                self.msg('Autosync is on. Auto-syncing FPC...')
                self.sync('fpc')
        else:
            self.msg('No new FPC data available.')

    def update_file_list(self, device):

        name_pattern_raw = '**/*.fit*' # for recursive indexing

        if device == 'all':

            for each_device in DEVICE_LIST:
                self.update_file_list(each_device)

        elif device == 'sti':

            self.msg('Globbing local directory tree {}...'
                                                    .format(DIR_LOCAL_STI))
            # glob using filename patterns
            paths_file_abs = glob.glob(
                os.path.join(DIR_LOCAL_STI, name_pattern_raw), recursive=True)
            # populate list raw
            self.ui.listRawSTI.clear()
            self.ui.listRawSTI.addItems(paths_file_abs)

        elif device == 'fvc':

            self.msg('Globbing local directory tree {}...'
                                                    .format(DIR_LOCAL_FVC))
            # glob using filename patterns
            paths_file_abs = glob.glob(
                os.path.join(DIR_LOCAL_FVC, name_pattern_raw), recursive=True)
            # populate list raw
            self.ui.listRawFVC.clear()
            self.ui.listRawFVC.addItems(paths_file_abs)

        elif device == 'fpc':

            self.msg('Globbing local directory tree {}...'
                                                    .format(DIR_LOCAL_FPC))
            # glob using filename patterns
            paths_file_abs = glob.glob(
                os.path.join(DIR_LOCAL_FPC, name_pattern_raw), recursive=True)
            # populate list raw
            self.ui.listRawFPC.clear()
            self.ui.listRawFPC.addItems(paths_file_abs)

            # Photometry tab
            # glob using filename patterns
            name_pattern_raw = '**/*master*.fit*'
            paths_file_abs = glob.glob(
                os.path.join(DIR_LOCAL_FPC, name_pattern_raw), recursive=True)
            # populate list master
            self.ui.listMaster.clear()
            self.ui.listMaster.addItems(paths_file_abs)

        # list updated for chosen device
        self.msg('File list updated for {}'.format(device))


    def update_raw(self, path_file_abs, device):

        self.msg(repr('Showing image:'+ path_file_abs))

        if device == 'sti':

            hdu = fits.open(path_file_abs)[0]

            # display header
            self.ui.headerViewSTI.setCurrentFont(QtGui.QFont(TEXT_BROWSER_FONT))
            self.ui.headerViewSTI.setFontPointSize(TEXT_BROWSER_FONT_SIZE)
            self.ui.headerViewSTI.setPlainText(repr(hdu.header))

            # display image selected
            minlevel = np.amin(hdu.data)
            maxlevel = 2000
            self.ui.rawViewSTI.show()
            self.ui.rawViewSTI.setImage(np.rot90(hdu.data, -1))
            self.ui.rawViewSTI.setLevels(minlevel, maxlevel)
            #self.ui.rawView.autoLevels()

        elif device == 'fvc':

            hdu = fits.open(path_file_abs)[0]

            # display header
            self.ui.headerViewFVC.setCurrentFont(QtGui.QFont(TEXT_BROWSER_FONT))
            self.ui.headerViewFVC.setFontPointSize(TEXT_BROWSER_FONT_SIZE)
            self.ui.headerViewFVC.setPlainText(repr(hdu.header))

            # display image selected
            minlevel = np.amin(hdu.data)
            maxlevel = 2000
            self.ui.rawViewFVC.show()
            self.ui.rawViewFVC.setImage(np.rot90(hdu.data, -1))
            self.ui.rawViewFVC.setLevels(minlevel, maxlevel)
            #self.ui.rawView.autoLevels()

        elif device == 'fpc':

            hdu = fits.open(path_file_abs)[0]

            # display header
            self.ui.headerViewFPC.setCurrentFont(QtGui.QFont(TEXT_BROWSER_FONT))
            self.ui.headerViewFPC.setFontPointSize(TEXT_BROWSER_FONT_SIZE)
            self.ui.headerViewFPC.setPlainText(repr(hdu.header))

            # display image selected
            minlevel = np.amin(hdu.data)
            maxlevel = 2000
            self.ui.rawViewFPC.show()
            self.ui.rawViewFPC.setImage(np.rot90(hdu.data, -1))
            self.ui.rawViewFPC.setLevels(minlevel, maxlevel)
            #self.ui.rawView.autoLevels()

    def update_photometry(self, path_file_abs):

        hdu = fits.open(path_file_abs)[0]
        data = hdu.data

        # display master preview
        self.msg('Showing master image:' + repr(path_file_abs))
        self.rmmpl()
        figure_photometry = Figure()
#       cmapRand = random_cmap(segm.max+1, random_state=12345)
        axes = figure_photometry.add_subplot(111)
        axes.imshow(data, origin='lower', cmap=plt.cm.gray)
        self.addmpl(figure_photometry)

        # perform photometry
        self.msg('Performing aperture photometry...')
        self.aperture_photometry(path_file_abs)

    def aperture_photometry(self, path_file_abs):

        """
        aperture photometry from source segmentation
        make_source_mask not yet available in photutils v0.2.1
        wait for v0.3 release

        aperture_photometry() assumes that the data have been
        background-subtracted.

        """

        # create preliminary mask
        #from photutils import make_source_mask
        #masterMask = make_source_mask(master, snr=2, npixels=5, dilate_size=11)

#        if LEDoff was used, get threshold from LEDoff/background
#        path_dataset = os.path.dirname(path_file_abs) + os.path.sep
#        filenameCombined = '\t'.join(
#            os.listdir(os.path.join(datasetDirLocal, 'master')))
#        if 'master_ledoff_subtracted' in filename:
#            self.msg('Using master_ledoff')
#            # path_file_abs = os.path.join(datasetDir, 'master', filename)
#            hdu = fits.open(path_file_abs)[0]
#            data_subtracted = hdu.data
#            # calculate threadhold
#            ledoff_pred = np.mean(data_subtracted) * \
#                np.ones(data_subtracted.shape)
#            mse = mean_squared_error(data_subtracted, ledoff_pred)
#            rmse = np.sqrt(mse)
#            threshold = 7.0 * rmse
#            threshold_value = threshold
#         if no LEDoff was used, background subtraction is needed
#         there should exist no file named "subtracted"
#         if 'master.fit' in filenameCombined \
#             or 'master_normalised.fit' in filenameCombined:

        if 'master.fit' in path_file_abs:
            self.msg('Photometry using un-normalised master image')
        elif 'master_normalised.fit' in path_file_abs:
            self.msg('Photometry using normalised master image')

        hdu = fits.open(path_file_abs)[0]
        data = hdu.data

        if 'EXPTIME' in hdu.header:
            exptime = hdu.header['EXPTIME']
        else:
            exptime = hdu.header['EXPREQ']

        # === background subtraction ===
        """
        if no LEDoff was used, background subtraction is needed.
        there should exist no file named "subtracted".

        create 2D image of background and background rms and
        apply sigma-clipping to each region in the low-res
        background map to get mean, median, and std/rms.
        sigma-clipping is the most widely used method though not as
        good as using mask; still superior to robust standard
        deviation using median absolute deviation (MAD-STD).

        """

        # create background
        # [mean, median, std] = sigma_clipped_stats(master, sigma=3.0, iters=5)
        # bkg = Background(master, (50, 50), filter_size=(3, 3), method='median')
        bkg = Background(data, (100, 100), filter_shape=(3, 3), method='median')

        # plot background image
        # plt.imshow(bkg.background, norm=normalisation, origin='lower', cmap=plt.cm.gray)
        plt.imshow(bkg.background, origin='lower', cmap=plt.cm.gray)
        [fig, ax] = plt.subplots(figsize=(8, 8))

        # make background-substracted image
        data_subtracted = data - bkg.background

        # plot subtracted image
        plt.imshow(data_subtracted, origin='lower', cmap=plt.cm.gray)

        # save background subtracted image
        path_dataset = os.path.dirname(path_file_abs)
        if 'master.fit' in path_file_abs:
            savepath = os.path.join(
                path_dataset, 'master_subtracted.fits')
        elif 'master_normalised.fit' in path_file_abs:
            savepath = os.path.join(
                path_dataset, 'master_normalised_subtracted.fits')
        hdu_subtracted = fits.PrimaryHDU(data_subtracted)
        hdu_subtracted.writeto(savepath, clobber = True)

        # === segmentation at a given sigma level ===
        # perform segmentation whether flat is available or not

        self.msg('Determining threshold for target detection...')
        # because data is background-subtracted
        threshold_array = 5.0 * bkg.background_rms
        # print out threshold value
        threshold_value = threshold_array.flat[0]
        self.msg('Threshold for target detection is: ' + repr(threshold_value))

        self.msg('Detecting sources and performing segmentation...')
        segm = detect_sources(data_subtracted, threshold_array, npixels=5)

        self.msg('Segmentation labels are:')
        self.msg((repr(segm.labels)))

        # === regional properties ===
        # measure regional source properties from segmentation
        # the centroid is from image moments, already intensity-weighted

        self.msg('Measuring source properties...')
        if 'bkg' in locals():
            # use the background determined from master_subtracted
            props = source_properties(data_subtracted, segm,
                error = bkg.background_rms, background = bkg.background)

#        elif 'master_ledoff_subtracted' in filenameCombined:
#            path_file_abs = os.path.join(
#                datasetDirLocal, 'master', 'master_ledoff_subtracted.fits')
#            hdu = fits.open(path_file_abs)[0]
#            master_ledoff_subtracted = hdu.data
#            props = source_properties(data_subtracted, segm,
#                error = master_ledoff_subtracted \
#                        - np.mean(master_ledoff_subtracted),
#                background = master_ledoff_subtracted)

        # add instrumental magnitude to properties
        # instrumental magnitude = -2.5 * log10(flux)
        for i in range(len(props)):
            # source_sum is by definition background-subtracted already
            props[i].mag_instr = -2.5 * np.log10(props[i].source_sum/exptime)

        # create table from props object
        # there are other properties available, see list of SourceProperties:
        # http://goo.gl/rkfQ9V

        propsTableColumns = ['id', 'xcentroid', 'ycentroid', 'area', 'max_value',
            'source_sum', 'mag_instr']
        propsTable = properties_table(props, columns = propsTableColumns)
        self.msg(repr(propsTable))

        # === update UI ===
        # plot segmentated image
        self.rmmpl()
        figure_photometry = Figure()
        cmap_rand = random_cmap(segm.max+1, random_state=12345)
        axes = figure_photometry.add_subplot(111)
        axes.imshow(segm, origin='lower', cmap=cmap_rand)
        axes.plot(
            propsTable['xcentroid'], propsTable['ycentroid'],
            ls='none', color='red', marker='+', ms=10, lw=1.5)
        self.addmpl(figure_photometry)

        # set properties table font and font size
        self.ui.tablePhot.setCurrentFont(QtGui.QFont(TEXT_BROWSER_FONT))
        self.ui.tablePhot.setFontPointSize(TEXT_BROWSER_FONT_SIZE)
        self.ui.tablePhot.setPlainText(repr(propsTable))

        self.msg('Photometry completed for {}.'.format(path_file_abs))

#        # old plots for visualisation
#
#        apertures = []
#        for prop in props:
#            position = (prop.xcentroid.value, prop.ycentroid.value)
#            a = prop.semimajor_axis_sigma.value * 3.0
#            b = prop.semiminor_axis_sigma.value * 3.0
#            theta = prop.orientation.value
#            apertures.append(EllipticalAperture(position, a, b, theta=theta))
#        norm = ImageNormalize(stretch=SqrtStretch())
#        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 18))
#
#        if 'bkg' in locals():
#            ax1.imshow(master_subtracted, origin='lower', cmap='Greys_r', norm=norm)
#        else:
#            ax1.imshow(master_subtracted_normalised, origin='lower', cmap='Greys_r', norm=norm)
#        ax2.imshow(segm, origin='lower', cmap=cmapRand)
#        for aperture in apertures:
#            aperture.plot(color='blue', lw=1.5, alpha=0.5, ax=ax1)
#            aperture.plot(color='white', lw=1.5, alpha=1.0, ax=ax2)

#%% Start Qt event loop unless running in interactive mode or using pyside

main = MainWindow()

if __name__ == '__main__':

    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):

        def list_item_changed_sti(curr, prev):
            path_file_abs = curr.text()
            main.update_raw(path_file_abs, 'sti')

        def list_item_changed_fvc(curr, prev):
            path_file_abs = curr.text()
            main.update_raw(path_file_abs, 'fvc')

        def list_item_changed_fpc(curr, prev):
            path_file_abs = curr.text()
            main.update_raw(path_file_abs, 'fpc')

        def master_item_changed(curr, prev):
            path_file_abs = curr.text()
            main.update_photometry(path_file_abs)

        # connect button trigger callbacks
        main.ui.actionSyncAll.triggered.connect(
                                        functools.partial(main.sync, 'all'))
        main.ui.actionSyncSTI.triggered.connect(
                                        functools.partial(main.sync, 'sti'))
        main.ui.actionSyncFVC.triggered.connect(
                                        functools.partial(main.sync, 'fvc'))
        main.ui.actionSyncFPC.triggered.connect(
                                        functools.partial(main.sync, 'fpc'))

        def set_path_to_flat():

            # prompt input parent directory
            main.ui.path_flat = input('Enter absolute path of directory \
                containing flat frames, no quotes: \n')

        def reduce_dataset_dir():

            main.ui.path_dataset = input('Enter absolute path of directory \
                containing object frames, no quotes: \n')
            reduce_dataset(main.ui.path_dataset,
                           path_flat = main.ui.path_flat)

        def reduce_current_dataset_dir():

            main.msg('Performing data reduction for currently selected \
                dataset...')
            # when a sync action is triggered, get active tab
            # and parse current selection
            active_tab_index = main.ui.tabWidget.currentIndex()
            if active_tab_index == 0: # ST-i tab
                path_file = main.ui.listRawSTI.currentItem.text()
            elif active_tab_index == 1: # FVC tab
                path_file = main.ui.listRawFVC.currentItem.text()
            elif active_tab_index == 2: # FPC tab
                path_file = main.ui.listRawFPC.currentItem.text()
            # get parent folder of file selected
            main.ui.path_dataset = os.path.dirname(path_file)
            reduce_dataset(main.ui.path_dataset, path_flat = main.ui.path_flat)

        def reduce_parent_dir():

            main.ui.path_parent = input('Enter absolute path of directory \
                containing folders, each of which is a dataset, no quotes: \n')
            main.msg('Performing data reduction for parent directory...')

            paths_dataset = glob.glob(
                os.path.abspath(
                    main.ui.path_parent + os.path.sep + '*')
                + os.path.sep)

            for path_dataset in paths_dataset:
                reduce_dataset(paths_dataset, path_flat = main.ui.path_flat)

        # -------

        main.ui.actionSetPathToFlat.triggered.connect(
            set_path_to_flat)
        main.ui.reduce_parent_dataset.triggered.connect(
            reduce_dataset_dir)
        main.ui.actionReduceCurrentDatasetDir.triggered.connect(
            reduce_current_dataset_dir)
        main.ui.actionReduceParentDir.triggered.connect(
            reduce_parent_dir)

        # connect list item change trigger callbacks
        main.ui.listRawSTI.currentItemChanged.connect(list_item_changed_sti)
        main.ui.listRawFVC.currentItemChanged.connect(list_item_changed_fvc)
        main.ui.listRawFPC.currentItemChanged.connect(list_item_changed_fpc)
        main.ui.listMaster.currentItemChanged.connect(master_item_changed)

        QtGui.QApplication.instance().exec_()
