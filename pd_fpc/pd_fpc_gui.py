#%% -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 22:49:36 2016

@author: Duan Yutong (dyt@physics.bu.edu, dyt@lbl.gov)

"""

import os
import sys
import functools
import glob
# import stat
import numpy as np
# import scipy
from astropy.io import fits
# import paramiko
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

os.chdir(os.path.dirname(os.path.abspath(__file__)))
from pd_fpc_gui_utils import open_sftp, sftp_get_all, ListeningHost
from data_reduction import reduce_dataset

#%% Global Settings

autosync = False
printEnabled = True
syncWriteEnabled = False

#%% Parameters

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

# display
TEXT_BROWSER_FONT = 'Courier'
TEXT_BROWSER_FONT_SIZE = 13

#%% Qt setup

pg.mkQApp()
uiPath = os.path.dirname(os.path.abspath(__file__))
uiFile = os.path.join(uiPath, 'pd_fpc_gui.ui')
[WindowTemplate, TemplateBaseClass] = uic.loadUiType(uiFile)

#%% Qt MainWindow Class

class MainWindow(TemplateBaseClass):

    def __init__(self):

        TemplateBaseClass.__init__(self)
        # super(MainWindow, self).__init__()

        # create the main window
        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        self.ui.path_flat = None
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

            if syncWriteEnabled:
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
            self.ui.listRawSTI.addItems(sorted(paths_file_abs))

        elif device == 'fvc':

            self.msg('Globbing local directory tree {}...'
                                                    .format(DIR_LOCAL_FVC))
            # glob using filename patterns
            paths_file_abs = glob.glob(
                os.path.join(DIR_LOCAL_FVC, name_pattern_raw), recursive=True)
            # populate list raw
            self.ui.listRawFVC.clear()
            self.ui.listRawFVC.addItems(sorted(paths_file_abs))

        elif device == 'fpc':

            self.msg('Globbing local directory tree {}...'
                                                    .format(DIR_LOCAL_FPC))
            # glob using filename patterns
            paths_file_abs = glob.glob(
                os.path.join(DIR_LOCAL_FPC, name_pattern_raw), recursive=True)
            # populate list raw
            self.ui.listRawFPC.clear()
            self.ui.listRawFPC.addItems(sorted(paths_file_abs))

            # Photometry tab
            # glob using filename patterns
            name_pattern_raw = '**/*master*.fit*'
            paths_file_abs = glob.glob(
                os.path.join(DIR_LOCAL_FPC, name_pattern_raw), recursive=True)
            print(paths_file_abs)
            # populate list master
            self.ui.listMaster.clear()
            self.ui.listMaster.addItems(sorted(paths_file_abs))

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
        if hasattr(self.ui, 'canvas'):
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

        Note that photutils.source_properties() assumes that the data have been
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

        props_table_columns = [
            'id', 'xcentroid', 'ycentroid',
            'area', 'max_value', 'source_sum', 'mag_instr']
        props_table = properties_table(
                                props, columns = props_table_columns)
        props_table_save = properties_table(props)
#        self.msg(repr(props_table_display))
        print(repr(props_table))

        # check and create analysis folder if it doesn't exist
        path_dataset = os.path.dirname(path_file_abs)
        path_analysis = path_dataset.replace(
            '/data/images/fpc/',
            '/data/images/fpc_analysis/')
        if not os.path.exists(path_analysis):
            os.makedirs(path_analysis)

        # save background subtracted image
        if 'master.fit' in path_file_abs:
            path_save = os.path.join(
                path_dataset, 'master_subtracted.fits')
        elif 'master_object.fit' in path_file_abs:
            path_save = os.path.join(
                path_dataset, 'master_object_subtracted.fits')
        elif 'master_normalised.fit' in path_file_abs:
            path_save = os.path.join(
                path_dataset, 'master_normalised_subtracted.fits')
        hdu_subtracted = fits.PrimaryHDU(data_subtracted)
        hdu_subtracted.writeto(path_save, clobber = True)

        # save properties to table file
        path_save = os.path.join(path_dataset, 'props_table.csv')
        ascii.write(props_table_save, path_save, format = 'csv')

        # === update UI ===
        # plot segmentated image
        self.rmmpl()
        figure_photometry = Figure()
        cmap_rand = random_cmap(segm.max+1, random_state=12345)
        axes = figure_photometry.add_subplot(111)
        axes.imshow(segm, origin='lower', cmap=cmap_rand)
        axes.plot(
            props_table['xcentroid'], props_table['ycentroid'],
            ls='none', color='red', marker='+', ms=10, lw=1.5)
        self.addmpl(figure_photometry)

        # set properties table font and font size
        self.ui.tablePhot.setCurrentFont(QtGui.QFont(TEXT_BROWSER_FONT))
        self.ui.tablePhot.setFontPointSize(TEXT_BROWSER_FONT_SIZE)
        self.ui.tablePhot.setPlainText(repr(props_table))

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
                reduce_dataset(paths_dataset,
                               path_flat = main.ui.path_flat,
                               save_dir_rel = '',
                               clobber = False)

        #%% triggers

        main.ui.actionSetPathToFlat.triggered.connect(
            set_path_to_flat)
        main.ui.actionReduceDatasetDir.triggered.connect(
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
