import pyqtgraph as pg
from PyQt4 import QtCore, QtGui, uic
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from astropy.io import fits
from sklearn.metrics import mean_squared_error
from photutils.background import Background
from photutils import detect_sources
from photutils.utils import random_cmap
from photutils import source_properties, properties_table
import paramiko

# set directory
#parentDir = r'K:\Google Drive\DESI\protoDESI\images\fpc_data'
dirRemote = r'/home/msdos/SBIG/'
dirLocal = os.path.join(os.path.expanduser("~"), 'Downloads', 'SBIG')
datasetID = ''
datasetDirRemote = os.path.join(dirRemote, datasetID)
datasetDirLocal = os.path.join(dirLocal, datasetID)
domain = 'desisti.kpno.noao.edu'
username = 'msdos'
password = 'MS-d0s'
autosync = True
printEnabled = True
pg.mkQApp()

#%% Define main window class from UI file

uiPath = os.path.dirname(os.path.abspath(__file__))
uiFile = os.path.join(uiPath, 'sti_gui.ui')
[WindowTemplate, TemplateBaseClass] = uic.loadUiType(uiFile)

def compare_sets(sftp):
    
    setRemote = set(sftp.listdir(datasetDirRemote))
    listRemote = list(setRemote)
    listRemoteConverted = list(listRemote)
    for i in range(len(listRemote)):
        listRemoteConverted[i] = listRemoteConverted[i].replace(':','_')
    setRemoteConverted = set(listRemoteConverted)
    setLocal = set(os.listdir(datasetDirLocal))
    setTransferConverted = setRemoteConverted - (setRemoteConverted & setLocal)
    setTransfer = set(setTransferConverted)
    for filename in setTransferConverted:
        i = listRemoteConverted.index(filename)
        setTransfer.remove(listRemoteConverted[i])
        setTransfer.add(listRemote[i])
    return setTransfer

class ListeningHost(QtCore.QThread):
    
    emitter = QtCore.pyqtSignal(object)
    
    def __init__(self, sftp):
        
        self.newDataAvailable = False
        self.sftp = sftp
        QtCore.QThread.__init__(self)
        
    def run(self):
        
        while self.newDataAvailable is False:
            
#            print('Listening host is actively checking for new data...')
            setTransfer = compare_sets(self.sftp)

            if len(setTransfer) is 0:
                self.newDataAvailable = False
            else:
                self.newDataAvailable = True
            self.emitter.emit(None)
            
#        while self.newDataAvailable is True:
#            self.emitter.emit(None)

class MainWindow(TemplateBaseClass):
    
    def __init__(self):
        
        TemplateBaseClass.__init__(self)
        # super(MainWindow, self).__init__()
    
        # create the main window
        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        #self.ui.plotBtn.clicked.connect(self.plot)
        self.show()
        self.setWindowTitle('ProtoDESI ST-i Focusing Viewer')
        self.ui.actionSync.triggered.connect(self.sync)
        
        # establish ssh
        self.open_ssh(domain, username, password)        

        self.updateFileList()
        segFig = Figure()
        self.addmpl(segFig)
        # show first image initially
        if len(self.ui.listRaw) > 0:
            filename = self.ui.listRaw.item(0).text()
            self.updateRaw(filename)
        
        # configure listener
        self.listeningHost = ListeningHost(self.sftp)
        self.listeningHost.emitter.connect(self.receiver)
        self.listeningHost.start()
        self.msg('Listening host started.')
        
    def msg(self, text):
        
        self.ui.statusbar.showMessage(text)
        if printEnabled:
            print(text)
            
    def paramiko_sftp_progress(self, transferred, total):
        
        percentage = transferred/total*100
        total_mb = total/(1024**2)
        self.msg('{0:.2f}%, Total {1:.2f} MB \r'.format(percentage, total_mb))
#        self.ui.statusbar.showMessage('{0:.2f}%, Total {1:.2f} MB \r'.format(percentage, total_mb))

        
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

    def open_ssh(self, domain, user, pw):
        
        self.msg('Establishing SSH...')
        self.ssh = paramiko.SSHClient()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        self.ssh.connect(domain, username=user, password=pw)
        
        self.msg('Establishing SFTP...')
        self.sftp = self.ssh.open_sftp()
        
    def sync(self):
        
        # respond to user click on sync button, and perform actual sync
        
        self.msg('Syncing new data to local drive...')
        
        if self.listeningHost.newDataAvailable:
            
            setTransfer = compare_sets(self.sftp)
            self.msg('Files to be transferred: ' + repr(setTransfer))
            
            for filename in setTransfer:
                self.msg('Transferring {} ...'.format(filename))
                filenameLocal = filename.replace(':', '_')
                try:
                    self.sftp.get(
                        os.path.join(datasetDirRemote, filename), 
                        os.path.join(datasetDirLocal, filenameLocal),
                        callback = self.paramiko_sftp_progress
                        )
                    self.msg('Transfer is complete: {}.'.format(filename))
                except:
                    self.msg('Skipping file: {}'.format(filename))
                    pass
            
            self.listeningHost.newDataAvailable = False
            self.listeningHost.start()
            self.updateFileList()
        
        else:
            self.msg('No new data available for sync.')
            
    def receiver(self):
        
        # receiver for listener emitter
        if self.listeningHost.newDataAvailable:
            self.msg('New data available. Sync Now.')
            if autosync:
                self.msg('Autosync is on. Auto-syncing...')
                self.sync()
        else:
            self.msg('No new data available.')
    
    def updateFileList(self):
        
        self.msg('Updating file lists...')
        
        # filename patterns
        os.chdir(dirLocal)
        namePatternRaw = '**/*.fit*'
        imgEntries = glob.glob(namePatternRaw, recursive=True)

        # populate list raw
        self.ui.listRaw.clear()
        self.ui.listRaw.addItems(imgEntries)
        
        # populate list master
        self.ui.listMaster.clear()
        self.ui.listMaster.addItems(imgEntries)
            
    def updateRaw(self, filename):
        
        self.msg(repr('Showing image:'+ filename))
        filepath = os.path.join(dirLocal, filename)
        hdu = fits.open(filepath)[0]
        
        # display header
        self.ui.headerView.setCurrentFont(QtGui.QFont('Courier'))
        self.ui.headerView.setFontPointSize(9)
        self.ui.headerView.setPlainText(repr(hdu.header))
        
        # display image selected
        minlevel = np.amin(hdu.data)
        maxlevel = 2000
        self.ui.rawView.show()
        self.ui.rawView.setImage(np.rot90(hdu.data, -1))
        self.ui.rawView.setLevels(minlevel, maxlevel)
        #self.ui.rawView.autoLevels()

    def updatePhot(self, filename):
        
        filepath = os.path.join(datasetDirLocal, 'master', filename)
        hdu = fits.open(filepath)[0]
        data = hdu.data
        
        # display master preview
        self.msg('Showing master image:' + repr(filename))
        self.rmmpl()
        segFig = Figure()
#       cmapRand = random_cmap(segm.max+1, random_state=12345)
        axes = segFig.add_subplot(111)
        axes.imshow(data, origin='lower', cmap=plt.cm.gray)
        self.addmpl(segFig)

        # perform photometry
        self.msg('Performing aperture photometry...')
        self.aperture_photometry(filename)
        
    def aperture_photometry(self, filename):

        # aperture photometry from source segmentation
        
        # determine threshold for background detection
        # if LEDoff was used, get threshold from LEDoff/background
        filepath = os.path.join(datasetDirLocal, 'master', filename)
        filenameCombined = '\t'.join(os.listdir(os.path.join(datasetDirLocal, 'master')))
        if 'master_ledoff_subtracted' in filename:
            self.msg('Using master_ledoff')
            # filepath = os.path.join(datasetDir, 'master', filename)
            hdu = fits.open(filepath)[0]
            data_subtracted = hdu.data
            # calculate threadhold
            ledoff_pred = np.mean(data_subtracted) * np.ones(data_subtracted.shape)
            mse = mean_squared_error(data_subtracted, ledoff_pred)  
            rmse = np.sqrt(mse)
            threshold = 7.0 * rmse
            threshold_value = threshold
            
        # if no LEDoff was used, background subtraction is needed
        # there should exist no file named "subtracted"
        elif 'master.fit' in filenameCombined \
             or 'master_normalised.fit' in filenameCombined:
            self.ui.statusbar.showMessage('Using master or master_normalised')
            
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
    
    
            hdu = fits.open(filepath)[0]
            data = hdu.data
            if 'EXPTIME' in hdu.header:
                exptime = hdu.header['EXPTIME']
            else:
                exptime = hdu.header['EXPREQ']
            
            self.msg('Determining threshold for target detection...')
            # calculate threashold
            # [mean, median, std] = sigma_clipped_stats(master, sigma=3.0, iters=5)
            bkg = Background(data, (100, 100), filter_shape=(3, 3), method='median')
            # bkg = Background(master, (50, 50), filter_size=(3, 3), method='median')
            # plt.imshow(bkg.background, norm=normalisation, origin='lower', cmap=plt.cm.gray)
            plt.imshow(bkg.background, origin='lower', cmap=plt.cm.gray)
            [fig, ax] = plt.subplots(figsize=(8, 8))
            # make background-substracted image
            data_subtracted = data - bkg.background
            # plot
            plt.imshow(data_subtracted, origin='lower', cmap=plt.cm.gray)
            
            # save background subtracted image
            if 'master.fit' in filename:
                hdu_subtracted = fits.PrimaryHDU(data_subtracted)
                hdu_subtracted.writeto('master_subtracted.fits', clobber = True)
            elif 'master_normalised.fit' in filename:
                hdu_normalised_subtracted = fits.PrimaryHDU(data_subtracted)
                hdu_normalised_subtracted.writeto('master_normalised_subtracted.fits', clobber = True)
    
            # segmentation at a given sigma level, for regional properties
            threshold = 5.0 * bkg.background_rms # since data is background-subtracted
            threshold_value = threshold.flat[0]
            
        self.msg('Threshold for target detection is: ' + repr(threshold_value))
        # perform segmentation whether flat was available or not
        self.msg('Performing segmentation...')
        segm = detect_sources(data_subtracted, threshold, npixels=5)
        
        self.msg('Segmentation labels are:')
        self.msg((str(segm.labels)))
        # measure regional source properties from segmentation
        # the centroid is from image moments, already intensity-weighted
        self.msg('Measuring source properties')
        if 'bkg' in locals():
            props = source_properties(data_subtracted, segm,
                error = bkg.background_rms, background = bkg.background)
        elif 'master_ledoff_subtracted' in filenameCombined:
            filepath = os.path.join(datasetDirLocal, 'master', 'master_ledoff_subtracted.fits')
            hdu = fits.open(filepath)[0]
            master_ledoff_subtracted = hdu.data
            props = source_properties(data_subtracted, segm,
                error = master_ledoff_subtracted - np.mean(master_ledoff_subtracted),
                background = master_ledoff_subtracted)
                
        # instrumental magnitude = -2.5 * log10(flux)
        for i in range(len(props)):
            props[i].mag_instr = -2.5 * np.log10(props[i].source_sum/exptime)
            # source_sum are by definition background-subtracted already
        propsTableColumns = ['id', 'xcentroid', 'ycentroid', 'area', 'max_value',
            'source_sum', 'mag_instr']
        # there are other properties available, see list of SourceProperties
        # http://photutils.readthedocs.io/en/latest/api/photutils.segmentation.SourceProperties.html#photutils.segmentation.SourceProperties
            
        propsTable = properties_table(props, columns = propsTableColumns)
        self.ui.statusbar.showMessage(repr(propsTable))
        
        # plot segmentated image
        self.rmmpl()
        segFig = Figure()
        cmapRand = random_cmap(segm.max+1, random_state=12345)
        axes = segFig.add_subplot(111)
        axes.imshow(segm, origin='lower', cmap=cmapRand)
        axes.plot(propsTable['xcentroid'], propsTable['ycentroid'], ls='none', color='red',
                 marker='+', ms=10, lw=1.5)
        self.addmpl(segFig)
        
        # set properties table font and font size
        self.ui.tablePhot.setCurrentFont(QtGui.QFont('Courier'))
        self.ui.tablePhot.setFontPointSize(9)
        self.ui.tablePhot.setPlainText(repr(propsTable))
        
        self.msg('Photometry completed')
#        # plots for visualisation
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
 
        def raw_item_changed(curr, prev):
            filename = curr.text()
            main.updateRaw(filename)
        
        def master_item_changed(curr, prev):
            filename = curr.text()
            main.updatePhot(filename)

        main.ui.listRaw.currentItemChanged.connect(raw_item_changed)
        main.ui.listMaster.currentItemChanged.connect(master_item_changed)
        
        QtGui.QApplication.instance().exec_()