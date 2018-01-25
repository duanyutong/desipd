# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 23:53:06 2016

@author: Duan Yutong
"""

import os
# import sys
# import functools
# import glob
import stat
# import numpy as np
# import scipy
# from astropy.io import fits
import paramiko
# import pyqtgraph as pg
from PyQt4 import QtCore #, QtGui, uic
# from sklearn.metrics import mean_squared_error

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
            if isdir(sftp, path):
                print('Checking subdirectory {}'.format(path))
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
