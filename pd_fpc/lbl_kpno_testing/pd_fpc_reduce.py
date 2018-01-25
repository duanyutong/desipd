# -*- coding: utf-8 -*-
"""
Created on Mon May 16 14:41:22 2016

@author: Duan Yutong (dyt@physics.bu.edu, dyt@lbl.gov)

Reduction mode: 1 for single dataset, 2 for a parent directory

This script works for both object and fake flat (white board calibration) datasets.

"""

#%% import modules

import os
import glob

os.chdir(os.path.dirname(os.path.abspath(__file__)))
from data_reduction_std import reduce_dataset

mode = 1    # set reduction mode here
dir_parent = '/data/images/fpc/20160827/linear'
dir_dataset = r'C:\Users\givoltage\Downloads\flats'
clobber = True
white_board_calibration = True

# reduce a single dataset: same exptime, same object repeatedly exposed
if mode == 1:

    dir_save = dir_dataset.replace('/data/images/fpc/',
                                   '/data/images/fpc_analysis/')
    reduce_dataset(dir_dataset,
                   dir_flat = None,
                   dir_save = dir_save,
                   clobber = clobber)

    if white_board_calibration:
        #rename output 'master_object.fits' to fpc flat
        try:
            os.rename(os.path.join(dir_save, 'master_object.fits'),
                      os.path.join(dir_save, 'master_object-fpc_flat.fits'))
        except:
            print('Cannot rename master_object.fits.')

# reduce a parent directory containing static datasets
elif mode == 2:

    dirs_dataset = glob.glob(
                        os.path.abspath(
                            dir_parent + os.path.sep + '*')
                        + os.path.sep)

    for dir_dataset in dirs_dataset:
        print('Path to child dataset is {}'.format(dir_dataset))
        dir_save = dir_dataset.replace('/data/images/fpc/',
                                       '/data/images/fpc_analysis/')
        reduce_dataset(dir_dataset,
                       dir_flat = None,
                       dir_save = dir_save,
                       clobber = clobber)

        if white_board_calibration:
        #rename output 'master_object.fits' to fpc flat
            try:
                os.rename(os.path.join(dir_save, 'master_object.fits'),
                          os.path.join(dir_save, 'master_object-fpc_flat.fits'))
            except:
                print('Cannot rename master_object.fits.')
