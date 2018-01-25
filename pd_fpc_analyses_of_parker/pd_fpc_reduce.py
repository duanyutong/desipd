# -*- coding: utf-8 -*-
"""
Created on Mon May 16 14:41:22 2016

@author: Duan Yutong (dyt@lbl.gov)
"""

#%% import modules


#import sys

#import matplotlib.pyplot as plt
# from skimage import measure
# import cv2

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
# from photutils.morphology import centroid_com, centroid_1dg, centroid_2dg

#%% function test run


import os
# import glob
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from fpc_functions import reduce_dataset

path_dataset = '/data/images/fpc/linear/'
path_flat = '/data/images/fpc/flats/'
save_dir_rel = ''
clobber = False

reduce_dataset(path_dataset,
               path_flat = path_flat,
               save_dir_rel = '',
               clobber = clobber)

#path_parent = '/data/images/fpc/linear/'

#paths_dataset = glob.glob(
#        os.path.abspath(
#            path_parent + os.path.sep + '*')
#        + os.path.sep)
#
#for path_dataset in paths_dataset:
#    print('Path to child dataset is {}'.format(path_dataset))
#    reduce_dataset(path_dataset,
#                   path_flat = path_flat,
#                   save_dir_rel = '',
#                   clobber = clobber)

