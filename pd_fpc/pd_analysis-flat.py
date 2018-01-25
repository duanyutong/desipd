# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 02:02:35 2016

@author: givoltage
"""

# import os
# os.chdir(r'K:\Google Drive\DESI\repository\repository_kpno\pd_fpc')
from pd_fpc_analysis import FPCAnalysis

analysis = FPCAnalysis('dome_flat')
analysis.CLOBBER = False

# specify dataset exptypes and expids
analysis.specify_dataset('bias stack', 5310, 5320,
                         'bias')
analysis.specify_dataset('dark stack', 6172, 6182,
                         'dark',
                         bias_name = 'bias')
analysis.specify_dataset('image stack', 1000, 1001,
                         '1s',
                         bias_name = 'bias',
                         dark_name = 'dark')
analysis.specify_dataset('image stack', 1002, 1003,
                         '2s',
                         bias_name = 'bias',
                         dark_name = 'dark')

analysis.read_datasets_header()
# reduce datasets, order is important
analysis.reduce_stack('bias')
analysis.reduce_stack('dark')
analysis.reduce_image_stack('1s')
analysis.reduce_image_stack('2s')
analysis.photometry('1s')
analysis.photometry('2s')
analysis.dump_dict()


analysis.load_dict()
analysis.plot_flats(['1s', '2s'])
