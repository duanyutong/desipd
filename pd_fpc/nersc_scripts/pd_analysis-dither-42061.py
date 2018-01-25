# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 23:52:15 2016

@author: givoltage
"""

from pd_fpc_analysis import FPCAnalysis

analysis = FPCAnalysis('dither-42061')

#analysis.CLOBBER = False
#analysis.DUMP_PICKLE = True
#
## specify dataset exptypes and expids
#analysis.specify_dataset('bias stack', 5517, 5527,
#                         'bias')
#analysis.specify_dataset('dark stack', 5852, 5860,
#                         'dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('image sequence', 3592, 3616,
#                         '42061',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '42061',
#                         focus = -10800,
#                         dither_mode = 'positioner',
#                         dither_grid = (5, 5),
#                         dither_step = 0.5)
#
#analysis.read_datasets_header()
#
## reduce datasets, order is important
#analysis.reduce_stack('bias')
#analysis.reduce_stack('dark')
#analysis.reduce_image_sequence('42061')
#analysis.photometry('42061')
#
#analysis.dump_dict()

analysis.load_dict()

analysis.plot_dither_2d('42061')
analysis.plot_dither_3d('42061')