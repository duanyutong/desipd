# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 23:52:15 2016

@author: givoltage
"""

from pd_fpc_analysis import FPCAnalysis

analysis = FPCAnalysis('dither-33012')
#
#analysis.CLOBBER = False
#analysis.DUMP_PICKLE = True
#
## specify dataset exptypes and expids
#analysis.specify_dataset('bias stack', 5517, 5527,
#                         'bias')
#analysis.specify_dataset('dark stack', 6172, 6182,
#                         'dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('image sequence', 3911, 3935,
#                         '33012_tel',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '33012',
#                         focus = -11167,
#                         dither_mode = 'telescope',
#                         dither_grid = (5, 5),
#                         dither_step = 1.0)
#analysis.specify_dataset('image sequence', 3936, 3960,
#                         '33012_pos',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '33012',
#                         focus = -11167,
#                         dither_mode = 'positioner',
#                         dither_grid = (5, 5),
#                         dither_step = 0.5)
#analysis.read_datasets_header()
#
## reduce datasets, order is important
#analysis.reduce_stack('bias')
#analysis.reduce_stack('dark')
#analysis.reduce_image_sequence('33012_tel')
#analysis.reduce_image_sequence('33012_pos')
#analysis.photometry('33012_tel')
#analysis.photometry('33012_pos')
#analysis.dump_dict()

analysis.load_dict()
analysis.plot_dither_2d('33012_tel')
analysis.plot_dither_3d('33012_tel')
analysis.plot_dither_2d('33012_pos')
analysis.plot_dither_3d('33012_pos')
