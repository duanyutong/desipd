# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 23:52:15 2016

@author: givoltage
"""

from pd_fpc_analysis import FPCAnalysis

analysis = FPCAnalysis('dither-hour3')

#analysis.CLOBBER = False
#analysis.DUMP_PICKLE = True
#
## specify dataset exptypes and expids
#analysis.specify_dataset('bias stack', 5517, 5527,
#                         'bias')
#analysis.specify_dataset('dark stack', 5852, 5860,
#                         'dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('image sequence', 5896, 5920,
#                         'hour3_1',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = 'hour3',
#                         focus = -10078,
#                         dither_mode = 'telescope',
#                         dither_grid = (5, 5),
#                         dither_step = 1.0)
#analysis.specify_dataset('image sequence', 5924, 5948,
#                         'hour3_2',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = 'hour3',
#                         focus = -10072,
#                         dither_mode = 'telescope',
#                         dither_grid = (5, 5),
#                         dither_step = 1.0)
#analysis.specify_dataset('image sequence', 5985, 6009,
#                         'hour3_3',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = 'hour3',
#                         focus = -10072,
#                         dither_mode = 'telescope',
#                         dither_grid = (5, 5),
#                         dither_step = 1.0)
#
#
#
#analysis.read_datasets_header()
#
## reduce datasets, order is important
#analysis.reduce_stack('bias')
#analysis.reduce_stack('dark')
#analysis.reduce_image_sequence('hour3_1')
#analysis.reduce_image_sequence('hour3_2')
#analysis.reduce_image_sequence('hour3_3')
#
#analysis.photometry('hour3_1')
#analysis.photometry('hour3_2')
#analysis.photometry('hour3_3')
#
#analysis.dump_dict()

analysis.load_dict()


analysis.plot_dither_2d('hour3_1')
analysis.plot_dither_2d('hour3_2')
analysis.plot_dither_2d('hour3_3')

analysis.plot_dither_3d('hour3_1')
analysis.plot_dither_3d('hour3_2')
analysis.plot_dither_3d('hour3_3')