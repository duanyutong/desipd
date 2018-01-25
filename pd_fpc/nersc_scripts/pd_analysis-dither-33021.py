# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 23:52:15 2016

@author: givoltage
"""

from pd_fpc_analysis import FPCAnalysis

analysis = FPCAnalysis('dither-33021')

#analysis.CLOBBER = False
#analysis.DUMP_PICKLE = True
#
## specify dataset exptypes and expids
#analysis.specify_dataset('bias stack', 5517, 5527,
#                         'bias')
#analysis.specify_dataset('dark stack', 6172, 6182,
#                         'dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('image sequence', 4102, 4126,
#                         '33021_11166',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '33021',
#                         focus = -11166,
#                         dither_mode = 'telescope',
#                         dither_grid = (5, 5),
#                         dither_step = 1.0)
#analysis.specify_dataset('image sequence', 6582, 6606,
#                         '33021_10283',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '33021',
#                         focus = -10283,
#                         dither_mode = 'positioner',
#                         dither_grid = (5, 5),
#                         dither_step = 1.0)
#analysis.specify_dataset('image sequence', 6677, 6701,
#                         '33021_10394',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '33021',
#                         focus = -10394,
#                         dither_mode = 'positioner',
#                         dither_grid = (5, 5),
#                         dither_step = 1.0)
#analysis.specify_dataset('image sequence', 6703, 6727,
#                         '33021_10474',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '33021',
#                         focus = -10474,
#                         dither_mode = 'positioner',
#                         dither_grid = (5, 5),
#                         dither_step = 1.0)
#analysis.read_datasets_header()
#
## reduce datasets, order is important
#analysis.reduce_stack('bias')
#analysis.reduce_stack('dark')
#analysis.reduce_image_sequence('33021_11166')
#analysis.reduce_image_sequence('33021_10283')
#analysis.reduce_image_sequence('33021_10394')
#analysis.reduce_image_sequence('33021_10474')
#analysis.photometry('33021_11166')
#analysis.photometry('33021_10283')
#analysis.photometry('33021_10394')
#analysis.photometry('33021_10474')
#analysis.dump_dict()

analysis.load_dict()

#analysis.plot_dither_2d('33021_11166')
#analysis.plot_dither_2d('33021_10283')
#analysis.plot_dither_2d('33021_10394')
#analysis.plot_dither_2d('33021_10474')

analysis.plot_dither_3d('33021_11166')
analysis.plot_dither_3d('33021_10283')
analysis.plot_dither_3d('33021_10394')
analysis.plot_dither_3d('33021_10474')
