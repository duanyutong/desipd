# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 23:52:15 2016

@author: givoltage
"""

from pd_fpc_analysis import FPCAnalysis

analysis = FPCAnalysis('dither-53001')

#analysis.CLOBBER = False
#analysis.DUMP_PICKLE = True
#
## specify dataset exptypes and expids
#analysis.specify_dataset('bias stack', 5517, 5527,
#                         'bias')
#analysis.specify_dataset('dark stack', 6172, 6182,
#                         'dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('image sequence', 6263, 6287,
#                         '53001_1',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '53001',
#                         focus = -10254,
#                         dither_mode = 'telescope',
#                         dither_grid = (5, 5),
#                         dither_step = 2.0)
#analysis.specify_dataset('image sequence', 6334, 6358,
#                         '53001_2',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '53001',
#                         focus = -10244,
#                         dither_mode = 'telescope',
#                         dither_grid = (5, 5),
#                         dither_step = 2.0)
#analysis.specify_dataset('image sequence', 6368, 6392,
#                         '53001_3',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '53001',
#                         focus = -10241,
#                         dither_mode = 'telescope',
#                         dither_grid = (5, 5),
#                         dither_step = 2.0)
#analysis.read_datasets_header()
#
## reduce datasets, order is important
#analysis.reduce_stack('bias')
#analysis.reduce_stack('dark')
#analysis.reduce_image_sequence('53001_1')
#analysis.reduce_image_sequence('53001_2')
#analysis.reduce_image_sequence('53001_3')
#analysis.photometry('53001_1')
#analysis.photometry('53001_2')
#analysis.photometry('53001_3')
#analysis.dump_dict()

analysis.load_dict()
analysis.plot_dither_2d('53001_1')
analysis.plot_dither_2d('53001_2')
analysis.plot_dither_2d('53001_3')

analysis.plot_dither_3d('53001_1')
analysis.plot_dither_3d('53001_2')
analysis.plot_dither_3d('53001_3')
