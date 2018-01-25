# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 23:52:15 2016

@author: givoltage
"""

from pd_fpc_analysis import FPCAnalysis

analysis = FPCAnalysis('dither_offsets-53002')

#analysis.CLOBBER = False
#analysis.DUMP_PICKLE = True
#
## specify dataset exptypes and expids
#analysis.specify_dataset('bias stack', 5517, 5527,
#                         'bias')
#analysis.specify_dataset('dark stack', 6172, 6182,
#                         'dark',
#                         bias_name = 'bias')
#
## 5 datasets with manual offsets on 20160927
#analysis.specify_dataset('image sequence', 7663, 7689,
#                         '53002_offset_1',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '53002',
#                         focus = -9931,
#                         dither_mode = 'positioner',
#                         dither_grid = (5, 5),
#                         dither_step = 0.5)
#analysis.specify_dataset('image sequence', 7735, 7759,
#                         '53002_offset_2',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '53002',
#                         focus = -9931,
#                         dither_mode = 'positioner',
#                         dither_grid = (5, 5),
#                         dither_step = 1.0)
#analysis.specify_dataset('image sequence', 7761, 7785,
#                         '53002_offset_3',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '53002',
#                         focus = -9931,
#                         dither_mode = 'positioner',
#                         dither_grid = (5, 5),
#                         dither_step = 1.0)
##analysis.specify_dataset('image sequence', 7787, 7821,
##                         '53002_offset_4',
##                         bias_name = 'bias',
##                         dark_name = 'dark',
##                         tile_id = '53002',
##                         focus = -10050,
##                         dither_mode = 'positioner',
##                         dither_grid = (5, 5),
##                         dither_step = 1.0)
#analysis.specify_dataset('image sequence', 7823, 7847,
#                         '53002_offset_5',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '53002',
#                         focus = -9950,
#                         dither_mode = 'positioner',
#                         dither_grid = (5, 5),
#                         dither_step = 1.0)
#
#analysis.read_datasets_header()
#
## reduce datasets, order is important
#analysis.reduce_stack('bias')
#analysis.reduce_stack('dark')
#
#analysis.reduce_image_sequence('53002_offset_1')
#analysis.reduce_image_sequence('53002_offset_2')
#analysis.reduce_image_sequence('53002_offset_3')
##analysis.reduce_image_sequence('53002_offset_4')
#analysis.reduce_image_sequence('53002_offset_5')
#
#analysis.photometry('53002_offset_1')
#analysis.photometry('53002_offset_2')
#analysis.photometry('53002_offset_3')
##analysis.photometry('53002_offset_4')
#analysis.photometry('53002_offset_5')
#
#analysis.dump_dict()

analysis.load_dict()

#analysis.plot_dither_2d('53002_offset_1')
#analysis.plot_dither_2d('53002_offset_2')
#analysis.plot_dither_2d('53002_offset_3')
##analysis.plot_dither_2d('53002_offset_4')
#analysis.plot_dither_2d('53002_offset_5')

analysis.plot_dither_3d('53002_offset_1')
analysis.plot_dither_3d('53002_offset_2')
analysis.plot_dither_3d('53002_offset_3')
#analysis.plot_dither_3d('53002_offset_4')
analysis.plot_dither_3d('53002_offset_5')
