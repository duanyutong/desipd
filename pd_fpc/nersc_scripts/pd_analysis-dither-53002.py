# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 23:52:15 2016

@author: givoltage
"""

from pd_fpc_analysis import FPCAnalysis

analysis = FPCAnalysis('dither-53002')

#analysis.CLOBBER = False
#analysis.DUMP_PICKLE = True
#
## specify dataset exptypes and expids
#analysis.specify_dataset('bias stack', 5517, 5527,
#                         'bias')
#analysis.specify_dataset('dark stack', 6172, 6182,
#                         'dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('image sequence', 6394, 6418,
#                         '53002_1',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '53002',
#                         focus = -10233,
#                         dither_mode = 'telescope',
#                         dither_grid = (5, 5),
#                         dither_step = 2.0)
#analysis.specify_dataset('image sequence', 6456, 6480,
#                         '53002_2',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '53002',
#                         focus = -10228,
#                         dither_mode = 'telescope',
#                         dither_grid = (5, 5),
#                         dither_step = 2.0)
#analysis.specify_dataset('image sequence', 6486, 6510,
#                         '53002_3',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '53002',
#                         focus = -10228,
#                         dither_mode = 'telescope',
#                         dither_grid = (5, 5),
#                         dither_step = 2.0)
#analysis.specify_dataset('image sequence', 7566, 7590,
#                         '53002_4',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '53002',
#                         focus = -9931,
#                         dither_mode = 'positioner',
#                         dither_grid = (5, 5),
#                         dither_step = 1.0)
#analysis.specify_dataset('image sequence', 7637, 7661,
#                         '53002_5',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = '53002',
#                         focus = -9931,
#                         dither_mode = 'telescope',
#                         dither_grid = (5, 5),
#                         dither_step = 2.0)
#
#analysis.read_datasets_header()
#
## reduce datasets, order is important
#analysis.reduce_stack('bias')
#analysis.reduce_stack('dark')
#
#analysis.reduce_image_sequence('53002_1')
#analysis.reduce_image_sequence('53002_2')
#analysis.reduce_image_sequence('53002_3')
#analysis.reduce_image_sequence('53002_4')
#analysis.reduce_image_sequence('53002_5')
#
#analysis.photometry('53002_1')
#analysis.photometry('53002_2')
#analysis.photometry('53002_3')
#analysis.photometry('53002_4')
#analysis.photometry('53002_5')
#
#analysis.dump_dict()

analysis.load_dict()

#analysis.plot_dither_2d('53002_1')
#analysis.plot_dither_2d('53002_2')
#analysis.plot_dither_2d('53002_3')
#analysis.plot_dither_2d('53002_4')
#analysis.plot_dither_2d('53002_5')

analysis.plot_dither_3d('53002_1')
analysis.plot_dither_3d('53002_2')
analysis.plot_dither_3d('53002_3')
analysis.plot_dither_3d('53002_4')
analysis.plot_dither_3d('53002_5')

