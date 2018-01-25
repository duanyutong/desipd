# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 23:52:15 2016

@author: givoltage
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from pd_fpc_analysis import FPCAnalysis

analysis = FPCAnalysis('dither-trumpler37')

#analysis.CLOBBER = False
#analysis.DUMP_PICKLE = True
#
## specify dataset exptypes and expids
#analysis.specify_dataset('bias stack', 5517, 5527,
#                         'bias')
#analysis.specify_dataset('dark stack', 6172, 6182,
#                         'dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('image sequence', 3232, 3256,
#                         'trumpler37_1_tel',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = 'Trumpler37',
#                         focus = -10800,
#                         dither_mode = 'telescope',
#                         dither_grid = (5, 5),
#                         dither_step = 1.0)
#analysis.specify_dataset('image sequence', 3257, 3281,
#                         'trumpler37_1_pos',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = 'Trumpler37',
#                         focus = -10800,
#                         dither_mode = 'positioner',
#                         dither_grid = (5, 5),
#                         dither_step = 0.5)
#analysis.specify_dataset('image sequence', 4044, 4068,
#                         'trumpler37_2_tel',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = 'Trumpler37',
#                         focus = -11098,
#                         dither_mode = 'telescope',
#                         dither_grid = (5, 5),
#                         dither_step = 1.0)
#analysis.specify_dataset('image sequence', 4069, 4093,
#                         'trumpler37_2_pos',
#                         bias_name = 'bias',
#                         dark_name = 'dark',
#                         tile_id = 'Trumpler37',
#                         focus = -11098,
#                         dither_mode = 'positioner',
#                         dither_grid = (5, 5),
#                         dither_step = 0.5)
#analysis.read_datasets_header()
#
## reduce datasets, order is important
#analysis.reduce_stack('bias')
#analysis.reduce_stack('dark')
#
#analysis.reduce_image_sequence('trumpler37_1_tel')
#analysis.reduce_image_sequence('trumpler37_1_pos')
#analysis.reduce_image_sequence('trumpler37_2_tel')
#analysis.reduce_image_sequence('trumpler37_2_pos')
#
#analysis.photometry('trumpler37_1_tel')
#analysis.photometry('trumpler37_1_pos')
#analysis.photometry('trumpler37_2_tel')
#analysis.photometry('trumpler37_2_pos')
#analysis.dump_dict()

analysis.load_dict()
analysis.plot_dither_2d('trumpler37_1_tel')
analysis.plot_dither_3d('trumpler37_1_tel')
analysis.plot_dither_2d('trumpler37_1_pos')
analysis.plot_dither_3d('trumpler37_1_pos')
analysis.plot_dither_2d('trumpler37_2_tel')
analysis.plot_dither_3d('trumpler37_2_tel')
analysis.plot_dither_2d('trumpler37_2_pos')
analysis.plot_dither_3d('trumpler37_2_pos')
