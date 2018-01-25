# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 02:02:35 2016

@author: givoltage
"""

from pd_fpc_analysis import FPCAnalysis

analysis = FPCAnalysis('stability_dither-33022')
#
#analysis.CLOBBER = False
#analysis.DUMP_PICKLE = True
#
## specify dataset exptypes and expids
#analysis.specify_dataset('bias stack', 5517, 5527,
#                         'bias')
#analysis.specify_dataset('dark stack', 5852, 5860,
#                         'dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('image sequence', 3448, 3472,
#                         '33022',
#                         tile_id = '33022',
#                         bias_name = 'bias',
#                         dark_name = 'dark',)
#
#analysis.read_datasets_header()
## reduce datasets, order is important
#analysis.reduce_stack('bias')
#analysis.reduce_stack('dark')
#analysis.reduce_image_sequence('33022')
#
#analysis.photometry('33022')
#
#analysis.dump_dict()

analysis.load_dict()
analysis.plot_stability('33022')

