# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 02:02:35 2016

@author: givoltage
"""

from pd_fpc_analysis import FPCAnalysis

analysis = FPCAnalysis('stability_dither-53002')
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
#analysis.specify_dataset('image sequence', 5541, 5565,
#                         '53002_tel_1',
#                         tile_id = '53002',
#                         bias_name = 'bias',
#                         dark_name = 'dark',)
#analysis.specify_dataset('image sequence', 5620, 5644,
#                         '53002_tel_2',
#                         tile_id = '53002',
#                         bias_name = 'bias',
#                         dark_name = 'dark',)
#
#analysis.read_datasets_header()
## reduce datasets, order is important
#analysis.reduce_stack('bias')
#analysis.reduce_stack('dark')
#analysis.reduce_image_sequence('53002_tel_1')
#analysis.reduce_image_sequence('53002_tel_2')
#analysis.photometry('53002_tel_1')
#analysis.photometry('53002_tel_2')
#analysis.dump_dict()

analysis.load_dict()
analysis.plot_stability('53002_tel_1')
analysis.plot_stability('53002_tel_2')
