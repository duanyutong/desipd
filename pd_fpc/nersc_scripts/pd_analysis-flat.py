# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 02:02:35 2016

@author: givoltage
"""

from pd_fpc_analysis import FPCAnalysis

analysis = FPCAnalysis('dome_flat')
#
#analysis.CLOBBER = False
#analysis.DUMP_PICKLE = True
#
## specify dataset exptypes and expids
#analysis.specify_dataset('bias stack', 5517, 5527,
#                         'bias')
#analysis.specify_dataset('dark stack', 7165, 7175,
#                         '0.2dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 7132, 7142,
#                         '0.4dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 7097, 7107,
#                         '0.6dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 7075, 7085,
#                         '0.8dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 7064, 7074,
#                         '1dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 7031, 7041,
#                         '2dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 7020, 7030,
#                         '3dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 6987, 6997,
#                         '4dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 6969, 6979,
#                         '5dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 6936, 6946,
#                         '6dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 6925, 6935,
#                         '7dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 6896, 6906,
#                         '8dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 6863, 6873,
#                         '9dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 6172, 6182,
#                         '10dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 6834, 6844,
#                         '11dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 6823, 6833,
#                         '12dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 6783, 6793,
#                         '13dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 6772, 6782,
#                         '14dark',
#                         bias_name = 'bias')
#analysis.specify_dataset('dark stack', 6750, 6760,
#                         '15dark',
#                         bias_name = 'bias')
#
## image exposures
#analysis.specify_dataset('image stack', 7154, 7164,
#                         '0.2image',
#                         bias_name = 'bias',
#                         dark_name = '0.2dark')
#analysis.specify_dataset('image stack', 7121, 7131,
#                         '0.4image',
#                         bias_name = 'bias',
#                         dark_name = '0.4dark')
#analysis.specify_dataset('image stack', 7143, 7153,
#                         '0.6image',
#                         bias_name = 'bias',
#                         dark_name = '0.6dark')
#analysis.specify_dataset('image stack', 7086, 7096,
#                         '0.8image',
#                         bias_name = 'bias',
#                         dark_name = '0.8dark')
#analysis.specify_dataset('image stack', 7053, 7063,
#                         '1image',
#                         bias_name = 'bias',
#                         dark_name = '1dark')
#analysis.specify_dataset('image stack', 7042, 7052,
#                         '2image',
#                         bias_name = 'bias',
#                         dark_name = '2dark')
#analysis.specify_dataset('image stack', 7009, 7019,
#                         '3image',
#                         bias_name = 'bias',
#                         dark_name = '3dark')
#analysis.specify_dataset('image stack', 6998, 7008,
#                         '4image',
#                         bias_name = 'bias',
#                         dark_name = '4dark')
#analysis.specify_dataset('image stack', 6958, 6968,
#                         '5image',
#                         bias_name = 'bias',
#                         dark_name = '5dark')
#analysis.specify_dataset('image stack', 6947, 6957,
#                         '6image',
#                         bias_name = 'bias',
#                         dark_name = '6dark')
#analysis.specify_dataset('image stack', 6907, 6917,
#                         '7image',
#                         bias_name = 'bias',
#                         dark_name = '7dark')
#analysis.specify_dataset('image stack', 6885, 6895,
#                         '8image',
#                         bias_name = 'bias',
#                         dark_name = '8dark')
#analysis.specify_dataset('image stack', 6874, 6884,
#                         '9image',
#                         bias_name = 'bias',
#                         dark_name = '9dark')
#analysis.specify_dataset('image stack', 6161, 6171,
#                         '10image',
#                         bias_name = 'bias',
#                         dark_name = '10dark')
#analysis.specify_dataset('image stack', 6845, 6855,
#                         '11image',
#                         bias_name = 'bias',
#                         dark_name = '11dark')
#analysis.specify_dataset('image stack', 6812, 6822,
#                         '12image',
#                         bias_name = 'bias',
#                         dark_name = '12dark')
#analysis.specify_dataset('image stack', 6801, 6811,
#                         '13image',
#                         bias_name = 'bias',
#                         dark_name = '13dark')
#analysis.specify_dataset('image stack', 6761, 6771,
#                         '14image',
#                         bias_name = 'bias',
#                         dark_name = '14dark')
#analysis.specify_dataset('image stack', 6739, 6749,
#                         '15image',
#                         bias_name = 'bias',
#                         dark_name = '15dark')
#
#analysis.read_datasets_header()
## reduce datasets, order is important
#analysis.reduce_stack('bias')
#
#for dataset in ['0.2dark', '0.4dark', '0.6dark', '0.8dark', '1dark', 
#                '2dark', '3dark', '4dark', '5dark', '6dark', 
#                '7dark', '8dark', '9dark', '10dark', '11dark', 
#                '12dark', '13dark', '14dark', '15dark',
#                '0.2image', '0.4image', '0.6image', '0.8image', '1image', 
#                '2image', '3image', '4image', '5image', '6image', 
#                '7image', '8image', '9image', '10image', '11image', 
#                '12image', '13image', '14image', '15image']:
#    analysis.reduce_stack(dataset)
#    
#for dataset_image in ['0.2image', '0.4image', '0.6image', '0.8image', '1image', 
#                      '2image', '3image', '4image', '5image', '6image', 
#                      '7image', '8image', '9image', '10image', '11image', 
#                      '12image', '13image', '14image', '15image']:
#    analysis.photometry(dataset_image)
#
#analysis.dump_dict()

analysis.load_dict()
analysis.plot_flats(['0.2image', '0.4image', '0.6image', '0.8image', '1image', 
                      '2image', '3image', '4image', '5image', '6image', 
                      '7image', '8image', '9image', '10image', '11image', 
                      '12image', '13image', '14image', '15image'])
