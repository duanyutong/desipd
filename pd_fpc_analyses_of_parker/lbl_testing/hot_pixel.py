# -*- coding: utf-8 -*-
"""
Hot pixel analysis
Created on Mon Jul 18 10:15:51 2016

@author: givoltage
"""

import os
from astropy.io import fits
from astropy.io import ascii
#from astropy.table import Table, Column
from photutils import detect_sources, source_properties, properties_table


# get labels of hot pixels from master bias and master dark
# pixel label, x coord, y coord, master bias, sd, master dark, sd, master object, sd, 

os.chdir(r'C:\Users\givoltage\Google Drive\DESI\protoDESI\ssl\20160524')

filepath_bias = r'C:\Users\givoltage\Downloads\20160718_-15c\1.0.fit'
filepath_dark = r'C:\Users\givoltage\Google Drive\DESI\protoDESI\ssl\20160524\1sec_dark_009.FIT'

m_bias = fits.open(filepath_bias)[0].data
m_dark = fits.open(filepath_dark)[0].data

segm_bias = detect_sources(m_bias, 1300, npixels=1)
segm_dark = detect_sources(m_dark, 1300, npixels=1)

props_bias = source_properties(m_bias, segm_bias)
props_dark = source_properties(m_dark, segm_dark)

# propsTableColumns = ['id', 'xcentroid', 'ycentroid', 'max_value', 'min_value','source_sum']
table_bias = properties_table(props_bias)
table_dark = properties_table(props_dark)

ascii.write(table_bias, '1.0.csv', format = 'csv')
ascii.write(table_dark, 'protodesi-FPC_00000009(1).csv', format = 'csv')
#log_bias = open('labels_bias.log', 'w')
#log_dark = open('labels_dark.log', 'w')
#
#log_bias.write(log_bias_text)
#log_dark.write(log_dark_text)
#
#log_bias.close()
#log_dark.close()