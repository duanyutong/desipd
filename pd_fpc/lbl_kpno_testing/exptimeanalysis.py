# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 13:47:07 2016

@author: givoltage
"""

import os
from astropy.io import fits
from astropy.io import ascii
import numpy as np
#from astropy.table import Table, Column
from photutils import detect_sources, source_properties, properties_table
import glob

# need dvipng, ghostscript installed and set in path (with latex distribution)
# need type1cm and zhmetrics packages


# get labels of hot pixels from master bias and master dark
# pixel label, x coord, y coord, master bias, sd, master dark, sd, master object, sd,

datasetDir = r'C:\Users\givoltage\Downloads\20160718'
path_pattern = os.path.join(datasetDir, '*.fit*')
paths = glob.glob(path_pattern)

exptime = np.arange(0.1,1.1,0.1)
pixels = [[1345, 1248],
          [ 345, 1735],
          [2079,  100],
          [ 286, 2422],
          [2248, 1354],
          [2843,  414],
          [ 476,  989]]
dict = {'exptime': exptime}
for i in range(len(pixels)):
    dict[i] = []

for path in paths:
    data = fits.open(path)[0].data
    for i in range(len(pixels)):
        x = pixels[i][0]
        y = pixels[i][1]
        value = data[y,x]
        print(path)
        print(pixels[i])
        print('value: ' + repr(value))
        dict[i].append(value)
print(dict)

#segm_bias = detect_sources(m_bias, 1300, npixels=1)
#segm_dark = detect_sources(m_dark, 520, npixels=1)
#
#props_bias = source_properties(m_bias, segm_bias)
#props_dark = source_properties(m_dark, segm_dark)
#
## propsTableColumns = ['id', 'xcentroid', 'ycentroid', 'max_value', 'min_value','source_sum']
#table_bias = properties_table(props_bias)
#table_dark = properties_table(props_dark)
#
#ascii.write(table_bias, 'table_bias.csv', format = 'csv')
#ascii.write(table_dark, 'table_dark.csv', format = 'csv')
##log_bias = open('labels_bias.log', 'w')
##log_dark = open('labels_dark.log', 'w')
##
##log_bias.write(log_bias_text)
##log_dark.write(log_dark_text)
##
##log_bias.close()
##log_dark.close()