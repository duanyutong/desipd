# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 12:50:50 2016

@author: givoltage
"""

import os
import glob

os.chdir(os.path.dirname(os.path.abspath(__file__)))
filelist = glob.glob('pd_analysis-*.py')
for file in filelist:
    print('Executing script: {}'.format(file))
    exec(open(file).read())