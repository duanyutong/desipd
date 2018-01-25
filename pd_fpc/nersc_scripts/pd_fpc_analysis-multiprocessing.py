# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 12:50:50 2016

@author: givoltage
"""

import os
import glob
import multiprocessing

def runfile(file):
    exec(open(file).read())

if __name__ == '__main__':
    
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    filelist = glob.glob('pd_analysis-*.py')
    print('Multiprocessing scripts {}'.format(repr(filelist)))
    with multiprocessing.Pool(8) as pool:
        pool.map(runfile, filelist)
