# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 12:50:50 2016

@author: givoltage
"""

import os
import multiprocessing

def runfile(file):
    exec(open(file).read())

if __name__ == '__main__':
    
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    filelist = ['pd_analysis-dither-33001.py',
                'pd_analysis-dither-33002.py',
                'pd_analysis-dither-33003.py',
                'pd_analysis-dither-33011.py',
                'pd_analysis-dither-33012.py',
                'pd_analysis-dither-33013.py']
    
    print('Multiprocessing scripts {}'.format(repr(filelist)))
    with multiprocessing.Pool(8) as pool:
        pool.map(runfile, filelist)
