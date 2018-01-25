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
    
    filelist = ['pd_analysis-dither-trumpler37.py',
                'pd_analysis-flat.py',
                'pd_analysis-stability_dither-33022.py',
                'pd_analysis-stability_dither-42061.py',
                'pd_analysis-stability_dither-53001.py',
                'pd_analysis-stability_dither-53002.py',
                'pd_analysis-stability-53002.py',
                'pd_analysis-stability-hour3.py']
    
    print('Multiprocessing scripts {}'.format(repr(filelist)))
    with multiprocessing.Pool(8) as pool:
        pool.map(runfile, filelist)
