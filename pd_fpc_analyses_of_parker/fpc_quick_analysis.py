import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
from pd_phot_tools import PhotTools

directory = '/Volumes/data/images/fpc' #'/Volumes/data/images/fits/' #'/Users/parkerf/Desktop/pd_fpc_analysis/3by3_3001'


#### INPUTS ####
FIRST_EXP_ID = 5047 
FINAL_EXP_ID = 5071
SEQ_TYPE = 'telescope' # 'positioner', 'stability', 'single'
STEP_SIZE = 2.5 #arcsec. No need to change if stability or single
FIBER = 3001 #3002, 3002, 'all'
FIELD = '33011'
DATE = 20160921

#### SCRIPT ####
list = np.arange(FIRST_EXP_ID, FINAL_EXP_ID+1,1)
name = "%s_%s_%s" % (FIELD,SEQ_TYPE,str(DATE))
if FIBER == 'all':
    match = 'location'
else: 
    match = 'max'

if SEQ_TYPE == 'telescope':
    PhotTools('teldither', list, False, match,STEP_SIZE, FIBER, directory,name).run()
elif SEQ_TYPE == 'positioner':
    PhotTools('first_posdither', list, False, match, STEP_SIZE, FIBER, directory,name).run()
elif SEQ_TYPE == 'stability':
    PhotTools('first_stability', list, False, match, STEP_SIZE, FIBER, directory,name).run()
elif SEQ_TYPE == 'single':
    PhotTools('one', list, False, match, STEP_SIZE, FIBER, directory,name).run()
else:
    print("THAT ISN'T A VALID SEQUENCE TYPE. PLEASE TRY AGAIN") 
