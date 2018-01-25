#!/usr/bin/env python

# argument is the cvs file that lists the exposures in a fvc_stability test
# display the first one and the output of all the red files

import numpy as np
import matplotlib.pyplot as plt
import sys,os,time
import subprocess
import re

data_dir = "/data/fvc/"

filelist = sys.argv[1]
regionfile = 'fvcstab.reg'

flist = []
with open(filelist, 'r') as fcsv:
    for aline in fcsv:
        if re.search('red', aline) != None:
            flist.append(aline)

clist = ['green','red','blue','cyan','magenta','yellow','black','white'] 
regionf = open(regionfile, 'w')

for idx in range(len(flist)):
    redfilename = flist[idx]
    redfilename = re.sub('\n', '', redfilename)
    if idx == 0:
        imfilename = re.sub('red', 'fits', redfilename, count=1)
        os.system("scp msdos@desifvc.kpno.noao.edu:%s ." % (data_dir+'/'+imfilename))
    os.system("scp msdos@desifvc.kpno.noao.edu:%s ." % (data_dir+'/'+redfilename))
    data = np.genfromtxt(redfilename, usecols=(0,1,2))
    cval = clist[idx%len(clist)]
    with open(redfilename, 'r') as fobj:
        for row in fobj:
            rowl = row.split()
            id = int(rowl[0])
            if (id > 100): # this is a fiducual or a fiber
                outrow = str("circle %0.2f %0.2f 15.0 # color = %s\n" % (float(rowl[1]), float(rowl[2]),cval))
                regionf.write(outrow)
                
    fobj.close

regionf.close
ds9str = 'ds9 ' + imfilename + ' -scale limits 1800 3000 -regions load ' + regionfile +' &'
subprocess.run(ds9str,shell=True)
time.sleep(2)
os.system('rm fvcstab.reg')
for fname in flist:
    os.system('rm '+fname)
