#!/usr/bin/env python

# argument is 'fvc.<expnum>' without the suffix
# read in a fvc red file, make regions and show in DS9
import numpy as np
import matplotlib.pyplot as plt
import sys,os,time
import subprocess

data_dir = "/data/fvc/"

redfilename =  sys.argv[1] + '.red'
imfilename = sys.argv[1] + '.fits'

os.system("scp msdos@desifvc.kpno.noao.edu:%s ." % (data_dir+'/'+redfilename))
os.system("scp msdos@desifvc.kpno.noao.edu:%s ." % (data_dir+'/'+imfilename))

data = np.genfromtxt(redfilename, usecols=(0,1,2))

regionfile = 'fvcred.reg'

outf = open(regionfile, 'w')
outf.write('global color=green\n')

with open(redfilename, 'r') as fobj:
    for row in fobj:
        rowl = row.split()
        id = int(rowl[0])
        if (id > 100): # this is a fiducual or a fiber
            outrow = str("circle %0.2f %0.2f 15.0\n" % (float(rowl[1]), float(rowl[2])))
            outf.write(outrow)

fobj.close
outf.close

ds9str = 'ds9 ' + imfilename + ' -scale limits 1800 3000 -regions load ' + regionfile +' &'
subprocess.run(ds9str,shell=True)
# Needs the & in order to actually load the region file. Don't understand.

# can't clean up or the region display wont work. Don't undertand. 
#cleanstr = "rm " + sys.argv[1] + ".*"
#os.system(cleanstr)
#os.system("rm fvcred.reg")


