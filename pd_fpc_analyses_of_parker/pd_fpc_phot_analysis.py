#%% -*- coding: utf-8 -*-

"""
Created on Sun Aug 28 00:26:37 2016

@author: Duan Yutong (dyt@physics.bu.edu, dyt@lbl.gov)

"""

import os
import numpy as np
# import glob
import matplotlib
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
# from astropy.io import fits

os.chdir(os.path.dirname(os.path.abspath(__file__)))
from pd_fpc_photometry import pd_fpc_photometry_dict

# Specify parent dir containing all datasets
# each dataset folder contains only one master image
dir_parent = r'/data/images/fpc/linear'
# photometry parameters
bkg_sigma = 3.0
source_snr = 3.0
fwhm_kernel = 20.0
x_size_kernel = 15
y_size_kernel = 15
clobber = True

# Fibre info
FIBRE_N = 3

#%% Example Analysis
# see how intensity and flux scales with exptime

photometry_dict = pd_fpc_photometry_dict(dir_parent,
                                         bkg_sigma = 3.0,
                                         source_snr = 3.0,
                                         fwhm_kernel = 20.0,
                                         x_size_kernel = 15,
                                         y_size_kernel = 15,
                                         clobber = True)

n_master = len(photometry_dict['id_master'])
expreq = np.array(photometry_dict['expreq'])

source_sum = np.empty([n_master, FIBRE_N])
flux = np.empty([n_master, FIBRE_N])
# fill in source_sum and flux lists
for i in range(len(n_master)):
    for j in range(FIBRE_N):
        # j runs in the order of fibre ids defined in pd_fpc_photometry
        source_sum[i,j] = photometry_dict['props_img'][i][j].source_sum
        flux[i,j] = source_sum[i,j]/expreq[i]

# make plots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12, 6))

[line1] = ax1.plot(expreq, source_sum[:,0],
                    c='r', marker='o', label = 'Fibre 1')
[line2] = ax1.plot(expreq, source_sum[:,1],
                    c='g', marker='x', label = 'Fibre 2')
[line3] = ax1.plot(expreq, source_sum[:,2],
                    c='b', marker='*', label = 'Fibre 3')
ax1.set_title('intensity - expreq')
ax1.set_xlabel('Exposure Time /s')
ax1.set_ylabel('Source Intensity Sum / ADU')
ax1.set_xlim(0, 6)
ax1.set_ylim(0, 1e7)
ax1.legend(handler_map={line1: HandlerLine2D(numpoints=3),
                        line2: HandlerLine2D(numpoints=3),
                        line3: HandlerLine2D(numpoints=3)},
           loc=4)

[line1] = ax2.plot(expreq, flux[:,0],
                    c='r', marker='o', label = 'Fibre 1')
[line2] = ax2.plot(expreq, flux[:,1],
                    c='g', marker='x', label = 'Fibre 2')
[line3] = ax2.plot(expreq, flux[:,2],
                    c='b', marker='*', label = 'Fibre 3')
ax2.set_title('flux - expreq')
ax2.set_xlabel('Exposure Time /s')
ax2.set_ylabel(r'Source Flux / ' +
               r'$\mathrm{} \cdot \mathrm{} ^{}$'.format(
                   '{ADU}', '{s}', '{-1}'))
ax2.set_xlim(0, 6)
ax2.set_ylim(0, 2e6)
ax2.legend(handler_map={line1: HandlerLine2D(numpoints=3),
                        line2: HandlerLine2D(numpoints=3),
                        line3: HandlerLine2D(numpoints=3)},
           loc=4)

# padding to avoid overlap between axes
plt.tight_layout()
# save
fig.savefig(os.path.join(dir_parent, 'source-exptime.png'), dpi=600)
pp = PdfPages(os.path.join(dir_parent, 'source-exptime.pdf'))
pp.savefig(fig)
pp.close()
