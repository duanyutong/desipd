# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits

os.chdir('/Volumes/data/images/gfa/')
filefull = np.fromfile('the_answer.dat', dtype=int)
filefull_new = np.fromfile('the_answer.dat', dtype=int)
file_new.where(file_new==913)
file = filefull[len(filefull)-2340576*2 : len(filefull)-2340576*1]
image = np.reshape(file, (516, 1134*4))
# plt.imshow(image)
image_cropped = image[:,0:1134]
hdu = fits.PrimaryHDU(image)
hdu_cropped = fits.PrimaryHDU(image_cropped)
hdu.writeto('/Users/msdos/the_answer3.fits', clobber=True)
hdu_cropped.writeto('/Users/msdos/the_answer_cropped3.fits', clobber=True)

image_cropped_shifted = np.array(image_cropped)
for rownumber in range(516):
    # find value > 60k
    rowdata = image_cropped[rownumber, :]
    indices = [i for i, v in enumerate(rowdata>60000) if v]
    print('number of pixels > 60k is: {}'.format(len(indices)))
    if len(indices)>0:
        index = indices[0]
        index_array_initial= np.arange(index, 1134, 1)
        
        index_array_append = np.arange(0, index, 1)
        index_array = np.append(index_array_initial, index_array_append)
        if len(index_array) == 1134:
            rowdata_new = rowdata[index_array]
            image_cropped_shifted[rownumber, :] = rowdata_new
            print('row {} shifted'.format(rownumber))
        else:
            print('wrong index array dimension')
        
hdu_cropped_shifted = fits.PrimaryHDU(image_cropped_shifted)
hdu_cropped_shifted.writeto('/Users/msdos/the_answer_cropped_shifted_q3.fits', clobber=True)
'''
exact 4 amps for 4 quadrants
test
5390  1281  1412  1622  1712  1790  1646  1620  1526  1793]
'''

search = np.where(filefull_new==65390)[0]
searchlist = list(search)
for index in searchlist:
    print(filefull_new[index:index+10])
index = searchlist[7]
print(filefull_new[index:index+10])
datafull = filefull_new[index:index+2340576]
datareshape = np.reshape(datafull, (516, 1134*4))
for ampind in range(4):
    ampno = ampind+1
    print('processing amplifier: {}'.format(ampno))
    ampimage = datareshape[:,1134*ampind:1134*(ampind+1)]
    # perform shift
    ampimage_shifted = np.array(ampimage)
    for rownumber in range(516):
        # find value > 60k
        rowdata = ampimage[rownumber, :]
        indices = [i for i, v in enumerate(rowdata>60000) if v]
        print('number of pixels > 60k is: {}'.format(len(indices)))
        if len(indices)>0:
            index = indices[0]
            index_array_initial= np.arange(index, 1134, 1)
            index_array_append = np.arange(0, index, 1)
            index_array = np.append(index_array_initial, index_array_append)
            if len(index_array) == 1134:
                rowdata_new = rowdata[index_array]
                ampimage_shifted[rownumber, :] = rowdata_new
                print('row {} shifted'.format(rownumber))
            else:
                print('wrong index array dimension')
    # save
    amphdu = fits.PrimaryHDU(ampimage)
    amphdu_shifted = fits.PrimaryHDU(ampimage_shifted)
    amphdu.writeto('/Users/msdos/PROTODESI_GFA0G_00001000_0011_AMP_'+str(ampno)+'.fits')
    amphdu_shifted.writeto('/Users/msdos/PROTODESI_GFA0G_00001000_0011_AMP_'+str(ampno)+'_shifted.fits')
    
'''
Reassemble 4 quadrants
'''
quadrant_mapping = np.array([[4, 3],
                             [1, 2]])
os.chdir(r'C:\Users\givoltage\Google Drive\DESI\protoDESI\images\kpno\gfa_shift')
amp1 = fits.open('PROTODESI_GFA0G_00001000_0011_AMP_1_SHIFTED.fits')[0].data
amp2 = fits.open('PROTODESI_GFA0G_00001000_0011_AMP_2_SHIFTED.fits')[0].data
amp3 = fits.open('PROTODESI_GFA0G_00001000_0011_AMP_3_SHIFTED.fits')[0].data
amp4 = fits.open('PROTODESI_GFA0G_00001000_0011_AMP_4_SHIFTED.fits')[0].data
ampdata = np.array([amp1, amp2, amp3, amp4])
q1 = ampdata[quadrant_mapping[0,1]-1]
q2 = ampdata[quadrant_mapping[0,0]-1]
q3 = ampdata[quadrant_mapping[1,0]-1]
q4 = ampdata[quadrant_mapping[1,1]-1]
q1 = np.fliplr(np.flipud(q1))
q2 = np.flipud(q2)
q4 = np.fliplr(q4)
image = np.concatenate([np.concatenate([q3, q4], axis=1),
                        np.concatenate([q2, q1], axis=1)])
hdu = fits.PrimaryHDU(image)
hdu.writeto('PROTODESI_GFA0G_00001000_0011_RECONSTRUCTED_SHIFTED.fits')


#======================

quadrant_mapping = np.array([[4, 3],
                             [1, 2]])
os.chdir(r'C:\Users\givoltage\Google Drive\DESI\protoDESI\images\kpno\gfa_shift')
amp1 = fits.open('PROTODESI_GFA0G_00001000_0011_AMP_1.fits')[0].data
amp2 = fits.open('PROTODESI_GFA0G_00001000_0011_AMP_2.fits')[0].data
amp3 = fits.open('PROTODESI_GFA0G_00001000_0011_AMP_3.fits')[0].data
amp4 = fits.open('PROTODESI_GFA0G_00001000_0011_AMP_4.fits')[0].data
ampdata = np.array([amp1, amp2, amp3, amp4])
q1 = ampdata[quadrant_mapping[0,1]-1]
q2 = ampdata[quadrant_mapping[0,0]-1]
q3 = ampdata[quadrant_mapping[1,0]-1]
q4 = ampdata[quadrant_mapping[1,1]-1]
q1 = np.fliplr(np.flipud(q1))
q2 = np.flipud(q2)
q4 = np.fliplr(q4)
image = np.concatenate([np.concatenate([q3, q4], axis=1),
                        np.concatenate([q2, q1], axis=1)])
hdu = fits.PrimaryHDU(image)
hdu.writeto('PROTODESI_GFA0G_00001000_0011_RECONSTRUCTED.fits')
    
''' from fits'''

import numpy as np
import os
from astropy.io import fits

os.chdir('/Volumes/data/images/gfa/')
hdu = fits.open('PROTODESI_GFA0G_00001000_0009.fits')[0]
image = hdu.data
# plt.imshow(image)
image_cropped = image[0:516, 0:1024]
hdu = fits.PrimaryHDU(image)
hdu_cropped = fits.PrimaryHDU(image_cropped)
# hdu.writeto('/Users/msdos/the_answer3.fits', clobber=True)
hdu_cropped.writeto('/Users/msdos/the_answer_cropped3.fits', clobber=True)


image_cropped_shifted = np.array(image_cropped)
for rownumber in range(516):
    # find value > 60k
    rowdata = image_cropped[rownumber, :]
    indices = [i for i, v in enumerate(rowdata>60000) if v]
    print('number of pixels > 60k is: {}'.format(len(indices)))
    if len(indices)>0:
        index = indices[0]
        index_array_initial= np.arange(index, 1134, 1)
        
        index_array_append = np.arange(0, index, 1)
        index_array = np.append(index_array_initial, index_array_append)
        if len(index_array) == 1134:
            rowdata_new = rowdata[index_array]
            image_cropped_shifted[rownumber, :] = rowdata_new
            print('row {} shifted'.format(rownumber))
        else:
            print('wrong index array dimension')
        
hdu_cropped_shifted = fits.PrimaryHDU(image_cropped_shifted)
hdu_cropped_shifted.writeto('/Users/msdos/the_answer_cropped_shifted3.fits', clobber=True)
