import os
import numpy as np
from pd_fpc_photometry_paf import pd_fpc_photometry_dict
from photometry_paf import segmentation_photometry
from photutils import properties_table
from astropy.io import fits
import pickle


directory = '/Volumes/data/images/fpc/' #'/Users/parkerf/Desktop/pd_fpc_analysis/3by3_3001'

def simple_photometry(lst_of_fits_files,dither):
    """
    Takes list of FPC exposure images, makes a folder to hold the segmentation results,    and returns a list of results from the segmentation photometry.

    INPUT: list of exposure IDs that identify the FPC images
    OUTPUT: python dictionary
        KEY: exposure ID
        ITEM: list of information about the photometry
              - path file
              - area:
              - Sum: sum of flux in pixels
              - exptime: from FITS header (requested exposure time)
              - background: mean background of image
              - x: x_centroid pixel value
              - y: x_centroid pixel value
              - flux: sum /exptime
    """


    d = {}
    print(lst_of_fits_files)
    list = np.array(lst_of_fits_files)
    print(len(list))
    for x in list:
        if dither == 'tel':
            id = x
            d[id] = {}
            filename = 'PROTODESI_FPC_%08d.fits' % x
        if dither == 'pos':
            filename = x
            id = os.path.splitext(x)[0][-3:]
            d[id] = {}
        print(filename)
        path_file_abs = directory+'/'+filename
        d[id]['path'] = filename
        
        [segm, props_list_return] = segmentation_photometry(path_file_abs,
                                                bkg_sigma = 3.0,
                                                source_snr = 3.0,
                                                fwhm_kernel = 2.0,
                                                x_size_kernel = 3,
                                                y_size_kernel = 3,
                                                clobber = False)

        FIBRE_IDS = {0:(1407, 1208), 1:(1209, 1072), 2: (1093, 1864)} #fibre names
        DIST_THRESHOLD = 20

        expreq = fits.open(path_file_abs)[0].header['EXPREQ']

        data = np.array(properties_table(props_list_return,columns=('area','source_sum','background_mean','xcentroid','ycentroid')))
        area = data['area']
        sume = data['source_sum']
        i = np.where(area == max(area))
        flux = sume[i]/expreq

        d[id]['area'] = area[i]
        d[id]['sum'] = sume[i]
        d[id]['exptime'] = expreq
        d[id]['background'] = data['background_mean'][i]
        d[id]['x'] = data['xcentroid'][i]
        d[id]['y'] = data['ycentroid'][i]
        d[id]['flux'] = flux
        

    return d

stability_list1 = (2745,2746,2747,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2763,2764)
stab_phot_dict = simple_photometry(stability_list1,'tel')
picle.dump(stab_phot_dict,("stability1.p","wb"))

#dither_list = os.listdir(directory)
#phot_dict = simple_photometry(dither_list,'pos')
#pickle.dump(phot_dict,open("positioner_dither.p","wb"))
#five_by_five_a = [2701,2692,2691,2713,2702,2700,2693,2690,2712,2703,2699,2694,2689,2711,2704,2698,2695,2709,2710,2705,2697,2696,2708,2707,2706]
#phot_dict = simple_photometry(five_by_five_a)
#pickle.dump( phot_dict, open( "five_by_five_a.p", "wb" ) )
#print("FIVE BY FIVE A IS READY TO GO")

#five_by_five_b = [2716,2717,2718,2719,2720,2725,2724,2723,2722,2721,2726,2727,2728,2729,2730,2735,2734,2715,2733,2732,2731,2736,2740,2741,2742,2743]
#phot_dict2 = simple_photometry(five_by_five_b)
#pickle.dump(phot_dict2,open( "five_by_five_b.p","wb"))
