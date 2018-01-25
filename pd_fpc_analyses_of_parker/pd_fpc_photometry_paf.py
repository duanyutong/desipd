# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 00:48:30 2016

@author: Duan Yutong (dyt@physics.bu.edu, dyt@lbl.gov)

"""

def pd_fpc_photometry_dict(dir_parent,
                           bkg_sigma = 3.0,
                           source_snr = 3.0,
                           fwhm_kernel = 2.0,
                           x_size_kernel = 3,
                           y_size_kernel = 3,
                           clobber = False):

    import os
    import numpy as np
    import glob
    #import matplotlib
    #matplotlib.rcParams['text.usetex'] = True
    #matplotlib.rcParams['text.latex.unicode'] = True
    #from matplotlib.backends.backend_pdf import PdfPages
    #import matplotlib.pyplot as plt
    #from matplotlib.legend_handler import HandlerLine2D
    from astropy.io import fits

    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    from photometry import segmentation_photometry

    #%% Define fibre identification parameters

    FIBRE_IDS = [0, 1, 2] #fibre names
    # Change fibre ID mapping here in coordinates
    # used for identifying known fibres
    FIBRE_COORDS = [(1208, 1407),
                    (1205, 1092),
                    (1864, 1093)]

    DIST_THRESHOLD = 20 # pixels to known centroid
    FIBRE_N = len(FIBRE_IDS)

    #%% start with globbing
    dirs_dataset = glob.glob(
                        os.path.abspath(
                            dir_parent + os.path.sep + '*')
                        + os.path.sep)
    dataset_ids = [os.path.basename(os.path.normpath(dir_dataset))
                    for dir_dataset in dirs_dataset]

    # build a dictionary of quantities shared by all fibres in master fits
    photometry_dict = {'id_master': dataset_ids,
                       'props_img': [None]*len(dataset_ids),
                       'expreq': [None]*len(dataset_ids),
                       'background': [None]*len(dataset_ids)}

    # photometry for each master image
    for i in range(len(dirs_dataset)):

        # processing only one master object image now
        path_file_abs = glob.glob(os.path.join(
                            dirs_dataset[i], 'master*object*.fit*'))[0]
        # add exptime
        try:
            expreq = fits.open(path_file_abs)[0].header['EXPREQ']
            photometry_dict['expreq'][i] = expreq
            print('EXPREQ = {} acquired from {}'
                    .format(expreq, os.path.dirname(path_file_abs)))
        except:
            print('EXPREQ not read from {}'
                    .format(expreq, os.path.dirname(path_file_abs)))
        # in props_list, each source produces a props object
        [segm, props_list_return] = segmentation_photometry(path_file_abs,
                                                bkg_sigma = bkg_sigma,
                                                source_snr = source_snr,
                                                fwhm_kernel = fwhm_kernel,
                                                x_size_kernel = x_size_kernel,
                                                y_size_kernel = y_size_kernel,
                                                clobber = clobber)

        # there may be more than 3 sources returned, depending on parameters
        # for each of the fibres, look for exactly one match in sources
        props_img = [None]*FIBRE_N

        for j in range(FIBRE_N):
                # for each segmentation label identified, remove it
            for props in props_list_return:
                # check distance to known fibre position
                # prop.centroid is in (y, x) format
                dist = np.linalg.norm(
                            np.fliplr([props.centroid.value])[0]
                            - FIBRE_COORDS[j])

                if dist < DIST_THRESHOLD:
                    # this is a fibre
                    # add props to dict
                    print('Fibre {} ID {} identified, distance to' \
                            'known centroid {} pixels'
                            .format(
                                j+1,
                                FIBRE_IDS[j],
                                dist
                                ))
                    props_img[j] = props
                    # remove the props of identified source from props list
                    props_list_return.remove(props)
                    # stop searching for current fibre
                    break

                else:
                    # this is not the fibre propwe are looking for
                    print('Cannot match {} to known fibre ID {},' \
                            ' distance to known centroid {} pixels'
                            .format(
                                np.fliplr([props.centroid.value])[0],
                                FIBRE_IDS[j],
                                dist
                        ))

        # verify in props_img each fibre has been identified correctly
        checklist = []
        for j in range(FIBRE_N):

            dist = np.linalg.norm(
                np.fliplr([props_img[j].centroid.value])[0]
                - FIBRE_COORDS[j])

            checklist.append(dist < DIST_THRESHOLD)

        if False in checklist:
            print('Fibre identification mismatch for dataset:', dirs_dataset[i])

        else:
            # ordering is good, centroid from props match known coordinates
            print('All fibres identified and well ordered for dataset:',
                      dirs_dataset[i])
            photometry_dict['props_img'][i] = list(props_img)
            photometry_dict['background'][i] = props_img[0].background_mean

    return photometry_dict
