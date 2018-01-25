# -*- coding: utf-8 -*-
"""
Created on Sun Aug 28 20:47:11 2016

@author: Duan Yutong (dyt@physics.bu.edu, dyt@lbl.gov)

"""

#%% function

def segmentation_photometry(path_file_abs,
                            bkg_sigma = 3.0,
                            source_snr = 3.0,
                            fwhm_kernel = 2.0,
                            x_size_kernel = 3,
                            y_size_kernel = 3,
                            save_fig_pdf = True,
                            clobber = False):

    """

    given a fits file (master image), this function calculates
    aperture photometry by source segmentation.

    make_source_mask not yet available in photutils v0.2.2, this version
    manually creates a source mask for determining background.

    """

    import os
    import copy
    import glob
    import pickle
    import numpy as np
    from scipy import ndimage
    import matplotlib
    import matplotlib.pyplot as plt
    from astropy.io import fits, ascii
    from astropy.convolution import Gaussian2DKernel
    from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
#    from astropy.table import Table
    from astropy.visualization import (LogStretch, mpl_normalize)
#    from astropy.extern.six.moves import StringIO
    from photutils import (detect_threshold, EllipticalAperture,
        source_properties, properties_table)
    from photutils.detection import detect_sources
    from photutils.utils import random_cmap

    # create preliminary mask
    #from photutils import make_source_mask
    #masterMask = make_source_mask(master, snr=2, npixels=5, dilate_size=11)

#        if LEDoff was used, get threshold from LEDoff/background
#        path_dataset = os.path.dirname(path_file_abs) + os.path.sep
#        filenameCombined = '\t'.join(
#            os.listdir(os.path.join(datasetDirLocal, 'master')))
#        if 'master_ledoff_subtracted' in filename:
#            print('Using master_ledoff')
#            # path_file_abs = os.path.join(datasetDir, 'master', filename)
#            hdu = fits.open(path_file_abs)[0]
#            data_subtracted = hdu.data
#            # calculate threadhold
#            ledoff_pred = np.mean(data_subtracted) * \
#                np.ones(data_subtracted.shape)
#            mse = mean_squared_error(data_subtracted, ledoff_pred)
#            rmse = np.sqrt(mse)
#            threshold = 7.0 * rmse
#            threshold_value = threshold
#         if no LEDoff was used, background subtraction is needed
#         there should exist no file named "subtracted"
#         if 'master.fit' in filenameCombined \
#             or 'master_normalised.fit' in filenameCombined:

    filename = os.path.basename(path_file_abs)
    dir_save = os.path.dirname(path_file_abs)
    filenames_combined = '\t'.join(os.listdir(os.path.dirname(path_file_abs)))

    if clobber == False \
        and 'segm.obj' in filenames_combined  \
        and 'props.obj' in filenames_combined \
        and 'props.csv' in filenames_combined\
        and 'props.ecsv' in filenames_combined:
        print('Photometry properties table already exists. Reading objects...')
        segm = pickle.load(open(glob.glob(os.path.join(
                                            dir_save, '*segm.obj*'))[0],
                           'rb'))
        props = pickle.load(open(glob.glob(os.path.join(
                                            dir_save, '*props.obj*'))[0],
                           'rb'))

        return [segm, props]

    if 'master' in path_file_abs:
        if 'normalised' in path_file_abs:
            print('Performing photometry to '
                + 'normalised master object image {}...'.format(path_file_abs))
        else:
            print('Performing photometry to '
                + 'un-normalised master image {}...'.format(path_file_abs))
    else:
        print('Warning: Photometry being performed to '
            + 'a single exposure {}...'.format(path_file_abs))

    hdu = fits.open(path_file_abs)[0]
    data = hdu.data
    header = hdu.header

    if 'EXPREQ' in header:
        exptime = header['EXPREQ']
    elif 'EXPTIME' in header:
        exptime = header['EXPTIME']
    else:
        print('Exposure time not found in header. Cannot determine magnitude.')
        exptime = np.nan

    # === Iteratively determine background level ===

    # assuming background is homogenous, estimate background by sigma clipping
    # if background noise varies across image, generate 2D background instead
    print('Determining background noise level...''')
    [mean, median, std] = sigma_clipped_stats(data, sigma=bkg_sigma, iters=5)
    threshold = median + (std * 2.0)
    segm = detect_sources(data, threshold, npixels=5)
    # turn segm into a mask
    mask = segm.data.astype(np.bool)
    # dilate the source mask to ensure complete masking of detected sources
    dilate_structure = np.ones((5, 5))
    mask_dilated = ndimage.binary_dilation(mask, structure=dilate_structure)
    # get sigma clipping stats of background, without sources that are masekd
    [bkg_mean, bkg_median, bkg_std] = sigma_clipped_stats(
        data, sigma=bkg_sigma, mask=mask_dilated, iters = 3)

    # === Detect sources by segmentation ===

    print('Determining threshold for source detection...')
    # determine threshold for source detection
    # in current implementation, if all inputs are present, the formula is
    # threshold = background + (background_error * snr)
    threshold = detect_threshold(data,
                                 background = bkg_median,
                                 error = bkg_std,
                                 snr = source_snr)
    print('Preparing 2D Gaussian kernal...')
    sigma_kernel = fwhm_kernel * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma_kernel,
                           x_size = x_size_kernel,
                           y_size = y_size_kernel)
    # normalise kernel
    # The kernel models are normalized per default, ∫∞−∞f(x)dx=1∫−∞∞f(x)dx=1.
    # But because of the limited kernel array size, the normalization
    # for kernels with an infinite response can differ from one.
    kernel.normalize()
    # obtain a  SegmentationImage object with the same shape as the data,
    # where sources are labeled by different positive integer values.
    # A value of zero is always reserved for the background.
    # if the threshold includes the background level as above, then the image
    # input into detect_sources() should not be background subtracted.
    print('Segmentation processing...')
    segm = detect_sources(data, threshold, npixels=5, filter_kernel=kernel)
    print('Segmentation labels are: ', repr(segm.labels))

    # === Measure regional source properties ===

    # source_properties() assumes that the data have been background-subtracted.
    # Background is the background level that was previously present
    # in the input data.
    # The input background does not get subtracted from the input data,
    # which should already be background-subtracted.
    print('Extracting source properties...')
    props = source_properties(data-bkg_median, segm, background=bkg_median)
    # add flux and instrumental magnitude to properties
    # flux = source_sum / exptime
    # instrumental magnitude = -2.5 * log10(flux)
    for i in range(len(props)):
        # source_sum is by definition background-subtracted already
        props[i].flux = props[i].source_sum/exptime
        props[i].mag_instr = -2.5 * np.log10(props[i].flux)
    # make plots and save to images
    # define approximate isophotal ellipses for each object
    apertures = []
    r = 2.8 # approximate isophotal extent
    for prop in props:
        position = (prop.xcentroid.value, prop.ycentroid.value)
        a = prop.semimajor_axis_sigma.value * r
        b = prop.semiminor_axis_sigma.value * r
        theta = prop.orientation.value
        apertures.append(EllipticalAperture(position, a, b, theta=theta))

    # create a table of properties
    try:
        props_table = properties_table(props)
    except:
        print('No source detected in {}'.format(path_file_abs))
        return [None, None]

    props_table['flux'] = [props[i].flux for i in range(len(props))]
    props_table['mag_instr'] = [props[i].mag_instr for i in range(len(props))]
    # add custom columns to the table: mag_instru and flux

    # plot centroid and segmentation using approximate elliptical apertures
    norm = mpl_normalize.ImageNormalize(stretch = LogStretch())
    rand_cmap = random_cmap(segm.max + 1, random_state=12345)
    [fig1, (ax1, ax2)] = plt.subplots(1, 2, figsize = (12, 6))
    ax1.imshow(data, origin='lower', cmap=plt.cm.gray, norm=norm)
    ax1.plot(
            props_table['xcentroid'], props_table['ycentroid'],
            ls='none', color='blue', marker='+', ms=10, lw=1.5)
    ax2.imshow(segm, origin='lower', cmap=rand_cmap)
    for aperture in apertures:
        aperture.plot(ax=ax1, lw=1.0, alpha=1.0, color='red')
        aperture.plot(ax=ax2, lw=1.0, alpha=1.0, color='red')
    # plot using actual segmentation outlines (to be improved)
    [fig2, ax3] = plt.subplots(figsize = (6, 6))
    ax3.imshow(data, origin='lower', cmap=plt.cm.gray, norm=norm)
    segm_outline = np.array(segm.outline_segments(), dtype=float)
    segm_outline[segm_outline<1] = np.nan
    # get a copy of the gray color map
    segm_outline_cmap = copy.copy(plt.cm.get_cmap('autumn'))
    # set how the colormap handles 'bad' values
    segm_outline_cmap.set_bad(alpha=0)
    ax3.imshow(segm_outline, origin='lower', cmap=segm_outline_cmap)

    # === save ===
    # Save segm, porps to object files, and also save props to table file.
    print('Saving segmentation and source properties to {}...'.format(dir_save))
    try:
        # if filename ends with fits, remove it in the  filename
        if filename[-5:] == '.fits':
            dir_save_prefix = os.path.join(dir_save, filename[0:-5])
        else:
            dir_save_prefix = os.path.join(dir_save, filename)
        # Enhanced CSV allows preserving table meta-data such as
        # column data types and units.
        # In this way a data table can be stored and read back as ASCII
        # with no loss of information.
        ascii.write(props_table, dir_save_prefix + '-phot_props.ecsv',
                    format = 'ecsv')
        # csv for readability in MS excel
        ascii.write(props_table,dir_save_prefix + '-phot_props.csv',
                    format = 'csv')

        # dump segmentation and properties to object files in binary mode
        file_segm = open(dir_save_prefix + '-phot_segm.obj', 'wb')
        pickle.dump(segm, file_segm)
        file_props = open(dir_save_prefix + '-phot_props.obj', 'wb')
        pickle.dump(props, file_props)

        # save figures
        fig1.savefig(dir_save_prefix + '-phot_segm_fig1.png', dpi=600)


        fig2.savefig(dir_save_prefix + '-phot_segm_fig2.png', dpi=600)

        if save_fig_pdf:

            matplotlib.rcParams['text.usetex'] = True
            matplotlib.rcParams['text.latex.unicode'] = True
            from matplotlib.backends.backend_pdf import PdfPages

            pp1 = PdfPages(dir_save_prefix + '-phot_segm_fig1.pdf')
            pp1.savefig(fig1)
            pp1.close()

            pp2 = PdfPages(dir_save_prefix + '-phot_segm_fig2.pdf')
            pp2.savefig(fig2)
            pp2.close()

        print('Segmentation, properties objects, tables, and images saved to',
                  dir_save)
    except:
        print('Unable to write to disk, check permissions.')

    return [segm, props]
