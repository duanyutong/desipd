# -*- coding: utf-8 -*-
"""
Created on Sun Aug 28 00:17:32 2016

@author: Duan Yutong (dyt@physics.bu.edu, dyt@lbl.gov)

"""

import os
import glob
import pickle
import numpy as np
from scipy import ndimage

import matplotlib
matplotlib.use('Agg') # in case no backend was automatically chosen
import matplotlib.pyplot as plt
plt.ioff() # turn off interative mode to suppress plot output in console
from matplotlib.backends.backend_pdf import PdfPages

from astropy.io import fits, ascii
from astropy.convolution import Gaussian2DKernel
from astropy.stats import (sigma_clipped_stats, gaussian_fwhm_to_sigma,
                           gaussian_sigma_to_fwhm)
from astropy.table import vstack, Column
from astropy.visualization import LogStretch, mpl_normalize
from photutils import (detect_threshold, EllipticalAperture,
    source_properties, properties_table, aperture_photometry)
try:
    from photutils import detect_sources
except:
    from photutils.detection import detect_sources
try:
    from photutils.utils import calc_total_error, random_cmap
except:
    from photutils.utils import random_cmap

#%%
def segmentation_photometry(path_image_abs,
                            path_error_abs = None,
                            logger = None,
                            bkg_sigma = 1.5,
                            source_snr = 1.05,
                            fwhm_kernel = 25,
                            x_size_kernel = 100,
                            y_size_kernel = 80,
                            dump_pickle = False,
                            clobber = True):

    """

    given a fits file (master image), this function calculates
    photometry by source segmentation.

    make_source_mask not yet available in photutils v0.2.2, this version
    manually creates a source mask for determining background.

    """
    
    def msg(string, msgtype = None):
        
        if logger == None:
            print(string)
        else:
            print(string)
            if msgtype == 'info':
                logger.info(string)
            if msgtype == 'error':
                logger.error(string)
            if msgtype == 'warning':
                logger.warning(string)   

    filename = os.path.basename(path_image_abs)
    dir_save = os.path.dirname(path_image_abs)
    filenames_combined = '\t'.join(os.listdir(dir_save))

    if clobber == False \
        and filename[0:-5]+'-segm.obj' in filenames_combined  \
        and filename[0:-5]+'-props.obj' in filenames_combined \
        and filename[0:-5]+'-centroid_outline.png' in filenames_combined \
        and filename[0:-5]+'-centroid_outline.pdf' in filenames_combined \
        and filename[0:-5]+'-segmentation.png' in filenames_combined \
        and filename[0:-5]+'-segmentation.pdf' in filenames_combined:
        msg('Photometry properties table already exists. '
            + 'Reading pickles...', 
            msgtype='info')

        try:
            segm = pickle.load(open(glob.glob(os.path.join(
                                dir_save, filename[0:-5]+'-segm.obj*'))[0],
                                'rb'))
            props_list = pickle.load(open(glob.glob(os.path.join(
                                dir_save, filename[0:-5]+'-props.obj*'))[0],
                                'rb'))
            return [segm, props_list]
        except:
            # pickle file corrupt or empty, proceed
            pass
    elif clobber == False \
        and filename[0:-5]+'-logstretch.png' in filenames_combined \
        and filename[0:-5]+'-logstretch.pdf' in filenames_combined:
        msg('Non-detection from previous results.',
            msgtype='info')
        return [None, []]
        
    # image type notifications
    if 'master' in path_image_abs:
        if 'normalised' in path_image_abs:
            msg('Performing photometry to '
                + 'normalised master object image {}...'.format(path_image_abs),
                msgtype='info')
        else:
            msg('Performing photometry to '
                + 'un-normalised master image {}...'.format(path_image_abs),
                msgtype='info')
    elif 'reduced' in path_image_abs:
        msg('Performing photometry to '
                + 'reduced image frame {}...'.format(path_image_abs),
                msgtype='info')
    else:
        msg('Warning: Photometry being performed to '
            + 'a single exposure {}...'.format(path_image_abs),
            msgtype='warning')
    
    # read in data
    try:
        hdu = fits.open(path_image_abs)['FPC']
        data = hdu.data
    except:
        hdu = fits.open(path_image_abs)[0]
        data = hdu.data
    # read in error in data
    msg('Reading master error image {}...'
        .format(path_error_abs))
    try:
        hdu_error = fits.open(path_error_abs)[0]
        data_error = hdu_error.data
    except:
        data_error = np.zeros(data.shape)
        msg('No master error image available for {}'
            .format(path_image_abs))

    header = hdu.header

    if 'EXPREQ' in header:
        exptime = header['EXPREQ']
    elif 'EXPTIME' in header:
        exptime = header['EXPTIME']
    else:
        msg('Exposure time not found in header. '
            + 'Cannot determine magnitude.',
            msgtype='error')
        exptime = np.nan

    # === Iteratively determine background level ===

    # assuming backcground is homogenous, estimate background by sigma clipping
    # if background noise varies across image, generate 2D background instead
    # using the Background function
    msg('Determining background noise level...',
        msgtype='info')
    [mean, median, std] = sigma_clipped_stats(data, sigma=bkg_sigma, iters=3)
    threshold = median + (std * 4)
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
    msg('Determining threshold for source detection...',
        msgtype='info')
    # determine threshold for source detection
    # in current implementation, if all inputs are present, the formula is
    # threshold = background + (background_error * snr)
    threshold = detect_threshold(data,
                                 background = bkg_median,
                                 error = data_error+bkg_std,
                                 snr = source_snr)
    # calculate total error including poisson statistics
    try:
        # this is for v0.3 and above
        msg('Calculating total errors including background and Poisson...',
            msgtype='info')
        err_tot = calc_total_error(data, 
                                   bkg_error= data_error+bkg_std, 
                                   effective_gain=0.37)
        gain = None
    # in version earlier than 0.3, this function is not available
    except:
        # error must be of the same shape as the data array
        # this is for v0.2.2
        err_tot = data_error + bkg_std
        gain = 0.37
    msg('Preparing 2D Gaussian kernal...',
        msgtype='info')
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
    msg('Segmentation processing...',
        msgtype='info')
    segm = detect_sources(data, threshold, npixels=5, filter_kernel=kernel)
    msg('Segmentation labels are: ' + repr(segm.labels),
        msgtype='info')

    # === Measure regional source properties ===

    # source_properties() assumes that the data have been background-subtracted.
    # Background is the background level that was previously present
    # in the input data.
    # The input background does not get subtracted from the input data,
    # which should already be background-subtracted.
    msg('Extracting source properties...',
        msgtype='info')
    if gain is None:
        # gain is no longer supported in v0.3 and included in total error array
        props_list = source_properties(data-bkg_median, segm, 
                               background = bkg_median,
                               error = err_tot)
    else: # still in v0.2.2
        props_list = source_properties(data-bkg_median, segm, 
                                       background = bkg_median,
                                       error = err_tot,
                                       effective_gain = gain)
    # add more properties that are not automatically calculated
    for i in range(len(props_list)):
        # source_sum is by definition background-subtracted already
        props_list[i].flux = props_list[i].source_sum/exptime
        props_list[i].flux_err = props_list[i].source_sum_err/exptime
        # flux = source_sum / exptime
        # instrumental magnitude = -2.5 * log10(flux)
        props_list[i].mag_instr = -2.5 * np.log10(props_list[i].flux)
        props_list[i].mag_instr_err = -2.5 / props_list[i].flux / np.log(10) \
                                        * props_list[i].flux_err
        # assuming fwhm of a circule gaussian of the same cross section area
        props_list[i].fwhm = gaussian_sigma_to_fwhm * np.sqrt(
                            props_list[i].semimajor_axis_sigma.value
                            * props_list[i].semiminor_axis_sigma.value)
    # make plots and save to images
    # define approximate isophotal ellipses for each object
    apertures = []
    r = 5 # approximate isophotal extent
    for props in props_list:
        position = (props.xcentroid.value, props.ycentroid.value)
        a = props.semimajor_axis_sigma.value * r
        b = props.semiminor_axis_sigma.value * r
        theta = props.orientation.value
        apertures.append(EllipticalAperture(position, a, b, theta=theta))
        
    # === plot and save ===
        
    # if filename ends with fits, remove it in the  filename
    if filename[-5:] == '.fits':
        path_save_prefix = os.path.join(dir_save, filename[0:-5])
    else:
        path_save_prefix = os.path.join(dir_save, filename)
        
    norm_log = mpl_normalize.ImageNormalize(vmin=0, vmax=2000,
                                            stretch = LogStretch())

    if len(props_list) > 0:
        
        # Save segm, porps to object files, and also save props to table file.
        msg('Saving segmentation and source properties to {}...'
            .format(dir_save),
            msgtype='info')
        
        # at least one source was detected
        # create a table of properties
        props_table = properties_table(props_list)
        # add custom columns to the table: mag_instru and flux
        props_table['flux'] = [props_list[i].flux
                                for i in range(len(props_list))]
        props_table['flux_err'] = [props_list[i].flux_err
                                for i in range(len(props_list))]
        props_table['mag_instr'] = [props_list[i].mag_instr
                                        for i in range(len(props_list))]
        props_table['mag_instr_err'] = [props_list[i].mag_instr_err
                                        for i in range(len(props_list))]
        props_table['fwhm'] = [props_list[i].fwhm
                                for i in range(len(props_list))]

        # plot centroid and segmentation outline
        [fig1, ax1] = plt.subplots(figsize=(4, 3))
        ax1.imshow(data, origin='lower', cmap=plt.cm.gray, norm=norm_log)
        ax1.plot(props_table['xcentroid'], props_table['ycentroid'],
                 linestyle='none', color='red',
                 marker='+', markersize=2, markeredgewidth=0.1, alpha=1)
        segm_outline = np.array(segm.outline_segments(), dtype=float)
        segm_outline[segm_outline<1] = np.nan
        # get a copy of the gray color map
        segm_outline_cmap = plt.cm.winter
        # set how the colormap handles 'bad' values
        segm_outline_cmap.set_bad(alpha=0)
        ax1.imshow(segm_outline, 
                   origin='lower', cmap=segm_outline_cmap, alpha=1)
        ax1.get_xaxis().set_visible(False)
        ax1.get_yaxis().set_visible(False)
        fig1.tight_layout()
    
        # segmentation image and aperture using approximate elliptical apertures
        [fig2, ax2] = plt.subplots(figsize=(4, 3))
        rand_cmap = random_cmap(segm.max + 1, random_state=8)
        ax2.imshow(segm, origin='lower', cmap=rand_cmap)
        ax2.plot(props_table['xcentroid'], props_table['ycentroid'],
                 linestyle='none', color='red',
                 marker='+', markersize=2, markeredgewidth=0.1, alpha=1)
        for aperture in apertures:
            aperture.plot(ax=ax2, lw=0.1, alpha=1, color='lime')
        ax2.axis('off')
        ax2.get_xaxis().set_visible(False)
        ax2.get_yaxis().set_visible(False)
        fig2.tight_layout()
        
        try:
            # Enhanced CSV allows preserving table meta-data such as
            # column data types and units.
            # In this way a data table can be stored and read back as ASCII
            # with no loss of information.
            ascii.write(props_table, path_save_prefix + '-props.ecsv',
                        format = 'ecsv')
            # csv for readability in MS excel
            ascii.write(props_table, path_save_prefix + '-props.csv',
                        format = 'csv')
    
            # save figures
            fig1.savefig(path_save_prefix + '-centroid_outline.png',
                         bbox_inches='tight', pad_inches=0, dpi=1200)
            fig2.savefig(path_save_prefix + '-segmentation.png',
                         bbox_inches='tight', pad_inches=0, dpi=2000)
                         
            pp1 = PdfPages(path_save_prefix + '-centroid_outline.pdf')
            pp1.savefig(fig1, dpi=1200)
            pp1.close()
    
            pp2 = PdfPages(path_save_prefix + '-segmentation.pdf')
            pp2.savefig(fig2, dpi=2000)
            pp2.close()

            if dump_pickle:
                # dump segmentation and properties to objects in binary mode
                file_segm = open(path_save_prefix + '-segm.obj', 'wb')
                pickle.dump(segm, file_segm)
                file_props = open(path_save_prefix + '-props.obj', 'wb')
                pickle.dump(props_list, file_props)
    
            msg('Segmentation, properties objects, tables, and images '
                + 'saved to {}'.format(dir_save),
                msgtype='info')
                
        except:
            msg('Unable to write to disk, check permissions.',
                msgtype='error')
            
        # memory leak?
        try:
            plt.close('all')
            del (hdu, hdu_error, data, data_error, header, mask, mask_dilated, 
                 err_tot, kernel, apertures, norm_log, props_table, 
                 segm_outline, segm_outline_cmap, rand_cmap, 
                 fig1, ax1, fig2, ax2, pp1, pp2, file_segm, file_props)
        except:
            pass
        return [segm, props_list]

    else:
        msg('No source detected in {}'.format(path_image_abs),
            msgtype='warning')
            
        # save log scale stretched image, if no source was detected
        [fig0, ax0] = plt.subplots(figsize=(4, 3))
        ax0.imshow(data, origin='lower', cmap=plt.cm.gray, norm=norm_log)
        ax0.get_xaxis().set_visible(False)
        ax0.get_yaxis().set_visible(False)
        
        try:
            fig0.savefig(path_save_prefix + '-logstretch.png',
                         bbox_inches='tight', pad_inches=0, dpi=1200)
            pp0 = PdfPages(path_save_prefix + '-logstretch.pdf')
            pp0.savefig(fig0, dpi=1200)
            pp0.close()
            
        except:
            msg('Unable to write to disk, check permissions.',
                msgtype='error')
            
        return [None, []]

#%%

def aper_photometry(path_image_abs, 
                        apertures, 
                        logger = None,
                        bkg_sigma = 1.5,
                        dump_pickle = True,
                        clobber = False):
    
    def msg(string, msgtype = None):
        
        if logger == None:
            print(string)
        else:
            print(string)
            if msgtype == 'info':
                logger.info(string)
            if msgtype == 'error':
                logger.error(string)
            if msgtype == 'warning':
                logger.warning(string)

    filename = os.path.basename(path_image_abs)
    dir_save = os.path.dirname(path_image_abs)
    filenames_combined = '\t'.join(os.listdir(dir_save))

    if clobber == False \
        and filename[0:-5]+'-aper_phot.png' in filenames_combined \
        and filename[0:-5]+'-aper_phot.pdf' in filenames_combined \
        and filename[0:-5]+'-aper_phot.csv' in filenames_combined \
        and filename[0:-5]+'-aper_phot.ecsv' in filenames_combined \
        and filename[0:-5]+'-aper_phot_table.obj' in filenames_combined:
        msg('Aperture photometry table already exists. '
            + 'Reading pickle...', 
            msgtype='info')
        try:
            phot_table = pickle.load(open(glob.glob(os.path.join(
                                    dir_save, 
                                    filename[0:-5]+'-aper_phot_table.obj'))[0],
                                    'rb'))
            return phot_table
        except:
            # pickle file corrupt or empty, proceed
            pass
                
    # image type notifications
    if 'master' in path_image_abs:
        if 'normalised' in path_image_abs:
            msg('Performing aperture photometry to'
                + 'normalised master object image {}...'.format(path_image_abs),
                msgtype='info')
        else:
            msg('Performing aperture photometry to'
                + 'un-normalised master image {}...'.format(path_image_abs),
                msgtype='info')
    elif 'reduced' in path_image_abs:
        msg('Performing aperture photometry to '
                + 'reduced image frame {}...'.format(path_image_abs),
                msgtype='info')
    else:
        msg('Warning: Aperture photometry being performed to '
            + 'a single exposure {}...'.format(path_image_abs),
            msgtype='warning')
    
    # read in data and header, including exptime
    try:
        hdu = fits.open(path_image_abs)['FPC']
    except:
        hdu = fits.open(path_image_abs)[0]
    data = hdu.data
    header = hdu.header

    if 'EXPREQ' in header:
        exptime = header['EXPREQ']
    elif 'EXPTIME' in header:
        exptime = header['EXPTIME']
    else:
        msg('Exposure time not found in header. '
            + 'Cannot determine magnitude.',
            msgtype='error')
        exptime = np.nan
        
    # iteratively determine background level
    # assuming backcground is homogenous, estimate background by sigma clipping
    # if background noise varies across image, generate 2D background instead
    msg('Determining background noise level...',
        msgtype='info')
    [mean, median, std] = sigma_clipped_stats(data, sigma=bkg_sigma, iters=3)
    threshold = median + (std * 4)
    segm = detect_sources(data, threshold, npixels=5)
    # turn segm into a mask
    mask = segm.data.astype(np.bool)
    # dilate the source mask to ensure complete masking of detected sources
    dilate_structure = np.ones((5, 5))
    mask_dilated = ndimage.binary_dilation(mask, structure=dilate_structure)
    # get sigma clipping stats of background, without sources that are masekd
    [bkg_mean, bkg_median, bkg_std] = sigma_clipped_stats(
        data, sigma=bkg_sigma, mask=mask_dilated, iters = 3)
    # perform aperture photometry and create a table
    msg('Extracting aperture photometry data...',
        msgtype='info')
    phot_tables = []
    for aperture in apertures:
        # aperture_photometry only supports one unique aperture at a time
        phot_table = aperture_photometry(data-bkg_median, aperture)
        phot_tables.append(phot_table)
    # vertically stack tables into one, given that they have the same columns
    phot_table = vstack(phot_tables)
    # clean up table - some entries may be a 1x1 array instead of a number
    # have to change datatype of column
    for column in ['xcenter', 'ycenter']:
        if phot_table[0][column].shape is not ():
            msg('Replacing arrays in column {} by float64...'
                .format(column))
            # get column data in a list of float64 numbers
            column_data = [phot_table[row][column][0] 
                           for row in range(len(phot_table))]
            phot_table.replace_column(column, column_data)

    # add columns
    flux_data = [phot_table[row]['aperture_sum']/exptime 
                                 for row in range(len(phot_table))]
    column_flux = Column(name='flux', 
                         data = flux_data)
    column_mag = Column(name='mag_instr',
                        data = [-2.5 * np.log10(flux) for flux in flux_data])
    column_bkg_mean = Column(name='background_mean', 
                             data = [bkg_mean]*len(phot_table))
    column_bkg_median = Column(name='background_median',
                               data = [bkg_median]*len(phot_table))
    phot_table.add_columns([column_flux, column_mag, 
                            column_bkg_mean, column_bkg_median])

    # plot and save
    msg('Saving aperture photometry results to {}...'
        .format(dir_save),
        msgtype='info')
    # if filename ends with fits, remove it in the  filename
    if filename[-5:] == '.fits':
        path_save_prefix = os.path.join(dir_save, filename[0:-5])
    else:
        path_save_prefix = os.path.join(dir_save, filename)
        
    norm_log = mpl_normalize.ImageNormalize(vmin=0, vmax=2000,
                                            stretch = LogStretch())
                                            
    [fig, ax] = plt.subplots(figsize=(4,3))
    ax.imshow(data, origin='lower', cmap=plt.cm.gray, norm=norm_log)
    for aperture in apertures:
        aperture.plot(ax=ax, lw=0.1, alpha=1, color='lime')
    ax.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    
    # save all files
    try:
        ascii.write(phot_table, path_save_prefix + '-aper_phot.ecsv',
                    format = 'ecsv')
        # csv for readability in MS excel
        ascii.write(phot_table, path_save_prefix + '-aper_phot.csv',
                    format = 'csv')
        # save figures
        fig.savefig(path_save_prefix + '-aper_phot.png',
                     bbox_inches='tight', pad_inches=0, dpi=1200)
                     
        pp = PdfPages(path_save_prefix + '-aper_phot.pdf')
        pp.savefig(fig, dpi=1200)
        pp.close()
        
        if dump_pickle:
            file_phot = open(path_save_prefix + '-aper_phot_table.obj', 'wb')
            pickle.dump(phot_table, file_phot)
            
        msg('Aperture photometry object, tables, and images '
            + 'saved to {}'.format(dir_save),
            msgtype='info')

    except:
        msg('Unable to write to disk, check permissions.',
            msgtype='error')
    
    # memory leak?
    try:
        plt.close('all')
        del (hdu, data, header, mask, mask_dilated, phot_tables, flux_data, 
         column_flux, column_mag, column_bkg_mean, column_bkg_median, 
         norm_log, fig, ax, apertures, pp, file_phot)
    except:
        pass
    
    return phot_table
