# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 18:25:16 2016

@author: givoltage
"""
import os
import glob
import logging
from astropy.io import fits

logger = logging.getLogger()

#%% worker function that gets wrapped

def reduce_dataset_worker(path_dataset,
                          save_dir_rel = '',
                          clobber = False,
                          dataset_is_flat = False,
                          suppress_return = True):

    logger.info('Reducing directory {}...'.format(path_dataset))

    # if any master file already exists, and clobber = False
    # return existing master file
    filenames_combined = '\t'.join(os.listdir(path_dataset))
    if 'master' in filenames_combined.lower() and not clobber:
        logger.info('Master files already exist. Overwriting not allowed.')
        # return data from existing master file
        if dataset_is_flat and not suppress_return:
            logger.info('Readeing {}'.format(
                os.path.join(
                    path_dataset, 'master_flat_normalised.fits')))
            return fits.open(
                    os.path.join(
                        path_dataset, 'master_flat_normalised.fits')
                        ) [0].data
        elif not dataset_is_flat and not suppress_return:
            logger.info('Reading {}'.format(
                os.path.join(
                    path_dataset, 'master_object.fits')))
            return fits.open(
                    os.path.join(
                        path_dataset, 'master_object.fits')
                        ) [0].data
        else:
            logger.info('No reduction done for {}'.format(path_dataset))
            return

    # check and create master folder
    dir_save = os.path.join(path_dataset, save_dir_rel)
    if not os.path.exists(dir_save):
        os.makedirs(dir_save)

    # read in fits file
    path_pattern = os.path.join(path_dataset, '*.fit*')
    paths = glob.glob(path_pattern)
    # exclude existing master files
    paths = [path for path in paths if 'master' not in path]
    # read obstype from all headers
    obstypes = []
    imagetyps = []
    for path in paths:
        header = fits.open(path)[0].header
        if 'OBSTYPE' in header:
            obstypes.append(header['OBSTYPE'])
        if 'IMAGETYP' in header:
            imagetyps.append(header['IMAGETYP'])
    frames_dict = {'path': paths,
                   'obstype': obstypes,
                   'imagetyp': imagetyps}
    # perform search based on obstype
    try:
        paths_bias = [frames_dict['path'][i]
            for i in range(len(paths))
            if frames_dict['obstype'][i].lower() in [
                'bias', 'zero'
                ]]
        paths_dark = [frames_dict['path'][i]
            for i in range(len(paths))
            if frames_dict['obstype'][i].lower() in [
                'dark'
                ]]
    # if no yield, search based on imagetyp
    except:
        paths_bias = [frames_dict['path'][i]
            for i in range(len(paths))
            if 'bias' in frames_dict['imagetyp'][i].lower()]
        paths_dark = [frames_dict['path'][i]
            for i in range(len(paths))
            if 'dark' in frames_dict['imagetyp'][i].lower()]

    if dataset_is_flat:
        try:
            paths_object = [frames_dict['path'][i]
                for i in range(len(paths))
                if frames_dict['obstype'][i].lower() in [
                   'flat' or 'dome flat' or 'sky flat'
                    ]]
        except:
            paths_object = [frames_dict['path'][i]
                for i in range(len(paths))
                if 'flat' in frames_dict['imagetyp'][i].lower()]
    else: # dataset is indeed object
        try:
            paths_object = [frames_dict['path'][i]
                for i in range(len(paths))
                if frames_dict['obstype'][i].lower() in [
                    'object' or 'light'
                    ]]
        except:
            paths_object = [frames_dict['path'][i]
                for i in range(len(paths))
                if 'light' in frames_dict['imagetyp'][i].lower()]

    # get info in order to preallocate arrays
    hdulist = fits.open(paths_bias[0])
    shape = hdulist[0].data.shape 	# get image dimension
    # dtype = hdulist[0].data.dtype 	# get data type # always float16
    datatype = np.float32   # let's just have the utmost precision
    # get numbers of files for each type
    num_bias = len(paths_bias)
    num_dark = len(paths_dark)
    num_object = len(paths_object)
    logger.info('Read {} bias frames'.format(num_bias))
    logger.info('Read {} dark frames'.format(num_dark))
    logger.info('Read {} object/flat frames'.format(num_object))

    # preallocate arrays
    bias3D = np.empty([shape[0], shape[1], num_bias], dtype = datatype)
    dark3D = np.empty([shape[0], shape[1], num_dark], dtype = datatype)
    object3D = np.empty([shape[0], shape[1], num_object], dtype = datatype)

    # fill arrays with fits data
    for i in range(num_bias):
        bias3D[:, :, i] = fits.open(paths_bias[i]) [0].data
    for i in range(num_dark):
        dark3D[:, :, i] = fits.open(paths_dark[i]) [0].data

    # master bias from median combining
    m_bias = np.median(bias3D, axis=2)
    # master dark from median combining, -m_bias
    m_dark = np.median(dark3D, axis=2) - m_bias
    # create hdu objects
    hdu_m_bias = fits.PrimaryHDU(m_bias)
    hdu_m_dark = fits.PrimaryHDU(m_dark)

    # # write exposure time in header for flux/magnitude calculation
    if 'EXPREQ' in fits.open(paths_dark[0]) [0].header:
        expreq_bias = fits.open(paths_bias[0]) [0].header['EXPREQ']
        expreq_dark = fits.open(paths_dark[0]) [0].header['EXPREQ']
        hdu_m_bias.header ['EXPREQ'] = expreq_bias
        hdu_m_dark.header ['EXPREQ'] = expreq_dark
    # if exptime is present, add it, too
    if 'EXPTIME' in fits.open(paths_dark[0]) [0].header:
        exptime_bias = fits.open(paths_bias[0]) [0].header['EXPTIME']
        exptime = fits.open(paths_dark[0]) [0].header['EXPTIME']
        hdu_m_bias.header ['EXPTIME'] = exptime_bias
        hdu_m_dark.header ['EXPTIME'] = exptime

    # write hdu to fits
    hdu_m_bias.writeto(
        os.path.join(dir_save, 'master_bias.fits'), clobber = clobber)
    logger.info('Master bias saved to {}'.format(
        os.path.join(dir_save, 'master_bias.fits')))
    hdu_m_dark.writeto(
        os.path.join(dir_save, 'master_dark.fits'), clobber = clobber)
    logger.info('Master dark saved to {}'.format(
        os.path.join(dir_save, 'master_dark.fits')))

    # master light which can be object or flat
    if dataset_is_flat:

        paths_flat = paths_object
        num_flat = len(paths_flat)
        flat3D = np.empty([shape[0], shape[1], num_flat], dtype = datatype)

        for i in range(len(paths_flat)):
            # subtract bias and dark
            flat3D[:,:,i] = fits.open(paths_flat[i]) [0].data \
                - m_bias - m_dark
            [mode, _ ] = stats.mode(flat3D[:, :, i], axis=None)
            # each flat is divide by mode
            flat3D[:,:,i] = flat3D[:,:,i]/mode
        # median combine all flat frames
        m_flat = np.median(flat3D, axis=2)
        # normalised by the mean of master flat
        m_flat_normalised = m_flat / np.mean(m_flat)
        # create hdu objects
        hdu_m_flat = fits.PrimaryHDU(m_flat)
        hdu_m_flat_normalised = fits.PrimaryHDU(m_flat_normalised)

        # write exposure time in header for flux/magnitude calculation
        if 'EXPREQ' in fits.open(paths_flat[0]) [0].header:
            expreq_flat = fits.open(paths_flat[0]) [0].header['EXPREQ']
            hdu_m_flat.header ['EXPREQ'] = expreq_flat
            hdu_m_flat_normalised.header ['EXPREQ'] = expreq_flat
        # if exptime is present, add it, too
        if 'EXPTIME' in fits.open(paths_flat[0]) [0].header:
            exptime = fits.open(paths_object[0]) [0].header['EXPTIME']
            hdu_m_flat.header ['EXPTIME'] = exptime
            hdu_m_flat_normalised.header ['EXPTIME'] = exptime

        # write hdu to fits
        hdu_m_flat.writeto(
            os.path.join(dir_save, 'master_flat.fits'),
            clobber = clobber)
        logger.info('Master flat saved to {}'.format(
            os.path.join(dir_save, 'master_flat.fits')))

        hdu_m_flat_normalised.writeto(
            os.path.join(dir_save, 'master_flat_normalised.fits'),
            clobber = clobber)
        logger.info('Master flat normalised saved to {}'.format(
            os.path.join(dir_save, 'master_flat_normalised.fits')))

        if not suppress_return:
            return m_flat_normalised

    else: # dataset is object, not flat

        for i in range(len(paths_object)):
            object3D[:, :, i] = fits.open(paths_object[i]) [0].data
        # master object from median combining, -m_bias, -m_dark
        m_object = np.median(object3D, axis=2) - m_bias - m_dark
        # create hdu objects
        hdu_m_object = fits.PrimaryHDU(m_object)

        # write exposure time in header for flux/magnitude calculation
        if 'EXPREQ' in fits.open(paths_object[0]) [0].header:
            expreq_object = fits.open(paths_object[0]) [0].header['EXPREQ']
            hdu_m_object.header ['EXPREQ'] = expreq_object
        # if exptime is present, add it, too
        if 'EXPTIME' in fits.open(paths_object[0]) [0].header:
            exptime = fits.open(paths_object[0]) [0].header['EXPTIME']
            hdu_m_object.header ['EXPTIME'] = exptime
        # write hdu to fits
        hdu_m_object.writeto(os.path.join(dir_save, 'master_object.fits'),
                             clobber = clobber)

        logger.info('Master object saved to {}'.format(
            os.path.join(dir_save, 'master_object.fits')))
        if not suppress_return:
            return m_object

#%%

# path_dataset = '/data/images/fpc/linear/'
path_flat = r'C:\Users\givoltage\Downloads\20160414\flat'
# path_flat = '/data/images/fpc/flats/'
save_dir_rel = ''
clobber = True

m_flat_normalised = reduce_dataset_worker(path_flat,
                                          save_dir_rel = save_dir_rel,
                                          clobber = clobber,
                                          dataset_is_flat = True,
                                          suppress_return = False)