# -*- coding: utf-8 -*-
"""
Created on Sun Aug 28 00:17:32 2016

@author: Duan Yutong (dyt@physics.bu.edu, dyt@lbl.gov)

"""

#%% worker function that gets wrapped

def reduce_dataset_worker(logger,
                          dir_dataset,
                          dir_save = '/data/images/fpc_analysis/',
                          clobber = False,
                          dataset_is_flat = False,
                          suppress_return = True):

    # import modules
    import os
    import glob
    import numpy as np
    from scipy import stats
    from astropy.io import fits

    # check and create master folder if it doesn't exist
    if not os.path.exists(dir_save):
        os.makedirs(dir_save)
    else:
        # if any master file already exists, and clobber = False
        # return existing master file
        filenames_combined = '\t'.join(os.listdir(dir_save))
        if 'master' in filenames_combined.lower() and not clobber:
            logger.info('Master files already exist. Overwriting not allowed.')
            # return data from existing master file
            if dataset_is_flat and not suppress_return:
                logger.info('Readeing {}'.format(
                    os.path.join(
                        dir_save, 'master_flat_normalised.fits')))
                return fits.open(
                        os.path.join(
                            dir_save, 'master_flat_normalised.fits')
                            ) [0].data
            elif not dataset_is_flat and not suppress_return:
                logger.info('Reading {}'.format(
                    os.path.join(
                        dir_save, 'master_object.fits')))
                return fits.open(
                        os.path.join(
                            dir_save, 'master_object.fits')
                            ) [0].data
            else:
                logger.info('No reduction done for {}'.format(dir_dataset))
                return

    # read in fits file
    path_pattern = os.path.join(dir_dataset, '*.fit*')
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
    logger.info('Read {} bias frames'.format(num_bias))
    logger.info('Read {} dark frames'.format(num_dark))

    # preallocate arrays
    bias3D = np.empty([shape[0], shape[1], num_bias], dtype = datatype)
    dark3D = np.empty([shape[0], shape[1], num_dark], dtype = datatype)

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
        logger.info('Read {} object/flat frames'.format(num_flat))
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
        num_object = len(paths_object)
        logger.info('Read {} object/flat frames'.format(num_object))
        object3D = np.empty([shape[0], shape[1], num_object], dtype = datatype)
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

#%% function for data reduction

def reduce_dataset(dir_dataset,
                   dir_flat = None,
                   dir_save = '/data/images/fpc_analysis',
                   clobber = False):

    '''
    To reduce a dataset with corresponding flat frames, specify flats in kw.
    To reduce a dataset without flat, leave path_flat as none.
    To reduce a dataset of flat images only, treat it as a dataset w/o flat.

    '''
    import os
    import logging
    from astropy.io import fits

    # check and create master folder if it doesn't exist
    if not os.path.exists(dir_save):
        os.makedirs(dir_save)

    # setup logger
    logger = logging.getLogger()
    logging.basicConfig()
    try:
        logger.removeHandler(logger.handlers[0])
    except:
        # this is the first logger instance with no handle attached
        pass
    logging.basicConfig(filename = os.path.join(dir_save, 'data_reduction.log'),
                        level = logging.INFO,
                        format = '%(asctime)s: %(message)s',
                        datefmt= '%Y/%m/%d %H:%M:%S')
    print('Log file is being saved to :', dir_save)
#    logging.debug('This message should go to the log file')
#    logging.info('So should this')
#    logging.warning('And this, too')

    logger.info('Reducing directory {}...'.format(dir_dataset))

    # determine whether there is flat available
    if dir_flat is None:

        logger.info('No flat provided, reducing given dataset only...')
        reduce_dataset_worker(logger,
                              dir_dataset,
                              dir_save = dir_save,
                              clobber = clobber,
                              dataset_is_flat = False,
                              suppress_return = True)
    else:
        # then we gotta reduce the dataset of flats
        # and then do flat field correction
        logger.info('Flat provided, using flat at {}'.format(dir_flat))
        logger.info('Reducing flat dataset...')
        m_flat_normalised = reduce_dataset_worker(logger,
                                                  dir_flat,
                                                  dir_save = dir_save,
                                                  clobber = clobber,
                                                  dataset_is_flat = True,
                                                  suppress_return = False)
        logger.info('Reducing object dataset...')
        m_object = reduce_dataset_worker(logger,
                                         dir_dataset,
                                         dir_save = dir_save,
                                         clobber = clobber,
                                         dataset_is_flat = False,
                                         suppress_return = False)

        logger.info('Flat field correction...')
        if 'master.fits' in os.listdir(dir_save) and not clobber:
            logger.info('master.fits already exists. Overwriting not allowed.')
        else:
            master = m_object / m_flat_normalised
            hdu_master = fits.PrimaryHDU(master)
            hdu_master.writeto(
                os.path.join(dir_save, 'master.fits'),
                clobber = clobber)
            logger.info('Master image with flat field correction saved to {}'
                .format(
                    os.path.join(dir_save, 'master.fits')))
