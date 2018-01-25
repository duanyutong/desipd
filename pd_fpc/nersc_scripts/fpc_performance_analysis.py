import os
from multiprocessing import Pool
import numpy as np
import pandas as pd
from astropy.io import fits
from photutils import EllipticalAperture, aperture_photometry
# metadata
DIR_SAVE = '/global/project/projectdirs/desi/protodesi/images/fpc_analysis/fpc_performance/'
DIR_FITS = '/global/project/projectdirs/desi/protodesi/images/fits/'
dataset_names = ['0.2s', '0.4s', '0.6s', '0.8s', '1.0s', '2.0s', '3.0s', '4.0s', '5.0s', '6.0s', '7.0s', '8.0s', '9.0s', '10.0s', '11.0s', '12.0s', '13.0s', '14.0s', '15.0s']
exptype = ['dark', 'light']
fibre_ids = ['3001', '3002', '3003']
ISOPHOTAL_EXT = 5
FIBRE_COORDS_FLAT =  {'3001': (1202, 1089), # when at white spot
                      '3002': (1860, 1090),
                      '3003': (1205, 1404)}
FIBRE_ELLIPTICAL_APERTURES_FLAT = {
            '3001':
                {'a': 12.4090*ISOPHOTAL_EXT,
                 'b':  9.1608*ISOPHOTAL_EXT,
                 'theta': 0.063960888},
            '3002':
                {'a': 11.3185*ISOPHOTAL_EXT,
                 'b':  8.7111*ISOPHOTAL_EXT,
                 'theta': -0.059023251},
            '3003':
                {'a': 16.6673*ISOPHOTAL_EXT,
                 'b': 11.8841*ISOPHOTAL_EXT,
                 'theta': -0.135632276}}
# initialise dataframe
iterable_names = ['dataset name', 'exptype', 'fibre id']
iterables = [dataset_names, exptype, fibre_ids]
index = pd.MultiIndex.from_product(iterables, names=iterable_names)
columns = list(np.arange(11))+['mean', 'std', 'percent']
df_aper_sum = pd.DataFrame(index = index, columns = columns).sort_index()
# generate three apertures
apertures = {}
for fibre_id in fibre_ids:
    position = FIBRE_COORDS_FLAT[fibre_id]
    a = FIBRE_ELLIPTICAL_APERTURES_FLAT[fibre_id]['a']
    b = FIBRE_ELLIPTICAL_APERTURES_FLAT[fibre_id]['b']
    theta = FIBRE_ELLIPTICAL_APERTURES_FLAT[fibre_id]['theta']
    apertures[fibre_id] = EllipticalAperture(position, a, b, theta=theta)

def expid_to_filename(expid_input):

    if isinstance(expid_input, (int, np.integer)):
        filename = 'PROTODESI_' + str(expid_input).zfill(8) + '.fits'
        return filename

    else:
        filenames = [None]*len(expid_input)
        for (i, expid) in enumerate(expid_input):
            filenames[i] = expid_to_filename(expid)
        return filenames

def expid_to_path(expid_input):

    if isinstance(expid_input, (int, np.integer)):
        path = os.path.join(DIR_FITS,
                            expid_to_filename(expid_input))
        return path
    else:
        paths_abs = [None]*len(expid_input)
        for (i, expid) in enumerate(expid_input):
            paths_abs[i] = expid_to_path(expid)
        return paths_abs
    
def fpc_performance_analysis(dataset_name, exptype, expid_initial):

    expids = list(np.arange(expid_initial, expid_initial+11, 1))
    paths = expid_to_path(expids)

    
    for i, path in enumerate(paths):
        data = fits.open(path)['FPC']
        for fibre_id in fibre_ids:
            aperture = apertures[fibre_id]
            phot_table = aperture_photometry(data, aperture)
            df_aper_sum.loc[dataset_name, exptype, fibre_id][i]=phot_table['aperture_sum'][0]
    # calculate std
    for fibre_id in fibre_ids:
        df_aper_sum.loc[dataset_name, exptype, fibre_id]['mean'] = \
            np.mean(df_aper_sum.loc[dataset_name, exptype, fibre_id][0:11])
        df_aper_sum.loc[dataset_name, exptype, fibre_id]['std'] = \
            np.std(df_aper_sum.loc[dataset_name, exptype, fibre_id][0:11])
        df_aper_sum.loc[dataset_name, exptype, fibre_id]['percent'] = \
            df_aper_sum.loc[dataset_name, exptype, fibre_id]['std'] \
            / df_aper_sum.loc[dataset_name, exptype, fibre_id]['mean']

if __name__ == '__main__':

    fpc_performance_analysis('0.2s', 'dark', 7165)
    fpc_performance_analysis('0.2s', 'light', 7154)
    fpc_performance_analysis('0.4s', 'dark', 7132)
    fpc_performance_analysis('0.4s', 'light', 7121)
    fpc_performance_analysis('0.6s', 'dark', 7097)
    fpc_performance_analysis('0.6s', 'light', 7143)
    fpc_performance_analysis('0.8s', 'dark', 7075)
    fpc_performance_analysis('0.8s', 'light', 7086)
    fpc_performance_analysis('1.0s', 'dark', 7064)
    fpc_performance_analysis('1.0s', 'light', 7053)
    fpc_performance_analysis('2.0s', 'dark', 7031)
    fpc_performance_analysis('2.0s', 'light', 7042)
    fpc_performance_analysis('3.0s', 'dark', 7020)
    fpc_performance_analysis('3.0s', 'light', 7009)
    fpc_performance_analysis('4.0s', 'dark', 6987)
    fpc_performance_analysis('4.0s', 'light', 6998)
    fpc_performance_analysis('5.0s', 'dark', 6969)
    fpc_performance_analysis('5.0s', 'light', 6958)
    fpc_performance_analysis('6.0s', 'dark', 6936)
    fpc_performance_analysis('6.0s', 'light', 6947)
    fpc_performance_analysis('7.0s', 'dark', 6925)
    fpc_performance_analysis('7.0s', 'light', 6907)
    fpc_performance_analysis('8.0s', 'dark', 6896)
    fpc_performance_analysis('8.0s', 'light', 6885)
    fpc_performance_analysis('9.0s', 'dark', 6863)
    fpc_performance_analysis('9.0s', 'light', 6874)
    fpc_performance_analysis('10.0s', 'dark', 6172)
    fpc_performance_analysis('10.0s', 'light', 6161)
    fpc_performance_analysis('11.0s', 'dark', 6834)
    fpc_performance_analysis('11.0s', 'light', 6845)
    fpc_performance_analysis('12.0s', 'dark', 6823)
    fpc_performance_analysis('12.0s', 'light', 6812)
    fpc_performance_analysis('13.0s', 'dark', 6783)
    fpc_performance_analysis('13.0s', 'light', 6801)
    fpc_performance_analysis('14.0s', 'dark', 6772)
    fpc_performance_analysis('14.0s', 'light', 6761)
    fpc_performance_analysis('15.0s', 'dark', 6750)
    fpc_performance_analysis('15.0s', 'light', 6739
                             )
    df_aper_sum.to_csv(os.path.join(DIR_SAVE, 'sum_table.csv'))