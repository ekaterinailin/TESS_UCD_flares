# Python3
# Some functions used for flare finding, data management


import pandas as pd
import numpy as np

from astropy.io import fits
from altaipony.flarelc import FlareLightCurve
import sys, os

def write_flares_to_file(flc, cluster):
    '''
    Write the resulting flare table to file, adding
    it  to the rest.
    '''
    cols = ['TIC', 'ampl_rec', 'cstart', 'cstop', 
            'ed_rec', 'ed_rec_err', 'istart', 'istop', 
            'tstart', 'tstop', 'Campaign']
            #'saturation_f10']
    try:
        df = pd.read_csv('{0}_flares.csv'.format(cluster),
            usecols=cols)
    except:
        print('Create a file named {}_flares.csv'.format(cluster))
        df = pd.DataFrame(columns=cols)
        
    lbefore = flc.flares[~flc.flares.ed_rec.isnull()].shape[0]
    if lbefore == 0:
        emptydf = pd.DataFrame(dict(zip(cols, [flc.targetid] + [np.nan] * (len(cols) - 1))), index=[0])
        flc.flares = flc.flares.append(emptydf, ignore_index=True)
        print('Added an empty row for TIC {}'.format(flc.targetid))
    flc.flares['TIC'] = flc.targetid
    flc.flares['Campaign'] = flc.campaign
    #flc.flares = flc.get_saturation(return_level=True) maybe later
    print(flc.flares)
    df = df.append(flc.flares, ignore_index=True)
    lafter = df[~df.ed_rec.isnull()].shape[0]
   # print('Added {} flares to the {} found so far in {}.'.format(lafter-lbefore, lbefore, cluster))
    df.to_csv('{0}_flares.csv'.format(cluster), index=False)
    return
    
def read_custom_aperture_lc(path):
    '''Read in custom aperture light curve
    from TESS. Needs specific naming convention.
    Applies pre-defined quality masks.
    
    Parameters:
    -----------
    path : str
        path to file
    
    Returns:
    --------
    FlareLightCurve
    '''
    hdu = fits.open(path)
    data = hdu[1].data

    sector = int(path.split("-")[1][-2:])
    TIC = int(path.split("-")[2])

    flc = FlareLightCurve(time=data["TIME"],
                        flux=data["FLUX"],
                        flux_err=data["FLUX_ERR"],
                        quality=data["QUALITY"],
                        cadenceno=data["CADENCENO"],
                        targetid=TIC,
                        campaign=sector)
    flc = fix_mask(flc)
    return flc

def fix_mask(flc):
    '''Here the masks for different TESS 
    sectors are defined and applied to 
    light curve fluxes.
    
    Parameters:
    ------------
    flc : FlareLightCurve
    
    Returns:
    ----------
    FlareLightCurve
    '''
    masks = {9: [(227352, 228550), (236220, 238250)],
             10: [(246227,247440),(255110,257370)],
             11:[(265912,268250),(275210,278500)],
             8: [(208722,209250)],
             6: [(179661,180680)],
             5: [(151586,151900),(160254,161353),(170000,170519)],
             4: [(140700,141161),(150652,150764)],
             12: [(286200,286300)]}

    if flc.campaign in masks.keys():
        for sta, fin in masks[flc.campaign]:
            flc.flux[np.where((flc.cadenceno >= sta) & (flc.cadenceno <= fin))] = np.nan
            
    flc.quality[:] = np.isnan(flc.flux)
    return flc
