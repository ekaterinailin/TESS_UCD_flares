import numpy as np
import pytest

from ..helper import (read_custom_aperture_lc,
                      fix_mask,
                      no_nan_inf,
                      #get_window_length_dict,
                      #write_flares_to_file
                      )


from astropy.io import fits

from altaipony.flarelc import FlareLightCurve

import os

CWD = "/".join(os.getcwd().split("/")[:-1])

print("CWDCWDCWDCWD: ", CWD) 

# -------------------------------- TESTING read_custom_aperture(path) ----------------------------

def test_read_custom_aperture_lc():
    
    # Create a light curve with a sinusoidal modulation
    hdr = fits.Header()
    hdr['OBSERVER'] = 'Mike'
    hdr['COMMENT'] = "Here's some commentary about this FITS file."
    start, stop, N = 1000, 1020, 1000
    c1 = fits.Column(name='TIME', array=np.linspace(start, stop, N), format='F10.4')
    c2 = fits.Column(name='FLUX', array=400 + 50*np.sin(np.linspace(start, stop, N)*2), format='F10.4')
    c3 = fits.Column(name='FLUX_ERR', array= 20*np.random.rand(N), format='F10.4')
    c4 = fits.Column(name='QUALITY', array= np.full(N,1), format='K')
    c5 = fits.Column(name='CADENCENO', array= np.arange(100,N+100), format='K')
    hdu = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5], header=hdr)
    PATH = f"{CWD}/notebook/funcs/tests/test.fits"
    hdu.writeto(PATH,overwrite=True)
    
    # Call the function
    flc = read_custom_aperture_lc(PATH, typ="custom", mission="TESS", mode="LC",
                                sector=10, TIC=1000)
    # Do some checks
    assert (flc.flux == 400 + 50*np.sin(np.linspace(start, stop, N)*2)).all()
    assert (flc.time == np.linspace(start, stop, N)).all()
    assert flc.campaign == 10
    assert flc.targetid == 1000
    
# -------------------------------- TESTING no_nan_inf(flc) ----------------------------

cases = [((0, 0, np.nan), False),
         ((0, 0, 0), True), 
         ((0, 3, 0, np.linspace(0, 1, 10 )), True),
         ((0, 3, 0, np.full(10,np.inf)), False),
         ((np.inf, 0, 0, 0), False),
         (np.array([9, 1, np.nan]), False),
         (np.array([9, 1, np.inf]), False),]

@pytest.mark.parametrize("l,expected", cases)
def test_no_nan_inf_succeed(l,expected):
    assert no_nan_inf(l) == expected

# -------------------------------- TESTING fix_mask(flc) ----------------------------

def test_fix_mask():
    # Select cadenceno range
    start, stop = int(1e5),int(3e5)

    # Define light curve
    c = np.arange(start, stop)
    t = np.linspace(1000,1030,stop-start)
    f = np.random.rand(stop-start)
    flc = FlareLightCurve(time=t, flux=f, cadenceno=c, campaign=10)

    # Call function
    flcc = fix_mask(flc)

    # Do some checks
    res = flcc.cadenceno[np.isnan(flcc.flux)] 
    assert (((res >= 246227) & (res <= 247440)) | ((res >= 255110) & (res <= 257370))).all()


    # A different case where the campaign has no custom mask
    c = np.arange(start, stop)
    t = np.linspace(1000,1030,stop-start)
    f = np.random.rand(stop-start)
    flc2 = FlareLightCurve(time=t, flux=f, cadenceno=c, campaign=18)
    flcc2 = fix_mask(flc2)
    assert (np.isnan(flcc2.flux) == False).all()
