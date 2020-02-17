# Inject and recover random synthetic flares.
# See FlareLightCurve.sample_flare_recovery for details
#
# Usage:
#
# python sample_injrec.py <path to light curve>
# 


from funcs.helper import read_custom_aperture_lc
from funcs.custom_detrending import custom_detrending

import sys
import numpy as np

def injrec(path, iterations=1000):
    
    # read the custom light curve
    flc = read_custom_aperture_lc(path, typ="custom")
    
    # alternatively, read a generic light curve
    # flc = read_custom_aperture_lc(path, typ="generic", mission="TESS", mode="LC")

    # De-trend the light curve
    flcd = flc.detrend("custom", func=custom_detrending)

    # alternatively, use a generic method, like a Savitzky-Golay filter
    # flcd = flc.detrend("savgol")
        
    # Find flares in the detrended light curve
    flcd = flcd.find_flares()

    # Constrain the paramater space so that it will cover the flare candidates
    max_ampl = 2 * flcd.flares.ampl_rec.max()
    min_dur = .5 * np.min(flcd.flares.tstop - flcd.flares.tstart)

    # if no flares are detected, set default values for minimum duration and maximum amplitude
    if np.isnan(max_ampl):
        max_ampl=.5

    if np.isnan(min_dur):
        min_dur=1/24./20.

    # Adjust mode and func if needed here:
    # See 
    # https://altaipony.readthedocs.io/en/latest/api/altaipony.flarelc.FlareLightCurve.html#altaipony.flarelc.FlareLightCurve.sample_flare_recovery 
    # for details
    flce, fake_flc = flcd.sample_flare_recovery(inject_before_detrending=True, mode="custom",
                                                func=custom_detrending,
                                                iterations=iterations, fakefreq=1e-3, ampl=[1e-2, max_ampl],
                                                 dur=[min_dur/6., 0.1], save=True,
                                                 path="{:012d}_s{:04d}.csv".format(flc.targetid,
                                                                                   flc.campaign))

if __name__ == "__main__":

    path = sys.argv[1]
    injrec(path)

