#python3
#run on TESS LCs

from altaipony.flarelc import FlareLightCurve
from astropy.io import fits

from funcs.funcs import read_custom_aperture_lc
from funcs.custom_detrending import custom_detrending

import sys

if __name__ == "__main__":
    path = sys.argv[1]
   # path = "103ca_tess2019090910582-s0008-0000000000926898-0000-s_lc.fits"
    fin =  "/home/ekaterina/Documents/001_Science/TESS_UCDs/TESS_UCD_flares/custom_aperture/"#"custom_aperture_lcs/"
    fout = ""
    iters = 2

    sector = int(path.split("-")[1][-2:])
    TIC = int(path.split("-")[2])

    print("\nStarted TIC {} ({})\n------------------------------\n".format(TIC, sector))


    flc = read_custom_aperture_lc(fin+path)
    
    
    #flcd = custom_detrending(flc)
    print(flc.flux[:10])

    flcd, fake_flc = flc.sample_flare_recovery(inject_before_detrending=True, mode="custom",
                                               func=custom_detrending,
                                              iterations=iters, fakefreq=1e-3, ampl=[1e-4, 0.5],
                                              dur=[.002/6., 0.12/6.], save=True,
                                              path="{}{}_{:012d}_s{:04d}.csv".format(fout,iters,
                                                                                     flc.targetid,
                                                                                     flc.campaign))
    print(flcd.flux[:10])
    
    print("\nFinished TIC {} ({})\n------------------------------\n".format(TIC, sector))