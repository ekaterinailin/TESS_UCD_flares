#python3
#run on TESS LCs

from altaipony.flarelc import FlareLightCurve
from astropy.io import fits

import sys

if __name__ == "__main__":
    path = sys.argv[1]
   # path = "103ca_tess2019090910582-s0008-0000000000926898-0000-s_lc.fits"
    fin = "custom_aperture_lcs/"
    fout = ""
    iters = 10
    hdu = fits.open("{}{}".format(fin,path))

    sector = int(path.split("-")[1][-2:])
    TIC = int(path.split("-")[2])

    print("\nStarted TIC {} ({})\n------------------------------\n".format(TIC, sector))

    flc = FlareLightCurve(time=data["TIME"],
                          flux=data["FLUX"],
                          flux_err=data["FLUX_ERR"],
                          quality=data["QUALITY"],
                          cadenceno=data["CADENCENO"],
                          targetid=TIC,
                          campaign=sector)


    flcd, fake_flc = flc.sample_flare_recovery(inject_before_detrending=True, mode="savgol",
                                              iterations=iters, fakefreq=1e-3, ampl=[1e-2, 0.5],
                                              dur=[.002/6., 0.1/6.], save=True,
                                              path="{}{}_{:012d}_s{:04d}.csv".format(fout,iters,
                                                                                     flc.targetid,
                                                                                     flc.campaign))

    flcc = flcd.characterize_flares( ampl_bins=80, dur_bins=160)

    flcc.flares.to_csv("{}characterized_flares_{:012d}_s{:04d}.csv".format(fout,flc.targetid, flc.campaign))

    print("\nFinished TIC {} ({})\n------------------------------\n".format(TIC, sector))