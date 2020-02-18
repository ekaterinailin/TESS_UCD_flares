import numpy as np
import pandas as pd

from funcs.helper import read_custom_aperture_lc
from funcs.custom_detrending import custom_detrending

    
import glob
import os

 
def get_customlc(TIC, c, clcs):
    
    for p in clcs:
        pp = p.split("/")[-1]
       # print(p, TIC, c)
        if (int(TIC) == int(pp.split("-")[2])) & (int(c) == int(pp.split("-")[1][1:])):
            pac = p
    return pac

def wrap_characterization(TIC, c, clcs, paths, tstamp):
    
    df = pd.DataFrame()
    for p in paths:
        if ((TIC in p) & (c in p)):
            df = df.append(pd.read_csv(p), ignore_index=True)

    # set number of synthetic flares per bin to 20 here:
    bins = int(np.rint(np.sqrt(df.shape[0]/20))) 
    
    # get path to custom LC
    pac = get_customlc(TIC, c, clcs)

    # read custom LC
    flc = read_custom_aperture_lc(pac)

    # detrend custom LC with custom prescription
    print("Started de-trending ...")
    flcd = flc.detrend("custom", func=custom_detrending)
    print("... Done.")
        
    # detect flare candidates
    flcd = flcd.find_flares()
    
    # add fake flares from previous sampling
    flcd.fake_flares = df
    
    if flcd.flares.shape[0] > 0:
        # characterize flare candidates
        flcd = flcd.characterize_flares(ampl_bins=bins, dur_bins=bins)
    
    return flcd.flares

   
if __name__ == "__main__":
    

    CWD = "/".join(os.getcwd().split("/")[:-1])
    clcs = glob.glob(CWD+"/custom_aperture/*.fits")
    paths = glob.glob(f"{CWD}/injrec/*.csv")
    tstamp = '20200218'
    
    ids = [p.split("/")[-1].split("_")[0]for p in paths]
    cs = [p.split("/")[-1].split("_")[1][1:-4]for p in paths]
    targs = list(zip(ids,cs))
    targs = list(set(targs))
    
    for TIC, c in targs:
        
        print(f"Analysing TIC {TIC} in sector {c}")
        # find and characterize all flares
        flares = wrap_characterization(TIC, c, clcs, paths, tstamp)
        
        # add results to file
        with open(f"/work1/eilin/TESS_UCDs/TESS_UCD_flares/flare_tables/{tstamp}_vetted_flares.csv", "a") as f:
            flares.to_csv(f, index=False, header=False)
        
