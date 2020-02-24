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
        if ((TIC in p) & ("s"+c in p)):
            print(p)
            df = df.append(pd.read_csv(p), ignore_index=True)
    print(f"{df.shape[0]} synthetic flares injected.")

    # set number of synthetic flares per bin to 20 here:
    bins = int(np.rint(np.sqrt(df.shape[0]/20))) 
    
    print("Get custom LC path")
    # get path to custom LC
    pac = get_customlc(TIC, c, clcs)

    print("Get custom LC")
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
    tstamp = '20200225'
    print(paths[4].split("/")[-1].split("_")[1][1:-4])
    ids = [p.split("/")[-1].split("_")[0] for p in paths]
    cs = [p.split("/")[-1].split("_")[1][1:5] for p in paths]
    print(cs)
    targs = list(zip(ids,cs))
    targs = list(set(targs))
    for i, t in enumerate(targs):
        print(i,t)
    for TIC, c in targs:
        try:    
            print(f"Analysing TIC {TIC} in sector {c}")
            # find and characterize all flares
            flares = wrap_characterization(TIC, c, clcs, paths, tstamp)
            
            # add results to file
            with open(f"/work1/eilin/TESS_UCDs/TESS_UCD_flares/flare_tables/{tstamp}_vetted_flares.csv", "a") as f:
                flares["ID"] = TIC
                flares["sector"] = c
                flares.to_csv(f, index=False, header=False)
        except:
            print(f"\n--------------------\nTIC {TIC}, sector {c} FAILED.\n---------------------\n")
            
