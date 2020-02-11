def injrec(path):

    path = f"lustre/eilin/custom_aperture/" + path
    cluster = f"/home/eilin/TESS_UCDs/flare_tables/2020_02_11_TESSUCDs"

    flc = read_custom_aperture_lc(path)

    flcd = flc.detrend("custom", func=custom_detrending)
    
    flcd = flcd.find_flares()
    
    #write_flares_to_file(flcd, cluster)

    flc.plot()
    plt.plot(flcd.time, flcd.detrended_flux+50)
    plt.ylim(500,700)
    
    max_ampl = 2*flcd.flares.ampl_rec.max()
    #print(max_ampl)
    min_dur = .5 * np.min(flcd.flares.tstop-flcd.flares.tstart)
    print(min_dur)
    flce, fake_flc = flcd.sample_flare_recovery(inject_before_detrending=True, mode="custom",
                                                    func=custom_detrending,
                                                  iterations=6, fakefreq=1e-3, ampl=[1e-2, max_ampl],
                                                  dur=[min_dur/6., 0.1/6.], save=True,
                                                  path="{}_{:012d}_s{:04d}.csv".format(5,
                                                                                       flc.targetid,
                                                                                       flc.campaign))
    print("\nFinished TIC {} ({})\n------------------------------\n".format(flc.targetid, flc.campaign))
    