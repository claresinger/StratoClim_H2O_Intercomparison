import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import seaborn
import datetime

def h2o_f2_f7(dat):
    fig,axes = plt.subplots(2,2,figsize=(16,12),gridspec_kw={"width_ratios":[1,3]})
    #c = plt.rcParams['axes.prop_cycle'].by_key()['color']
    #c = ["#FF1F58","#009ADE","#FFC61E"]
    c = ["#009ADE","#FF1F58","k"]
    ms = 4.0
    
    flno = [2,3,4,6,7,8]
    maxlag = [0,0,5,10,10,20]
    
    # add cloudy flag
    dat['CLOUDY'] = ((dat['NICE'] > 0) | (dat['MASBR'] >= 1.2)).astype(int)
    for i,f in enumerate(flno):
        for lag in np.arange(1,maxlag[i]):
            dat.loc[(dat['FLIGHT'] == f),'CLOUDY'] = np.maximum(dat.loc[(dat['FLIGHT'] == f),'CLOUDY'], 
                                                  dat[(dat['FLIGHT'] == f)].shift(periods=lag, fill_value=0.0)['CLOUDY'])
    
    # add ascent/descent flag
    dz = (dat['ALT'] - dat.shift(periods=1)['ALT'])*1e3
    dt = dat['TIME'] - dat.shift(periods=1)['TIME']
    vert = np.abs(dz / dt)
    vert_avg = vert.rolling(window=20).mean()
    dat['ASCENT_FLAG'] = ((vert_avg > 10) | (dat['ALT'] < 12)).astype(int)
    
    # add chiwis flag
    dat['CELL_FLAG'] = ((dat['PRES_CELL'] < 30.0) | (dat['PRES_CELL'] > 45.0) | (dat['FLAG'] == 1)).astype(int)
    dat['CELL_OOR'] = ((dat['PRES_CELL'] < 20.0) | (dat['PRES_CELL'] > 30.0) | (dat['FLAG'] == 1)).astype(int)
    
    # FL7 dive flag
    dat['F7_DIVE'] = ((dat['FLIGHT'] == 7) & (dat['TIME'] > 19.9e3) & (dat['TIME'] < 20.2e3)).astype('int')
    
    for f in [2,7]:
        if f == 2:
            ax0 = axes[0,0]
            ax1 = axes[0,1]
            sub0 = "a"
            sub1 = "b"
            t1 = 9.15
            t2 = 10.9
            t2 = 12.5
            ptlim = [360,480]
            wvlim = [3,14]
            wvtks = [4,8,12]
            maxlag = 0
        if f == 7:
            ax0 = axes[1,0]
            ax1 = axes[1,1]
            sub0 = "c"
            sub1 = "d"
            t1 = 10.58
            t2 = 12.2
            ptlim = [360,440]
            wvlim = [3,10]
            wvtks = [4,6,8]
            maxlag = 10
            
        dati = dat[(dat['FLIGHT'] == f) & (dat['ASCENT_FLAG'] == 0)]
        dat_clr_fish = dat[(dat['FLIGHT'] == f) & (dat['CLOUDY'] == 0) & (dat['ASCENT_FLAG'] == 0)]
        dat_flash = dat[(dat['FLIGHT'] == f) & (dat['ASCENT_FLAG'] == 0) & (dat['F7_DIVE'] == 0)
                        & ((dat['TIME']/3600.+5.75 > t1) | (dat['TIME']/3600.+5.75 < t2))]
        dat_chiwis = dat[(dat['FLIGHT'] == f) & (dat['CELL_FLAG'] == 0) & (dat['ASCENT_FLAG'] == 0)]
        dat_chiwis_oor = dat[(dat['FLIGHT'] == f) & (dat['CELL_OOR'] == 0) & (dat['ASCENT_FLAG'] == 0)]

        ax0.plot(dat_clr_fish['FIH2O'],dat_clr_fish['PT'],'.', ms=ms, color=c[1], label="clear-sky FISH")
        ax0.plot(dat_flash['FLH2O'],dat_flash['PT'],'.', ms=ms, color=c[0], label="FLASH")
        ax0.plot(dat_chiwis['H2O'],dat_chiwis['PT'],'.', ms=ms, color=c[2], label="ChiWIS")
        ax0.plot(dat_chiwis_oor['H2O'], dat_chiwis_oor['PT'],'.', ms=ms, color='grey')
        
        ax0.grid(which='major',linestyle=':')
        ax0.set_xlim(wvlim)
        ax0.set_xticks(wvtks)
        if f == 7:
            ax0.set_xlabel(r"H$_2$O (ppmv)")
        ax0.set_ylim(ptlim)
        ax0.set_ylabel("Potential Temperature (K)")
        ax0.set_title(sub0,weight="bold",loc="left")
        if f == 2:
            ax0.legend(loc=1, markerscale=3.0, labelspacing=0.4, handletextpad=0.1, fontsize=15)
        ax0.set_title("Flight "+str(f))
        
        ax1.plot(dat_clr_fish['TIME']/3600.+5.75,dat_clr_fish['FIH2O'],'.', ms=ms, color=c[1],label="clear-sky FISH")
        ax1.plot(dat_flash['TIME']/3600.+5.75,dat_flash['FLH2O'],'.', ms=ms, color=c[0],label="FLASH")
        ax1.plot(dat_chiwis['TIME']/3600.+5.75,dat_chiwis['H2O'],'.', ms=ms, color=c[2],label="ChiWIS")
        ax1.plot(dat_chiwis_oor['TIME']/3600.+5.75, dat_chiwis_oor['H2O'], '.', ms=ms, color="grey")
        
        if f == 7:
            ax1.set_xlabel("Kathmandu Local Time")
        ax1.set_xlim([t1,t2])
        locs = ax1.get_xticks()        
        labels = [str(datetime.timedelta(hours=x)).rsplit(':',1)[0] for x in locs]
        ax1.set_xticklabels(labels)
        ax1.set_ylabel(r"H$_2$O (ppmv)")
        ax1.set_ylim(wvlim)
        ax1.grid(which='major',linestyle=':')
        ax1.set_title(sub1,weight="bold",loc="left")
        ax1.set_title("Flight "+str(f))

        ax2 = ax1.twinx()
        dat_alt = dat[(dat['FLIGHT'] == f) & (dat['ALT'] > 10)]
        ax2.plot(dat_alt['TIME']/3600.+5.75,dat_alt['ALT'],"-",color="green",lw=2)
        ax2.set_ylabel("Altitude (km)",color="green")
        ax2.tick_params(axis='y', colors="green")
        ax2.set_yticks([10,12,14,16,18,20])

    plt.rcParams.update({'font.size': 20})
    fig.tight_layout()
    plt.savefig("./Paper-Figures/fig5-f2f7.png",dpi=300)
    plt.show()