import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import seaborn
import datetime

flno = [2,3,4,6,7,8]
#colors = ['black','red','orange','lime','black','green','blue','purple']
colors = ["k","#045275","#0C7BDC","#7CCBA2","k","#FED976","#F0746E","#7C1D6F"]
c = ["#009ADE","#FF1F58","k"]

def h2o_timeseries_all(dat):
    fig,axes = plt.subplots(6,1,figsize=(12,20))
    plt.rcParams.update({'font.size': 20})

    t1 = [9.12, 9.5, 14.5, 13.55, 10.55, 15.1]
    t2 = [12.5, 12,  17.85, 16.4, 12.2, 17.7]
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
    
    for i,f in enumerate(flno):
        ax1 = axes[i]
        
        time = dat['TIME']/3600.+5.75
        dati = dat[(dat['FLIGHT'] == f) & (dat['ASCENT_FLAG'] == 0) & (time > t1[i]) & (time < t2[i])]
        dat_chiwis = dati[(dati['CELL_FLAG'] == 0)]
        dat_chiwis_oor = dati[(dati['CELL_OOR'] == 0)]
        dat_flash = dati[(dati['F7_DIVE'] == 0)]
        dat_clr_fish = dati[(dati['CLOUDY'] == 0)]
        ax1.plot(dat_flash['TIME']/3600.+5.75,dat_flash['FLH2O'],'.',color=c[0],label="FLASH")
        ax1.plot(dat_clr_fish['TIME']/3600.+5.75,dat_clr_fish['FIH2O'],'.',color=c[1],label="clear-sky FISH")
        ax1.plot(dat_chiwis['TIME']/3600.+5.75,dat_chiwis['H2O'],'.',color=c[2],label="ChiWIS")
        ax1.plot(dat_chiwis_oor['TIME']/3600.+5.75,dat_chiwis_oor['H2O'],'.',color='grey')

        if i == len(flno)-1:
            ax1.set_xlabel("Kathmandu Local Time")
        if i == 0:
            ax1.legend(loc=(0.7,0.45),markerscale=3.0, labelspacing=0.4, handletextpad=0.1, fontsize=15)
        ax1.set_ylabel(r"H$_2$O (ppmv)")
        
        ax1.set_title("Flight {0}".format(str(f)))
        ax1.grid(which='major',linestyle=':')

        ax1.set_xlim([t1[i],t2[i]])
        locs = ax1.get_xticks()        
        labels = [str(datetime.timedelta(hours=x)).rsplit(':',1)[0] for x in locs]
        ax1.set_xticklabels(labels)
        
        ax1.set_ylim([2,14])
        ax1.set_yticks([2,8,14])

        ax2 = ax1.twinx()
        dat_alt = dat[(dat['FLIGHT'] == f) & (dat['ALT'] > 10)]
        ax2.plot(dat_alt['TIME']/3600.+5.75,dat_alt['ALT'],"-",color="green",lw=2)
        ax2.set_ylabel("Altitude (km)",color="green")
        ax2.tick_params(axis='y', colors="green")
        ax2.set_yticks([10,12,14,16,18,20])

    fig.tight_layout()
    plt.savefig("./Paper-Figures/supp-h2o-timeseries.png",dpi=300)
    plt.show()