import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

flno = [2,3,4,6,7,8]
maxlag = [0,0,5,10,10,20]
cmap = "YlGnBu"

def satppm(p,T):
    Psat = np.exp(9.550426 - 5723.265/T + 3.53068 * np.log(T) - 0.00728332 * T) / 100.
    return 1e6 * Psat / p
    
def rhi_temp_freq_oor(dat):
    fig = plt.figure(figsize=(16.3,6))
    plt.rcParams.update({"font.size":25})
    gs0 = gridspec.GridSpec(1,3,figure=fig,width_ratios=[10,10,1])
    cax = fig.add_subplot(gs0[-1])
    inst = "ChiWIS"
    key = "H2O"
    
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
    dat['CELL_GOOD'] = ((dat['PRES_CELL'] > 30.0) & (dat['PRES_CELL'] < 45.0) & (dat['FLAG'] == 0)).astype(int)
    dat['CELL_ALL'] = ((dat['PRES_CELL'] > 20.0) & (dat['PRES_CELL'] < 45.0) & (dat['FLAG'] == 0)).astype(int)
    
    # FL7 dive flag
    dat['F7_DIVE'] = ((dat['FLIGHT'] == 7) & (dat['TIME'] > 19.9e3) & (dat['TIME'] < 20.2e3)).astype('int')
    
    h2othresh = 100
    
    for i in [0,1]:
        ax = fig.add_subplot(gs0[i])
        if i == 0:
            dati = dat[(dat[key] < h2othresh) & (dat['CLOUDY'] == 0) & (dat['ASCENT_FLAG'] == 0) & (dat['CELL_GOOD'] == 1)]
            ax.set_ylabel(r'RH$_{ice}$')
        if i == 1:
            dati = dat[(dat[key] < h2othresh) & (dat['CLOUDY'] == 0) & (dat['ASCENT_FLAG'] == 0) & (dat['CELL_ALL'] == 1)]
            plt.setp(ax.get_yticklabels(), visible=False)
        
        hi = np.sum(~np.isnan(dati[key])) / 3600.
        ax.text(186, 0.05, "{:.1f} hrs".format(hi), fontsize=20)
        ax.set_xlabel('Temperature (K)')
        ax.text(207,1.85,"clear-sky "+inst,fontsize=20)
        
        temp = dati['TEMP']+273.15
        rhi = dati[key]/dati['SATPPM']
        print(inst, np.nanmean(rhi), np.nanmedian(rhi))

        m = ax.hist2d(temp,rhi,bins=[120,40],range=[[180, 240],[0, 2]],
            density=1,norm=mcolors.PowerNorm(gamma=0.5),vmin=0,vmax=1,
            cmin=1e-10,cmap=cmap)

        ax.plot([175,240],[1.67,1.42],'k:', lw=3)
        ax.plot([175,240],[1,1],'k-')
        
        ax.grid(which='major',linestyle=':')
        ax.set_xlim([185,225])
        ax.set_ylim([0,2])
        
        let = ["a","b"]
        ax.set_title(let[i],loc="left",weight="bold")
        cp = [30,20]
        ax.set_title("cell pressure$\geq${:.0f} mbar".format(cp[i]), fontsize=25)
        
    plt.colorbar(m[3], cax=cax, ticks=[0,0.1,0.3,0.6,1.0], label="Frequency")
    plt.rcParams.update({"font.size":25})
    plt.savefig("./Paper-Figures/supp-rhi-freq-oor.png",bbox_inches="tight",dpi=200)
    plt.show()
