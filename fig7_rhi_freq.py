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

def rhi_temp_freq(dat):
    fig = plt.figure(figsize=(24,12))
    plt.rcParams.update({"font.size":25})
    gs0 = gridspec.GridSpec(1,4,figure=fig,width_ratios=[10,10,10,1])
    cax = fig.add_subplot(gs0[-1])
    insts = ["ChiWIS","FLASH","FISH"]
    
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
    
    # FL7 dive flag
    dat['F7_DIVE'] = ((dat['FLIGHT'] == 7) & (dat['TIME'] > 19.9e3) & (dat['TIME'] < 20.2e3)).astype('int')
    
    h2othresh = 100
    
    for i, inst in enumerate(insts):
        gs00 = gridspec.GridSpecFromSubplotSpec(2,1,subplot_spec=gs0[i])
        
        if i == 0:
            ax1 = fig.add_subplot(gs00[0])
            plt.setp(ax1.get_xticklabels(), visible=False)
            axlt = ax1
            ax2 = fig.add_subplot(gs00[1], sharex=ax1)
            axlb = ax2
            
            inst = "ChiWIS"
            key = "H2O"
            # chiwis clear-sky
            dat0 = dat[(dat[key] < h2othresh) & (dat['CLOUDY'] == 0) & (dat['ASCENT_FLAG'] == 0) & (dat['CELL_FLAG'] == 0)]
            # chiwis cloudy
            dat1 = dat[(dat[key] < h2othresh) & (dat['CLOUDY'] == 1) & (dat['ASCENT_FLAG'] == 0) & (dat['CELL_FLAG'] == 0)]
            
            h0 = np.sum(~np.isnan(dat0[key])) / 3600.
            h1 = np.sum(~np.isnan(dat1[key])) / 3600.
            print(inst)
            print(h0, h1, " hrs (clear, cloudy)")
            
        if i == 1:
            ax1 = fig.add_subplot(gs00[0], sharey=axlt)
            ax2 = fig.add_subplot(gs00[1], sharex=ax1, sharey=axlb)
            plt.setp(ax1.get_xticklabels(), visible=False)
            plt.setp(ax1.get_yticklabels(), visible=False)
            plt.setp(ax2.get_yticklabels(), visible=False)
            
            inst = "FLASH"
            key = "FLH2O"
            # flash clear-sky
            dat0 = dat[(dat[key] < h2othresh) & (dat['FLIGHT'] != 1) & (dat['FLIGHT'] != 5) 
                       & (dat['CLOUDY'] == 0) & (dat['ASCENT_FLAG'] == 0) & (dat['F7_DIVE'] == 0)]
            # flash cloudy
            dat1 = dat[(dat[key] < h2othresh) & (dat['FLIGHT'] != 1) & (dat['FLIGHT'] != 5) 
                       & (dat['CLOUDY'] == 1) & (dat['ASCENT_FLAG'] == 0) & (dat['F7_DIVE'] == 0)]
            
            h0 = np.sum(~np.isnan(dat0[key])) / 3600.
            h1 = np.sum(~np.isnan(dat1[key])) / 3600.
            print(inst)
            print(h0, h1, " hrs (clear, cloudy)")
            
        if i == 2:
            ax1 = fig.add_subplot(gs00[0], sharey=axlt)
            plt.setp(ax1.get_yticklabels(), visible=False)
            
            inst = "FISH"
            key = "FIH2O"
            # fish clear-sky
            dat0 = dat[(dat[key] < h2othresh) & (dat['CLOUDY'] == 0) & (dat['ASCENT_FLAG'] == 0)]
            
            h0 = np.sum(~np.isnan(dat0[key])) / 3600.
            print(inst)
            print(h0, " hrs (clear)")
        
        temp = dat0['TEMP']+273.15
        rhi = dat0[key]/dat0['SATPPM']
        print("clear ", inst, "mean = ", np.nanmean(rhi), "median = ", np.nanmedian(rhi))

        m = ax1.hist2d(temp,rhi,bins=[120,40],range=[[180, 240],[0, 2]],
            density=1,norm=mcolors.PowerNorm(gamma=0.5),vmin=0,vmax=1,
            cmin=1e-10,cmap=cmap)

        ax1.plot([175,240],[1.67,1.42],'k:', lw=3)
        ax1.plot([175,240],[1,1],'k-', lw=2)
        ax1.text(186,0.05,"{:.1f} hrs".format(h0))

        if i < 2:
            temp = dat1['TEMP']+273.15
            rhi = dat1[key]/dat1['SATPPM']
            print("cloudy ", inst, "mean = ", np.nanmean(rhi), "median = ", np.nanmedian(rhi))

            m = ax2.hist2d(temp,rhi,bins=[120,40],range=[[180, 240],[0, 2]],
                density=1,norm=mcolors.PowerNorm(gamma=0.5),vmin=0,vmax=1,
                cmin=1e-10,cmap=cmap)

            ax2.plot([175,240],[1.67,1.42],'k:', lw=3, label="Koop homogenous\nice nucleation threshold")
            ax2.plot([175,240],[1,1],'k-',lw=2,label=r"saturation line, RH$_{ice} = 1$")
            ax2.text(186,0.05,"{:.1f} hrs".format(h1))

            # plot attributes
            ax2.set_xlabel('Temperature (K)')
            if i == 0:
                ax1.set_ylabel(r'RH$_{ice}$')
                ax2.set_ylabel(r'RH$_{ice}$')
        
            ax1.set_title("clear-sky "+inst)
            ax2.set_title("in-cloud "+inst)
        
            ax1.grid(which='major',linestyle=':')
            ax2.grid(which='major',linestyle=':')
            ax1.set_xlim([185,225])
            ax1.set_ylim([0,2])
            ax2.set_xlim([185,225])
            ax2.set_ylim([0,2])
            
            if i == 1:
                ax2.legend(loc=1,bbox_to_anchor=(2.25, 0.9))
            
            let = ["a","b","c"]
            ax1.set_title(let[i],loc="left",weight="bold")
            let = ["d","e"]
            ax2.set_title(let[i],loc="left",weight="bold")
        
        else:
            ax1.set_xlabel('Temperature (K)')
            ax1.set_title("clear-sky "+inst)
            ax1.grid(which='major',linestyle=':')
            ax1.set_xlim([185,225])
            ax1.set_ylim([0,2])
            
            let = ["a","b","c"]
            ax1.set_title(let[i],loc="left",weight="bold")
        
    plt.colorbar(m[3], cax=cax, ticks=[0,0.1,0.3,0.6,1.0], label="Frequency")
    plt.rcParams.update({"font.size":25})
    plt.savefig("./Paper-Figures/fig7-rhi-freq.png",bbox_inches="tight",dpi=200)
    plt.show()
