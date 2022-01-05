import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

flno = [2,3,4,6,7,8]
maxlag = [0,0,5,10,10,20]
cmap = "YlGnBu"
colors = np.array(["k","#045275","#0C7BDC","#7CCBA2","k","#FED976","#F0746E","#7C1D6F"])

def h2o_diff_profs_clr_cloud(dat):
    fig = plt.figure(figsize=(16,6))
    plt.rcParams.update({"font.size":20})
    gs0 = gridspec.GridSpec(1,4,figure=fig,width_ratios=[10,10,10,1])
    
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
    
    xrange=[-60,60]
    xbins = int((xrange[-1]-xrange[0])/2)
    yrange=[350,490]
    ybins = int((yrange[-1]-yrange[0])/2)
    
    ### clear-sky FISH
    ax = fig.add_subplot(gs0[0])
    ax0 = ax
    dati = dat[(dat['CLOUDY'] == 0) & (dat['ASCENT_FLAG'] == 0) & (dat['F7_DIVE'] == 0)]
    pt = dati["PT"]
    diff = (dati["FIH2O"] - dati["FLH2O"]) / dati["FLH2O"] * 100
    m = ax.hist2d(diff,pt,bins=[xbins,ybins],range=[xrange,yrange],
            norm=mcolors.PowerNorm(gamma=0.3),vmin=5,vmax=500,
            cmin=1e-10,cmap=cmap)
    ax.set_title("a",weight="bold",loc="left")
    ax.set_title("clear-sky FISH")
    ax.set_xlabel("% diff from FLASH")
    ax.set_ylabel("Potential Temperature (K)")
    ax.set_xticks(np.linspace(-60,60,7))
    ax.grid()
    ax.errorbar(np.mean(diff), 480, xerr=np.std(diff), fmt="ro", elinewidth=3, ecolor="r", capthick=3, capsize=5)
    ax.text(9.5,481,"{:.1f} $\pm$ {:.1f}%".format(np.mean(diff),np.std(diff)), color="r", fontsize=18)
    
    ### clear-sky ChiWIS
    ax = fig.add_subplot(gs0[1], sharex=ax0, sharey=ax0)
    plt.setp(ax.get_yticklabels(), visible=False)
    dati = dat[(dat['CLOUDY'] == 0) & (dat['CELL_FLAG'] == 0) & (dat['ASCENT_FLAG'] == 0) & (dat['F7_DIVE'] == 0)]
    pt = dati["PT"]
    diff = (dati["H2O"] - dati["FLH2O"]) / dati["FLH2O"] * 100
    m = ax.hist2d(diff,pt,bins=[xbins,ybins],range=[xrange,yrange],
            norm=mcolors.PowerNorm(gamma=0.3),vmin=5,vmax=500,
            cmin=1e-10,cmap=cmap)
    ax.set_title("b",weight="bold",loc="left")
    ax.set_title("clear-sky ChiWIS")
    ax.set_xlabel("% diff from FLASH")
    ax.grid()
    ax.errorbar(np.mean(diff), 480, xerr=np.std(diff), fmt="ro", elinewidth=3, ecolor="r", capthick=3, capsize=5)
    ax.text(9.5,481,"{:.1f} $\pm$ {:.1f}%".format(np.mean(diff),np.std(diff)), color="r", fontsize=18)
    
    ### in-cloud ChiWIS
    ax = fig.add_subplot(gs0[2], sharex=ax0, sharey=ax0)
    plt.setp(ax.get_yticklabels(), visible=False)
    dati = dat[(dat['CLOUDY'] == 1) & (dat['CELL_FLAG'] == 0) & (dat['ASCENT_FLAG'] == 0) & (dat['F7_DIVE'] == 0)]
    pt = dati["PT"]
    diff = (dati["H2O"] - dati["FLH2O"]) / dati["FLH2O"] * 100
    m = ax.hist2d(diff,pt,bins=[xbins,ybins],range=[xrange,yrange],
            norm=mcolors.PowerNorm(gamma=0.3),vmin=5,vmax=500,
            cmin=1e-10,cmap=cmap)
    ax.set_title("c",weight="bold",loc="left")
    ax.set_title("in-cloud ChiWIS")
    ax.set_xlabel("% diff from FLASH")
    ax.grid()
    ax.errorbar(np.mean(diff), 480, xerr=np.std(diff), fmt="ro", elinewidth=3, ecolor="r", capthick=3, capsize=5)
    ax.text(9.5,481,"{:.1f} $\pm$ {:.1f}%".format(np.mean(diff),np.std(diff)), color="r", fontsize=18)

    cax = fig.add_subplot(gs0[-1])
    plt.colorbar(m[3], cax=cax, ticks=np.geomspace(5,500,3), label="counts")
    plt.tight_layout()
    plt.savefig("Paper-Figures/fig4-diff-profiles.png",dpi=200)
    plt.show()