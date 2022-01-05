import matplotlib.pyplot as plt
import numpy as np

flno = [2,3,4,6,7,8]
colors = ["k","#045275","#0C7BDC","#7CCBA2","k","#FED976","#F0746E","#7C1D6F"]
maxlag = [0,0,0,10,10,20]
    
def dry_errors_2rows(dat):
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
    
    fig,axes = plt.subplots(2,len(flno),figsize=(len(flno)*5,10),sharey=True)
    plt.rcParams.update({"font.size":25})
    for i,f in enumerate(flno):
        ax = axes[0,i]
        ax.axvline(0, color="k",linestyle=":")

        H2OCUT = 10.0 # ppmv
        dati = dat[(dat['FLIGHT'] == f) & (dat['ASCENT_FLAG'] == 0) & (dat['F7_DIVE'] == 0) & (dat['FLH2O'] <= H2OCUT)]

        dat1 = dati[(dati['CLOUDY'] == 0)]
        e1 = (dat1['FIH2O'] - dat1['FLH2O']).dropna()
        dat2 = dati[(dati['CLOUDY'] == 0) & (dati['CELL_FLAG'] == 0)]
        e2 = (dat2['H2O'] - dat2['FLH2O']).dropna()
        bins = np.linspace(-2.5,2.5,21)
        
        print(len(e1), len(e2))
        print(np.std(e1), np.std(e2))
                    
        ax.hist(e2, color=colors[f-1], bins=bins, alpha=0.3, histtype="bar", 
                log=True, density=True, label="ChiWIS")
        ax.hist(e2, color=colors[f-1], bins=bins, linewidth=2, histtype="step", density=True)
        ax.axvline(np.median(e2), color="k", linestyle="-")
        ax.text(0.9,0.5,"median\n={:.2f}".format(np.median(e2)))
        ax.set_xlabel("ChiWIS$-$FLASH (ppmv)")
        ax.set_title("F"+str(f)+", $\\tau_c=$"+str(maxlag[i]))
        
        if i == 0:
            ax.set_ylabel("Normalized PDF")
            
        ax = axes[1,i]
        ax.axvline(0, color="k",linestyle=":")
        
        ax.hist(e1, color=colors[f-1], bins=bins, alpha=0.3, histtype="bar", 
                log=True, density=True, label="FISH")
        ax.hist(e1, color=colors[f-1], bins=bins, linewidth=2, histtype="step", density=True)
        ax.axvline(np.median(e1), color="k", linestyle="-")
        ax.text(0.9,0.5,"median\n={:.2f}".format(np.median(e1)))
        ax.set_ylim(1e-5,1e1)
        ax.set_xlabel("FISH$-$FLASH (ppmv)")

        if i == 0:
            ax.set_ylabel("Normalized PDF")
            
    fig.delaxes(axes[1,1])
                
    plt.rcParams.update({"font.size":25})
    plt.tight_layout()
    plt.savefig("Paper-Figures/fig3-diff-hist.png",dpi=200)
    plt.show()