import matplotlib.pyplot as plt
import numpy as np

flno = [2,3,4,6,7,8]
colors = ["k","#045275","#0C7BDC","#7CCBA2","k","#FED976","#F0746E","#7C1D6F"]

def frac_error_pdfs(dat):
    # add cloudy flag
    dat['CLOUDY'] = ((dat['NICE'] > 0) | (dat['MASBR'] >= 1.2)).astype(int)
    
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

    lags = [0,5,10,20]
    cols = len(lags) + 1
    fig,axes = plt.subplots(len(flno),cols,figsize=(5*cols,4*len(flno)),sharey=True,sharex=False)
    plt.rcParams.update({"font.size":25})

    for i,f in enumerate(flno):
        dati = dat[(dat['FLIGHT'] == f) & (dat['ASCENT_FLAG'] == 0) & (dat['F7_DIVE'] == 0) & (dat['FLH2O'] <= 100.0)]
        
        for j in np.arange(cols):
            ax = axes[i,j]
            ax.axvline(color="k",linestyle="--")
            
            if j == 0:
                dat1 = dati
                e1 = ((dat1['FIH2O'] - dat1['FLH2O']) / dat1['FLH2O'] * 100).dropna()
                dat2 = dati[(dati['CELL_FLAG'] == 0)]
                e2 = ((dat2['H2O'] - dat2['FLH2O']) / dat2['FLH2O'] * 100).dropna()
            
            if j > 0:
                MAXLAG = lags[j-1]
                for lag in np.arange(1,MAXLAG):
                    dati['CLOUDY'] = np.maximum(dati['CLOUDY'], dati.shift(periods=lag, fill_value=0.0)['CLOUDY'])

                dat1 = dati[(dati['CLOUDY'] == 0)]
                e1 = ((dat1['FIH2O'] - dat1['FLH2O']) / dat1['FLH2O'] * 100).dropna()
                dat2 = dati[(dati['CLOUDY'] == 0) & (dati['CELL_FLAG'] == 0)]
                e2 = ((dat2['H2O'] - dat2['FLH2O']) / dat2['FLH2O'] * 100).dropna()
            
            if j == cols - 1:
                bins = np.linspace(-50,62.5,10)
                binlabels = [-50,0,50]
                ax.set_xlim(-65,65)
            else:
                bins = np.linspace(-50,212.5,22)
                binlabels = [-50,0,50,100,150,200]
                ax.set_xlim(-65,215)
            
            ax.hist(e2, color="k", bins=bins, alpha=0.1, histtype="bar", 
                    log=True, density=True, label="F"+str(f)+", ChiWIS")
            ax.hist(e2, color="k", bins=bins, linewidth=2, histtype="step", density=True)
            
            ax.hist(e1, color=colors[f-1], bins=bins, alpha=0.4, histtype="bar", 
                    log=True, density=True, label="F"+str(f)+", FISH")
            ax.hist(e1, color=colors[f-1], bins=bins, linewidth=2, histtype="step", density=True)
            
            labels = ["{:d}".format(x) for x in binlabels]
            ax.set_xscale('linear')
            ax.set_xticks(binlabels)
            ax.set_xticklabels(labels)
            ax.set_ylim(1e-5,5e-1)

            if i == 0:
                if j == 0:
                    ax.set_title("all-sky")
                if j > 0:
                    ax.set_title("clear-sky, $\\tau_c=${:d}".format(MAXLAG))
            if j == 0:
                ax.legend(loc=1,frameon=False,borderpad=0.2,borderaxespad=0.2,handletextpad=0.4,handlelength=1.0)
                
    fig.text(0.5, -0.02, "% diff from FLASH", ha='center', fontsize=35)
    fig.text(-0.02, 0.5, "Normalized PDF", va='center', rotation='vertical', fontsize=35)

    plt.rcParams.update({"font.size":25})
    plt.tight_layout()
    plt.savefig('Paper-Figures/supp-frac-error-pdfs.png',dpi=200,bbox_inches="tight")
    plt.show()
