import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.stats as stats

flno = [2,3,4,6,7,8]
colors = np.array(["k","#045275","#0C7BDC","#7CCBA2","k","#FED976","#F0746E","#7C1D6F"])
maxlag = [0,0,5,10,10,20]
cmap = 'YlGnBu'
    
def rhi_pt_by_pt_whist_oor(dat):
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
    dat['CELL_LOW'] = ((dat['PRES_CELL'] > 20.0) & (dat['PRES_CELL'] < 30.0) & (dat['FLAG'] == 0)).astype(int)
    
    # FL7 dive flag
    dat['F7_DIVE'] = ((dat['FLIGHT'] == 7) & (dat['TIME'] > 19.9e3) & (dat['TIME'] < 20.2e3)).astype('int')
    
    fig,axes = plt.subplots(figsize=(13,9),ncols=2,nrows=2,constrained_layout=True)
    plt.rcParams.update({"font.size":22})
    
    axused = axes.flatten()
    for a,ax in enumerate(axused):        
        ax.plot([0,2],[0,2],"k-")
    
        # regression
        datx = dat[(dat['ASCENT_FLAG'] == 0)]
        if a == 0:
            dat1 = datx[(datx['CLOUDY'] == 0) & (datx['CELL_GOOD'] == 1)]
            x = dat1['H2O'] / dat1['SATPPM']
            y = dat1['FLH2O'] / dat1['SATPPM']
            title = "a"
        if a == 1:
            dat1 = datx[(datx['CLOUDY'] == 0) & (datx['CELL_LOW'] == 1)]
            x = dat1['H2O'] / dat1['SATPPM']
            y = dat1['FLH2O'] / dat1['SATPPM']
            title = "b"
        if a == 2:
            dat1 = datx[(datx['CLOUDY'] == 0) & (datx['CELL_GOOD'] == 1)]
            x = dat1['H2O'] / dat1['SATPPM']
            y = dat1['FLH2O'] / dat1['SATPPM']
            title = "c"
        if a == 3:
            dat1 = datx[(datx['CLOUDY'] == 0) & (datx['CELL_LOW'] == 1)]
            x = dat1['H2O'] / dat1['SATPPM']
            y = dat1['FLH2O'] / dat1['SATPPM']
            title = "d"
                
        mask = ~np.isnan(x) & ~np.isnan(y)
        slope, intercept, rvalue, pvalue, se = stats.linregress(x[mask],y[mask])
            
        bias = (x[mask] - y[mask]) / y[mask] * 100.0
        meanbias = np.mean(bias)
        stdbias = np.std(bias)
            
        if a < 2:
            print(title)
            print(a, "r2=",rvalue**2)
            print("mean bias = ", meanbias, "%")
            print("std bias = ", stdbias, "%")
            print()
        
        ax.set_title(title,weight="bold",loc="left")
        if a < 2:
            ax.set_title("diff={:.1f}%,  $r^2=${:.3f}".format(meanbias, rvalue**2),loc="right",fontsize=20)
        else:
            ax.text(1.5,0.05,"{:.1f} hrs".format(len(x[mask]) / 3600))
    
        # plot
        if a == 0:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0) & (dat['CELL_GOOD'] == 1)]
            x = np.array(dat1['H2O'] / dat1['SATPPM'])
            y = np.array(dat1['FLH2O'] / dat1['SATPPM'])
            fi = np.array(dat1['FLIGHT'])
            p = np.random.permutation(len(x))
            x, y, fi = x[p], y[p], fi[p]
            ax.text(0.2, 2.35, "cell pressure$\geq 30$mbar")
            
        if a == 1:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0) & (dat['CELL_LOW'] == 1)]
            x = np.array(dat1['H2O'] / dat1['SATPPM'])
            y = np.array(dat1['FLH2O'] / dat1['SATPPM'])
            fi = np.array(dat1['FLIGHT'])
            p = np.random.permutation(len(x))
            x, y, fi = x[p], y[p], fi[p]
            ax.text(0.1, 2.35, "$20 \leq$cell pressure$\leq 30$mbar")
            
        if a < 2:
            ax.scatter(x,y,20,c=colors[fi-1])
        
        if a == 2:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0) & (dat['CELL_GOOD'] == 1)]
            x = dat1['H2O'] / dat1['SATPPM']
            y = dat1['FLH2O'] / dat1['SATPPM']
            
            vmin, vmax = 1, 400
            bins = [60,60]
            r = [[0,2],[0,2]]
            cmin = 1e-5
            m = ax.hist2d(x,y,bins=bins,range=r, 
                          cmap=cmap,norm=mcolors.PowerNorm(gamma=0.3),
                          vmin=vmin,vmax=vmax,cmin=cmin)
                    
        if a == 3:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0) & (dat['CELL_LOW'] == 1)]
            x = dat1['H2O'] / dat1['SATPPM']
            y = dat1['FLH2O'] / dat1['SATPPM']
            
            m = ax.hist2d(x,y,bins=bins,range=r, 
                          cmap=cmap,norm=mcolors.PowerNorm(gamma=0.3),
                          vmin=vmin,vmax=vmax,cmin=cmin)
            plt.colorbar(m[3], ax=ax, ticks=[vmin, 10, 100, vmax], label="counts")
                          
        if a == 0:
            for fi in flno:
                ax.scatter([-1],[-1],20,c=colors[fi-1], label="F"+str(fi))

        ax.set_xlim([0,2])
        ax.set_ylim([0,2])
        ax.grid()
    
    axused[0].legend(loc=4, ncol=3, frameon=True,
                      labelspacing=0.1, handletextpad=0.1, columnspacing=0.1, 
                      borderpad = 0.2, borderaxespad = 0.4,
                      markerscale=2.0, fontsize=20, title_fontsize=20)
    
    fig.text(0.48, -0.05, r"clear-sky ChiWIS RH$_{ice}$", ha='center')
    fig.text(-0.05, 0.5, r"clear-sky FLASH RH$_{ice}$", va='center', rotation='vertical')
    
    plt.rcParams.update({"font.size":22})
    plt.savefig("./Paper-Figures/supp-scatter-rhi-hist-oor.png",dpi=300,bbox_inches="tight")
    plt.show()
