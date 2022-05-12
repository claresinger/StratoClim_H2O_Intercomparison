import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.stats as stats

flno = [2,3,4,6,7,8]
colors = np.array(["k","#045275","#0C7BDC","#7CCBA2","k","#FED976","#F0746E","#7C1D6F"])
maxlag = [0,0,5,10,10,20]
cmap = 'YlGnBu'

def rhi_pt_by_pt_whist6(dat):
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
    
    fig,axes = plt.subplots(figsize=(20,10),ncols=3,nrows=2,constrained_layout=True)
    plt.rcParams.update({"font.size":22})
    
    axused = axes.flatten()
    for a,ax in enumerate(axused):
        ax.plot([0,2],[0,2],"k-")
    
        # regression
        datx = dat[(dat['ASCENT_FLAG'] == 0)]
        if a == 0:
            dat1 = datx[(datx['CLOUDY'] == 0)]
            x = dat1['FIH2O'] / dat1['SATPPM']
            y = dat1['FLH2O'] / dat1['SATPPM']
            title = "a"
        if a == 1:
            dat1 = datx[(datx['CLOUDY'] == 0) & (datx['CELL_FLAG'] == 0)]
            x = dat1['H2O'] / dat1['SATPPM']
            y = dat1['FLH2O'] / dat1['SATPPM']
            title = "b"
        if a == 2:
            dat1 = datx[(datx['CLOUDY'] == 1) & (datx['CELL_FLAG'] == 0) & (datx['F7_DIVE'] == 0)]
            x = dat1['H2O'] / dat1['SATPPM']
            y = dat1['FLH2O'] / dat1['SATPPM']
            title = "c"
        if a == 3:
            dat1 = datx[(datx['CLOUDY'] == 0)]
            x = dat1['FIH2O'] / dat1['SATPPM']
            y = dat1['FLH2O'] / dat1['SATPPM']
            title = "d"
        if a == 4:
            dat1 = datx[(datx['CLOUDY'] == 0) & (datx['CELL_FLAG'] == 0)]
            x = dat1['H2O'] / dat1['SATPPM']
            y = dat1['FLH2O'] / dat1['SATPPM']
            title = "e"
        if a == 5:
            dat1 = datx[(datx['CLOUDY'] == 1) & (datx['CELL_FLAG'] == 0) & (datx['F7_DIVE'] == 0)]
            x = dat1['H2O'] / dat1['SATPPM']
            y = dat1['FLH2O'] / dat1['SATPPM']
            title = "f"

        mask = ~np.isnan(x) & ~np.isnan(y)
        slope, intercept, rvalue, pvalue, se = stats.linregress(x[mask],y[mask])

        bias = (x[mask] - y[mask]) / y[mask] * 100.0
        meanbias = np.mean(bias)
        stdbias = np.std(bias)
        
        if a < 3:
            print(a, "r2=",rvalue**2)
            print("mean bias = ", meanbias, "%")
            print("std bias = ", stdbias, "%")
            print()

        ax.set_title(title,weight="bold",loc="left")
        if a < 3:
            ax.set_title("diff={:.1f}%,  $r^2=${:.3f}".format(meanbias, rvalue**2),loc="right",fontsize=20)
        else:
            ax.text(1.5,0.05,"{:.1f} hrs".format(len(x[mask]) / 3600))

        # plot
        if a == 0:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0)]
            x = np.array(dat1['FIH2O'] / dat1['SATPPM'])
            y = np.array(dat1['FLH2O'] / dat1['SATPPM'])
            fi = np.array(dat1['FLIGHT'])
            p = np.random.permutation(len(x))
            x, y, fi = x[p], y[p], fi[p]
            ylabel = r"FLASH RH$_{ice}$"

        if a == 1:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0) & (dat['CELL_FLAG'] == 0)]
            x = np.array(dat1['H2O'] / dat1['SATPPM'])
            y = np.array(dat1['FLH2O'] / dat1['SATPPM'])
            fi = np.array(dat1['FLIGHT'])
            p = np.random.permutation(len(x))
            x, y, fi = x[p], y[p], fi[p]

        if a == 2:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 1) & (dat['CELL_FLAG'] == 0)  & (dat['F7_DIVE'] == 0)]
            x = np.array(dat1['H2O'] / dat1['SATPPM'])
            y = np.array(dat1['FLH2O'] / dat1['SATPPM'])
            fi = np.array(dat1['FLIGHT'])
            p = np.random.permutation(len(x))
            x, y, fi = x[p], y[p], fi[p]
            
            dat3a = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 1) & (dat['CELL_FLAG'] == 0) & (dat['F7_DIVE'] == 1)]
            xa = dat3a['H2O'] / dat3a['SATPPM']
            ya = dat3a['FLH2O'] / dat3a['SATPPM']
            fia = np.array(dat3a['FLIGHT'])
            
        if a < 3:
            ax.scatter(x,y,20,c=colors[fi-1])
        if a == 2:
            ax.scatter(xa,ya,50,facecolors='none',edgecolors=colors[fia-1])
        
        if a == 3:
            dat3 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0)]
            x = dat3['FIH2O'] / dat3['SATPPM']
            y = dat3['FLH2O'] / dat3['SATPPM']
            
            vmin, vmax = 1, 400
            m = ax.hist2d(x,y,bins=[60,60],range=[[0,2],[0,2]], 
                          cmap=cmap,norm=mcolors.PowerNorm(gamma=0.3),
                          vmin=vmin,vmax=vmax,cmin=1e-10)
            xlabel = r"FISH RH$_{ice}$"
            ylabel = r"FLASH RH$_{ice}$"
            
        if a == 4:
            dat3 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0) & (dat['CELL_FLAG'] == 0) & (dat['F7_DIVE'] == 0)]
            x = dat3['H2O'] / dat3['SATPPM']
            y = dat3['FLH2O'] / dat3['SATPPM']
            
            m = ax.hist2d(x,y,bins=[60,60],range=[[0,2],[0,2]], 
                          cmap=cmap,norm=mcolors.PowerNorm(gamma=0.3),
                          vmin=vmin,vmax=vmax,cmin=1e-10)
            xlabel = r"ChiWIS RH$_{ice}$"
            
        if a == 5:
            dat3 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 1) & (dat['CELL_FLAG'] == 0) & (dat['F7_DIVE'] == 0)]
            x = dat3['H2O'] / dat3['SATPPM']
            y = dat3['FLH2O'] / dat3['SATPPM']
            
            m = ax.hist2d(x,y,bins=[60,60],range=[[0,2],[0,2]], 
                          cmap=cmap,norm=mcolors.PowerNorm(gamma=0.3),
                          vmin=vmin,vmax=vmax,cmin=1e-10)
            plt.colorbar(m[3], ax=ax, ticks=[vmin, 10, 100, vmax], label="counts")
            xlabel = r"ChiWIS RH$_{ice}$"
              
        if a == 0:
            for fi in flno:
                ax.scatter([-1],[-1],20,c=colors[fi-1], label="F"+str(fi))
            
        if a > 2:
            ax.set_xlabel(xlabel)
        if a == 0 or a == 3:
            ax.set_ylabel(ylabel)
        ax.set_xlim([0,2])
        ax.set_ylim([0,2])
        ax.set_xticks([0,0.5,1,1.5,2])
        ax.set_yticks([0,0.5,1,1.5,2])
        ax.grid()
    
    axused[0].legend(loc=4, ncol=3, frameon=True,
                      labelspacing=0.1, handletextpad=0.1, columnspacing=0.1, 
                      borderpad = 0.2, borderaxespad = 0.4,
                      markerscale=2.0, fontsize=20, title_fontsize=20)
    
    plt.figtext(0.33,1.06,"Clear-sky", va="center", ha="center", size=25, weight="bold")
    plt.figtext(0.78,1.06,"In-cloud", va="center", ha="center", size=25, weight="bold")
    
    plt.figtext(0.175,1.02,"FISH vs. FLASH", va="center", ha="center", size=25, weight="bold")
    plt.figtext(0.48,1.02,"ChiWIS vs. FLASH", va="center", ha="center", size=25, weight="bold")
    plt.figtext(0.78,1.02,"ChiWIS vs. FLASH", va="center", ha="center", size=25, weight="bold")
    
    plt.rcParams.update({"font.size":22})
    plt.savefig("./Paper-Figures/fig8-scatter-rhi-hist6.png",dpi=300,bbox_inches="tight")
    plt.show()
