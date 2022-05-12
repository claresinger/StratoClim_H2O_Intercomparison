import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.stats as stats

flno = [2,3,4,6,7,8]
colors = np.array(["k","#045275","#0C7BDC","#7CCBA2","k","#FED976","#F0746E","#7C1D6F"])
maxlag = [0,0,5,10,10,20]
cmap = 'YlGnBu'
    
def h2o_pt_by_pt_whist_oor(dat):
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
        if a < 2:
            axin = ax.inset_axes([2,7,3,3], transform=ax.transData)
            axin.yaxis.set_label_position("right")
            axin.yaxis.tick_right()
            axin.plot([2,100],[2,100],"k-")
        
        # plot diagonal lines
        ax.plot([2,10],[2,10],"k-")
        if a > -1:
            ax.plot([0,12],[0,12*1.1],"k--")
            ax.plot([0,12],[0,12*0.9],"k--")
            ax.plot([0,12],[0,12*1.2],"k:")
            ax.plot([0,12],[0,12*0.8],"k:")
    
        # regression
        for i,h2ocut in enumerate([100,10]):
            datx = dat[(dat['ASCENT_FLAG'] == 0) & (dat['FLH2O'] <= h2ocut)]
            if a == 0:
                dat1 = datx[(datx['CLOUDY'] == 0) & (datx['CELL_GOOD'] == 1) & (datx['H2O'] <= h2ocut)]
                x = dat1['H2O']
                y = dat1['FLH2O']
                title = "a"
            if a == 1:
                dat1 = datx[(datx['CLOUDY'] == 0) & (datx['CELL_LOW'] == 1) & (datx['H2O'] <= h2ocut)]
                x = dat1['H2O']
                y = dat1['FLH2O']
                title = "b"
            if a == 2:
                dat1 = datx[(datx['CLOUDY'] == 0) & (datx['CELL_GOOD'] == 1) & (datx['H2O'] <= h2ocut)]
                x = dat1['H2O']
                y = dat1['FLH2O']
                title = "c"
            if a == 3:
                dat1 = datx[(datx['CLOUDY'] == 0) & (datx['CELL_LOW'] == 1) & (datx['H2O'] <= h2ocut)]
                x = dat1['H2O']
                y = dat1['FLH2O']
                title = "d"
                
            mask = ~np.isnan(x) & ~np.isnan(y)
            slope, intercept, rvalue, pvalue, se = stats.linregress(x[mask],y[mask])
            
            bias = (x[mask] - y[mask]) / y[mask] * 100.0
            meanbias = np.mean(bias)
            stdbias = np.std(bias)
            
            if a < 2:
                print(title)
                print(h2ocut, a, "r2=",rvalue**2)
                print("mean bias = ", meanbias, "%")
                print("std bias = ", stdbias, "%")
        
        ax.set_title(title,weight="bold",loc="left")
        if a < 2:
            ax.set_title("diff={:.1f}%,  $r^2=${:.3f}".format(meanbias, rvalue**2),loc="right",fontsize=20)
        else:
            ax.text(8.2,2.2,"{:.1f} hrs".format(len(x[mask]) / 3600))
    
        # plot
        if a == 0:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0) & (dat['CELL_GOOD'] == 1)]
            x = np.array(dat1['H2O'])
            y = np.array(dat1['FLH2O'])
            fi = np.array(dat1['FLIGHT'])
            p = np.random.permutation(len(x))
            x, y, fi = x[p], y[p], fi[p]
            ylabel = r"clear-sky FLASH H$_2$O"
            ax.text(3, 11.2, "cell pressure$\geq 30$mbar")
            
        if a == 1:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0) & (dat['CELL_LOW'] == 1)]
            x = np.array(dat1['H2O'])
            y = np.array(dat1['FLH2O'])
            fi = np.array(dat1['FLIGHT'])
            p = np.random.permutation(len(x))
            x, y, fi = x[p], y[p], fi[p]
            ax.text(2.5, 11.2, "$20 \leq$cell pressure$\leq 30$mbar")
            
        if a < 2:
            ax.scatter(x,y,20,c=colors[fi-1])
            axin.scatter(x,y,5,c=colors[fi-1])
        
        if a == 2:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0) & (dat['CELL_GOOD'] == 1)]
            x = dat1['H2O']
            y = dat1['FLH2O']
            
            vmin, vmax = 1, 100
            bins = [80,80]
            r = [[2,10],[2,10]]
            cmin = 1e-5
            m = ax.hist2d(x,y,bins=bins,range=r, 
                          cmap=cmap,norm=mcolors.PowerNorm(gamma=0.3),
                          vmin=vmin,vmax=vmax,cmin=cmin)
            xlabel = r"clear-sky ChiWIS H$_2$O"
            ylabel = r"clear-sky FLASH H$_2$O"
                    
        if a == 3:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0) & (dat['CELL_LOW'] == 1)]
            x = dat1['H2O']
            y = dat1['FLH2O']
            
            m = ax.hist2d(x,y,bins=bins,range=r, 
                          cmap=cmap,norm=mcolors.PowerNorm(gamma=0.3),
                          vmin=vmin,vmax=vmax,cmin=cmin)
            xlabel = r"clear-sky ChiWIS H$_2$O"  
            plt.colorbar(m[3], ax=ax, ticks=[vmin, 3, 30, vmax], label="counts")
                          
        if a == 0:
            for fi in flno:
                ax.scatter([-1],[-1],20,c=colors[fi-1], label="F"+str(fi))

        ax.set_xlim([2,10])
        ax.set_ylim([2,10])
        ax.grid()
    
        if a < 3:
            axin.set_xticks([25,50,75]); axin.set_yticks([25,50,75])
            axin.set_xlim(2,100), axin.set_ylim([2,100])
            axin.grid(which='both',linestyle=':')
    
    axused[0].legend(loc=4, ncol=3, frameon=True,
                      labelspacing=0.1, handletextpad=0.1, columnspacing=0.1,
                      borderpad = 0.2, borderaxespad = 0.4,
                      markerscale=2.0, fontsize=20, title_fontsize=20)
    
    fig.text(0.48, -0.05, r"clear-sky ChiWIS H$_2$O (ppmv)", ha='center')
    fig.text(-0.05, 0.5, r"clear-sky FLASH H$_2$O (ppmv)", va='center', rotation='vertical')
    
    plt.rcParams.update({"font.size":22})
    plt.savefig("./Paper-Figures/supp-scatter-h2o-hist-oor.png",dpi=300,bbox_inches="tight")
    plt.show()
