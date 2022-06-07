import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.stats as stats

flno = [2,3,4,6,7,8]
colors = np.array(["k","#045275","#0C7BDC","#7CCBA2","k","#FED976","#F0746E","#7C1D6F"])
maxlag = [0,0,5,10,10,20]
cmap = 'YlGnBu'
    
def h2o_pt_by_pt_whist6(dat):
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
        if a < 3:
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
                dat1 = datx[(datx['CLOUDY'] == 0) & (datx['FIH2O'] <= h2ocut)]
                x = dat1['FIH2O']
                y = dat1['FLH2O']
                title = "a"
            if a == 1:
                dat1 = datx[(datx['CLOUDY'] == 0) & (datx['CELL_FLAG'] == 0) & (datx['H2O'] <= h2ocut)]
                x = dat1['H2O']
                y = dat1['FLH2O']
                title = "b"
            if a == 2:
                dat1 = datx[(datx['CLOUDY'] == 1) & (datx['CELL_FLAG'] == 0) & (datx['F7_DIVE'] == 0) & (datx['H2O'] <= h2ocut)]
                x = dat1['H2O']
                y = dat1['FLH2O']
                title = "c"
            if a == 3:
                dat1 = datx[(datx['CLOUDY'] == 0) & (datx['FIH2O'] <= h2ocut)]
                x = dat1['FIH2O']
                y = dat1['FLH2O']
                title = "d"
            if a == 4:
                dat1 = datx[(datx['CLOUDY'] == 0) & (datx['CELL_FLAG'] == 0) & (datx['H2O'] <= h2ocut)]
                x = dat1['H2O']
                y = dat1['FLH2O']
                title = "e"
            if a == 5:
                dat1 = datx[(datx['CLOUDY'] == 1) & (datx['CELL_FLAG'] == 0) & (datx['F7_DIVE'] == 0) & (datx['H2O'] <= h2ocut)]
                x = dat1['H2O']
                y = dat1['FLH2O']
                title = "f"
                
            mask = ~np.isnan(x) & ~np.isnan(y)
            slope, intercept, rvalue, pvalue, se = stats.linregress(x[mask],y[mask])
            
            bias = (x[mask] - y[mask]) / y[mask] * 100.0
            absbias = np.abs(bias)
            meanbias = np.mean(bias)
            stdbias = np.std(bias)
            
            if a < 3:
                print(title)
                print(h2ocut, a, "r2=",rvalue**2)
                print("mean bias = ", meanbias, "%")
                print("std bias = ", stdbias, "%")
                if i == 1:
                    w = np.where(absbias <= 10.0)[0]
                    print(np.round(len(w)/len(absbias) * 100.0), "< 10% diff")
                    w = np.where(absbias <= 20.0)[0]
                    print(np.round(len(w)/len(absbias) * 100.0), "< 20% diff")
                print()
            
            if a == 1:
                dat1 = datx[(datx['CELL_FLAG'] == 0) & (datx['F7_DIVE'] == 0) & (datx['H2O'] <= h2ocut)]
                allx = dat1['H2O']
                ally = dat1['FLH2O']
                mask = ~np.isnan(allx) & ~np.isnan(ally)
                bias = (allx[mask] - ally[mask]) / ally[mask] * 100
                print("mean bias (all sky) = ", np.mean(bias), "%")
                print("std bias (all sky) = ", np.std(bias), "%")
                
        ax.set_title(title,weight="bold",loc="left")
        if a < 3:
            ax.set_title("diff={:.1f}%,  $r^2=${:.3f}".format(meanbias, rvalue**2),loc="right",fontsize=20)
        else:
            ax.text(8.2,2.2,"{:.1f} hrs".format(len(x[mask]) / 3600))

        # plot
        if a == 0:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0)]
            x = np.array(dat1['FIH2O'])
            y = np.array(dat1['FLH2O'])
            fi = np.array(dat1['FLIGHT'])
            p = np.random.permutation(len(x))
            x, y, fi = x[p], y[p], fi[p]
            ylabel = r"FLASH H$_2$O (ppmv)"

        if a == 1:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0) & (dat['CELL_FLAG'] == 0)]
            x = np.array(dat1['H2O'])
            y = np.array(dat1['FLH2O'])
            fi = np.array(dat1['FLIGHT'])
            p = np.random.permutation(len(x))
            x, y, fi = x[p], y[p], fi[p]
            ylabel = r"FLASH H$_2$O (ppmv)"

        if a == 2:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 1) & (dat['CELL_FLAG'] == 0)  & (dat['F7_DIVE'] == 0)]
            x = np.array(dat1['H2O'])
            y = np.array(dat1['FLH2O'])
            fi = np.array(dat1['FLIGHT'])
            p = np.random.permutation(len(x))
            x, y, fi = x[p], y[p], fi[p]
            
            dat3a = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 1) & (dat['CELL_FLAG'] == 0) & (dat['F7_DIVE'] == 1)]
            xa = dat3a['H2O'] 
            ya = dat3a['FLH2O']
            fia = np.array(dat3a['FLIGHT'])
            ylabel = r"FLASH H$_2$O (ppmv)"
            
        if a < 3:
            ax.scatter(x,y,20,c=colors[fi-1])
            axin.scatter(x,y,5,c=colors[fi-1])
        if a == 2:
            ax.scatter(xa,ya,50,facecolors='none',edgecolors=colors[fia-1])
            axin.scatter(xa,ya,10,facecolors='none',edgecolors=colors[fia-1])
        
        if a == 3:
            dat3 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0)]
            x = dat3['FIH2O']
            y = dat3['FLH2O']
            
            vmin, vmax = 1, 100
            bins = [80,80]
            r = [[2,10],[2,10]]
            cmin = 1e-5
            m = ax.hist2d(x,y,bins=bins,range=r, 
                          cmap=cmap,norm=mcolors.PowerNorm(gamma=0.3),
                          vmin=vmin,vmax=vmax,cmin=cmin)
            xlabel = r"FISH H$_2$O (ppmv)"
            ylabel = r"FLASH H$_2$O (ppmv)"
            
        if a == 4:
            dat3 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0) & (dat['CELL_FLAG'] == 0) & (dat['F7_DIVE'] == 0)]
            x = dat3['H2O']
            y = dat3['FLH2O']
            
            m = ax.hist2d(x,y,bins=bins,range=r, 
                          cmap=cmap,norm=mcolors.PowerNorm(gamma=0.3),
                          vmin=vmin,vmax=vmax,cmin=cmin)
            xlabel = r"ChiWIS H$_2$O (ppmv)"
            ylabel = r"FLASH H$_2$O (ppmv)"
            
        if a == 5:
            dat3 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 1) & (dat['CELL_FLAG'] == 0) & (dat['F7_DIVE'] == 0)]
            x = dat3['H2O']
            y = dat3['FLH2O']
            
            m = ax.hist2d(x,y,bins=bins,range=r, 
                          cmap=cmap,norm=mcolors.PowerNorm(gamma=0.3),
                          vmin=vmin,vmax=vmax,cmin=cmin)
            plt.colorbar(m[3], ax=ax, ticks=[vmin, 3, 30, vmax], label="counts")
            xlabel = r"ChiWIS H$_2$O (ppmv)"
            ylabel = r"FLASH H$_2$O (ppmv)"
              
        if a == 0:
            for fi in flno:
                ax.scatter([-1],[-1],20,c=colors[fi-1], label="F"+str(fi))
            
        if a > 2:
            ax.set_xlabel(xlabel)
        if a == 0 or a == 3:
            ax.set_ylabel(ylabel)
        ax.set_xlim([2,10])
        ax.set_ylim([2,10])
        ax.grid()
    
        if a < 3:
            axin.set_xticks([25,50,75]); axin.set_yticks([25,50,75])
            axin.set_xlim(2,100), axin.set_ylim([2,100])
            axin.grid(which='both',linestyle=':')
    
    plt.figtext(0.33,1.06,"Clear-sky", va="center", ha="center", size=25, weight="bold")
    plt.figtext(0.78,1.06,"In-cloud", va="center", ha="center", size=25, weight="bold")
    
    plt.figtext(0.175,1.02,"FISH vs. FLASH", va="center", ha="center", size=25, weight="bold")
    plt.figtext(0.48,1.02,"ChiWIS vs. FLASH", va="center", ha="center", size=25, weight="bold")
    plt.figtext(0.78,1.02,"ChiWIS vs. FLASH", va="center", ha="center", size=25, weight="bold")
    
    axused[0].legend(loc=4, ncol=3, frameon=True,
                      labelspacing=0.1, handletextpad=0.1, columnspacing=0.1, 
                      borderpad = 0.2, borderaxespad = 0.4,
                      markerscale=2.0, fontsize=20, title_fontsize=20)
    
    plt.savefig("./Paper-Figures/fig3-scatter-h2o-hist6.png",dpi=300,bbox_inches="tight")
    plt.show()

def key_figure(dat):
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
    
    fig,axes = plt.subplots(figsize=(20,5.5),ncols=3,nrows=1,constrained_layout=True)
    plt.rcParams.update({"font.size":22})
    
    axused = axes.flatten()
    for a,ax in enumerate(axused):
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
                dat1 = datx[(datx['CLOUDY'] == 0) & (datx['FIH2O'] <= h2ocut)]
                x = dat1['FIH2O']
                y = dat1['FLH2O']
                title = "a"
            if a == 1:
                dat1 = datx[(datx['CLOUDY'] == 0) & (datx['CELL_FLAG'] == 0) & (datx['H2O'] <= h2ocut)]
                x = dat1['H2O']
                y = dat1['FLH2O']
                title = "b"
            if a == 2:
                dat1 = datx[(datx['CLOUDY'] == 1) & (datx['CELL_FLAG'] == 0) & (datx['F7_DIVE'] == 0) & (datx['H2O'] <= h2ocut)]
                x = dat1['H2O']
                y = dat1['FLH2O']
                title = "c"

            mask = ~np.isnan(x) & ~np.isnan(y)
            slope, intercept, rvalue, pvalue, se = stats.linregress(x[mask],y[mask])
            
            bias = (x[mask] - y[mask]) / y[mask] * 100.0
            absbias = np.abs(bias)
            meanbias = np.mean(bias)
            stdbias = np.std(bias)

        ax.set_title(title,weight="bold",loc="left")
        ax.set_title("diff={:.1f}%,  $r^2=${:.3f}".format(meanbias, rvalue**2),loc="right",fontsize=20)
        ax.text(8.2,2.2,"{:.1f} hrs".format(len(x[mask]) / 3600))

        # plot
        if a == 0:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0)]
            x = np.array(dat1['FIH2O'])
            y = np.array(dat1['FLH2O'])
            fi = np.array(dat1['FLIGHT'])
            p = np.random.permutation(len(x))
            x, y, fi = x[p], y[p], fi[p]
            xlabel = r"FISH H$_2$O (ppmv)"
            ylabel = r"FLASH H$_2$O (ppmv)"

        if a == 1:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 0) & (dat['CELL_FLAG'] == 0)]
            x = np.array(dat1['H2O'])
            y = np.array(dat1['FLH2O'])
            fi = np.array(dat1['FLIGHT'])
            p = np.random.permutation(len(x))
            x, y, fi = x[p], y[p], fi[p]
            xlabel = r"ChiWIS H$_2$O (ppmv)"
            ylabel = r"FLASH H$_2$O (ppmv)"

        if a == 2:
            dat1 = dat[(dat['ASCENT_FLAG'] == 0) & (dat['CLOUDY'] == 1) & (dat['CELL_FLAG'] == 0)  & (dat['F7_DIVE'] == 0)]
            x = np.array(dat1['H2O'])
            y = np.array(dat1['FLH2O'])
            fi = np.array(dat1['FLIGHT'])
            p = np.random.permutation(len(x))
            x, y, fi = x[p], y[p], fi[p]
            xlabel = r"ChiWIS H$_2$O (ppmv)"
            ylabel = r"FLASH H$_2$O (ppmv)"

        ax.scatter(x,y,20,c=colors[fi-1])

        if a == 0:
            ax.set_ylabel(ylabel)
            for fi in flno:
                ax.scatter([-1],[-1],20,c=colors[fi-1], label="F"+str(fi))

        ax.set_xlabel(xlabel)  
        ax.set_xlim([2,10])
        ax.set_ylim([2,10])
        ax.grid()
    
    plt.figtext(0.35,1.05,"Clear-sky", va="center", ha="center", size=25, weight="bold")
    plt.figtext(0.83,1.05,"In-cloud", va="center", ha="center", size=25, weight="bold")
    
    axused[0].legend(loc=2, ncol=3, frameon=True,
                      labelspacing=0.1, handletextpad=0.1, columnspacing=0.1, 
                      borderpad = 0.2, borderaxespad = 0.4,
                      markerscale=2.0, fontsize=20, title_fontsize=20)
    
    plt.savefig("./Paper-Figures/key_figure.png",dpi=300,bbox_inches="tight")
    plt.show()