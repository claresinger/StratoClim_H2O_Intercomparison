import numpy as np
import matplotlib.pyplot as plt

flno = [2,3,4,6,7,8]
maxlag = [0,0,5,10,10,20]
c = ["#009ADE","#FF1F58","k","green","orange"]

def make_means(dat,pt):
    chiA, chiB = np.zeros((len(pt),3)), np.zeros((len(pt),3))
    flaA, flaB = np.zeros((len(pt),3)), np.zeros((len(pt),3))
    fishA, fishB = np.zeros((len(pt),3)), np.zeros((len(pt),3))
    chiA[:,0], flaA[:,0], fishA[:,0] = pt, pt, pt
    chiB[:,0], flaB[:,0], fishB[:,0] = pt, pt, pt
    
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
    
    for i,pti in enumerate(pt):
        datA = dat[(dat['ASCENT_FLAG'] == 0) & (dat['PT'] >= pti-2.0) & (dat['PT'] < pti+2.0) & (dat['FLIGHT'] < 5)]
        dat_chiwisA = datA[(datA['CELL_FLAG'] == 0)]
        dat_flashA = datA[(datA['F7_DIVE'] == 0)]
        dat_clr_fishA = datA[(datA['CLOUDY'] == 0)]
        chiA[i,1], chiA[i,2] = np.mean(dat_chiwisA['H2O']), np.std(dat_chiwisA['H2O'])
        flaA[i,1], flaA[i,2] = np.mean(dat_flashA['FLH2O']), np.std(dat_flashA['FLH2O'])
        fishA[i,1], fishA[i,2] = np.mean(dat_clr_fishA['FIH2O']), np.std(dat_clr_fishA['FIH2O'])
        
        datB = dat[(dat['ASCENT_FLAG'] == 0) & (dat['PT'] >= pti-2.0) & (dat['PT'] < pti+2.0) & (dat['FLIGHT'] > 5)]
        dat_chiwisB = datB[(datB['CELL_FLAG'] == 0)]
        dat_flashB = datB[(datB['F7_DIVE'] == 0)]
        dat_clr_fishB = datB[(datB['CLOUDY'] == 0)]
        chiB[i,1], chiB[i,2] = np.mean(dat_chiwisB['H2O']), np.std(dat_chiwisB['H2O'])
        flaB[i,1], flaB[i,2] = np.mean(dat_flashB['FLH2O']), np.std(dat_flashB['FLH2O'])
        fishB[i,1], fishB[i,2] = np.mean(dat_clr_fishB['FIH2O']), np.std(dat_clr_fishB['FIH2O'])

    return chiA, chiB, flaA, flaB, fishA, fishB

def plot_profs_v4v5(ax, chi, fla, fish, mls_v5, mls_v4, byesno=False, b=False):
    ax.plot(chi[:,1], chi[:,0], color=c[2], label="ChiWIS")
    ax.fill_betweenx(chi[:,0], chi[:,1]-chi[:,2], chi[:,1]+chi[:,2], color=c[2], alpha=0.2)
    ax.plot(fla[:,1], fla[:,0], color=c[0], label="FLASH")
    ax.fill_betweenx(fla[:,0], fla[:,1]-fla[:,2], fla[:,1]+fla[:,2], color=c[0], alpha=0.2)
    if byesno != False:
        ax.plot(b[:,1], b[:,0], color=c[3], label="balloon CFH")
        ax.fill_betweenx(b[:,0], b[:,1]-b[:,2], b[:,1]+b[:,2], color=c[3], alpha=0.2)
    ax.plot(mls_v4[:,1], mls_v4[:,0], color=c[4], label="MLS v4")
    ax.fill_betweenx(mls_v4[:,0], mls_v4[:,1]-mls_v4[:,2], mls_v4[:,1]+mls_v4[:,2], color=c[4], alpha=0.2)
    ax.plot(mls_v5[:,1], mls_v5[:,0], color=c[4], linestyle=":", lw=2, label="MLS v5")
    ax.fill_betweenx(mls_v5[:,0], mls_v5[:,1]-mls_v5[:,2], mls_v5[:,1]+mls_v5[:,2], color=c[4], alpha=0.2)
    
    ax.plot([1,100],[382,382],"k--")
    ax.plot([1,100],[405,405],"k:")
    ax.set_ylim([370,480])
    ax.set_xlim([2.5,14])

    axin = ax.inset_axes([8,420,6,60], transform=ax.transData)
    axin.set_xscale("log", nonposx='clip')
    axin.plot(chi[:,1], chi[:,0], color=c[2])
    axin.plot(fla[:,1], fla[:,0], color=c[0])
    if byesno!= False:
        axin.plot(b[:,1], b[:,0], color=c[3])
    axin.plot(mls_v4[:,1], mls_v4[:,0], color=c[4])
    axin.plot(mls_v5[:,1], mls_v5[:,0], color=c[4], linestyle=":", lw=2)
    axin.plot([1,300],[382,382],"k--")
    axin.plot([1,300],[405,405],"k:")

    axin.set_ylim([360,480])
    ytk = [360,380,400,420,440,460]
    axin.set_yticklabels(list(map(str,ytk)))
    axin.set_xlim([2.5,100])
    xtk = [4,6,10,20,50]
    axin.set_xticks(xtk)
    axin.set_xticklabels(list(map(str,xtk)))
    axin.grid(linestyle=':')
    
    return ax

def mean_profile_compare(dat):
    pt = np.arange(362,502,4)
    chiA, chiB, flaA, flaB, fishA, fishB = make_means(dat,pt)

    balloon = np.loadtxt("Data/balloon_avg_prof.csv",delimiter=',',skiprows=1)
    mlsA = np.loadtxt("Data/mls_perA_prof.csv",delimiter=',',skiprows=1)
    mlsB = np.loadtxt("Data/mls_perB_prof.csv",delimiter=',',skiprows=1)
    mlsA_v4 = np.loadtxt("Data/mls_perA_prof_v4.csv",delimiter=',',skiprows=1)
    mlsB_v4 = np.loadtxt("Data/mls_perB_prof_v4.csv",delimiter=',',skiprows=1)
    
    s = 3
    print("PT compared: ", mlsA[s:,0])
    print("v5 diff in warm/wet: ", np.nanmean((mlsA[s:,1] - mlsA_v4[s:,1]) / mlsA_v4[s:,1]) * 100)
    print("v5 diff in cold/dry: ", np.nanmean((mlsB[s:,1] - mlsB_v4[s:,1]) / mlsB_v4[s:,1]) * 100)
    
    # bias calc
    print()
    print("% bias from MLS v5")
    
    i,j=3,12
    x, y = mlsA, chiA
    print("PT compared: ", y[i:j,0])
    print("chiwis warm/wet: {:.1f} \pm {:.1f}".format(
        np.nanmean((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0),
        np.nanstd((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0)
    ))
    y = flaA
    print("flash warm/wet: {:.1f} \pm {:.1f}".format(
        np.nanmean((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0),
        np.nanstd((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0)
    ))
    
    x, y = mlsB, chiB
    print("chiwis cold/dry: {:.1f} \pm {:.1f}".format(
        np.nanmean((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0),
        np.nanstd((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0)
    ))
    y = flaB
    print("flash cold/dry: {:.1f} \pm {:.1f}".format(
        np.nanmean((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0),
        np.nanstd((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0)
    ))
    y = balloon
    print("balloon cold/dry: {:.1f} \pm {:.1f}".format(
        np.nanmean((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0),
        np.nanstd((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0)
    ))
    
    # bias calc
    print()
    print("% bias from MLS v4")
    
    i,j=3,12
    x, y = mlsA_v4, chiA
    print("PT compared: ", y[i:j,0])
    print("chiwis warm/wet: {:.1f} \pm {:.1f}".format(
        np.nanmean((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0),
        np.nanstd((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0)
    ))
    y = flaA
    print("flash warm/wet: {:.1f} \pm {:.1f}".format(
        np.nanmean((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0),
        np.nanstd((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0)
    ))
    x, y = mlsB_v4, chiB
    print("chiwis cold/dry: {:.1f} \pm {:.1f}".format(
        np.nanmean((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0),
        np.nanstd((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0)
    ))
    y = flaB
    print("flash cold/dry: {:.1f} \pm {:.1f}".format(
        np.nanmean((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0),
        np.nanstd((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0)
    ))
    y = balloon
    print("balloon cold/dry: {:.1f} \pm {:.1f}".format(
        np.nanmean((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0),
        np.nanstd((np.interp(y[i:j,0], x[:,0], x[:,1])-y[i:j,1])/y[i:j,1]*100.0)
    ))

    # plot
    plt.rcParams.update({"font.size":15})
    fig,axes = plt.subplots(1,2,figsize=(12,8),sharey=True)
    # plt.rcParams.update({'font.size': 15})
    # plt.rcParams.update({'font.size': 15})
    
    ax0 = axes[0]
    ax0 = plot_profs_v4v5(ax0, chiA, flaA, fishA, mlsA, mlsA_v4, byesno=False, b=False)
    ax0.text(10.4,378.2,"approx. CPT")
    ax0.text(9.5,397,"max. altitude\nconv. influence")
    ax0.set_title("a",loc="left",weight="bold")
    ax0.set_title("F2-F4, warm/wet period")

    ax1 = axes[1]
    ax1 = plot_profs_v4v5(ax1, chiB, flaB, fishB, mlsB, mlsB_v4, byesno=True, b=balloon)
    ax1.set_title("b",loc="left",weight="bold")
    ax1.set_title("F6-F8, cold/dry period")

    # plot annotations
    ax0.set_xlabel(r'H$_2$O (ppmv)')#, fontsize=15)
    ax1.set_xlabel(r'H$_2$O (ppmv)')#, fontsize=15)
    ax0.set_ylabel('Potential Temperature (K)')#, fontsize=15)
    ax0.grid(which='both',linestyle=':')
    ax1.grid(which='both',linestyle=':')
    plt.rcParams.update({'font.size': 15})
    lgd = ax1.legend(loc=1,bbox_to_anchor=(1.6, 1))
    plt.savefig("Paper-Figures/fig9-mean-profiles-v4v5.png",bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)
    plt.show()
