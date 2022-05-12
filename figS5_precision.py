import numpy as np
import matplotlib.pyplot as plt
import datetime

c = ["#009ADE","#FF1F58","k"]

def moving_average(x, w):
    xp = np.convolve(x, np.ones(w), 'valid') / w
    return xp

def precision_f4(dat):
    datx = dat[(dat['FLIGHT'] == 4) & (dat['TIME']/3600.+5.75 > 16.0) & (dat['TIME']/3600.+5.75 < 16.14)]
    time = datx['TIME']/3600.+5.75
    
    plt.rcParams.update({"font.size":15})
    plt.figure(figsize=(10,4))
    ax = plt.gca()

    ax.plot(time, datx['FLH2O'],'.', color=c[0], label="FLASH")
    ax.plot(time, datx['FIH2O'],'.', color=c[1], label="FISH")
    ax.plot(time, datx['H2O'],'.', color=c[2], label="ChiWIS")
    
    ax.errorbar([16.146], np.nanmean(datx['FLH2O']), yerr=np.nanstd(datx['FLH2O']), 
                color=c[0], capsize=5)
    ax.text(16.148, np.nanmean(datx['FLH2O'])+0.05, "{:.1f}".format(np.nanstd(datx['FLH2O'])), color=c[0])
    
    ax.errorbar([16.146], np.nanmean(datx['FIH2O']), yerr=np.nanstd(datx['FIH2O']), 
                color=c[1], capsize=5)
    ax.text(16.148, np.nanmean(datx['FIH2O'])+0.05, "{:.1f}".format(np.nanstd(datx['FIH2O'])), color=c[1])
    
    ax.errorbar([16.146], np.nanmean(datx['H2O']), yerr=np.nanstd(datx['H2O']), 
                color=c[2], capsize=5)
    ax.text(16.148, np.nanmean(datx['H2O'])+0.05, "{:.2f}".format(np.nanstd(datx['H2O'])), color=c[2])
    
    ax.plot(16.151, 0, 'o')
    
    xlim = ax.get_xlim()
    ax.plot(xlim, np.ones(2)*np.nanmean(datx['FLH2O']), '-', color=c[0], alpha=0.5)
    ax.plot(xlim, np.ones(2)*np.nanmean(datx['FIH2O']), '-', color=c[1], alpha=0.5)
    ax.plot(xlim, np.ones(2)*np.nanmean(datx['H2O']), '-', color=c[2], alpha=0.5)
    ax.set_xlim(xlim)
    
    locs = ax.get_xticks()        
    labels = [str(datetime.timedelta(hours=x)).rsplit(':',1)[0] for x in locs]
    ax.set_xticklabels(labels)
    
    plt.xlabel("Kathmandu Local Time")
    plt.ylabel("H$_2$O (ppmv)")
    plt.title("Constant H$_2$O Segment of Flight 4")
    
    plt.ylim([4.75,6.75])
    plt.legend(loc=2,ncol=3,markerscale=3,handletextpad=0.2,columnspacing=0.5)

    plt.grid(axis='x')
    plt.savefig('./Paper-Figures/supp-precision-f4.png',dpi=300,bbox_inches='tight')
    plt.show()
