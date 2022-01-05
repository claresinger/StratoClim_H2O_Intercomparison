import numpy as np
import matplotlib.pyplot as plt
import datetime

def cloud_mask(dat):
    flno = [4,6,7,8]
    t1 = [14.5, 13.55, 10.55, 15.1]
    t2 = [17.85, 16.4, 12.2, 17.7]
    
    fig = plt.figure(figsize=(14,14))
    ax_com = fig.add_subplot(111)    # The big subplot
    ax = [fig.add_subplot(411), fig.add_subplot(412), fig.add_subplot(413), fig.add_subplot(414)]
    plt.rcParams.update({'font.size':20})
    
    # Turn off axis lines and ticks of the big subplot
    ax_com.spines['top'].set_color('none')
    ax_com.spines['bottom'].set_color('none')
    ax_com.spines['left'].set_color('none')
    ax_com.spines['right'].set_color('none')
    ax_com.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

    for i,f in enumerate(flno):
        dati = dat[(dat['FLIGHT'] == f) & (dat['ALT'] > 12.0)]
        xvar = 'TIME'; xfac = 1/3600.; xoff = 5.75

        ax[i].plot(dati[xvar]*xfac+xoff,dati['NICE'],'k.-')
        dati_mc = dati[(dati['MASBR'] >= 1.2)]
        dati_nc = dati[(dati['NICE'] > 0.0)]
        dati_ic = dati[(dati['IWC'] > 0.0)]
        #dati_particles = dat[(dat['FLIGHT'] == f) & (dat['ALT'] > 15.0) & (dat['PARTICLES'] == 1)]

        ax[i].set_yscale('log')
        ax[i].set_ylim([1e-3,1e2])
        ax[i].set_ylabel('N$_{ice}$ (cm$^{-3}$)')
        
        ax[i].set_xlim([t1[i],t2[i]])
        locs = ax[i].get_xticks()        
        labels = [str(datetime.timedelta(hours=x)).rsplit(':',1)[0] for x in locs]
        ax[i].set_xticklabels(labels)
        
        ax[i].set_title("Flight {0}".format(str(f)))
        ax[i].grid(which='major',linestyle=':')

        mas, = ax[i].plot(dati_mc[xvar]*xfac+xoff,np.ones(len(dati_mc))*(70),'r.')
        nice, = ax[i].plot(dati_nc[xvar]*xfac+xoff,np.ones(len(dati_nc))*(45),'.',color="#ebc415")
        iwc, = ax[i].plot(dati_ic[xvar]*xfac+xoff,np.ones(len(dati_ic))*(28), 'b.')
        
        ax2 = ax[i].twinx()
        dat_alt = dat[(dat['FLIGHT'] == f) & (dat['ALT'] > 10)]
        ax2.plot(dat_alt['TIME']/3600.+5.75,dat_alt['ALT'],"-",color="green",lw=2)
        ax2.set_ylabel("Altitude (km)",color="green")
        ax2.tick_params(axis='y', colors="green")
        ax2.set_yticks([10,12,14,16,18,20])

        if i == len(flno)-1:
            ax[i].set_xlabel('Kathmandu Local Time')

    lgd = ax[0].legend([mas,nice,iwc],['BR $\geq$ 1.2','N$_{ice}$ > 0','IWC > 0'], markerscale=3,
                       frameon=True, loc=(0.77,0.05), fontsize=15)
    plt.tight_layout()
    plt.savefig('Paper-Figures/supp-cldmask-nice.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)
    plt.show()