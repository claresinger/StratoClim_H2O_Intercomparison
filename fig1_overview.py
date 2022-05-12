import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import seaborn
from mpl_toolkits.basemap import Basemap
import datetime

flno = [2,3,4,6,7,8]
date = ["27 Jul", "29 Jul", "31 Jul", "2 Aug", "4 Aug", "6 Aug", "8 Aug", "10 Aug"]
colors = ["k","#045275","#0C7BDC","#7CCBA2","k","#FED976","#F0746E","#7C1D6F"]

def map_alt_h2o(dat):
    plt.rcParams.update({'font.size': 15})
    fig = plt.figure(figsize=(11.5,8), constrained_layout=True)
    gs = fig.add_gridspec(6, 19)
    ax1 = fig.add_subplot(gs[0:-2, 0:-9])
    ax2 = fig.add_subplot(gs[-2:, 0:-9])
    ax3 = fig.add_subplot(gs[:,-9:])
    
    latmin, latmax = 18, 32
    lonmin, lonmax = 75, 95
    
    # setup basemap.
    m = Basemap(ax=ax1, projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax)    
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(fill_color='lightsteelblue')
    m.fillcontinents(color='white',lake_color='lightsteelblue')
    parallels = np.linspace(latmin,latmax,num=5,endpoint=True)
    m.drawparallels(parallels,labels=[True,False,True,False])
    meridians = np.linspace(lonmin,lonmax,num=5,endpoint=True)
    m.drawmeridians(meridians,labels=[True,False,False,True])
    
    # MLS flyovers
    coords = np.loadtxt("./Data/MLS4SClim_V5/profilesCoord.txt", skiprows=1)
    lon, lat = coords[:,0], coords[:,1]
    print(np.shape(lon))
    m.plot(lon,lat,'.',color='grey',latlon=True)

    # draw flight pathes on top
    for i in flno:
        dati = dat[(dat['FLIGHT'] == i)]
        lon = np.array(dati['LON'])
        lat = np.array(dati['LAT'])
        m.plot(lon,lat,'.',latlon=True,color=colors[i-1])
    
    # balloon launch location
    balloon_lon, balloon_lat = 85.5, 27.619
    m.plot(balloon_lon, balloon_lat, latlon=True, 
           marker="*", markersize=20, color="white", markeredgewidth=2, markeredgecolor="k")
    
    ax1.set_xlabel("Longitude",labelpad=25)
    ax1.set_ylabel("Latitutde",labelpad=55)
    ax1.set_title("a",loc="left",weight="bold")
    
    x, y = m(76, 24); ax1.text(x, y, 'India',fontsize=20,fontweight='bold')
    x, y = m(88, 30); ax1.text(x, y, 'China',fontsize=20,fontweight='bold')
    x, y = m(81.5, 29.5); ax1.text(x, y, 'Nepal',fontsize=20,fontweight='bold')
    x, y = m(86, 19); ax1.text(x, y, 'Bay of Bengal',fontsize=20,fontweight='bold')
      
    # plot altitude timeseries
    for i in flno:
        dati = dat[(dat['FLIGHT'] == i)]
        ax2.plot(dati['TIME']/3600.+5.75,dati['ALT'],'.',color=colors[i-1],label="F"+str(i)+", "+date[i-1])
    
    ax2.set_xlabel("Kathmandu Local Time")
    ax2.set_ylabel("Altitude (km)")
    ax2.set_ylim([1,22])
    ax2.set_yticks([5,10,15,20])
    
    ticks = [8,10,12,14,16,18]
    ticks_name = [str(datetime.timedelta(hours=x)).rsplit(':',1)[0] for x in ticks]
    ax2.set_xlim([8, 19]) 
    ax2.set_xticks(ticks)
    ax2.set_xticklabels(ticks_name)
    
    ax2.grid(which='both',linestyle=':')
    #lgd = ax2.legend(loc=1, markerscale=3.0, frameon=True, bbox_to_anchor=(1.4, 1))
    ax2.set_title("b",loc="left",weight="bold")
    
    # plot h2o profiles
    for i in flno:
        dati = dat[(dat['FLIGHT'] == i)]
        ax3.semilogx(dati['FLH2O'], dati['PT'], '.', color=colors[i-1], label="F"+str(i)+", "+date[i-1])
        
    ax3.set_xlabel("H$_2$O (ppmv)")
    ax3.set_xlim([2,200])
    ax3.set_xticks([2,5,10,20,50,100,200,200])
    ax3.set_xticklabels(np.array([2,5,10,20,50,100,200,200]).astype(str))
    ax3.set_ylabel("Potential Temperature (K)")
    ax3.grid(which='both',linestyle=':')
    ax3.set_title("c",loc="left",weight="bold")
    lgd = ax3.legend(loc=1, markerscale=3, frameon=True)


    plt.rcParams.update({'font.size': 15})
    plt.savefig("./Paper-Figures/fig1-flight-overview.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)
    plt.show()