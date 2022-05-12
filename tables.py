import numpy as np
import datetime
import scipy.stats as stats

flno = [2,3,4,6,7,8]
timeper = [2,3,4,6,7,8,51,52,0]
maxlag = [0,0,5,10,10,20]

# table 2 calculate flight hours
def table2(dat):
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

    h2ocut = 10
    # flight
    for fi in timeper:
        if fi == 51:
            datx = dat[(dat['FLIGHT'] < 5) & (dat['ASCENT_FLAG'] == 0)]
        elif fi == 52:
            datx = dat[(dat['FLIGHT'] > 5) & (dat['ASCENT_FLAG'] == 0)]
        elif fi == 0:
            datx = dat[((dat['FLIGHT'] < 5) | (dat['FLIGHT'] > 5)) & (dat['ASCENT_FLAG'] == 0)]
        else:
            datx = dat[(dat['FLIGHT'] == fi) & (dat['ASCENT_FLAG'] == 0)]
        
        print()
        print(fi)
        
        # conditions
        for cond in ["all", "cloudy", "clear"]:
            if cond == "all":
                dat1 = datx
            if cond == "cloudy":
                dat1 = datx[(datx['CLOUDY'] == 1)]
            if cond == "clear":
                dat1 = datx[(datx['CLOUDY'] == 0)]
            
            # instrument
            for inst in ['H2O', 'FLH2O', 'FIH2O']:
                if inst == 'H2O':
                    x = dat1[(dat1['CELL_FLAG'] == 0) & (dat1[inst] <= h2ocut)][inst]
                    #x = dat1[(dat1[inst] <= h2ocut)][inst]
                    if cond != "all":
                        print("& ", end="")
                if inst == 'FLH2O':
                    x = dat1[(dat1['F7_DIVE'] == 0) & (dat1[inst] <= h2ocut)][inst]
                    print("/ ", end="")
                if inst == 'FIH2O':
                    x = dat1[(dat1['CLOUDY'] == 0) & (dat1[inst] <= h2ocut)][inst]
                    print("/ ", end="")

                N = len(x[~np.isnan(x)])
                if (cond == "cloudy") & (inst == "FIH2O"):
                    print("- ", end="")
                else:
                    print("{:.2f} ".format(N / 3600), end="")
                
import numpy as np
import datetime
import scipy.stats as stats

flno = [2,3,4,6,7,8]
timeper = [2,3,4,6,7,8,51,52,0]
maxlag = [0,0,5,10,10,20]






# table 3 calculate mean and stdev differences between pairs for each flight
def table3(dat):
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

    h2ocut = 10
    # pairs
    for p in (1,2,3,4,5):
        #fish vs. flash (clear sky)
        if p == 1:
            instx = 'FIH2O'
            insty = 'FLH2O'
            datx = dat[(dat['CLOUDY'] == 0) & (dat[instx] <= h2ocut) & (dat['ASCENT_FLAG'] == 0)]
            daty = dat[(dat['CLOUDY'] == 0) & (dat[insty] <= h2ocut) & (dat['ASCENT_FLAG'] == 0) & (dat['F7_DIVE'] == 0)]
        #fish vs. chiwis (clear sky)
        if p == 2:
            instx = 'FIH2O'
            insty = 'H2O'
            datx = dat[(dat['CLOUDY'] == 0) & (dat[instx] <= h2ocut) & (dat['ASCENT_FLAG'] == 0)]
            daty = dat[(dat['CLOUDY'] == 0) & (dat[insty] <= h2ocut) & (dat['ASCENT_FLAG'] == 0) & (dat['CELL_FLAG'] == 0)]
        # chiwis vs. flash (clear sky)
        if p == 3:
            instx = 'H2O'
            insty = 'FLH2O'
            datx = dat[(dat['CLOUDY'] == 0) & (dat[instx] <= h2ocut) & (dat['ASCENT_FLAG'] == 0) & (dat['CELL_FLAG'] == 0)]
            daty = dat[(dat['CLOUDY'] == 0) & (dat[insty] <= h2ocut) & (dat['ASCENT_FLAG'] == 0) & (dat['F7_DIVE'] == 0)]
        # chiwis vs. flash (cloudy sky)
        if p == 4:
            instx = 'H2O'
            insty = 'FLH2O'
            datx = dat[(dat['CLOUDY'] == 1) & (dat[instx] <= h2ocut) & (dat['ASCENT_FLAG'] == 0) & (dat['CELL_FLAG'] == 0)]
            daty = dat[(dat['CLOUDY'] == 1) & (dat[insty] <= h2ocut) & (dat['ASCENT_FLAG'] == 0) & (dat['F7_DIVE'] == 0)]
        # chiwis vs. flash (all sky)
        if p == 5:
            instx = 'H2O'
            insty = 'FLH2O'
            datx = dat[(dat[instx] <= h2ocut) & (dat['ASCENT_FLAG'] == 0) & (dat['CELL_FLAG'] == 0)]
            daty = dat[(dat[insty] <= h2ocut) & (dat['ASCENT_FLAG'] == 0) & (dat['F7_DIVE'] == 0)]
        
        print()
        print(p)

        # flight
        for fi in timeper:
            if fi == 51:
                x = datx[(datx['FLIGHT'] < 5)][instx]
                y = daty[(daty['FLIGHT'] < 5)][insty]
            elif fi == 52:
                x = datx[(datx['FLIGHT'] > 5)][instx]
                y = daty[(daty['FLIGHT'] > 5)][insty]
            elif fi == 0:
                x = datx[(datx['FLIGHT'] < 5) | (datx['FLIGHT'] > 5)][instx]
                y = daty[(daty['FLIGHT'] < 5) | (datx['FLIGHT'] > 5)][insty]
            else:
                x = datx[(datx['FLIGHT'] == fi)][instx]
                y = daty[(daty['FLIGHT'] == fi)][insty]

            mask = ~np.isnan(x) & ~np.isnan(y)
            if len(x[mask]) > 0:
                slope, intercept, rvalue, pvalue, se = stats.linregress(x[mask],y[mask])

                bias = (x[mask] - y[mask]) / y[mask] * 100.0
                absbias = np.abs(bias)
                meanbias = np.mean(bias)
                stdbias = np.std(bias)

                print("{:.1f} ({:.1f}) & ".format(meanbias, stdbias), end="")
            else:
                print("-- & ", end="")
