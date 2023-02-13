import sys
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import kstest
from numpy import linspace
import netCDF4 as nc
import numpy as np

def get_discrete_cdf(values):
    values = (values - np.min(values)) / (np.max(values) - np.min(values))

    values_sort = np.sort(values)
    values_sum = np.sum(values)

    values_sums = []
    cur_sum = 0
    for it in values_sort:
        cur_sum += it
        values_sums.append(cur_sum)

    cdf = [values_sums[np.searchsorted(values_sort, it)]/values_sum for it in values]
    return cdf

    
filen = 'Figerprint_APHRO_Intensity_1958-2015.slope.nc'
fpath = '/home/suyeon/front/data/Figure_data/'+filen
opath = '/home/suyeon/front/pic/'

ncfile = nc.Dataset(fpath,mode='r')
jra = ncfile.variables['jra'][:]
nat = ncfile.variables['nat'][:]
ant = ncfile.variables['ant'][:]
jra_all = 1.430519
jra_ext = 1.887262
bootmin, bootmed, bootmax = -2.631987,-1.909958, -1.309236
bootminH, bootmedH, bootmaxH = -0.2253365, 	0.3775114, 1.035205

kde_x = gaussian_kde(nat)
kde_h = gaussian_kde(ant)


## Cumulative Distribution Function (CDF)
cdf = get_discrete_cdf(nat)
x_p = list(zip(nat,cdf))
x_p.sort(key=lambda it: it[0])
x = [it[0] for it in x_p]
y = [it[1] for it in x_p]

nbin = 200
dist_space_x = linspace( min(nat), max(nat),nbin)
dist_space_h = linspace( min(ant), max(ant),nbin)
print(min(nat),max(nat),np.mean(nat),np.median(nat))


## Draw figure ======================================================
fig = plt.figure(figsize = (10,6), dpi = 100)

lw = 1.7 
color_b = 'mediumblue'
color_r = 'crimson'
plt.plot( dist_space_h, kde_h(dist_space_h), color=color_r, linewidth = lw)
plt.plot( dist_space_x, kde_x(dist_space_x), color=color_b, linewidth = lw)

plt.vlines(jra, 0,1, color='k', linestyles='-', linewidth=lw)
jra_c = 'limegreen'
plt.vlines(jra_all, 0,1, color= jra_c, linestyles='-', linewidth= lw)
plt.vlines(jra_ext, 0,1, color= jra_c, linestyles='--', linewidth= lw)

plt.vlines(bootmed, 0, 1, color=color_b, linestyles=':', linewidth = lw) 
plt.vlines(bootmedH, 0, 1, color=color_r, linestyles=':', linewidth = lw) 
plt.fill_betweenx(np.arange(0., 0.4, 0.1), bootmin, bootmax, color=color_b, alpha=0.15 ,edgecolor=None)
plt.fill_betweenx(np.arange(0., 0.4, 0.1), bootminH, bootmaxH, color=color_r, alpha=0.15, edgecolor=None )


# plot setting
plt.title('Frontal rainfall intensity', fontsize = 19, loc='left')
plt.xlabel("Linear slope of $\it{s}$", fontsize=18, labelpad=16)
plt.ylabel("Probability density", fontsize=18, labelpad=16)
plt.xticks(np.arange(-10,10,1), fontsize=16)
plt.yticks(np.arange(0,0.35,0.05), fontsize=16)
plt.ylim((0,0.3))
plt.xlim((-6,6))
plt.legend(("HIST","XGHG","APHRODITE","JRA55","JRA55_expan"),fontsize=15, loc='upper right')
plt.text(-5.5, 0.27, r'P < 0.01', fontsize=17, color=color_b) #fontdict=font)

plt.tight_layout()
plt.show()
sys.exit()

plt.savefig(opath+'Fig_fingerprint_FRI.png',dpi=300)



