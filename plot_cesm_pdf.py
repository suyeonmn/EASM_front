import sys
from scipy.stats import norm
import scipy.stats as ss 
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from numpy import linspace
import netCDF4 as nc
import numpy as np
import seaborn as sns

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

def get_risk_index(x0, y0, x1, y1, cri):
    data_ext = np.sum([jj < cri for jj in x0])
    data_extH = np.sum([ss < cri for ss in x1])
    P0 = 1-y0[data_ext-1]
    P1 = 1-y1[data_extH-1] 
    RR = P1/P0
    FAR = 1-(P0/P1)
    return RR, FAR

## Input file
filen = 'Fig1c_LE_MEAN_Bootstrap_PDF_input.1991-2015.5yr_running.nc'
fpath = '/home/suyeon/front/data/Figure_data/'+filen
opath = '/home/suyeon/front/pic/'
ncfile = nc.Dataset(fpath,mode='r')
hist = ncfile.variables['hist_raw'][:]
xghg = ncfile.variables['xghg_raw'][:]
hist_bt = ncfile.variables['hist'][:]
xghg_bt = ncfile.variables['xghg'][:]
hist_hi = ncfile.variables['hist_High'][:]
xghg_hi = ncfile.variables['xghg_High'][:]
hist_lw = ncfile.variables['hist_Low'][:]
xghg_lw = ncfile.variables['xghg_Low'][:]

## Changing rate
a = np.mean(hist_bt)
b = np.mean(xghg_bt)
c = (a-b)/b*100.
print("Change rate",np.mean(hist_bt),np.mean(xghg_bt),c)

min_a = np.min(hist_bt)
max_a = np.max(hist_bt)
min_b = np.min(xghg_bt)
max_b = np.max(xghg_bt)
c_max = (max_a - min_b)/min_b*100.
c_min = (min_a - max_b)/max_b*100.
print("Min/Max/Mean= ", c_min, c_max, c)


## Find the threshold
mean = np.mean(xghg)
std = np.std(xghg)
pct_99 = norm.ppf(0.99, mean, std)
pct_95 = norm.ppf(0.95, mean, std)
print("XGHG", mean, std, pct_95)

## Cumulative Distribution Function (CDF) 
cdf = get_discrete_cdf(xghg)
x_p = list(zip(xghg,cdf))
x_p.sort(key=lambda it: it[0])
x = [it[0] for it in x_p]
y = [it[1] for it in x_p]

cdfH = get_discrete_cdf(hist)
h_p = list(zip(hist,cdfH))
h_p.sort(key=lambda it: it[0])
xh = [it[0] for it in h_p]
yh = [it[1] for it in h_p]  #plt.plot(xh, yh)

rr95, far95 =  get_risk_index(x, y, xh, yh, pct_95)
rr99, far99 =  get_risk_index(x, y, xh, yh, pct_99)
print("RR",rr95, rr99)
print("FAR",far95,far99)

## PDF
# this create the kernel, given an array it will estimate the probability over that values
kde_h = gaussian_kde(hist)
kde_x = gaussian_kde(xghg)

# these are the values over wich your kernel will be evaluated
nbin = 50
dist_space_x = linspace( 8,22, nbin )
dist_space_h = linspace( 8,22, nbin )


## Draw figure ======================================================
fig = plt.figure(figsize = (10,7), dpi = 100)

# PDF plot
plt.plot( dist_space_x, kde_x(dist_space_x), color='blue', linewidth = 1.5)
plt.plot( dist_space_h, kde_h(dist_space_h), color='red', linewidth = 1.5)
plt.vlines(np.mean(xghg_bt), 0, 3.5, color='blue', linestyles='--', linewidth=1.8)
plt.vlines(np.mean(hist_bt), 0, 3.5, color='red', linestyles='--', linewidth=1.8)

plt.fill_betweenx(np.arange(0., 0.4, 0.1), xghg_lw, xghg_hi, color='blue',alpha=0.15)
plt.fill_betweenx(np.arange(0., 0.4, 0.1), hist_lw, hist_hi, color='red',alpha=0.15)

# Extreme
plt.fill_between(dist_space_h, kde_h(dist_space_h),  where=( pct_95 < dist_space_h) , color='red',alpha=0.2)
plt.fill_between(dist_space_x, kde_x(dist_space_x),  where=( pct_95 < dist_space_x) , color='blue',alpha=0.2)
plt.fill_between(dist_space_h, kde_h(dist_space_h),  where=( pct_99 < dist_space_h) , color='red',alpha=0.4)
plt.fill_between(dist_space_x, kde_x(dist_space_x),  where=( pct_99 < dist_space_x) , color='blue',alpha=0.4)

# plot setting
plt.xlabel("Intensity of frontal rainfall (mm/day)", fontsize=25, labelpad=16)
plt.ylabel("Probability density", fontsize=25, labelpad=16)
plt.xticks(np.arange(8,24,2), fontsize=20) 
plt.yticks(np.arange(0,0.4,0.05), fontsize=20)
plt.ylim((0,.25))
plt.xlim((8,22))
plt.legend(("XGHG","HIST"),fontsize=20)

loc1 = 16.3
plt.annotate(r'$RR_{95}=$'+'{:.2f}'.format(rr95), xy=(loc1, 0.142), ha='left', fontsize=20)
plt.annotate(r'$FAR_{95}=$'+'{:.2f}'.format(far95), xy=(loc1, 0.125), ha='left', fontsize=20)

loc2 = 17.3
plt.annotate(r'$RR_{99}=$'+'{:.2f}'.format(rr99), xy=(loc2, 0.082), ha='left', fontsize=20)
plt.annotate(r'$FAR_{99}=$'+'{:.2f}'.format(far99), xy=(loc2, 0.065), ha='left', fontsize=20)

plt.title('1991-2015', fontsize=25, loc='right')
plt.tight_layout()
plt.show()
sys.exit()

plt.savefig(opath+'Fig1_PDF_LE1_MEAN.Boot_RR.png',dpi=300)

