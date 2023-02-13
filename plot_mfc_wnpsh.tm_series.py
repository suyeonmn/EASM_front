import numpy as np
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerLine2D
import netCDF4 as nc

filen = 'Fig3a_MFC_PSL_Diff2bx.137-153E_15-25N_minus_117-133E_26-36N.1958-2020.nc'
fpath = '/home/suyeon/front/data/'+filen
opath = '/home/suyeon/front/pic/'
ncfile = nc.Dataset(fpath,mode='r')

hist_up = ncfile.variables['hist_up'][:]
hist_low = ncfile.variables['hist_low'][:]
xghg_up = ncfile.variables['xghg_up'][:]
xghg_low = ncfile.variables['xghg_low'][:]
hist_mme = ncfile.variables['hist_mme'][:]
xghg_mme = ncfile.variables['xghg_mme'][:]
prob = ncfile.variables['prob'][:]

diff_mme = hist_mme - xghg_mme
diff_mme_sig = hist_mme - xghg_mme
diff_up = hist_up - xghg_up
diff_low = hist_low - xghg_low
diff_range = (abs(diff_up)+abs(diff_low))/2./2.

diff_mme_sig_mask = prob <0.95 
diff_mme_sig[diff_mme_sig_mask] = np.nan


## Draw figure ======================================================
fig,ax = plt.subplots(figsize=(8,4))
color1 = "darkorange"
color2 = "seagreen"
ms = 5
alphaN = 0.15

ax.fill_between(np.linspace(1958,2020,63), diff_mme[0,:]+diff_range[0,:], diff_mme[0,:]-diff_range[0,:], alpha=alphaN, color=color2, edgecolor=None,zorder=3)
ax.plot(np.linspace(1958,2020,63), diff_mme[0,:], color=color2, label='MFC', marker='o', markerfacecolor='white', markersize=ms, lw=1.2, zorder=3)
ax.plot(np.linspace(1958,2020,63), diff_mme_sig[0,:], color=color2, marker='o',  markersize=ms, lw=1.2, zorder=3)
ax.set_ylabel(ylabel='\u0394MFC (g/kg/s/10\N{SUPERSCRIPT SIX})',size=17, color=color2)
ax.set_ylim((-1,5))
ax.tick_params(labelsize=13, labelcolor=color2)
ax.set_xlim((1958,2015))
ax.set_xlabel(xlabel='Year',size=17)
ax.tick_params(labelsize=13, labelcolor=color2)
ax.tick_params(axis='x', labelcolor='k',labelsize=13)

ax2 = ax.twinx()
ax2.fill_between(np.linspace(1958,2020,63), (diff_mme[1,:]+diff_range[1,:]), (diff_mme[1,:]-diff_range[1,:]), alpha=alphaN, color=color1, edgecolor=None, zorder=1)
ax2.plot(np.linspace(1958,2020,63), diff_mme[1,:], color=color1, marker='o', markersize=ms, markerfacecolor='white', lw=1.2 ,zorder=1)
ax2.plot(np.linspace(1958,2020,63), diff_mme_sig[1,:], color=color1, marker='o', markersize=ms, lw=1.2, zorder=1)
ax2.set_ylabel(ylabel='\u0394WNPSH (hPa)',size=17, labelpad=10, color=color1)
ax2.set_yticks(np.arange(-0.3, 1.8, 0.3))
ax2.set_ylim((-0.3 ,1.5))
ax2.tick_params(labelsize=13, labelcolor=color1)

plt.title('HIST-XGHG', fontsize=17, loc='left')
plt.tight_layout()
plt.show()
sys.exit()

plt.savefig(opath+"Fig3b_DIFF_MFC_WNPSH2bx.137-153E_15-25N_minus_117-133E_26-36N.1958-2015.png",dpi=400)

