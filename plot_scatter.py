import sys
from scipy.stats import norm
import scipy.stats as ss
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from numpy import linspace
import netCDF4 as nc
import numpy as np
import seaborn as sns

## Input file
#filen = 'Fig2_Trend_Variable_NOV17.nc'

#filen = 'Fig2_Trend_Variable_DEC15_med.nc'
#filen = 'Fig4_Trend_MFC_WNPSH2bx_FEB6_med.nc'
filen = 'Fig4_Trend_2bx.137-153E_15-25N_minus_117-133E_26-36N.nc'

fpath = '/home/suyeon/front/data/Figure_data/'+filen
ncfile = nc.Dataset(fpath,mode='r')
print(ncfile.variables['hist_p1'][:].shape)
hist_p1 = ncfile.variables['hist_p1'][:]
hist_p2 = ncfile.variables['hist_p2'][:]
xghg_p1 = ncfile.variables['xghg_p1'][:]
xghg_p2 = ncfile.variables['xghg_p2'][:]

xNo = 0     #MFC
yNo = 1     #WNPSH

hist_p1[1,:] = hist_p1[1,:]*1.
hist_p2[1,:] = hist_p2[1,:]*1.
xghg_p1[1,:] = xghg_p1[1,:]*1.
xghg_p2[1,:] = xghg_p2[1,:]*1.        #Unit change Pa -> hPa

#print(xghg_p1[1,:])
#print(hist_p2[1,:])
#sys.exit()
color1 = 'crimson'
color0 = 'mediumblue'

## Scatter plot 
# definitions for the axes
bottom, height = 0.1, 0.48
left, width = 0.13, 0.51
spacing = 0.005
box_size = 0.14
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, box_size]
rect_histx2 = [left, bottom + height + spacing + box_size + spacing, width, box_size]
rect_histy = [left + width + spacing, bottom, box_size, height]
rect_histy2 = [left + width + spacing+ box_size + spacing, bottom, box_size, height]

# start with a rectangular Figure
plt.figure(figsize=(10,10))
ax_scatter = plt.axes(rect_scatter)
ax_scatter.tick_params(direction='in', top=True, right=True)
ax_histx = plt.axes(rect_histx)
ax_histx2 = plt.axes(rect_histx2)
ax_histx.tick_params(direction='in', labelbottom=False)
ax_histx2.tick_params(direction='in', labelbottom=False)
ax_histy = plt.axes(rect_histy)
ax_histy2 = plt.axes(rect_histy2)
ax_histy.tick_params(direction='in', labelleft=False)
ax_histy2.tick_params(direction='in', labelleft=False)

# the scatter plot
lw1 = 1.0
lw2 = 1.0
trans = 0.3 #0.4 #;0.6 #0.7
ms = 15
ms2 = 22
ax_scatter.set_xlabel(xlabel='Trend of MFC (g/kg/s/10\N{SUPERSCRIPT SIX}/yr)', size=20,labelpad=10)
ax_scatter.set_ylabel(ylabel='Trend of WNPSH (Pa/yr)', size=20, labelpad=10)

## Scatter_plot
ax_scatter.scatter(hist_p1[xNo,:],hist_p1[yNo,:], facecolor='none',  edgecolors=color1, alpha= trans, linewidths=lw1, s=ms, zorder = 4)
#ax_scatter.scatter(hist_p2[xNo,:],hist_p2[yNo,:], mark='o'facecolor='none',  edgecolors=color1 , alpha= trans, linewidths=lw2,  zorder =4)
ax_scatter.scatter(hist_p2[xNo,:],hist_p2[yNo,:], c=color1,facecolor=color1, edgecolors='none',  s=ms2,  alpha= trans, linewidths=lw2, zorder =4)
ax_scatter.scatter(xghg_p1[xNo,:],xghg_p1[yNo,:], facecolors='none',s = ms,  edgecolors=color0, alpha=trans, linewidths=lw1,  zorder =4 )
ax_scatter.scatter(xghg_p2[xNo,:],xghg_p2[yNo,:], facecolors=color0, edgecolors='none', c=color0, s=ms2, alpha=trans, linewidths=lw2,  zorder = 4)
ax_scatter.tick_params(labelsize=16)

## Shading
xdata2 = xghg_p2[xNo,:]
ydata2 = xghg_p2[yNo,:]
k2 = gaussian_kde([xdata2,ydata2])
xi2, yi2 = np.mgrid[xdata2.min()-2:xdata2.max()+2:300j,ydata2.min()-2:ydata2.max()+2:300j]
zi2=k2(np.vstack([xi2.ravel(),yi2.ravel()]))
level0 = np.arange(5,16,2)
#ax_scatter.contourf(xi2,yi2,zi2.reshape(xi2.shape), levels=level0, cmap='Blues', linestyles='solid', alpha=0.8,  zorder = 3)
#ax_scatter.contour(xi2,yi2,zi2.reshape(xi2.shape), levels=level0, colors=color0, linewidths=1., linestyles='solid', alpha=0.4,  zorder = 3)

xdata = hist_p2[xNo,:]
ydata = hist_p2[yNo,:]
k = gaussian_kde([xdata,ydata])
xi, yi = np.mgrid[xdata.min()-2:xdata.max()+2:300j,ydata.min()-2:ydata.max()+2:300j]
zi=k(np.vstack([xi.ravel(),yi.ravel()]))
level1 = np.arange(5,16,2)
#level1 = np.arange(10,50,5)
ax_scatter.contourf(xi,yi,zi.reshape(xi.shape), levels=level1, cmap='Reds', linestyles='solid', alpha=0.8,  zorder = 3)

## Contour for P1
xdata_x1 = xghg_p1[xNo,:]
ydata_x1 = xghg_p1[yNo,:]
k_x1 = gaussian_kde([xdata_x1,ydata_x1])
xi_x1, yi_x1 = np.mgrid[xdata_x1.min()-2:xdata_x1.max()+2:300j,ydata_x1.min()-2:ydata_x1.max()+2:300j]
zi_x1=k_x1(np.vstack([xi_x1.ravel(),yi_x1.ravel()]))
#level0 = np.arange(0.1,0.6,0.05)
#level0 = np.arange(0.8,5.,0.3)
level0 = np.arange(5,16,2)
#ax_scatter.contour(xi_x1,yi_x1,zi_x1.reshape(xi_x1.shape), levels=level0, colors=color0,  linewidths=1. ,linestyles='dashed', alpha=0.4,  zorder = 2)

xdata_h1 = hist_p1[xNo,:]
ydata_h1 = hist_p1[yNo,:]
k_h1 = gaussian_kde([xdata_h1,ydata_h1])
xi_h1, yi_h1 = np.mgrid[xdata_h1.min()-2:xdata_h1.max()+2:300j,ydata_h1.min()-2:ydata_h1.max()+2:300j]
zi_h1=k_h1(np.vstack([xi_h1.ravel(),yi_h1.ravel()]))
level1 = np.arange(5,16,2)
#level1 = np.arange(0.8,5,0.3)
#level1 = np.arange(10,50,5)
ax_scatter.contour(xi_h1,yi_h1,zi_h1.reshape(xi_h1.shape), levels=level1, colors=color1, linewidths=1., linestyles='dashed', alpha=0.8,  zorder = 2)



# now determine nice limits by hand:
binwidth = 0.25
lim = np.ceil(np.abs([hist_p1[xNo,:],hist_p1[yNo,:]]).max() / binwidth) * binwidth
#ax_scatter.set_xlim((-0.75,1.2))
#ax_scatter.set_xlim((-0.9,1.2))
#ax_scatter.set_ylim((-0.7, 0.7))
ax_scatter.set_ylim((-4.5, 4.5))
#ax_scatter.set_ylim((-0.2, 0.2))

# Reference line
line_config = dict()
line_config['linestyle'] = '--'
line_config['color'] = 'gray'
line_config['linewidth'] = 0.6
ax_scatter.axhline(0,**line_config)
ax_scatter.axvline(0,**line_config)

# legend with custom labels
#lines, labels = ax_scatter.legend_elements()
labels = ["HIST_P1","HIST_P2","XGHG_P1","XGHG_P2"]
legend = ax_scatter.legend(labels, loc="upper left", fontsize=15)
ax_scatter.add_artist(legend)

bins = np.arange(-lim, lim + binwidth, binwidth)
nbin = 50

#y_min, y_max = -0.075, 0.075
#y_min, y_max = -4.5, 4.5
y_min, y_max = -0.12,0.18
x_min, x_max = -0.85, 0.85
#x_min, x_max = -0.55, 0.7
dist_space_x = linspace( x_min, x_max, nbin )
kde_h1 = gaussian_kde(hist_p1[xNo,:])
kde_h2 = gaussian_kde(hist_p2[xNo,:])
kde_x1 = gaussian_kde(xghg_p1[xNo,:])
kde_x2 = gaussian_kde(xghg_p2[xNo,:])
#ax_histx2.plot(dist_space_x, kde_h1(dist_space_x), color='r', linewidth = 1.3,  linestyle='-', alpha=0.5)
ax_histx.plot(dist_space_x, kde_h1(dist_space_x), color=color1, linewidth = lw1,  linestyle='--') 
ax_histx2.plot(dist_space_x, kde_h2(dist_space_x), color= color1, linewidth = lw2,  linestyle='-')
ax_histx.plot(dist_space_x, kde_x1(dist_space_x), color=color0, linewidth = lw1,  linestyle='--')
ax_histx2.plot(dist_space_x, kde_x2(dist_space_x), color=color0, linewidth = lw2,  linestyle='-')
ax_histx.axvline(0,**line_config)
ax_histx2.axvline(0,**line_config)
ax_histx.set_yticks(np.arange(0, 3, 1))
ax_histx2.set_yticks(np.arange(0, 3, 1))
ax_histx.set_ylim((0,2.3))
ax_histx2.set_ylim((0,2.3))
ax_histx.set_xticks(np.arange(-1.6, 1.6, 0.4))
ax_histx2.set_xticks(np.arange(-1.6, 1.6, 0.4))
ax_histx.set_xlim((x_min, x_max))
ax_histx2.set_xlim((x_min, x_max))

labels2 = ["HIST_P2","XGHG_P2"] #"XGHG_P1","XGHG_P2"]
legend2 = ax_histx2.legend(labels2, bbox_to_anchor=(1.38, 1.073), fontsize=15)
ax_histx2.add_artist(legend2)

labels3 = ["HIST_P1","XGHG_P1"] #"XGHG_P1","XGHG_P2"]
legend3 = ax_histx.legend(labels3, bbox_to_anchor=(1.38, 1.073), fontsize=15)
ax_histx.add_artist(legend3)

ax_scatter.set_xticks(np.arange(-1.5, 1.5, 0.3))
ax_scatter.set_xlim((x_min, x_max))

dist_space_y = linspace( y_min, y_max, nbin )
vkde_h1 = gaussian_kde(hist_p1[yNo,:])
vkde_h2 = gaussian_kde(hist_p2[yNo,:])
vkde_x1 = gaussian_kde(xghg_p1[yNo,:])
vkde_x2 = gaussian_kde(xghg_p2[yNo,:])
ax_histy.plot(vkde_h1(dist_space_y), dist_space_y,  color=color1, linewidth = lw1,  linestyle='--', alpha=1.)
ax_histy2.plot(vkde_h2(dist_space_y), dist_space_y,  color=color1, linewidth = lw2,  linestyle='-')
ax_histy.plot(vkde_x1(dist_space_y), dist_space_y,  color=color0, linewidth = lw1,  linestyle='--', alpha=1.)
ax_histy2.plot(vkde_x2(dist_space_y), dist_space_y,  color=color0, linewidth = lw2,  linestyle='-')
ax_histx.tick_params(labelsize=16)
ax_histx2.tick_params(labelsize=16)
ax_histy.tick_params(labelsize=16)
ax_histy2.tick_params(labelsize=16)
ax_histy.axhline(0,**line_config)
ax_histy2.axhline(0,**line_config)
ax_histy.set_xticks(np.arange(0, 12,3))
ax_histy2.set_xticks(np.arange(0, 12, 3))
ax_histy.set_xlim((0, 8.9))
ax_histy2.set_xlim((0,8.9))
ax_histy.set_ylim((y_min, y_max))
ax_histy2.set_ylim((y_min, y_max))
ax_scatter.set_ylim((y_min, y_max))

ax_histx.text(-1.22, 1.0, 'Density Dist.', fontsize=20, rotation=90.)#, labelpad=10)
ax_histy.text( 3., -0.16 , 'Density Dist.', fontsize=20)#, labelpad=10)

# P1 & P2 info 
ax_histx.text(-0.8, 1.8, 'P1', fontsize=16)
ax_histx2.text(-0.8, 1.8, 'P2', fontsize=16)
ax_histy.text(6, 0.16, 'P1', fontsize=16) #0.9
ax_histy2.text(6, 0.16, 'P2', fontsize=16)


ax_histx.set_xlim(ax_scatter.get_xlim())
ax_histy.set_ylim(ax_scatter.get_ylim())
#plt.tight_layout()
plt.show()
sys.exit()

plt.savefig('../../pic/Fig4_Scatter_MFC_gradWNPSH.Feb7.png',dpi=300)
#plt.savefig('../../pic/Fig3_Scatter_WNPSH_MFC_rev2_all.png',dpi=300)



sys.exit()
