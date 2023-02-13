import numpy as np
import sys
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
import netCDF4 as nc
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER,LATITUDE_FORMATTER

minlon, maxlon = 113, 140+0.01 
minlat, maxlat = 24.95, 40.01 

filen = 'Fig1_APHRO_DelI_Diff.P2-P1.NOV22_10per.nc'
fpath = '/home/suyeon/front/data/Figure_data/'+filen
opath = '/home/suyeon/front/pic/'
ncfile = nc.Dataset(fpath,mode='r')
lon,lat = ncfile.variables['lon'][:], ncfile.variables['lat'][:]

latid = np.where((lat>=minlat) * (lat<=maxlat))[0]
lonid = np.where((lon>minlon) * (lon<=maxlon))[0]

yy0, yy1 = latid[0],latid[-1]+1
xx0, xx1 = lonid[0],lonid[-1]+1     

lat, lon = lat[latid], lon[lonid]
lon2d, lat2d = np.meshgrid(lon,lat) 

diff = ncfile.variables['diff'][yy0:yy1,xx0:xx1]
prob = ncfile.variables['prob'][yy0:yy1,xx0:xx1]
ncfile.close()


## Draw figure ======================================================
#  map setting 
fig = plt.figure(figsize = (10,7), dpi = 100)
gs1 = gridspec.GridSpec(1, 1, left=0.1, right=0.85, top=0.95, bottom=0.1)
cgs = gridspec.GridSpec(1, 1, left=0.86, right=0.88, top=0.88, bottom=0.17)

# Shading 
cax = plt.subplot(cgs[0])
projection_type = ccrs.Mercator()
ax= plt.subplot(gs1[0],projection = projection_type)
ax.set_extent([minlon, maxlon, minlat, maxlat], crs = ccrs.PlateCarree())
gl = ax.gridlines(crs = ccrs.PlateCarree(), color='k', linestyle = ':', zorder = 7, draw_labels = True, \
    xlocs = np.arange(120,150,10), ylocs=np.arange(25,45,5)) 
gl.xlabels_top = False
gl.ylabels_right = False

gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 20, 'color': 'black'}
gl.ylabel_style = {'size': 20, 'color': 'black'}
ax.coastlines(linewidth = 0.7, zorder = 7)


# Colorbar setting 
vmin, vmax = -8., 8.
clevs = np.linspace(vmin, vmax, 9)     
cs = ax.pcolormesh(lon2d, lat2d, diff, cmap = 'BrBG', \
    transform=ccrs.PlateCarree(), vmin=vmin, vmax= vmax, zorder=1)
cbar = plt.colorbar(cs, extend='both', cax=cax, ticks=clevs, orientation='vertical')
cbar.ax.tick_params(labelsize=20)
cbar.set_label(label='$\Delta$Frontal rainfall (mm/day)', size=25, rotation=270, labelpad=40)

plev = 0.95
cont1 = ax.contourf(lon2d,lat2d, prob, colors = 'r', levels = [plev, 1.], zorder=8, vmin=plev, hatches=['//'],  alpha=0, transform=ccrs.PlateCarree())


# Draw rectangle
x = [117, 135, 135, 117, 117]
y = [28, 28, 36, 36, 28]
ax.plot(x, y, transform=ccrs.PlateCarree(), color ='k', zorder = 9)


# Title setting 
plt.title('APHRODITE', fontsize=25, loc='left')
plt.title('+19.8%', fontsize=25, loc='right')

plt.show()
sys.exit()
plt.savefig(opath+'Fig1_APHRO_DelI_Mean.P2-P1.BOX.png',dpi=300)
