import numpy as np
import sys
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import netCDF4 as nc
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER,LATITUDE_FORMATTER
import matplotlib.colors as mcolors

minlon, maxlon = 105-0.01, 160+0.1
minlat, maxlat = 8, 45+0.01 

filen = 'Reg_Delta_MFC_H500_SLP_TS_PR.P2.Nor.nc'
fpath = '/home/suyeon/front/data/Figure_data/'+filen
filen2 = 'Reg_Delta_Variable_uvq.two_periods.nc'
fpath2 = '/home/suyeon/front/data/Figure_data/'+filen2
opath = '/home/suyeon/front/pic/'
ncfile = nc.Dataset(fpath,mode='r')
ncfile2 = nc.Dataset(fpath2,mode='r')
print('variables in nc file:',ncfile.variables.keys())
lon,lat = ncfile.variables['lon'][:], ncfile.variables['lat'][:]

latid = np.where((lat>=minlat) * (lat<=maxlat))[0]
lonid = np.where((lon>minlon) * (lon<=maxlon))[0]

yy0, yy1 = latid[0],latid[-1]+1
xx0, xx1 = lonid[0],lonid[-1]+1     

lat, lon = lat[latid], lon[lonid]   
lon2d, lat2d = np.meshgrid(lon,lat) 

#iexp = 0
#expn = '1958-1982'
iexp = 1
expn = '1991-2015'
reg = ncfile.variables['reg'][:,yy0:yy1,xx0:xx1]  
prob= ncfile.variables['prob'][:,yy0:yy1,xx0:xx1]
uq = ncfile2.variables['uq'][iexp,yy0:yy1,xx0:xx1]     
vq = ncfile2.variables['vq'][iexp,yy0:yy1,xx0:xx1]
ncfile.close()
ncfile2.close()
uq = uq*3
vq = vq*3


## Draw figure ======================================================
fig = plt.figure(figsize = (12,8), dpi = 100)
gs1 = gridspec.GridSpec(1, 1, left=0.1, right=0.85, top=0.95, bottom=0.1)
cgs = gridspec.GridSpec(1, 1, left=0.865, right=0.885, top=0.9, bottom=0.15)

# Shading 
cax = plt.subplot(cgs[0])
projection_type = ccrs.PlateCarree() 
ax= plt.subplot(gs1[0],projection = projection_type)   
ax.set_extent([minlon, maxlon, minlat, maxlat], crs = ccrs.PlateCarree())
gl = ax.gridlines(crs = ccrs.PlateCarree(), color='k', linestyle = ':', zorder = 7, draw_labels = True, \
    xlocs = np.arange(60,180,20), ylocs=np.arange(0,90,10)) 
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 22, 'color': 'black'}
gl.ylabel_style = {'size': 22, 'color': 'black'}
ax.coastlines(linewidth = 0.7, zorder = 7)

# Colorbar setting 
vmin, vmax = -0.6, 0.6   
clevs = np.linspace(vmin, vmax, 7)    
cs = ax.pcolormesh(lon2d, lat2d, reg[2,:,:], cmap='seismic', \
    transform=ccrs.PlateCarree(), vmin=vmin, vmax= vmax, zorder=1) 
cbar = plt.colorbar(cs, extend='both', cax=cax, ticks=clevs, orientation='vertical')
cbar.ax.tick_params(labelsize=20)
cbar.set_label(label='\u0394SLP (hPa)', size=26, rotation=270, labelpad=35)

# significant SLP
plev = 0.95
cont2 = ax.contourf(lon2d, lat2d, prob[2,:,:], levels=[plev, 1.], zorder=9, vmin=plev, hatches=['//'], alpha=0, transform=ccrs.PlateCarree(), facecolor='red', color="none", edgecolors="b",)

# Draw vector 
skipp = 3 
windscale =1 
quiv = ax.quiver(lon2d[::skipp,::skipp], lat2d[::skipp,::skipp], uq[::skipp,::skipp]*1.2, vq[::skipp,::skipp]*1.2, units='xy', angles='uv', zorder=11, scale_units='inches', scale = 40, headwidth=3, linewidths=3.5,transform=ccrs.PlateCarree(), color='k')


# Draw rectangle
lwbox = 3.5
c_box1 = 'lawngreen'
x_sb = [115, 130, 130, 115, 115]   #MFC
y_sb = [22 , 22, 32, 32, 22]
ax.plot(x_sb, y_sb, transform=ccrs.PlateCarree(), color = c_box1, linewidth=lwbox, zorder = 10)

c_box = 'darkorange'            #WNPSH1
x_wnp = [137, 153, 153, 137, 137]
y_wnp = [15, 15, 25, 25, 15]        
ax.plot(x_wnp, y_wnp, transform=ccrs.PlateCarree(), color = c_box, linewidth=lwbox, zorder = 10)

x_wnp2 = [117, 133, 133, 117, 117]      #WNPSH2
y_wnp2 = [26, 26, 36, 36, 26]
ax.plot(x_wnp2, y_wnp2, transform=ccrs.PlateCarree(), color = c_box, linewidth=lwbox, zorder = 10)


# Title setting 
plt.title('Reg(\u0394FRI, \u0394SLP & \u0394MF850)', fontsize=26, loc='left')
plt.title(expn, fontsize=24, loc='right')
plt.show()
sys.exit()

plt.savefig(opath+'Fig3a_Reg_SLP_MF850.'+expn+'.png',dpi=300)
