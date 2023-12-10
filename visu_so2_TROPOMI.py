import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

ifileS ='/home/mhaghigh/solforia/sandra_so2/TROPOMI_SO2_mass_loading_LaSoufriere_SZA.nc'
nc = Dataset(ifileS,'r')
fig, ((axs0,axs1),(axs2,axs3))= plt.subplots(2,2,figsize=(25,20))

def read_data(n):
    so2= nc.variables['SO2_mass_loading'][:,:,n]
    lon= nc.variables['longitude'][:]
    lat= nc.variables['latitude'][:]
    return so2,lon,lat

def plot(axs,day):
    cmap = [(0.0,0.0,0.0)] + [(cm.jet(i)) for i in range(1,256)]
    cmap = mpl.colors.ListedColormap(cmap)
    bounds = np.arange(0,1,5e-2)
    bounds_2= np.arange(0,1.1,1e-1)
    bounds_2=np.around(bounds_2,2)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    norm_2 = mpl.colors.BoundaryNorm(bounds_2, cmap.N)
    cmap_2 = LinearSegmentedColormap.from_list("", ["white","lightskyblue","steelblue","green","yellowgreen","yellow","gold","red","firebrick","darkred"])
    plt.subplot(axs)
    m = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=30,\
           llcrnrlon=-80,urcrnrlon=-20,resolution='c')
    m.drawcoastlines()
    n=day-9
    so2_m,lon_m,lat_m=read_data(n)
    x, y = m(lon_m, lat_m)
    day=str(day)
    m.drawparallels(np.arange( -90.,90.,10.),labels=[1,0,0,0],fontsize=15)
    m.drawmeridians(np.arange(-180.,180,10.),labels=[0,0,0,1],fontsize=15)
    axs.set_title(day+' April 2021', fontsize=25)
    axs=m.pcolormesh(x, y, so2_m, cmap=cmap_2,  norm=norm)
    cax = fig.add_axes([0.15,0.05,0.7,0.02])
    cbar_bounds = bounds_2
    cbar_ticks =  bounds_2
    cbar = fig.colorbar(axs,cax=cax, norm=norm_2, boundaries=cbar_bounds, ticks=cbar_ticks, orientation='horizontal')
    cbar.set_label('so2 mass loading(g/m2)', fontsize=35)
    cbar.ax.tick_params(labelsize=30)
    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.05)
    return
plot(axs0,9)
plot(axs1,10)
plot(axs2,11)
plot(axs3,12)
plt.savefig('so2_profile.png')
plt.show()

plt.close