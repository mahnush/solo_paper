import scipy.interpolate as sci
import numpy as np
from netCDF4 import Dataset
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from global_land_mask import globe

ipath = '/home/mhaghigh/solforia/pdf/'
nc_icon = ipath + 'per_run_10_2data.nc'
mdata=nc_icon
nc_m=Dataset(mdata,'r')
lat_fine=nc_m.variables['lat'][:]
lon_fine=nc_m.variables['lon'][:]
lwp=nc_m.variables['nd_dw'][:,:]

nc_control= ipath + 'control_run_10_2data.nc'
mdata=nc_control
nc_c=Dataset(mdata,'r')
lwp_c=nc_c.variables['nd_dw'][:,:]

#--------------------------------
ipath='/home/mhaghigh/solforia/sandra_so2/'
nc_model='TROPOMI_SO2_mass_loading_LaSoufriere_SZA.nc'
mdata=ipath+nc_model
nc_m=Dataset(mdata,'r')
lat_m=nc_m.variables['latitude'][:]
lon_m=nc_m.variables['longitude'][:]
so2=nc_m.variables['SO2_mass_loading'][:,:,1]
so2[so2<0]=0
so2[np.isnan(so2)]=0

#----------------------------------
#print(so2_mask)
print('gg')
f = sci.RectBivariateSpline(lat_m, lon_m,so2 )
scale_interp = f(lat_fine, lon_fine)
print(np.shape(scale_interp))

lon_mesh,lat_mesh=np.meshgrid(lon_fine,lat_fine)
print(np.shape(lon_mesh))
globe_land_mask = globe.is_land(lat_mesh,lon_mesh)
lwp_ocean=ma.masked_where(globe_land_mask==True, lwp)
lwp_ocean_c=ma.masked_where(globe_land_mask==True, lwp_c)

lon_mask_in=ma.masked_where(scale_interp<0.02,lon_mesh)
lwp_in=ma.masked_where(lon_mask_in==True,lwp_ocean)
lon_mask_out=ma.masked_where(scale_interp>0.02,lon_mesh)
lwp_out=ma.masked_where(lon_mask_out==True,lwp_ocean)


lwp_in_c=ma.masked_where(lon_mask_in==True,lwp_ocean_c)
lwp_out_c=ma.masked_where(lon_mask_out==True,lwp_ocean_c)
#lwp_in= ma.masked_where(lwp_in>2000,lwp_in)
#lwp_in_c= ma.masked_where(lwp_in_c>1000,lwp_in_c)
#lwp_out = ma.masked_where(lwp_out>2000,lwp_out)
#lwp_out_c = ma.masked_where(lwp_out_c>500,lwp_out_c)

#===================================================================

lwp_in_nomask=ma.compressed(lwp_in)
lwp_in_nomask=lwp_in_nomask.flatten()

lwp_out_nomask=ma.compressed(lwp_out)
lwp_out_nomask=lwp_out_nomask.flatten()

lwp_in_nomask_c=ma.compressed(lwp_in_c)
lwp_in_nomask_c=lwp_in_nomask_c.flatten()

lwp_out_nomask_c=ma.compressed(lwp_out_c)
lwp_out_nomask_c=lwp_out_nomask_c.flatten()


#print(np.shape(re_in_nomask))
#print(np.shape(re_out_nomask))
print('in')
mean_in= np.mean(lwp_in_nomask)
mean_in=str(round(mean_in,2))
medi_in= np.median(lwp_in_nomask)
medi_in=str(round(medi_in,2))
print(mean_in)

print('out')
mean_out=np.mean(lwp_out_nomask)
mean_out=str(round(mean_out,2))
medi_out= np.median(lwp_out_nomask)
medi_out=str(round(medi_out,2))

print(mean_out)
print('in')
mean_in_c= np.mean(lwp_in_nomask_c)
mean_in_c=str(round(mean_in_c,2))
medi_in_c=np.median(lwp_in_nomask_c)
medi_in_c=str(round(medi_in_c,2))
print(mean_in_c)

print('out')
mean_out_c=np.mean(lwp_out_nomask_c)
mean_out_c=str(round(mean_out_c,2))
medi_out_c= np.median(lwp_out_nomask_c)
medi_out_c=str(round(medi_out_c,2))

print(mean_out_c)
weight_in =(1+np.zeros(len(lwp_in_nomask)))/len(lwp_in_nomask)
weight_out= (1+np.zeros(len(lwp_out_nomask)))/len(lwp_out_nomask)
weight_in_c =(1+np.zeros(len(lwp_in_nomask_c)))/len(lwp_in_nomask_c)
weight_out_c =(1+np.zeros(len(lwp_out_nomask_c)))/len(lwp_out_nomask_c)

range_r=(1,1000)

fig, (axs0,axs1) = plt.subplots(1,2,figsize=(30,20))
axs0.hist(lwp_in_nomask,   bins=50 , range=range_r, weights=weight_in, color="red", histtype='step', linewidth=4, label='Volcano')
axs0.hist(lwp_in_nomask_c,   bins=50 , range=range_r, weights=weight_in_c, color="blue", histtype='step', linewidth=4, label='No-Volcano')
axs0.legend(loc='upper right',fontsize=40)
ticks = np.arange(0,0.20,0.02)
axs0.set_yticks(ticks)
axs0.tick_params(axis='x', labelsize=40 ) #to Set Matplotlib Tick Labels Font Size
axs0.tick_params(axis='y', labelsize=40 )
axs0.set_xlabel('$ \mathrm{N_d}$ ($\mathrm{cm^{-3}}$)',fontsize=40)
axs0.set_ylabel('Relative Frequency',fontsize=40)



#------------------------------------
axs1.hist(lwp_out_nomask,   bins=50 , range=range_r, weights=weight_out, color="red", histtype='step', linewidth=4, label='Volcano')
axs1.hist(lwp_out_nomask_c,   bins=50 , range=range_r, weights=weight_out_c, color="blue", histtype='step', linewidth=4, label='No-Volcano')
axs1.legend(loc='upper right',fontsize=40)
ticks = np.arange(0,0.20,0.02)
axs1.set_yticks(ticks)
axs1.tick_params(axis='x', labelsize=40 ) #to Set Matplotlib Tick Labels Font Size
axs1.tick_params(axis='y', labelsize=40 )
axs1.set_xlabel('$ \mathrm{N_d}$ ($\mathrm{cm^{-3}}$)',fontsize=40)
axs1.set_ylabel('Relative Frequency',fontsize=40)





axs0.set_title('Inside Plume',fontsize=40)
axs1.set_title('Outside Plume',fontsize=40)

axs0.annotate('(a)',xy=(-30,0.17),size=40)
axs1.annotate('(b)',xy=(-30,0.17),size=40)
plt.tight_layout()
#axs0.grid(True)
#axs1.grid(True)

plt.xticks(fontsize = 40)
plt.yticks(fontsize = 40)
plt.savefig('Nd_No_label.png')
plt.show()
#====================================================================
#fig, axs = plt.subplots(2,2,figsize=(30,20))
#cmap = [(0.0,0.0,0.0)] + [(cm.jet(i)) for i in range(1,256)]
#cmap = mpl.colors.ListedColormap(cmap)
#bounds=np.arange(0,35,5)
#norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
#cmap_2 = LinearSegmentedColormap.from_list("", ["white","lightskyblue","steelblue","green","yellowgreen","yellow","gold","red","firebrick","darkred"])
#m = Basemap(ax=axs[0,0],projection='cyl',llcrnrlat=50,urcrnrlat=80,\
#           llcrnrlon=-60,urcrnrlon=20,resolution='c')
#m.drawcoastlines()
#x, y = m(lon_fine,lat_fine)
#axs[0,0]=m.pcolormesh(x, y,re_out, cmap=cmap_2, norm=norm)
#cax = fig.add_axes([0.15,0.05,0.7,0.02])
#cbar_bounds = bounds
#cbar_ticks =  bounds
#cbar = fig.colorbar(axs[0,0],cax=cax, norm=norm, boundaries=cbar_bounds, ticks=cbar_ticks, orientation='horizontal')

#axs[1,0]=m.pcolormesh(x, y,re_out, cmap=cmap_2, norm=norm)
#cbar.set_label('The total vertical column amount SO2 (in DU) in lower Troposphere', fontsize=14)
#plt.savefig('out.png')