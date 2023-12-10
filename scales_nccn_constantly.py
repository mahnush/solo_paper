import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import postpro_func
import numpy.ma as ma

# read in the control NCCN file
data_path = '/home/mhaghigh/solforia/nccn/run_new_version_icon/icon_sel_2021-04-10_control_test_input.nc'

# open the SO2 file
omps_path_10Apr = '/home/mhaghigh/solforia/omps/grid_files/omps_10Apr_test_domain.nc'

def so2_scale_factor(omps_path) :
    so2 = postpro_func.read_nc(omps_path, 'so2_TRL')
    lat_so2 = postpro_func.read_nc(omps_path, 'lat')
    lon_so2 = postpro_func.read_nc(omps_path, 'lon')
    so2 = np.transpose(so2)
    lon_so2 = np.transpose(lon_so2)
    lat_so2 = np.transpose(lat_so2)
    so2 = np.ma.filled(so2, fill_value = 0)
    so2_back = so2[so2 <= 1.]
    so2_plume = so2[so2 > 1.0]
    mean_plume = np.ma.mean(so2_plume)
    print(mean_plume)
    mean_back = np.ma.mean(so2_back)
    print(mean_back)
    scale_fac = so2/4.
    scale_fac[scale_fac <= 1.0] = 1.0
    scale_fac[so2 < 1.0] = 1.0
    return scale_fac, lat_so2, lon_so2

def get_control_ccn():
    nccn = postpro_func.read_nc(data_path, 'N_CCN')
    height = postpro_func.read_nc(data_path, 'LH')
    w = postpro_func.read_nc(data_path, 'w')
    lon = postpro_func.read_nc(data_path, 'lon')
    lat = postpro_func.read_nc(data_path, 'lat')
    lev = postpro_func.read_nc(data_path, 'lev')
    time = postpro_func.read_nc(data_path, 'time')
    print(np.shape(nccn))
    return nccn, w, height, lat, lon, lev, time

def interpolate_scal_factor(omps_path):
    import scipy.interpolate as sci
    ori_scal_factor = so2_scale_factor(omps_path)[0]
    ori_lat_so2 = so2_scale_factor(omps_path)[1]
    ori_lon_so2 = so2_scale_factor(omps_path)[2]
    ori_lat_so2 = ori_lat_so2[:, 0]
    ori_lon_so2 = ori_lon_so2[0, :]
    lat_fine = get_control_ccn()[3]
    lon_fine = get_control_ccn()[4]
    f = sci.RectBivariateSpline(ori_lat_so2, ori_lon_so2, ori_scal_factor )
    scale_interp = f(lat_fine, lon_fine)
    return scale_interp

def visualize(row, j, updraft, panels, avg_ccn):
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap

    font_size = 30
    cmap = [(0.0, 0.0, 0.0)] + [(cm.jet(i)) for i in range( 1 , 256 )]
    cmap = mpl.colors.ListedColormap( cmap )
    bounds = np.arange(0, 1000, 100)
    norm = mpl.colors.BoundaryNorm( bounds , cmap.N )
    cmap_2 = LinearSegmentedColormap.from_list( "" , ["lightskyblue" , "steelblue" , "green" , "yellowgreen" ,
                                                      "yellow" , "gold" , "red" , "firebrick" , "darkred"] )

    m = Basemap( ax = axs[row , j] , projection = 'cyl' , llcrnrlat = -5 , urcrnrlat = 35 , llcrnrlon = -90 ,
                 urcrnrlon = -30 , resolution = 'c' )
    m.drawcoastlines( )
    m.drawmeridians( np.arange( -180. , 180. , 5. ) , linewidth = 1.2 , labels = [0 , 0 , 0 , 1] , fontsize = 20 ,
                     color = 'black' , zorder = 3 , latmax = 90 ,
                     )
    m.drawparallels( np.arange( 0. , 85. , 5. ) , linewidth = 1.2 , labels = [1 , 0 , 0 , 0] , fontsize = 20 ,
                     color = 'black' , zorder = 3 , latmax = 90 ,
                     )

    axs[row , j].set_title( updraft , fontsize = font_size , loc = 'center' )
    axs[row , j].set_title( panels , loc = 'left' , fontsize = font_size )

    axs[row , j] = m.imshow( avg_ccn , cmap = cmap_2 , norm = norm )
    cax = fig.add_axes( [0.15 , 0.05 , 0.69 , 0.02] )  # left, bottom, width,height
    cbar_bounds = bounds
    cbar_ticks = cbar_bounds
    cbar = fig.colorbar( axs[row , j] , cax = cax , norm = norm , boundaries = cbar_bounds ,
                         ticks = cbar_ticks , orientation = 'horizontal' )
    cbar.ax.tick_params( labelsize = font_size )
    cbar.set_label( 'CCN ($\mathrm{cm^{-3}}$)' , fontsize = font_size )

    return
first_ts = int(0)
a = int(8)
b = int(8)
nt = np.size(get_control_ccn()[6])
nw = np.size(get_control_ccn()[1])
nlev = np.size(get_control_ccn()[5])
nlat = np.size(get_control_ccn()[3])
nlon = np.size(get_control_ccn()[4])
nccn_perturbed = np.zeros((nt, nw, nlev, nlat, nlon))

for t in range (nt):
    nccn_control = get_control_ccn()[0]
    nccn_control = nccn_control[t, :, :, :, :]
    scale_factor = so2_scale_factor(omps_path_10Apr)[0]
    for k in range(60):
        for w in range(10) :
            nccn_perturbed[t, w, k, :, :] = nccn_control[w, k, :, :]*scale_factor[:, :]


w = get_control_ccn()[1]
print(np.max(nccn_perturbed)*1e-6)
print(np.max(get_control_ccn()[0]*1e-6))
z_avg = get_control_ccn()[2]
lat = get_control_ccn()[3]
lon = get_control_ccn()[4]
lev = get_control_ccn()[5]
time = get_control_ccn()[6]


ncout = Dataset('/home/mhaghigh/solforia/nccn/run_new_version_icon/icon_perturbed_10Apr.nc', mode="w",
                format='NETCDF4_CLASSIC')
ncout.description = 'ICON_CCN_file'
ncout.createDimension('lon', nlon)
ncout.createDimension('lat', nlat)
ncout.createDimension('lev', nlev)
ncout.createDimension('w', nw)
ncout.createDimension('time', nt)
lon_o = ncout.createVariable('lon', np.float32, ('lon',))
lat_o = ncout.createVariable('lat', np.float32, ('lat',))
lev_o = ncout.createVariable('lev', np.float32, ('lev',))
W = ncout.createVariable('w', np.float32, ('w',))
time_o = ncout.createVariable('time', np.float32, ('time',))
N_CCN = ncout.createVariable('N_CCN', np.float32, ('time', 'w', 'lev', 'lat', 'lon'))
LH = ncout.createVariable('LH', np.float32, ('lev', 'lat', 'lon'))

lon_o.long_name = 'longitude'
lat_o.long_name = 'lattitude'
time_o.long_name = 'time'
N_CCN.long_name = 'CCN number concentration'
W.long_name = 'vertical velocity'
LH.long_name = 'top height of grid cell'
lev_o.long_name = 'level number'
lat_o.units = 'degrees_north'
lon_o.units = 'degrees_east'
time_o.units = 'hours since 2021-04-10T00:00:00'
LH.units = 'meters'
W.units = 'meters per second'
N_CCN.units = 'particles per cubic meter'


lon_o[:] = lon[:]
lat_o[:] = lat[:]
lev_o[:] = lev[:]
W[:] = w[:]
time_o[:] = time[:]
LH[:] = z_avg[:]
N_CCN[:] = nccn_perturbed[:]


panels_1 = ['perturbed', 'control', 'c', 'd']
panel_1 = ['w = 0.215m/s' , 'w = 4.64 m/s']

nccn_perturbed_lev_mean_22dec = np.mean(nccn_perturbed[0, 3, :, :, :],axis= 0)*1e-6
nccn_control = get_control_ccn()[0]
nccn_control_lev_mean_22dec = np.mean(nccn_control[0, 3, :, :, :],axis= 0)*1e-6
print(np.shape(nccn_perturbed_lev_mean_22dec))
fig , axs = plt.subplots(2, 2, figsize = (30, 20))
visualize(0, 0, panel_1[0], panels_1[0], nccn_perturbed_lev_mean_22dec)
visualize(0, 1, panel_1[0], panels_1[1], nccn_control_lev_mean_22dec)
nccn_perturbed_lev_mean_22dec = np.mean(nccn_perturbed[0, 6, :, :, :],axis= 0)*1e-6
nccn_control = get_control_ccn()[0]
nccn_control_lev_mean_22dec = np.mean(nccn_control[0, 6, :, :, :],axis= 0)*1e-6
visualize(1, 0, panel_1[1], panels_1[0], nccn_perturbed_lev_mean_22dec)
visualize(1, 1, panel_1[1], panels_1[1], nccn_control_lev_mean_22dec)
plt.savefig('per_con_nccn_10Apr.png')
plt.show()