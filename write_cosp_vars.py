import numpy as np
from netCDF4 import Dataset

ncout = Dataset('/home/mhaghigh/solforia/new_version_run/new_ccn_all_data_per_second_run/allvars/pwwerturbed_second_run_cosp.nc'
                , mode="w", format='NETCDF4_CLASSIC')
FILE_NAME = '/home/mhaghigh/solforia/new_version_run/new_ccn_all_data_per_second_run/allvars/perturbed_second_run.nc'
nlon = 2501
nlat = 1501
ns = 23
nc_dw = Dataset(FILE_NAME, 'r')
re_dw = nc_dw.variables['reff'][:, :, :]
tau_dw = nc_dw.variables['tau'][:, :, :]
lwp_dw = nc_dw.variables['lwp'][:, :, :]

#re_mask = np.ma.masked_where(re_dw == 0, re_dw)
nd_dw = 1.37e-5*(tau_dw**0.5)*(re_dw**-2.5)*1e-6

lon = nc_dw.variables['lon'][:]
lat = nc_dw.variables['lat'][:]

ncout.createDimension('lon', nlon)
ncout.createDimension('lat', nlat)
ncout.createDimension('time', ns)
lon_o = ncout.createVariable('lon', np.float32, ('lon',))
lat_o = ncout.createVariable('lat', np.float32, ('lat',))
re_dw_mean_o = ncout.createVariable('reff', np.float32, ('time','lat','lon'))
tau_dw_mean_o = ncout.createVariable('tau', np.float32, ('time','lat','lon'))
lwp_dw_mean_o = ncout.createVariable('lwp', np.float32, ('time', 'lat','lon'))
qnc_dw_mean_o = ncout.createVariable('nd', np.float32, ('time', 'lat','lon'))
lon_o[:] = lon[:]
lat_o[:] = lat[:]
re_dw_mean_o[:] = re_dw[:]
tau_dw_mean_o[:] = tau_dw[:]
lwp_dw_mean_o[:] = lwp_dw[:]
qnc_dw_mean_o[:] = nd_dw[:]

print('done')