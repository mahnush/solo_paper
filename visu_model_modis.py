import matplotlib.pyplot as plt
import numpy as np
import postpro_func
file_name_modis = '/home/mhaghigh/solforia/modis/MODIS_grid_data/terra_modis_solo_2km_icon_domain_10Apr.nc'
#file_name_model = '/home/mhaghigh/solforia/model_run/nc/control_run_10_2data.nc'
file_name_model = '/home/mhaghigh/solforia/pdf/con_half_mean.nc'
var_name_modis = ['nd','lwp', 'tau', 're']
var_name_model = ['nd_dw','lwp_dw', 'tau', 'reff']
titel_var = ['Nd', 'LWP', 'tau', 'reff']
titel_kind = [' (MODIS)', ' (ICON)']
cbar_label = ['$\mathrm{cm^{-3}}$', '$\mathrm{g\,m^{-2}}$', '','$\mathrm{{\mu}m}$']
fs_label = 20
fs_titel = 20
nd_bound = np.arange(0, 201, 1)
re_bound = np.arange(0, 32, 2)
lwp_bound = np.arange(0, 202, 2)
tau_bound = np.arange(0, 110,10)
limit = np.zeros(4)
limit[0] = -1
limit[1] = 29
limit[2] = -85
limit[3] = -35
titel = '10 April 2021'
import glob
files = glob.glob('/home/mhaghigh/solforia/new_version_run/old_nccn_all_data/lwp/*.nc', recursive=True)
def read_nc_mulitple(files, var_name):
    from netCDF4 import Dataset
    import netCDF4
    #nc = Dataset(file, 'r')
    nc = netCDF4.MFDataset(files)
    var = nc.variables[var_name][:, :, :]
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    var_mean = np.mean(var, axis= 0)
    var_mean = var_mean *1e3
    return var_mean
re_modis = postpro_func.read_nc(file_name_modis, var_name_modis[0])
#re_model = read_nc_mulitple(files, var_name_model[1])
re_model = postpro_func.read_nc(file_name_model, var_name_model[0])
lat_model = postpro_func.read_nc(file_name_model, 'lat')[:]
lon_model = postpro_func.read_nc(file_name_model, 'lon')[:]
lat_modis = postpro_func.read_nc(file_name_modis, 'lat')
lon_modis = postpro_func.read_nc(file_name_modis, 'lon')

#print(lat_modis)
#print(lon_modis)
fig, ((ax0,ax1)) = plt.subplots(2, 1, figsize=(17, 12))
fig.suptitle(titel, fontsize=fs_titel)
postpro_func.visulize_sat(ax0, re_modis, lat_modis , lon_modis , lwp_bound, cbar_label[0],
 titel_kind[0], limit)
#2end subplot
postpro_func.visulize_model(ax1, re_model, lwp_bound, cbar_label[0], titel_kind[1], limit)
#lat = postpro_func.read_nc(file_name_model,'lat')
#print(lat)
#3rd subplot
#postpro_func.visulize_sat(ax2, postpro_func.read_nc(file_name, var_name[0]),
#postpro_func.read_nc(file_name,'lat') , postpro_func.read_nc(file_name,'lon') , nd_bound, cbar_label[0],
#titel_var[0] + titel_kind[0], limit)
#4th subplot
#postpro_func.visulize_sat(ax3, postpro_func.read_nc(file_name, var_name[3]),
#postpro_func.read_nc(file_name,'lat') , postpro_func.read_nc(file_name,'lon') , re_bound, cbar_label[3],
#titel_figure[3], limit)
#postpro_func.visulize_model(ax3, postpro_func.read_nc(file_name_model, var_name_model[0])
#, nd_bound, cbar_label[0],titel_var[0] + titel_kind[1], limit)
plt.savefig('lwp_old_icon_version_17&18UTC.png')
#plt.savefig('28Dec.pdf')
plt.show()