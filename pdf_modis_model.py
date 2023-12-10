import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('TkAgg')
import glob
files = glob.glob('/home/mhaghigh/solforia/new_version_run/old_nccn_all_data/re/*.nc', recursive=True)

def read_nc_mulitple(files, var_name):
    from netCDF4 import Dataset
    import netCDF4
    #nc = Dataset(file, 'r')
    nc = netCDF4.MFDataset(files)
    var = nc.variables[var_name][:, :, :]
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    #var_mean = np.ma.mean(var, axis= 0)
    #print(np.shape(var_mean))
    return var, lat, lon
def read_nc(file_name, var_name):
    from netCDF4 import Dataset
    nc = Dataset(file_name, 'r')
    var = nc.variables[var_name][0, :, :]
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    return var, lat, lon
def read_nc_modis(file_name, var_name):
    from netCDF4 import Dataset
    file = file_name
    nc = Dataset(file, 'r')
    var = nc.variables[var_name][:, :]
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    return var, lat, lon
#def get_mean_model_var(files, var_name) :
#    var = read_nc_mulitple(files, var_name )

files_per = '/home/mhaghigh/solforia/new_version_run/new_ccn_all_data_per/all_var/perturbed.nc'
files_con = '/home/mhaghigh/solforia/new_version_run/new_ccn_all_data/all_var/control.nc'
files_modis = '/home/mhaghigh/solforia/modis/MODIS_grid_data/aqua_modis_solo_1km_10Apr.nc'
var_name = 'reff'
var_name_modis = 're'
#data = [files_per, files_con]

#var_con, lat, lon = read_nc_mulitple(files, var_name)
var_con, lat, lon = read_nc(files_con, var_name)
var_per, lat, lon = read_nc(files_per, var_name)
var_modis, lat, lon = read_nc_modis(files_modis, var_name_modis)
var_modis = np.transpose(var_modis)
var_con = var_con[:, :]*1e6
var_per = var_per[:, :]*1e6
var_modis = var_modis[:, :]


print(np.shape(var_modis))
print(np.shape(var_con))
#var_con = var_con[250:1151, 450:1651]
var_con_test = var_con[var_con > 4]
var_per_test = var_per[var_per > 4]
var_modis_test = var_modis[var_modis > 4.]
var_modis_test = var_modis_test.compressed()
print(np.min(var_modis_test))
print(np.shape(var_modis_test))
print(np.shape(var_con_test))
def weight(var):
    #weight = (1 + np.zeros(len(var))) / len(var)
    weight = np.zeros_like(var) + 1. / (var.size)
    return weight
def lable_hist(var):
    median = str(np.median(var))
    mean = str((np.mean(var)))

    std = str(np.std(var))
    lable = '('+'median = ' + mean +')'
    return lable
fig, axs0 = plt.subplots(figsize= (20, 15))
numbin = np.arange(0, 31, 1)
#numbin = np.arange(0, 500, 10)
#numbin = np.arange(0,200,5)
font_tick = 30
font_legend = 30
font_lable = 30
line_width = 4
font_tick = 20
#name = '$ \mathrm{N_d}$ ($\mathrm{cm^{-3}}$)'
#name =  '$ \mathrm{LWP}$ ($\mathrm{g m^{-2}}$)'
name = '\u03BC'+'m'
#name = ''
axs0.hist(var_per_test,bins=numbin, weights=weight(var_per_test) , histtype='step',
         linewidth=line_width, color='red', label='perturbed '+lable_hist(var_per_test))

axs0.hist(var_con_test,  bins=numbin, weights=weight(var_con_test),  histtype='step',
         linewidth=line_width, color='blue', label='control '+lable_hist(var_con_test))
axs0.hist(var_modis_test, bins=numbin, weights=weight(var_modis_test),  histtype='step',
         linewidth=line_width, color='black',  label='MODIS '+lable_hist(var_modis_test))
axs0.legend(loc='upper right', fontsize=font_legend, frameon=True)
#ticks = np.arange(0, 0.12, 0.02)
#ticks_x = np.arange(0,1200,200)
#axs0.set_yticks(ticks)
axs0.tick_params(axis='x', labelsize=font_tick)  # to Set Matplotlib Tick Labels Font Size
axs0.tick_params(axis='y', labelsize=font_tick)
axs0.set_xlabel(name, fontsize=font_lable)
axs0.set_ylabel('Relative Frequency', fontsize=font_lable)

axs0.set_title('Reff', fontsize= font_lable)
# axs0.annotate('(a)',xy=(-30,0.079),size=font_lable)
# plt.tight_layout()
axs0.grid(True)
plt.show()
plt.savefig('levante_model_con_per_modis_re_10Apr.png')


