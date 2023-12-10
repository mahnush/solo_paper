import matplotlib.pyplot as plt
import postpro_func
import numpy as np
import numpy.ma as ma

modis_10Apr = '/home/mhaghigh/solforia/modis/MODIS_grid_data/aqua_temperature_modis_solo_2km_icon_domain_10Apr.nc'

model_con_10Apr = '/home/mhaghigh/solforia/new_version_run/new_ccn_all_data/all_var/control.nc'


var_name_model = ['nd', 'lwp', 'tau', 'reff']
var_name_modis = ['nd', 'lwp', 'tau', 're']



control_nd = postpro_func.read_nc(model_con_10Apr, var_name_model[1])[:,:,:]*1000
control_nd = np.ma.mean(control_nd, axis=0)
re_modis = postpro_func.read_nc(modis_10Apr, var_name_modis[1])

lat_modis = postpro_func.read_nc(modis_10Apr, 'lat')
lon_modis = postpro_func.read_nc(modis_10Apr, 'lon')

#print(lat_modis)
#print(lon_modis)
itel_var = ['Nd', 'LWP', 'tau', 'reff']
titel_kind = [' (MODIS)', ' (ICON)']
cbar_label = ['$\mathrm{cm^{-3}}$', '$\mathrm{g\,m^{-2}}$', '','$\mathrm{{\mu}m}$']
fs_label = 20
fs_titel = 20
nd_bound = np.arange(0, 5000, 10)
re_bound = np.arange(0, 32, 2)
lwp_bound = np.arange(0, 202, 2)
tau_bound = np.arange(0, 110,10)
limit = np.zeros(4)
limit[0] = -1
limit[1] = 29
limit[2] = -85
limit[3] = -35
titel = '10 April 2021'
fig, ((ax0,ax1)) = plt.subplots(2, 1, figsize=(17, 12))
fig.suptitle(titel, fontsize=fs_titel)
postpro_func.visulize_sat(ax0, re_modis, lat_modis , lon_modis , lwp_bound, cbar_label[1],
 titel_kind[0], limit)
#2end subplot
postpro_func.visulize_model(ax1, control_nd, lwp_bound, cbar_label[1], titel_kind[1], limit)
plt.show()

re_modis_1dim = re_modis.flatten()
re_modis_1dim_n = re_modis_1dim[re_modis_1dim >0.]
re_modis_1dim_c = re_modis_1dim_n.compressed()
control_nd_mean_1dim = control_nd.flatten()
control_nd_mean_1dim_n = control_nd_mean_1dim[control_nd_mean_1dim > 5.]
control_nd_mean_1dim_c = control_nd_mean_1dim_n.compressed()

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
numbin = np.arange(0, 400, 1)
#numbin = np.arange(0, 500, 10)
#numbin = np.arange(0,200,5)
font_tick = 30
font_legend = 30
font_lable = 30
line_width = 4
font_tick = 20
#name = '$ \mathrm{N_d}$ ($\mathrm{cm^{-3}}$)'
name =  '$ \mathrm{LWP}$ ($\mathrm{g m^{-2}}$)'
#name = '\u03BC'+'m'
#name = ''
#axs0.hist(,bins=numbin, weights=weight(var_per_test) , histtype='step',
#         linewidth=line_width, color='red', label='perturbed '+lable_hist(var_per_test))

axs0.hist(control_nd_mean_1dim_c,  bins=numbin, weights=weight(control_nd_mean_1dim_c),  histtype='step',
         linewidth=line_width, color='blue', label='control '+lable_hist(control_nd_mean_1dim_c))
axs0.hist(re_modis_1dim_c, bins=numbin, weights=weight(re_modis_1dim_c),  histtype='step',
         linewidth=line_width, color='black',  label='MODIS '+lable_hist(re_modis_1dim_c))
axs0.legend(loc='upper right', fontsize=font_legend, frameon=True)
#ticks = np.arange(0, 0.12, 0.02)
#ticks_x = np.arange(0,1200,200)
#axs0.set_yticks(ticks)
axs0.tick_params(axis='x', labelsize=font_tick)  # to Set Matplotlib Tick Labels Font Size
axs0.tick_params(axis='y', labelsize=font_tick)
axs0.set_xlabel(name, fontsize=font_lable)
axs0.set_ylabel('Relative Frequency', fontsize=font_lable)

axs0.set_title('LWP', fontsize= font_lable)
# axs0.annotate('(a)',xy=(-30,0.079),size=font_lable)
# plt.tight_layout()
axs0.grid(True)
plt.show()
plt.savefig('tqc_version_modis_geo.png')


