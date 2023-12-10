#this code creates a geogrophical distribusion  MODIS
#different variables can be chossen but here the Nd and LWP are plotted.
import matplotlib.pyplot as plt
import numpy as np
import postpro_func
#necessary inputs
file_name  = '/home/mhaghigh/solforia/modis/MODIS_grid_data/aqua_modis_solo_2km_icon_domain_10Apr.nc'
var_name = ['nd','lwp', 'tau', 're']
var_name_model = ['nd_dw','lwp_dw', 'tau_dw', 're_dw']
titel_var = ['Nd', 'LWP', 'tau', 'reff']
titel_kind = [' (MODIS)', ' (ICON)']
cbar_label = ['$\mathrm{cm^{-3}}$', '$\mathrm{g\,m^{-2}}$', '','$\mathrm{{\mu}m}$']
fs_label = 20
fs_titel = 20
nd_bound = np.arange(0,201,1)
re_bound = np.arange(0,62,2)
lwp_bound = np.arange(0,201,2)
tau_bound = np.arange(0,100,1)
limit = np.zeros(4)
limit[0] = -1
limit[1] = 29
limit[2] = -85
limit[3] = -35
gridSize = 0.02
#limit[0] = 5
#limit[1] = 22
#limit[2] = -70
#limit[3] = -35
#limit[0] = 0
#limit[1] = 20
#limit[2] = -60
#limit[3] = -20
#gridSize = 0.02
#gridSize = 0.01
titel = '10 APR 2021'

lat = postpro_func.read_nc_sat(file_name,'lat')[:,0]
lon = postpro_func.read_nc_sat(file_name,'lon')[0,:]
re = postpro_func.read_nc_sat(file_name, var_name[3])
print(lon)
print(lat)
print(np.shape(lat))
re_o = postpro_func.mask_land(lat,lon,re)
fig, ((ax0,ax1), (ax2,ax3)) = plt.subplots(2, 2, figsize=(17, 12))
fig.suptitle(titel, fontsize=fs_titel)
postpro_func.visulize_sat(ax0, re_o, lat , lon, re_bound, cbar_label[3],
titel_var[3] + titel_kind[0], limit)
#2end subplot
#postpro_func.visulize_sat(ax1, postpro_func.read_nc(file_name, var_name[2]),
#postpro_func.read_nc(file_name,'lat') , postpro_func.read_nc(file_name,'lon') , tau_bound, cbar_label[2],
#titel_var[2] + titel_kind[0], limit)

#3rd subplot
#postpro_func.visulize_sat(ax2, postpro_func.read_nc(file_name, var_name[0]),
#postpro_func.read_nc(file_name,'lat') , postpro_func.read_nc(file_name,'lon') , nd_bound, cbar_label[0],
#titel_var[0] + titel_kind[0], limit)
#4th subplot
#postpro_func.visulize_sat(ax3, postpro_func.read_nc(file_name, var_name[1]),
#postpro_func.read_nc(file_name,'lat') , postpro_func.read_nc(file_name,'lon') , lwp_bound, cbar_label[1],
#titel_var[1] + titel_kind[0], limit)
#plt.tight_layout()
#plt.savefig('aqua_MODIS_VARS_10APR_2km_icon_domain.png')
#plt.savefig('28Dec.pdf')
plt.show()