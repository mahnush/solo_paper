import numpy as np
from netCDF4 import Dataset
import h5py
import MODIS_functions
import postpro_func
import matplotlib.pyplot as plt
#change the domain 8.8.2022
limit = np.zeros(4)
#limit[0] = 0
#limit[1] = 35
#limit[2] = -80
#limit[3] = 20
#gridSize = 1.0
limit[0] = -1
limit[1] = 29
limit[2] = -85
limit[3] = -35
gridSize = 1
scientific_data = 'SCIENCE_DATA'
so2_PBL = 'ColumnAmountSO2_PBL'
so2_TRL = 'ColumnAmountSO2_TRL'
so2_TRM = 'ColumnAmountSO2_TRM'
so2_STL = 'ColumnAmountSO2_STL'
geolocation_data = 'GEOLOCATION_DATA'
latitute = 'Latitude'
longitude = 'Longitude'


def read_data(first_data_set , second_data_set):
     hf_in = h5py.File( FILE_NAME , 'r' )
     #print(list(hf_in.keys()))
     #print(list(hf_in['SCIENCE_DATA']))
     read_var = hf_in[first_data_set][second_data_set]
     temp_var = read_var[:]
     temp_var = temp_var.flatten( )
     return temp_var


file_path = '/home/mhaghigh/solforia/omps/10Apr/'
fileList = open( file_path + 'file_list.txt' , 'r+' )
allLat = []
allLon = []
allso2_PBL = []
allso2_TRL = []
allso2_TRM = []
allso2_STL = []

for FILE_NAME in fileList :
    FILE_NAME = FILE_NAME.strip( )
    FILE_NAME = file_path + FILE_NAME
    print( FILE_NAME )

    if len( allso2_PBL ) == 0 :
        allso2_PBL = read_data( scientific_data, so2_PBL)
        allso2_TRL = read_data( scientific_data, so2_TRL)
        allso2_TRM = read_data( scientific_data, so2_TRM)
        allso2_STL = read_data( scientific_data, so2_STL)
        allLat = read_data( geolocation_data, latitute)
        allLon = read_data( geolocation_data, longitude)
    elif len( allso2_PBL ) > 0 :
        allso2_PBL = np.concatenate( (allso2_PBL , read_data( scientific_data , so2_PBL )) , axis = 0 )
        allso2_TRL = np.concatenate( (allso2_TRL , read_data( scientific_data , so2_TRL )) , axis = 0 )
        allso2_TRM = np.concatenate( (allso2_TRM , read_data( scientific_data , so2_TRM )) , axis = 0 )
        allso2_STL = np.concatenate( (allso2_STL , read_data( scientific_data , so2_STL )) , axis = 0 )
        allLat = np.concatenate( (allLat , read_data( geolocation_data , latitute )) , axis = 0 )
        allLon = np.concatenate( (allLon , read_data( geolocation_data , longitude )) , axis = 0 )

so2_PBL_grid = MODIS_functions.grid( limit , gridSize , allso2_PBL , allLat , allLon )
so2_TRL_grid = MODIS_functions.grid( limit , gridSize , allso2_TRL , allLat , allLon )
so2_TRM_grid = MODIS_functions.grid( limit , gridSize , allso2_TRM , allLat , allLon )
so2_STL_grid = MODIS_functions.grid( limit , gridSize , allso2_STL , allLat , allLon )
lat_grid , lon_grid, nlon, nlat = MODIS_functions.grid_coordinate( limit , gridSize )
print(np.shape(so2_PBL_grid))
#so2_TRL_grid[0:150,0:15] = 0.0
#so2_PBL_grid[:, 0:30] = 0.0

fs_titel = 12
titel = '10 April 2021'
bounds = np.arange(0, 30, 1)
cbar_label = 'So2 amount (DU)'
title = ['Boundary Layer', 'Lower Troposphere(3 km)', 'Middle Troposphere (8 km)', 'Stratosphere(18 km)']
fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2, figsize = (8, 6), sharey = True, sharex = True)
#fig.suptitle( titel , fontsize = fs_titel )
cax = fig.add_axes([0.25, 0.08, 0.5, 0.02])  # left, bottom, width,height
postpro_func.visulize_sat(ax0, so2_PBL_grid, lat_grid, lon_grid, bounds, cbar_label, title[0], limit, cax)
postpro_func.visulize_sat(ax1, so2_TRL_grid, lat_grid, lon_grid, bounds, cbar_label, title[1], limit, cax)
postpro_func.visulize_sat(ax2, so2_TRM_grid, lat_grid, lon_grid, bounds, cbar_label, title[2], limit, cax)
postpro_func.visulize_sat(ax3, so2_STL_grid, lat_grid, lon_grid, bounds, cbar_label, title[3], limit, cax)
#plt.tight_layout()
plt.savefig('so2_concentration_10pr_icondomain_figure1.png')
plt.savefig('so2_concentration_10pr_icondomain_figure1.pdf')
plt.show()

#to plot it in g/m2
fig , ((ax0 , ax1) , (ax2 , ax3)) = plt.subplots( 2 , 2 , figsize = (17 , 12) )
fig.suptitle( titel , fontsize = fs_titel )
so2_PBL_grid_geram = 0.0285*so2_PBL_grid
so2_TRL_grid_geram = 0.0285*so2_TRL_grid
so2_TRM_grid_geram = 0.0285*so2_TRM_grid
so2_STL_grid_geram = 0.0285*so2_STL_grid
bounds = np.arange(0,1.1,1e-1)
cbar_label = 'So2 amount (g m-2)'
postpro_func.visulize_sat(ax0, so2_PBL_grid_geram, lat_grid, lon_grid, bounds, cbar_label, title[0], limit, cax)
postpro_func.visulize_sat(ax1, so2_TRL_grid_geram, lat_grid, lon_grid, bounds, cbar_label, title[1], limit, cax)
postpro_func.visulize_sat(ax2, so2_TRM_grid_geram, lat_grid, lon_grid, bounds, cbar_label, title[2], limit, cax)
postpro_func.visulize_sat(ax3, so2_STL_grid_geram, lat_grid, lon_grid, bounds, cbar_label, title[3], limit, cax )
#plt.savefig( 'so2_concentration_10Apr_newdomain_DU.png' )
plt.show()
#plt.savefig('./outputs/so2_concentration_21DEC_new_domain.pdf')
# to write OMPS data on netcdf file
ncout = Dataset( '/home/mhaghigh/solforia/omps/grid_files/paper_omps_10Apr_icon_domain.nc' , mode = "w" ,
                 format = 'NETCDF4_CLASSIC' )
#nlon = 52
#nlat = 41
print(nlon, nlat)

ncout.createDimension('lat', nlat )
ncout.createDimension('lon', nlon )
lat_o = ncout.createVariable('lat', np.float32 , ('lon' , 'lat') )
lon_o = ncout.createVariable( 'lon' , np.float32 , ('lon' , 'lat') )
so2_PBL_o = ncout.createVariable( 'so2_PBL' , np.float32 , ('lon' , 'lat') )
so2_TRL_o = ncout.createVariable( 'so2_TRL' , np.float32 , ('lon' , 'lat') )
so2_TRM_o = ncout.createVariable( 'so2_TRM' , np.float32 , ('lon' , 'lat') )
so2_STL_o = ncout.createVariable( 'so2_STL' , np.float32 , ('lon' , 'lat') )

so2_PBL_o[:] = so2_PBL_grid[:]
so2_TRL_o[:] = so2_TRL_grid[:]
so2_TRM_o[:] = so2_TRM_grid[:]
so2_STL_o[:] = so2_STL_grid[:]
lat_o[:] = lat_grid[:]
lon_o[:] = lon_grid[:]

print('do')