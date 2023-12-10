from pyhdf.SD import SD, SDC
import numpy as np
from netCDF4 import Dataset
import MODIS_functions
import numpy.ma as ma

# import modis_postprocess

# --set the limit for Lat and Lon and grid size
limit = np.zeros(4)
#limit[0] = 6
#limit[1] = 34
limit[0] = -1
limit[1] = 29
limit[2] = -85
limit[3] = -35
gridSize = 0.02
#gridSize = 0.01
# --opening L2 and geo file list
fileList = open('/home/mhaghigh/solforia/modis/modis_data/10apr/terra_file_list.txt', 'r+')
fileList_geo = open('/home/mhaghigh/solforia/modis/modis_data/10apr/terra_g_file_list.txt', 'r+')
geo_file = [line for line in fileList_geo.readlines()]

latgrid, longrid,nlon, nlat= MODIS_functions.grid_coordinate(limit, gridSize)
# --defining the ar nray to put all lat and lon and data
allLat = []
allLon = []
allreff = []
alltau = []
allctp = []
alllwp = []
allphase = []

ipath = '/home/mhaghigh/solforia/modis/modis_data/10apr/'
# name of variables
reff = 'Cloud_Effective_Radius'
tau = 'Cloud_Optical_Thickness'
ctp = 'cloud_top_temperature_1km'
optical_phase = 'Cloud_Phase_Optical_Properties'
lwp = 'Cloud_Water_Path'
# --reading the data on data set
k = 0

for FILE_NAME in fileList:
    FILE_NAME = FILE_NAME.strip()
    file_gn = geo_file[k]
    file_gn = file_gn.strip()
    # print(k)
    FILE_NAME = ipath + FILE_NAME
    print(FILE_NAME)
    file_gn = ipath + file_gn
    print(file_gn)
    # set the modis file L2 & geo
    file_i = SD(FILE_NAME, SDC.READ)
    file_g = SD(file_gn, SDC.READ)
    if len(MODIS_functions.read_cloud_var(file_i, reff)) == 0:
        allLat, allLon = MODIS_functions.read_coordinate(file_g)
        allreff = MODIS_functions.read_cloud_var(file_i, reff)
        alltau = MODIS_functions.read_cloud_var(file_i, tau)
        allctp = MODIS_functions.read_cloud_var(file_i, ctp)
        allphase = MODIS_functions.read_cloud_var(file_i, optical_phase)
        alllwp = MODIS_functions.read_cloud_var(file_i, lwp)
    elif len(MODIS_functions.read_cloud_var(file_i, reff)) > 0:
        allLat = np.concatenate((allLat, MODIS_functions.read_coordinate(file_g)[0]), axis=0)
        allLon = np.concatenate((allLon, MODIS_functions.read_coordinate(file_g)[1]), axis=0)
        allreff = np.concatenate((allreff, MODIS_functions.read_cloud_var(file_i, reff)), axis=0)
        alltau = np.concatenate((alltau, MODIS_functions.read_cloud_var(file_i, tau)), axis=0)
        allctp = np.concatenate((allctp, MODIS_functions.read_cloud_var(file_i, ctp)), axis=0)
        allphase = np.concatenate((allphase, MODIS_functions.read_cloud_var(file_i, optical_phase)), axis=0)
        alllwp = np.concatenate((alllwp, MODIS_functions.read_cloud_var(file_i, lwp)), axis=0)
    k = k + 1

reff_grid = MODIS_functions.grid(limit, gridSize, allreff, allLat, allLon)
lwp_grid = MODIS_functions.grid(limit, gridSize, alllwp, allLat, allLon)
phase_grid = MODIS_functions.grid(limit, gridSize, allphase, allLat, allLon)
tau_grid = MODIS_functions.grid(limit, gridSize, alltau, allLat, allLon)
ctp_grid = MODIS_functions.grid(limit, gridSize, allctp, allLat, allLon)
# mask ice phase clouds
reff_grid = ma.masked_where(phase_grid == 3, reff_grid)
lwp_grid = ma.masked_where(phase_grid == 3, lwp_grid)
tau_grid = ma.masked_where(phase_grid == 3, tau_grid)
#reff_grid_test = reff_grid[ctp_grid>273]
#print(np.max(reff_grid_test))
# to compute Nd
nd_grid = (1.37e-5 * (tau_grid ** 0.5) * ((reff_grid * 1e-6) ** -2.5)) * 1e-6

#plt.savefig('nd_warm.pdf')
#plt.show()

ncout = Dataset('/home/mhaghigh/solforia/modis/MODIS_grid_data/terra_modis_solo_2km_icon_domain_10Apr.nc', mode="w", format='NETCDF4_CLASSIC')
#nlon = 1801
#nlat = 1401

ncout.createDimension('lat', nlat)
ncout.createDimension('lon', nlon)
lat_o = ncout.createVariable('lat', np.float32, ('lon', 'lat'))
lon_o = ncout.createVariable('lon', np.float32, ('lon', 'lat'))
re_o = ncout.createVariable('re', np.float32, ('lon', 'lat'))
tau_o = ncout.createVariable('tau', np.float32, ('lon', 'lat'))
lwp_o = ncout.createVariable('lwp', np.float32, ('lon', 'lat'))
nd_o = ncout.createVariable('nd', np.float32, ('lon', 'lat'))

re_o[:] = reff_grid[:]
tau_o[:] = tau_grid[:]
lwp_o[:] = lwp_grid[:]
nd_o[:] = nd_grid[:]
lat_o[:] = latgrid[:]
lon_o[:] = longrid[:]

print('done')