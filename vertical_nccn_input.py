import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import postpro_func
import numpy.ma as ma
from scipy.signal import savgol_filter

data_path = '/home/mhaghigh/solforia/nccn/run_new_version_icon/icon_per_10Apr.nc'
data_path_check_input = '/home/mhaghigh/solforia/nccn/run_new_version_icon/icon_sel_2021-04-10_control_test_input.nc'

def get_pressure():
    import math
    nlev = 60
    p = np.zeros((nlev))
    scale_lev = np.zeros((nlev))
    hyam = [10, 29.2126712799072, 51.0365734100342, 79.6423835754395,
        115.060134887695, 157.533828735352, 207.681701660156, 266.637420654297,
        336.233856201172, 419.295028686523, 520.134567260742, 644.434539794922,
        798.439300537109, 989.247619628906, 1225.65466308594, 1518.55743408203,
        1881.45709228516, 2331.08129882812, 2888.15515136719, 3578.35656738281,
        4433.5, 5462.36401367188, 6662.32543945312, 8035.84252929688,
        9570.59033203125, 11226.7866210938, 12926.3857421875, 14577.5654296875,
        16099.6401367188, 17432.3291015625, 18536.439453125, 19391.40234375,
        19988.6572265625, 20326.0341796875, 20407.171875, 20240.94140625,
        19840.8662109375, 19224.5400390625, 18413.0537109375, 17430.4130859375,
        16302.9580078125, 15058.7856445312, 13727.1713867188, 12337.9887695312,
        10921.1298828125, 9505.9287109375, 8120.57983398438, 6791.55908203125,
        5543.04663085938, 4396.34582519531, 3369.30493164062, 2475.73815917969,
        1724.84619140625, 1120.63717651367, 661.347671508789, 338.863739013672,
        138.141567230225, 36.6284935474396, 3.68387150764465, 0 ]
    hybm = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        3.79117482225411e-05, 0.000268609241175, 0.00113827553286683,
        0.00344813754782081, 0.0081120147369802, 0.0159103930927813,
        0.0273995194584131, 0.0429057851433754, 0.0626121200621128,
        0.0866042636334896, 0.114848602563143, 0.147203415632248,
        0.183430127799511, 0.223204538226128, 0.266128048300743,
        0.311738923192024, 0.359523519873619, 0.408927544951439,
        0.459367260336876, 0.510240748524666, 0.560939162969589,
        0.610857933759689, 0.659408032894135, 0.706027209758759,
        0.750191211700439, 0.79142501950264, 0.829314172267914,
        0.863515913486481, 0.893770396709442, 0.919912099838257,
        0.941880911588669, 0.959733366966248, 0.973653972148895,
        0.983966410160065, 0.991144776344299, 0.995824784040451, 0.998815059661865]
    p_s = np.zeros((nlev))+1.e5
    p[:] = hybm[:] * p_s[:] +  hyam[:]
    p = p[::-1]
    return (p)


def get_control_ccn(data_path):
    nccn = postpro_func.read_nc(data_path, 'N_CCN')
    height = postpro_func.read_nc(data_path, 'LH')
    w = postpro_func.read_nc(data_path, 'w')
    lon = postpro_func.read_nc(data_path, 'lon')
    lat = postpro_func.read_nc(data_path, 'lat')
    lev = postpro_func.read_nc(data_path, 'lev')
    time = postpro_func.read_nc(data_path, 'time')
    return nccn, w, height, lat, lon, lev, time


w = get_control_ccn(data_path)[1]
z_avg = get_control_ccn(data_path)[2]
lat = get_control_ccn(data_path)[3]
lon = get_control_ccn(data_path)[4]
lev = get_control_ccn(data_path)[5]
time = get_control_ccn(data_path)[6]
nt = np.size(time)
nlev = np.size(lev)
nlon = np.size(lon)
nlat = np.size(lat)
nw = np.size(w)

def get_specific_time_w(data_path,  date, w ) :
    nccn = get_control_ccn(data_path)[0]
    ts_1 = date*8
    ts_2 = ts_1 + 8
    nccn_daily = np.mean(nccn[ts_1:ts_2, w, :, :, :], axis = 0)
    w_values = get_control_ccn(data_path)[1]
    nccn_daily_lat = np.mean(nccn_daily, axis = 1)
    nccn_daily_lon = np.mean(nccn_daily_lat, axis = 1)
    nccn_daily_lon = nccn_daily_lon*1e-6
    return nccn_daily_lon, w_values


def get_height(data_path):
    height = get_control_ccn(data_path)[2]
    height_lat_mean = np.mean(height, axis = 1)
    height_mean = np.mean(height_lat_mean, axis = 1)
    height_axis = np.round(height_mean) * 1e-3
    return height_axis

def vertical_profile_CCN( data_path, date, w, ax):
    ccn_mean_cor = get_specific_time_w(data_path, date, w)[0]
    ccn_mean_check_input = get_specific_time_w(data_path_check_input, date, w)[0]
    height = get_height(data_path)
    w_values = get_specific_time_w(data_path, date, w)[1]
    ccn_mean_cor = ccn_mean_cor[0:30]
    ccn_mean_check_input = ccn_mean_check_input[0:30]
    height = height[0:30]
    line_width = 2
    fs = 10
    w = 'w = ' + str(w_values[w]) + ' m/s'
    ax.plot(ccn_mean_cor, height, linewidth = line_width, label = 'corect_input')
    ax.plot(ccn_mean_check_input, height, linewidth = line_width, label = 'not_corect_input')
    ax.set_xlabel("CCN ($\mathrm{cm^{-3}}$)", fontsize = fs)
    ax.set_ylabel("Height (km)", fontsize = fs)
    #ax.annotate(w, xy= (10, 10), fontsize = fs)
    ax.set_title(w, fontsize = fs, loc = 'left')
    ax.legend(fontsize = fs)
    return

fig, axes = plt.subplots(nrows = 5, ncols = 2, figsize = (15, 12))
i = 0
num_step = 0
first_timestep = num_step*8
second_timestep = first_timestep + 8
for ax in axes.flatten():
    # the number 0,1,2,.. indicates the index of omps_data which encounters for the day of simulation
    vertical_profile_CCN(data_path, num_step, i, ax)
    i = i + 1

plt.savefig('ccn_vertical_profile.png')
plt.show()