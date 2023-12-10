import postpro_func
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
old_con_run = '/work/bb1036/b380900/work_mistral//output_solo/10Apr_newrun/NWP_LAM_DOM01_20210410T100000Z_0010.nc'
lbc_con_run = '/work/bb1093/b380900/experiments/icon_lam_lim_solo_control_boundary/NWP_LAM_DOM01_20210410T100000Z_0010.nc'
ini_con_run = '/work/bb1093/b380900/experiments/icon_lam_lim_solo_control/NWP_LAM_DOM01_20210410T100000Z_0010.nc'

def weight(var):
    weight = np.zeros_like(var) + 1. / (var.size)
    return weight

def plot_hist(data_path , var_name, label):
    var = postpro_func.read_nc(data_path, var_name)
    var = var.flatten()
    var_nozero = np.ma.masked_where(var == 0.0, var)
    var_nozero = var_nozero.flatten()
    var_nozero = var_nozero * 1000.0
    numbin = np.arange(0, 500, 10)
    ax.hist(var_nozero, numbin, weights = weight(var_nozero), label = label, histtype = 'step')
    ax.legend()
    return
#fig, ax = plt.subplots(1,1, figsize = (15,12))
fig, (ax) = plt.subplots(1, 1, figsize = (15,12))

var_name = 'tqc'

plot_hist(old_con_run, var_name, 'old_run')
plot_hist(lbc_con_run, var_name, 'lbc_new_run')
plot_hist(ini_con_run, var_name, 'ini_new_run')
plt.show()
plt.savefig('test.png')