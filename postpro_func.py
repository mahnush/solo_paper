import numpy as np


def read_nc(file_name, var_name):
    from netCDF4 import Dataset
    file = file_name
    print(file)
    nc = Dataset(file,'r')
    var = nc.variables[var_name][:]
    return var

def read_nc_sat(file_name, var_name):
    from netCDF4 import Dataset
    file = file_name
    nc = Dataset(file,'r')
    var = nc.variables[var_name][:,:]
    var = np.transpose(var)
    return var
def visulize_sat(ax, var,latgrid, longrid, bounds, cbar_label,
             titel_figure, map_limit, cax):
    fs_titel = 10
    fs_label = 10
    fs_coardinate = 6
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
    cmap = [(0.0,0.0,0.0)] + [(cm.jet(i)) for i in range(1,256)]
    cmap = mpl.colors.ListedColormap(cmap)
    bounds = bounds
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cmap_2 = LinearSegmentedColormap.from_list("", ["white","lightskyblue","steelblue","green","yellowgreen","yellow","gold","red","firebrick","darkred"])
    #cmap_2 = LinearSegmentedColormap.from_list("",["lightskyblue", "steelblue", "green", "yellowgreen", "yellow",
     #                                           "gold", "red", "firebrick", "darkred"])
    m = Basemap(ax=ax, projection='cyl', llcrnrlat= map_limit[0], urcrnrlat=map_limit[1],\
           llcrnrlon=map_limit[2], urcrnrlon=map_limit[3], resolution='c')
    m.drawmeridians(np.arange(-180., 180., 10.), linewidth=1.2, labels=[0, 0, 0, 1], color='grey', zorder=3, latmax=90,
               fontsize = fs_coardinate )
    m.drawparallels(np.arange(-85., 85., 10), linewidth=1.2, labels=[1, 0, 0, 0], color='grey', zorder=2, latmax=90,
                 fontsize = fs_coardinate)
    m.drawcoastlines()
    m.drawcountries()
    x, y = m(longrid, latgrid)
    ax.set_title(titel_figure, fontsize=fs_titel)
    l1=m.pcolormesh(x, y, var, cmap=cmap_2, norm=norm)
    cbar = plt.colorbar(l1, ax=ax, cax = cax, orientation = 'horizontal')
    #cbar = plt.colorbar(l1, shrink=0.50, ax=ax)
    cbar_bounds = bounds
    cbar_ticks =  bounds
    cbar.set_label(cbar_label, fontsize= fs_label)
    cbar.ax.tick_params(labelsize= 8)

def visulize_model(ax, var, bounds, cbar_label, titel_figure, map_limit):
    fs_titel = 20
    fs_label = 20
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
    cmap = [(0.0,0.0,0.0)] + [(cm.jet(i)) for i in range(1,256)]
    cmap = mpl.colors.ListedColormap(cmap)
    bounds = bounds
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    #cmap_2 = LinearSegmentedColormap.from_list("", ["white","lightskyblue","steelblue","green","yellowgreen","yellow","gold","red","firebrick","darkred"])
    cmap_2 = LinearSegmentedColormap.from_list("",["white", "lightskyblue", "steelblue", "green", "yellowgreen", "yellow",
                                                "gold", "red", "firebrick", "darkred"])
    m = Basemap(ax=ax, projection='cyl', llcrnrlat= map_limit[0], urcrnrlat=map_limit[1],\
           llcrnrlon=map_limit[2], urcrnrlon=map_limit[3], resolution='c')
    m.drawmeridians(np.arange(-180., 180., 10.), linewidth=1.2, labels=[0, 0, 0, 1], color='grey', zorder=3, latmax=90)
    m.drawparallels(np.arange(0., 85., 10.), linewidth=1.2, labels=[1, 0, 0, 0], color='grey', zorder=2, latmax=90)

    m.drawcoastlines()
    ax.set_title(titel_figure, fontsize=fs_titel)
    #l1 =m.imshow(var, cmap = cmap_2, norm = norm)
    lat = np.linspace(map_limit[0], map_limit[1], 1501)
    lon = np.linspace(map_limit[2], map_limit[3], 2501)
    xx, yy = np.meshgrid(lon, lat)
    print(lat)
    print(lon)
    varMin, varMax, varInt = 0, 202, 2
    levels = np.arange(varMin, varMax + varInt, varInt)
    # l1 =m.imshow(var,vmin=varMin, vmax=varMax)
    l1 = m.contourf(xx, yy, var, vmin = varMin,
                    vmax = varMax, levels = levels, cmap = cmap_2)
    cbar = plt.colorbar(l1, shrink=0.50, ax=ax)
    cbar_bounds = bounds
    cbar_ticks =  bounds
    cbar.set_label(cbar_label, fontsize= fs_label)
    cbar.ax.tick_params(labelsize='xx-large')

def visulize_model_diff(ax, var, bounds, cbar_label, titel_figure, map_limit):
    fs_titel = 30
    fs_label = 30
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
    cmap = [(0.0,0.0,0.0)] + [(cm.jet(i)) for i in range(1,256)]
    cmap = mpl.colors.ListedColormap(cmap)
    bounds = bounds
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    m = Basemap(ax=ax, projection='cyl', llcrnrlat= map_limit[0], urcrnrlat=map_limit[1],\
           llcrnrlon=map_limit[2], urcrnrlon=map_limit[3], resolution='c')
    m.drawmeridians(np.arange(-180., 180.,5.), linewidth=1.2, labels=[0, 0, 0, 1], color='grey', zorder=3, latmax=90,
                 fontsize = 20)
    m.drawparallels(np.arange(0., 85., 5.), linewidth=1.2, labels=[1, 0, 0, 0], color='grey', zorder=2, latmax=90,
                 fontsize = 20 )
    m.drawcoastlines()
    ax.set_title(titel_figure, fontsize=fs_titel)
    l1 =m.imshow(var, cmap = 'seismic', norm = norm)
    cbar = plt.colorbar(l1, ax=ax)
    cbar_bounds = bounds
    cbar_ticks =  bounds
    cbar.set_label(cbar_label, fontsize= fs_label)
    cbar.ax.tick_params(labelsize='xx-large')


def mask_land(var, lat, lon):
    from global_land_mask import globe
    import numpy as np
    import  numpy.ma as ma
    lon_mesh, lat_mesh = np.meshgrid(lon, lat)
    globe_land_mask = globe.is_land(lat_mesh, lon_mesh)
    var_ocean = ma.masked_where(globe_land_mask == True, var)
    return var_ocean

