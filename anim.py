from google.cloud import storage
import matplotlib.pylab as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.io.img_tiles as tiles
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import goes_air_temperature.goes_16_air_temperature as g16 # Revised script
import math
import numpy as np
from matplotlib import rcParams
import cartopy.io.img_tiles as cimgt
from matplotlib.animation import FuncAnimation
     
import animatplot as amp

### General plot settings    
font_sz = 14
rcParams['font.family'] = 'FreeSans'
rcParams['font.size'] = font_sz
sat_img = cimgt.GoogleTiles(style='satellite')
plt.rc('text', usetex=False)

######################################################################################################
### Basemapper Plot (to be replaced by plotter.py)
# Objective: Grab LST data from GOES-16 ABI L2+ on Google Cloud

def basemapper(crd,res):
    
    bbox = g16.site_info(crd,res)

    layout_map = tiles.QuadtreeTiles()
    fig = plt.figure(figsize=(12,9))
    # Create a GeoAxes in the tile's projection.
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    # Limit the extent of the map to a small longitude/latitude range.
    ax.set_extent([bbox[0],bbox[2],bbox[1],bbox[3]], crs=ccrs.PlateCarree())
    
    # # Add the satellite background data at zoom level 12.
    ax.add_image(sat_img, 12)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0, color='gray', alpha=1, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    # gl.xformatter = LONGITUDE_FORMATTER
    # gl.yformatter = LATITUDE_FORMATTER

    return fig, ax
    
## ################## END Basemapper #########################
# def animated_plot(data, crd, goes_grid_info, prods, data_type)
    
# '''
res=0.25
data = t_arr
crd = nyc_crd
[lats, lons, nc_indx_corners, nc_indx_bounds, lon_plot, lat_plot] = goes_grid_info
[prod_vals, prod_set, prod_names] = prods[0]
data_type = "t_air"

z = data.shape[2]

lon_adjust = 0.5*np.diff(lons[nc_indx_bounds[0]:nc_indx_bounds[1],
                          nc_indx_bounds[2]-1:nc_indx_bounds[3]],axis=1) # this shifts lons just for pcolormesh (lower-left corner plot shift)
lat_adjust = 0*np.diff(lats[nc_indx_bounds[0]:nc_indx_bounds[1]+1,
                          nc_indx_bounds[2]:nc_indx_bounds[3]],axis=0) # this shifts lats just for pcolormesh (lower-left corner plot shift)

fig, ax = basemapper(crd,res) # plot the air temp basemap

# Set colorbar limits, resolution, and number of increments
if data_type == "q_h":
    norm_t = [-50, 150] # Limits chosen based on min and max temperatures recorded in CONUS
    norm_res = 10
elif data_type == "t_air":
    norm_t = [290, 310]
    norm_res = 2
else:
    norm_t = [0, 10]
    norm_res = 1
norm_incs = math.ceil((norm_t[1] - norm_t[0])/norm_res) # Get ceiling to force integer
cmap = plt.get_cmap('RdYlGn_r')
norm = colors.BoundaryNorm(boundaries=np.linspace(norm_t[0], norm_t[1], norm_incs), ncolors=cmap.N)

cax = [ax.pcolormesh(lon_plot-lon_adjust, 
                    lat_plot-lat_adjust,
                    data[:,:,0],
                    cmap=cmap,
                    zorder=999,
                    alpha=0.5,
                    transform=ccrs.PlateCarree(), 
                    norm=norm, 
                    edgecolor='k', 
                    linewidth=0)] # plotting the air temp data

cb = fig.colorbar(cax[0], pad=0.01)

if data_type == 'q_h':
    plt.suptitle('Sensible Heat Flux, $Q_H$ \n', fontsize=font_sz*1.5, y=0.95)
    plt.title('{1} {2} - {3} UTC'.format(prod_set.variables[prod_names[0]].long_name,
                          prod_set.time_coverage_start.split('T')[0],
                          prod_set.time_coverage_start.split('T')[1].replace('Z',''),
                          prod_set.time_coverage_end.split('T')[1].replace('Z','')),
                  fontsize=font_sz*1.2, pad=10)
    cb.set_label('$Q_H$ [$W/m^2$]')
elif data_type == 't_air':
        plt.suptitle('Air Temperature, $T_{air}$ \n', fontsize=font_sz*1.5, y=0.95)
        plt.title('{1} {2} - {3} UTC'.format(prod_set.variables[prod_names[0]].long_name,
                              prod_set.time_coverage_start.split('T')[0],
                              prod_set.time_coverage_start.split('T')[1].replace('Z',''),
                              prod_set.time_coverage_end.split('T')[1].replace('Z','')),
                      fontsize=font_sz*1.2, pad=10)
        cb.set_label('$T_{air}$ [$K$]')
else:
    plt.title('{1} {2} - {3} UTC'.format(prod_set.variables[prod_names[0]].long_name,
                          prod_set.time_coverage_start.split('T')[0],
                          prod_set.time_coverage_start.split('T')[1].replace('Z',''),
                          prod_set.time_coverage_end.split('T')[1].replace('Z','')),
                  fontsize=font_sz*1.2, pad=10)
def animate(t, data, cax):
    cax[0].remove()
    cax[0] = ax.pcolormesh(lon_plot-lon_adjust, 
                        lat_plot-lat_adjust,
                        data[:,:,t],
                        cmap=cmap,
                        zorder=999,
                        alpha=0.5,
                        transform=ccrs.PlateCarree(), 
                        norm=norm, 
                        edgecolor='k', 
                        linewidth=0) 
    [prod_vals, prod_set, prod_names] = prods[t]
        
    plt.title('{1} {2} - {3} UTC'.format(prod_set.variables[prod_names[0]].long_name,
                          prod_set.time_coverage_start.split('T')[0],
                          prod_set.time_coverage_start.split('T')[1].replace('Z',''),
                          prod_set.time_coverage_end.split('T')[1].replace('Z','')),
                  fontsize=font_sz*1.2, pad=10)

ani = FuncAnimation(fig, animate, frames=z, interval=750, repeat_delay=0, blit=False, fargs=(data, cax))
plt.show()
# '''

