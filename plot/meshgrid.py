'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        Meshgrid Plotter
Package file path:  ~/plot/meshgrid.py
Objective:          Plot gridded data onto corresponding 2D map with animation option.
Author:             Gabriel Rios
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

import cartopy.crs as ccrs, cartopy.io.img_tiles as cimgt, datetime, numpy as np, matplotlib, matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

##############################################################################################
# END IMPORTS
##############################################################################################

##############################################################################################
# Method name:      base
# Method objective: Provide general formatting and projection settings to the 2D plots.
# Input(s):         domain [list], domain_center [list], zoom [float]
# Outputs(s):       matplotlib objects [fig, ax, projection]
##############################################################################################

def base(domain, domain_center, zoom, font_size):
    
    # Define projections for the base map
    proj_pc = ccrs.PlateCarree()
    proj_ortho = ccrs.Orthographic(central_latitude=domain_center[0], central_longitude=domain_center[1])
    fig, ax = plt.subplots(subplot_kw={'projection': proj_ortho}, dpi=144)
    # Limit the extent of the map to a small longitude/latitude range
    ax.set_extent([domain[2] - zoom, domain[3] + zoom, domain[0] - zoom, domain[1] + zoom])
    ax.coastlines()
    
    # Add the satellite background data at zoom level 12.
    # sat_img = cimgt.GoogleTiles(style='satellite')
    # ax.add_image(sat_img, 12)

    # Draw and format gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.1, linestyle='--')
    gl.xlabel_style, gl.ylabel_style = {'fontsize': font_size['plot_ticks']}, {'fontsize': font_size['plot_ticks']}
    gl.top_labels = False
    gl.right_labels = False
    gl.xpadding = 10
    gl.ypadding = 10

    return fig, ax, proj_pc

##############################################################################################
# Method name:      gridded_data
# Method objective: Define data to be plotted along with corresponding formatting.
# Input(s):         domain [list], domain_center [list], date [datetime or list], 
#                   lats [2D NumPy array], lons [2D NumPy array], data [2D or 3D array], 
#                   var [str], animated [bool]
# Outputs(s):       plot or .gif file
##############################################################################################

def gridded_data(domain, domain_center, date, lats, lons, data, data_limits, var, zoom=0.1, animated=False, savefig=False):
    
    matplotlib.rcParams['font.family'] = 'FreeSans'
    font_size = {'title': 13,
                 'subtitle': 10,
                 'plot_ticks': 10,
                 'axes_label': 10}
    
    # Define padding around data for map
    fig, ax, proj_pc = base(domain, domain_center, zoom, font_size)
    
    # Colormap definition
    # Logic: Color negative fluxes blue (to represent cooling), positive fluxes red (heating), net zero white
    if data_limits[0] < 0 and data_limits[1] > 0:
        norm = matplotlib.colors.TwoSlopeNorm(vmin=data_limits[0], vcenter=0, vmax=data_limits[1])
        cmap = plt.cm.get_cmap('RdBu_r')
    elif data_limits[0] < 0 and data_limits[1] < 0:
        norm = matplotlib.colors.Normalize(vmin=data_limits[0], vmax=0)
        cmap = plt.cm.get_cmap('Blues_r')
    else:
        norm = matplotlib.colors.Normalize(vmin=0, vmax=data_limits[1])
        cmap = plt.cm.get_cmap('Reds_r')
    
    # Plot the data
    # Colormap placed in list to enable animation function
    if animated:    
        # Note that pcolormesh data dimensions must 1 smaller along each axis than the X and Y arrays.
        im = [ax.pcolormesh(lons, lats, data[:-1, :-1, 0], cmap=cmap, norm=norm, transform=proj_pc)]
        colorbar = fig.colorbar(im[0], pad=0.01)
        subtitle_str = '{0} UTC'.format(date[0].strftime('%Y-%m-%d %H:%M:%S'))
    else: 
        # Note that pcolormesh data dimensions must 1 smaller along each axis than the X and Y arrays.
        im = ax.pcolormesh(lons, lats, data[:-1, :-1], cmap=cmap, norm=norm, transform=proj_pc)
        colorbar = fig.colorbar(im, pad=0.01)
        colorbar.set_label('${0}$ [${1}$]'.format(var['name'], var['units']))
        subtitle_str = '{0} UTC'.format(date)
    
    colorbar.set_label('${0}$ [${1}$]'.format(var['name'], var['units']), fontsize=font_size['axes_label'], rotation=270, labelpad=20)
    colorbar.ax.tick_params(labelsize=font_size['axes_label']) 
    
    # Figure title formatting
    title_str = '{0}'.format(var['full_name'])
    # Note: value used for x-position centers the figure title relative to the subplot
    fig_title = fig.suptitle(title_str, 
                             fontsize=font_size['title'], 
                             x=np.nanmean([ax.get_position().extents[0], ax.get_position().extents[2]]),
                             y=0.99)
    ax_title = ax.set_title(subtitle_str, fontsize=font_size['subtitle'])
    
    # Call animation function and save figure if animation is enabled.
    if animated:
        z = data.shape[2]
        
        ############################################################################
        # Method name:      animate
        # Method objective: Enable plot animation.
        # Input(s):         t [int], data [2D ndarray], im [matplotlib object]
        # Outputs(s):       none
        ############################################################################
        
        def animate(t, data, im):
            im[0].remove()
            # Note that pcolormesh data dimensions must 1 smaller along each axis than the X and Y arrays.
            im[0] = ax.pcolormesh(lons, lats, data[:-1, :-1, t], cmap=cmap, alpha=0.5, norm=norm, transform=ccrs.PlateCarree(), zorder=999, linewidth=0) 
            subtitle_str = '{0} UTC'.format(date[t].strftime('%Y-%m-%d %H:%M:%S'))
            ax_title = ax.set_title(subtitle_str, fontsize=font_size['subtitle'])
        
        # Call animation function given the appropriate arguments.
        anim = FuncAnimation(fig, animate, frames=z, interval=500, repeat_delay=2000, blit=False, fargs=(data, im))
        plt.show()
        # Save the figure
        if savefig:
            # File name format: 
            # (NW latitude)_(NW longitude)_(SE latitude)_(SE longitude)_(start date)_(end_date).gif
            file_str = '{0}_{1}_{2}_{3}_s{4}_e{5}.gif'.format(domain[1], domain[2], domain[0], domain[3], date[0], date[-1])
            anim.save(file_str, writer='pillow', fps=2, dpi=300, bitrate=-1, codec="libx264")
        
    return fig, ax

##############################################################################################
# Method name:      main
# Method objective: Control sequencing and relaying of plot data.
# Input(s):         domain [list], domain_center [list], date_range [datetime or list], 
#                   lats [2D NumPy array], lons [2D NumPy array], var [str], metadata [dict], 
#                   zoom [float], animated [bool]
# Outputs(s):       plot or .gif file
##############################################################################################

def main(domain, domain_center, date_range, lats, lons, var, metadata, zoom=0.0, animated=False, savefig=False):
    # Grab variable of interest from global variables collection.
    grid_data = var
    # Define absolute data limits over timeseries for coloring purposes
    data_limits = [np.nanmin(grid_data), np.nanmax(grid_data)]
    # Define depth of the gridded data matrix
    z = grid_data.shape[2]
    # If animation chosen, plot using this function. Else, use iterative prodecure below.
    if animated:
        datamap = gridded_data(domain, domain_center, date_range, lats, lons, grid_data, data_limits, metadata, zoom=zoom, animated=animated, savefig=savefig)
    else:
        for i in range(0, z):
            datamap = gridded_data(domain, domain_center, date_range[i], lats, lons, grid_data[:, :, i], data_limits, metadata, zoom=0.00)