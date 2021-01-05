'''
North American Mesoscale (NAM) Model Data Plotter
Objective: Store and plot NAM data for model verification purposes
Inputs: date, latitude, longitude
Outputs: surface plot with colormap
'''

import pygrib # Used to read grib files
from datetime import datetime
from datetime import timedelta
import numpy as np
import matplotlib.pyplot as plt
import time
import glob
import sys
import cartopy.crs as ccrs
from matplotlib import rcParams

def nam_data_access(date):
    # Set up time range of an hour
    time_interval = timedelta(hours=2).total_seconds()
    # Get all grib files in 'nam' directory
    nam_file_list = glob.glob('nam/*.grb2')
    # Define smallest time difference between NAM file and date argument.
    # Smallest value corresponds to the NAM file to be used
    min_timedelta = time_interval
    file = ''
    # Iterate through list of NAM files and select the one with the shortest timedelta from the argument date, if smaller than the given time interval of +/- 2 hours
    for f in nam_file_list:
        file_time = datetime.strptime(''.join(f.split('.')[-2].split('_')[-3:-1]), '%Y%m%d%H%M')
        if abs((file_time - date).total_seconds()) < min_timedelta:
            file = f
            
    if not file:
        print('NAM data needed closer to requested time for accurate error analysis. Exiting...')
        sys.exit()
                
    return file

def nam_data_processing(nam_file, central_lat, central_lon, bound_sz):
    
    ### Read in grib file(s)
    gr = pygrib.open(nam_file)
    nam_raw_data = gr.select(name='Sensible heat net flux')[0] # Only select Q_H
    gr.close()
    
    nam_data_time = ''.join(nam_file.split('.')[-2].split('_')[-3:-1])
    nam_data_time = datetime.strptime(nam_data_time, '%Y%m%d%H%M').strftime('%Y-%m-%d %H:%M')
    
    # NAM Grid Projection (see https://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID218)
    lats = nam_raw_data.latlons()[0]
    lons = nam_raw_data.latlons()[1]
    nam_vals = nam_raw_data.values
        
    ### Filter data by coordinates of interest using masking
    lats_filter = np.where((lats > central_lat - bound_sz) & (lats < central_lat + bound_sz), 
                    lats, np.nan) 
    lons_filter = np.where((lons > central_lon - bound_sz) & (lons < central_lon + bound_sz), 
                    lons, np.nan) 
    nam_vals = np.where(~np.isnan(lats_filter) & ~np.isnan(lons_filter), nam_vals, np.nan)
    
    return lats, lons, nam_vals, nam_data_time

def nam_plot(lats, lons, nam_vals, nam_data_time, central_lat, central_lon, bound_sz):
    
    # Figure setup
    font_sz = 14
    rcParams['font.family'] = 'FreeSans'
    rcParams['font.size'] = font_sz
    
    fig, ax = plt.subplots(figsize=[8, 8])
    ax = plt.axes(projection=ccrs.LambertConformal(central_longitude=central_lon, central_latitude=central_lat))
    ax.set_extent([central_lon - bound_sz, 
                   central_lon + bound_sz, 
                   central_lat - bound_sz, 
                   central_lat + bound_sz])
    
    # Plot data formatting
    pc = plt.pcolormesh(lons, lats, nam_vals, cmap=plt.get_cmap('RdYlGn_r'), transform=ccrs.PlateCarree())
    colorbar = fig.colorbar(pc, pad=0.01)
    
    # Plot axis formatting
    gl = ax.gridlines(linestyle=":", draw_labels=True)
    gl.top_labels, gl.left_labels, gl.bottom_labels, gl.right_labels = [False, True, False, False]
    
    # Figure formatting
    title = 'Sensible heat net flux, ' + nam_data_time
    plt.title(title)
    plt.show()

def main(date, lat, lon):
    start_time = time.time() # Start program timer
    
    ### Define location and time of interest
    bound_sz = 5 # Window size outward from central point, in degrees
    date = datetime.strptime(date, '%Y%m%d%H%M')
    
    nam_file = nam_data_access(date)
    lats, lons, nam_vals, nam_data_time = nam_data_processing(nam_file, lat, lon, bound_sz)
    nam_plot(lats, lons, nam_vals, nam_data_time, lat, lon, bound_sz)
    
    end_time = time.time() - start_time
    print('Program runtime: %.3f s' % end_time)
    
# Call function with date (YYYYMMDDHHMM format), latitude, longitude (degrees E)
main('201907281700', 40.631762, -73.953678)