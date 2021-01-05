'''
High-Resolution Rapid Refresgh (HRRR) Model Data Plotter
Objective: Store and plot HRRR data for model verification purposes
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
import os
import sys
import re
import cartopy.crs as ccrs
from matplotlib import rcParams
import urllib.request

def hrrr_data_download(date):
    # Reformat the URL and define the local target file strings
    date_path = date.strftime('%Y%m%d')
    url = 'https://pando-rgw01.chpc.utah.edu/hrrr/sfc/' + date.strftime('%Y%m%d') + '/hrrr.t' + date.strftime('%H') + 'z.wrfsfcf00.grib2'
    target = 'hrrr/' + date.strftime('%Y%m%d') + '_hrrr.t' + date.strftime('%H') + 'z.wrfsfcf00.grib2'
    
    # Check if file exists. If so, download. Else, exit.
    if not os.path.isfile(target):
        try:
            print('Downloading file...')
            urllib.request.urlretrieve(url, target)
            print('Downloaded file!')
        except:
            print('Could not download file, exiting...')
            sys.exit()
    else:
        print('File already exists, moving onto data access and processing...')    
            
def hrrr_data_access(date):    
    
    # Set up time range of an hour
    time_interval = timedelta(hours=2).total_seconds()
    # Get all grib files in 'hrrr' directory
    hrrr_file_list = glob.glob('hrrr/*.grib2')
    # Define smallest time difference between hrrr file and date argument.
    # Smallest value corresponds to the hrrr file to be used
    min_timedelta = time_interval
    file = ''
    
    # Iterate through list of hrrr files and select the one with the shortest timedelta from the argument date, if smaller than the given time interval of +/- 2 hours
    for f in hrrr_file_list:
        time_str = f.split('/')[-1][0:8] + f.split('/')[-1].split('.')[1][1:3]
        file_time = datetime.strptime((time_str + '00'), '%Y%m%d%H%M')
        if abs((file_time - date).total_seconds()) < min_timedelta:
            file = f
            
    if not file:
        print('HRRR data needed closer to requested time for accurate error analysis. Exiting...')
        sys.exit()
                
    return file

def hrrr_data_processing(hrrr_file, central_lat, central_lon, bound_sz):
    
    # List of variables: https://www.nco.ncep.noaa.gov/pmb/docs/on388/table2.html
    print(hrrr_file)
    
    ### Read in grib file(s)
    gr = pygrib.open(hrrr_file)
    hrrr_raw_data = gr.select(name='Sensible heat net flux')[0] # Only select Q_H
    gr.close()
    
    hrrr_data_time = hrrr_file.split('/')[-1][0:8] + hrrr_file.split('/')[-1].split('.')[1][1:3] + '00'
    hrrr_data_time = datetime.strptime(hrrr_data_time, '%Y%m%d%H%M').strftime('%Y-%m-%d %H:%M')
    
    # hrrr Grid Projection (see https://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID218)
    lats = hrrr_raw_data.latlons()[0]
    lons = hrrr_raw_data.latlons()[1]
    hrrr_vals = hrrr_raw_data.values
        
    ### Filter data by coordinates of interest using masking
    lats_filter = np.where((lats > central_lat - bound_sz) & (lats < central_lat + bound_sz), 
                    lats, np.nan) 
    lons_filter = np.where((lons > central_lon - bound_sz) & (lons < central_lon + bound_sz), 
                    lons, np.nan) 
    hrrr_vals = np.where(~np.isnan(lats_filter) & ~np.isnan(lons_filter), hrrr_vals, np.nan)
    
    return lats, lons, hrrr_vals, hrrr_data_time


def hrrr_plot(lats, lons, hrrr_vals, hrrr_data_time, central_lat, central_lon, bound_sz):
    
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
    pc = plt.pcolormesh(lons, lats, hrrr_vals, cmap=plt.get_cmap('RdYlGn_r'), transform=ccrs.PlateCarree())
    colorbar = fig.colorbar(pc, pad=0.01)
    
    # Plot axis formatting
    gl = ax.gridlines(linestyle=":", draw_labels=True)
    gl.top_labels, gl.left_labels, gl.bottom_labels, gl.right_labels = [False, True, False, False]
    
    # Figure formatting
    title = 'Sensible heat net flux, ' + hrrr_data_time
    plt.title(title)
    plt.show()

def main(date, lat, lon):
    start_time = time.time() # Start program timer
    
    ### Define location and time of interest
    bound_sz = 2 # Window size outward from central point, in degrees
    date = datetime.strptime(date, '%Y%m%d%H%M')
    
    hrrr_data_download(date)
    hrrr_file = hrrr_data_access(date)
    lats, lons, hrrr_vals, hrrr_data_time = hrrr_data_processing(hrrr_file, lat, lon, bound_sz)
    hrrr_plot(lats, lons, hrrr_vals, hrrr_data_time, lat, lon, bound_sz)
    
    end_time = time.time() - start_time
    print('Program runtime: %.3f s' % end_time)
    
# Call function with date (YYYYMMDDHHMM format), latitude, longitude (degrees E)
main('201907281800', 40.631762, -73.953678)