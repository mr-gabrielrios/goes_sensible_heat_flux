'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        HRRR Data Processing
Package file path:  ~/vdtn/wrf/hrrr/reader.py
Objective:          Collect HRRR data and use for model input or validation.
Author:             Gabriel Rios
Note:               Credit to the University of Utah (Brian Blaylock) and the Climate Corporation 
                    for hosting an HRRR archive.
                    Link: http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/cgi-bin/hrrr_download.cgi
                    DOI: https://doi.org/10.7278/S5JQ0Z5B
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

import pygrib # Used to read grib files
from datetime import datetime
from datetime import timedelta
import numpy as np
import matplotlib.pyplot as plt
import time
import glob
import os
import sys
import cartopy.crs as ccrs
from matplotlib import rcParams
import urllib.request
import re
import requests

##############################################################################################
# END IMPORTS
##############################################################################################

##############################################################################################
# Method name:      download_HRRR_subset
# Method objective: Download HRRR data for specified parameter. Assume that forecast hour is 0.
# Input(s):         date [datetime], searchString [str], savedir [str], dryrun [bool]
# Outputs(s):       none
# Reference:        https://doi.org/10.1016/j.cageo.2017.08.005
##############################################################################################

def download_HRRR_subset(date, searchString, savedir='data/', dryrun=False):
    
    url = 'https://pando-rgw01.chpc.utah.edu/hrrr/sfc/' + date.strftime('%Y%m%d') + '/hrrr.t' + date.strftime('%H') + 'z.wrfsfcf00.grib2'
    target = savedir + 'subset_' + searchString.split(':')[1] + '_' + date.strftime('%Y%m%d') + '_hrrr.t' + date.strftime('%H') + 'z.wrfsfcf00.grib2'
    target = os.path.join(os.path.dirname(__file__), target)
    
    if os.path.isfile(target):
        print('File already exists, moving onto data access and processing...')    
    else:
        # Make SAVEDIR if path doesn't exist
        if not os.path.exists(savedir):
            os.makedirs(savedir)
            print(f'Created directory: {savedir}')
        
        # Make a request for the .idx file for the above URL
        idx = url + '.idx'
        r = requests.get(idx)
    
        # Check that the file exists. If there isn't an index, you will get a 404 error.
        if not r.ok: 
            print('Failure tp access file. Status Code:', r.status_code, r.reason)
            print(f'It does not look like the index file exists: {idx}')
    
        # Read the text lines of the request
        lines = r.text.split('\n')
        
        # Search expression
        expr = re.compile(searchString)
    
        # Store the byte ranges in a dictionary
        #     {byte-range-as-string: line}
        byte_ranges = {}
        for n, line in enumerate(lines, start=1):    
            # Use the compiled regular expression to search the line
            if expr.search(line):   
                # Get the beginning byte in the line we found
                parts = line.split(':')
                rangestart = int(parts[1])
                # Get the beginning byte in the next line...
                if n+1 < len(lines):
                    parts = lines[n].split(':')
                    rangeend = int(parts[1])
                else:
                    rangeend = ''
    
                # Store the byte-range string in our dictionary, 
                # and keep the line information too so we can refer back to it.
                byte_ranges[f'{rangestart}-{rangeend}'] = line
                
        for i, (byteRange, line) in enumerate(byte_ranges.items()):
            if i == 0:
                # If we are working on the first item, overwrite the existing file.
                curl = f'curl -s --range {byteRange} {url} > {target}'
            else:
                # If we are working on not the first item, append the existing file.
                curl = f'curl -s --range {byteRange} {url} >> {target}'
                
            num, byte, date, var, level, forecast, _ = line.split(':')
            
            if dryrun:
                print(f'\t Dry Run: Found GRIB line [{num:>3}]: variable={var}, level={level}, forecast={forecast}')
            else:
                print(f'\t Downloading GRIB line [{num:>3}]: variable={var}, level={level}, forecast={forecast}')    
                os.system(curl)
        
        if dryrun:
            print(f'Dry Run: Success! Searched for [{searchString}] and found [{len(byte_ranges)}] GRIB fields. Would save as {target}')
        else:
            print(f'\t Success! Searched for [{searchString}] and got [{len(byte_ranges)}] GRIB fields and saved as {target}')
        
            return target
        
        # Reformat the URL and define the local target file strings
        url = 'https://pando-rgw01.chpc.utah.edu/hrrr/sfc/' + date.strftime('%Y%m%d') + '/hrrr.t' + date.strftime('%H') + 'z.wrfsfcf00.grib2'
        target = savedir + 'subset_' + searchString.split(':')[1] + '_' + date.strftime('%Y%m%d') + '_hrrr.t' + date.strftime('%H') + 'z.wrfsfcf00.grib2'
        target = os.path.join(os.path.dirname(__file__), target)
        
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
     
##############################################################################################
# Method name:      hrrr_data_access
# Method objective: Specify file to be used for a given time.
# Input(s):         date [datetime], savedir [str]
# Outputs(s):       file [str]
##############################################################################################
     
def hrrr_data_access(date, savedir='data/'):    
    
    # Set up time range of an hour
    time_interval = timedelta(hours=2).total_seconds()
    # Get all grib files in 'hrrr' directory
    hrrr_file_list = glob.glob(os.path.join(os.path.dirname(__file__), savedir) + '/*.grib2')
    # Define smallest time difference between hrrr file and date argument.
    # Smallest value corresponds to the hrrr file to be used
    min_timedelta = time_interval
    filename = ''
    
    # Iterate through list of hrrr files and
    # select the one with the shortest timedelta from the argument date, 
    # if smaller than the given time interval of +/- 2 hours
    for f in hrrr_file_list:
        time_str = f.split('_')[-2] + f.split('.')[-3][1:3]
        file_time = datetime.strptime(time_str, '%Y%m%d%H')
        if abs((file_time - date).total_seconds()) < min_timedelta:
            filename = f
            min_timedelta = abs((file_time - date).total_seconds())
            
    if not filename:
        print('HRRR data needed closer to requested time for accurate error analysis. Exiting...')
        sys.exit()
                
    return filename

##############################################################################################
# Method name:      hrrr_data_processing
# Method objective: Filter HRRR data for a given file and time.
# Input(s):         file [str], domain [list] 
# Outputs(s):       lats [2D ndarray], lons [2D ndarray], hrrr_vals [2D ndarray], 
#                   hrrr_data_time [str], hrrr_dataset_name [str]
##############################################################################################
     
def hrrr_data_processing(file, domain):
    
    # List of variables:
    # http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/HRRR_archive/hrrr_sfc_table_f00-f01.html
    
    # Read in grib file(s)
    gr = pygrib.open(file)
    # Grab name of dataset desired. Assume only one parameter's data is downloaded.
    for g in gr:
        hrrr_dataset_name = str(g).split(':')[1]
    # Only select the data for the chosen parameter
    hrrr_raw_data = gr.select(name=hrrr_dataset_name)[0] 
    gr.close()
    
    hrrr_data_time = file.split('_')[-2] + file.split('.')[-3][1:3]
    hrrr_data_time = datetime.strptime(hrrr_data_time, '%Y%m%d%H').strftime('%Y-%m-%d %H:%M')    

    hrrr_vals, lats, lons = hrrr_raw_data.data(lat1=domain[0], lat2=domain[1], lon1=domain[2], lon2=domain[3])
    
    return lats, lons, hrrr_vals, hrrr_data_time, hrrr_dataset_name

##############################################################################################
# Method name:      hrrr_plot
# Method objective: Filter HRRR data for a given file and time.
# Input(s):         file [str], domain [list] 
# Outputs(s):       lats [2D ndarray], lons [2D ndarray], hrrr_vals [2D ndarray], 
#                   hrrr_data_time [str], hrrr_dataset_name [str]
##############################################################################################
    
def hrrr_plot(lats, lons, hrrr_vals, domain, hrrr_data_time, hrrr_dataset_name):
    
    # Figure setup
    font_sz = 14
    rcParams['font.family'] = 'Segoe UI'
    rcParams['font.size'] = font_sz
    
    proj = ccrs.LambertConformal(central_longitude=np.nanmean([domain[2], domain[3]]), central_latitude=np.nanmean([domain[0], domain[1]]))
    fig, ax = plt.subplots(figsize=[8, 8], dpi=144)
    ax = plt.axes(projection=proj)
    ax.set_extent([domain[2], domain[3], domain[0], domain[1]])
    ax.coastlines()
    
    # Plot data formatting
    pc = plt.pcolormesh(lons, lats, hrrr_vals, cmap=plt.get_cmap('RdYlGn_r'), transform=ccrs.PlateCarree())
    colorbar = fig.colorbar(pc, pad=0.01)
    
    # Plot axis formatting
    gl = ax.gridlines(linestyle=":", draw_labels=True)
    gl.top_labels, gl.left_labels, gl.bottom_labels, gl.right_labels = [False, True, False, False]
    
    # Figure formatting
    title = hrrr_dataset_name + '\n' + hrrr_data_time
    plt.title(title)
    plt.show()

if __name__ == '__main__':
    start_time = time.time() # Start program timer
    
    ### Define location and time of interest
    bound_sz = 0.25 # Window size outward from central point, in degrees
    dates = [datetime(year=2019, month=7, day=28, hour=18)]
    
    for date in dates:
        print(date)
        download_HRRR_subset(date, ':WIND:')
    
    for date in dates:
        print(date)
        
        hrrr_file = hrrr_data_access(date) 
        domain = [40.60, 40.80, -74.05, -73.85]  
        lats, lons, hrrr_vals, hrrr_data_time, hrrr_dataset_name = hrrr_data_processing(hrrr_file, domain)        
        t = time.time()
        # hrrr_plot(lats, lons, hrrr_vals, domain, hrrr_data_time, hrrr_dataset_name)
    
    end_time = time.time() - start_time
    print('Program runtime: %.3f s' % end_time)