##############################################################################################
# MIT License
# 
# Copyright (c) [2020] [Gabriel Rios]
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
##############################################################################################

'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        Main
Package file path:  ~
Objective:          Estimate sensible heat flux (Q_H) (known as 'surface upward heat flux in air' per the Climate          
                    and Forecasting Conventions, v77)
Author:             Gabriel Rios
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

# External imports
from google.cloud import storage
import datetime, numpy as np, os, pandas, time
# Internal imports
from grid import canopy_height, elevation, goes_grid_setup, land_cover 
# from grid import land_cover
from misc import asos_data_reader, time_adjust
from phys import qh, t_lst, t_air_2m
from plot import meshgrid, timeseries, varnames
from user import inputs
from vdtn.obs import functions
from vdtn.obs.mesonet import reader as mesonet
from vdtn.wrf.hrrr import reader as hrrr

##############################################################################################
# END IMPORTS
##############################################################################################
     
##############################################################################################
# Method name:      gcp_bucket
# Method objective: Set up Google Cloud credentials and access publicly-hosted bucket.
# Input(s):         none
# Outputs(s):       bucket [N/A]
##############################################################################################

def gcp_bucket():
    # Reference local Google Credentials
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "GOOGLE_AUTH_CREDS.json"
    
    # Start the storage client
    client = storage.Client() 
    
    # Call the GOES-16 storage bucket
    bucket = client.get_bucket('gcp-public-data-goes-16') 
    
    return bucket

##############################################################################################
# Method name:      input_setup
# Method objective: Grab user-provided data and package it for model use.
# Input(s):         domain [bool or list], date_range [bool or list], plot_var [bool or str],
#                   timeseries [bool], gridmap [bool]
# Outputs(s):       domain [list], date_range [list], plot_var [str]
##############################################################################################
    
def input_setup(domain=False, date_range=False, plot_var=False, time_plot=False, grid_plot=False):
    # Get user input for spatial and temporal domains
    domain, date_range, plot_var = inputs.user_input(domain=domain, date_range=date_range, plot_var=plot_var)
    
    # Get central point of spatial domain for timezone calculation
    domain_center = [domain[0]+(domain[1] - domain[0])/2, domain[2] + (domain[3] - domain[2])/2]
    
    # Define the UTC offset in hours and adjust the date range for UTC time
    utc_offset = time_adjust.time_adjuster(domain_center[0], domain_center[1])
    local_time = [date + utc_offset for date in date_range]
    return domain, domain_center, date_range, utc_offset, local_time, plot_var, time_plot, grid_plot

##############################################################################################
# Method name:      grid_setup
# Method objective: Set up grid that serves as basis for model data structure and geographic info.
# Input(s):         domain [list]
# Outputs(s):       lats [2D ndarray], lons [2D ndarray], goes_idx [1D ndarray], 
#                   lcd [2D ndarray], water_pixel [2D ndarray], h_0 [2D ndarray], 
#                   elevation_grid [2D ndarray]
##############################################################################################

def grid_setup(domain):
    # Define coordinate grid based on defined spatial domain
    lats, lons, goes_idx = goes_grid_setup.goes_grid(domain)
    print('Grid set up...')
    
    # Generate land cover data grid based on defined domain
    lcd, water_pixel = land_cover.land_cover_grid(goes_idx)
    print('Land cover grid generated...')
    
    # Generate canopy height data grid based on land cover types
    h_0 = canopy_height.canopy_height_grid(lcd)
    print('Element height grid generated...')
    
    # Generate elevation data grid
    elevation_grid = elevation.elevation_generator(lats, lons, domain)
    print('Elevation grid generated...')
    return lats, lons, goes_idx, lcd, water_pixel, h_0, elevation_grid

##############################################################################################
# Method name:      wrf_setup
# Method objective: Collect WRF data for the defined spatial and temporal domains.
# Input(s):         lats [2D ndarray], local_time [list], wrf [bool]
# Outputs(s):       wrf_data [2D ndarray]
##############################################################################################

def wrf_setup(lats, local_time, wrf=False):
    # WRF wind speed, HRRR (u_wrf)
    # Generate 3D data matrix, with a different time at each z-index of the grid
    wrf_data = np.empty([np.shape(lats)[0], np.shape(lats)[1], len(local_time)])
    if wrf:
        wrf = wrf.lower()
        for (z, date) in enumerate(local_time):
            if wrf == 'hrrr':
                hrrr.download_HRRR_subset(date, ':WIND:')
                hrrr_file = hrrr.hrrr_data_access(date)
                hrrr_lats, hrrr_lons, _hrrr_data, _, _ = hrrr.hrrr_data_processing(hrrr_file, domain)  
                wrf_data[:, :, z] = goes_grid_setup.grid_alignment([lons, lats], [hrrr_lons, hrrr_lats, _hrrr_data])
        return wrf_data
    else:
        return None

##############################################################################################
# Method name:      physics_setup
# Method objective: Perform physics-based calculations or algorithms.
# Input(s):         See above for explanations of each variable.
# Outputs(s):       T_s [2D ndarray], T_a [2D ndarray], Q_H [2D ndarray]
# Note:             The outputs of this function can easily be manipulated to be whatever
#                   the use case is. The output must be in a 3D ndarray, however.
##############################################################################################

def physics_setup(lats, lons, date_range, local_time, domain_center, goes_idx, bucket, h_0, p_air, u_r, T_dew):
    # Land surface temperature (T_s)
    # Generate 3D data matrix, with a different time at each z-index of the grid
    T_s = np.empty([np.shape(lats)[0], np.shape(lats)[1], len(date_range)])
    # Note: z is used as variable name to maintain intuition for the depth of the grid
    for z, date in enumerate(date_range): 
        goes_product, goes_varnames = t_lst.data_download(domain_center, date, bucket)
        for i in range(0, lats.shape[0]):
            for j in range(0, lats.shape[1]):
                # If the current pixel is over water, set data to nan
                if water_pixel[i, j]:
                    T_s[i, j, z] = np.nan
                else:
                    pixel = [lats[i, j], lons[i, j]]
                    idx = [goes_idx[0] + i, goes_idx[0] + i + 1, goes_idx[2] + j, goes_idx[2] + j + 1]
                    T_s[i, j, z] = t_lst.lst(goes_product, goes_varnames, date, idx)
    print('Land surface temperature calculation complete...')
        
    # Air temperature at 2m AGL (T_a)
    # Issue: LST and air temperature pixels don't always line up.
    # Cause: Discrepancy between NLCD and elevation land classification styles.  
    
    # Generate 3D data matrix, with a different time at each z-index of the grid
    T_a = np.empty([np.shape(lats)[0], np.shape(lats)[1], len(date_range)])
    for z, date in enumerate(local_time): # Note: z is used to maintain intuition for the depth of the grid
        print(date.strftime('%Y-%m-%d %H:%M:%S'))
        for i in range(0, lats.shape[0]):
            for j in range(0, lats.shape[1]):
                # If the current pixel is over water, set data to nan
                if water_pixel[i, j]:
                    T_a[i, j, z] = np.nan
                else:
                    pixel = [lats[i, j], lons[i, j]]
                    T_a[i, j, z] = t_air_2m.air_temp(domain, pixel, date, goes_idx, [i, j], lcd[i, j], T_s[i, j, z])
    print('Air temperature calculation complete...')
    
    # Sensible heat flux (Q_H)
    # Generate 3D data matrix, with a different time at each z-index of the grid
    Q_H = np.empty([np.shape(lats)[0], np.shape(lats)[1], len(date_range)])
    for z, date in enumerate(date_range): # Note: z is used to maintain intuition for the depth of the grid
        for i in range(0, lats.shape[0]):
            for j in range(0, lats.shape[1]):
                Q_H[i, j, z] = qh.hfx(5*h_0[i, j], h_0[i, j], p_air[z], u_r[z], T_s[i, j, z], T_a[i, j, z], T_dew[z])
                
    return T_s, T_a, Q_H

##############################################################################################
# Method name:      vdtn_setup
# Method objective: Get data from observation station nearest to spatial domain's center.
# Input(s):         lats [2D ndarray], date_range [list], domain_center [list], plot_var [str]
# Outputs(s):       station_data [3D ndarray]
##############################################################################################

def vdtn_setup(lats, date_range, domain_center, plot_var):
    # Generate 3D data matrix, with a different time at each z-index of the grid
    station_data = np.empty([np.shape(lats)[0], np.shape(lats)[1], len(date_range)])
    
    # Choose observation station for validation based on domain center distance to station
    station_path = functions.station_finder(domain_center)
    
    # Get data for selected station
    obs = mesonet.csv_reader(date_range, station_path)
    for z, date in enumerate(date_range): # Note: z is used to maintain intuition for the depth of the grid
        for i in range(0, lats.shape[0]):
            for j in range(0, lats.shape[1]):
                station_data[i, j, z] = obs[plot_var][z]
                
    return station_data

def timeseries_plot(date_range, lats, lons, crd, var_list, plot_var):
    # Retrieve grid point for selected coordinate
    i, j = goes_grid_setup.crd_idx(crd, lats, lons)
    
    # List one variable listed in each row of var_list for metadata gathering purposes
    timeseries_vars = []
    for var in var_list:
        timeseries_vars.append(var[0])
        
    # Get metadata for each variable
    timeseries_metadata = [] 
    for var in timeseries_vars:
        timeseries_metadata.append(varnames.var_metadata(var, plot_var))
        
    # Get point-based data for corresponding coordinate
    data_list = [[None for _ in range(len(var_list[0]))] for _ in range(len(var_list))]

    for ii in range(0, len(var_list)):
        for jj in range(0, len(var_list[ii])):
            var_data = globals()[var_list[ii][jj]]
            data_list[ii][jj] = var_data[i, j]
            
    # Plot the timeseries
    coord_timeseries = timeseries.time_plot(date_range, data_list, timeseries_metadata)
    return coord_timeseries

def gridmap_plot():
    
    return None

if __name__ == '__main__':
    t = time.time()
    
    # Grab inputs for spatial and temporal domains and variable of interest
    # To be removed when published
    domain = [40.50, 40.90, -74.15, -73.75]
    sample_coordinate = [40.6305, -73.9521]
    date_range = [datetime.datetime(year=2019, month=7, day=28, hour=17),
                  datetime.datetime(year=2019, month=7, day=28, hour=20)-datetime.timedelta(hours=1)]
    date_range = pandas.date_range(start=date_range[0], end=date_range[1], freq='H') 
    plot_var = 'Q_H'
    
    # Set up model inputs (either from the user or from pre-defined inputs)
    domain, domain_center, date_range, utc_offset, local_time, plot_var, time_plot, grid_plot = input_setup(domain=domain, date_range=date_range, plot_var=plot_var)
    
    # Set up grid base for the model
    lats, lons, goes_idx, lcd, water_pixel, h_0, elevation_grid = grid_setup(domain)
    
    # Pull ASOS data for defined spatial domain and central point for domain
    u_r, T_dew, p_air = asos_data_reader.data_read(local_time, domain_center, utc_offset)
    
    # Pull GOES-16 data for defined spatial and temporal domains
    bucket = gcp_bucket()
    
    # Grab WRF data if the user decides to incorporate it
    wrf_data = wrf_setup(lats, local_time, wrf=False)
    
    # Calculate physics-based variables
    T_s, T_a, Q_H = physics_setup(lats, lons, date_range, local_time, domain_center, goes_idx, bucket, h_0, p_air, u_r, T_dew)
    
    # Grab nearest observation data for a given variable
    station_data = vdtn_setup(lats, date_range, domain_center, plot_var)
        
    ''' POST-PROCESSING '''
    # To-do:
    # 1. Create single formatting file for all text
    
    # Timeseries plotting, point-based
    if (date_range[-1] - date_range[0]).total_seconds() >= 3600:
        var_list = [['Q_H', 'station_data'],
                    ['T_s', 'T_a']]
        coord_timeseries = timeseries_plot(date_range, lats, lons, sample_coordinate, var_list, plot_var)

    # Spatial gridded data plotting, domain-based    
    grid_var = globals()[plot_var]
    grid_metadata = varnames.var_metadata(plot_var, plot_var)    
    datamap = meshgrid.main(domain, domain_center, date_range, lats, lons, grid_var, grid_metadata, zoom=0.0, animated=False, savefig=False)
  
    runtime = time.time() - t
    print('Program runtime is: {0:.3f} sec'.format(runtime))
    