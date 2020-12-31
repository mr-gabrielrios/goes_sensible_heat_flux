### Objective
# The objective of this script is to generate the sensible heat flux based on input from GOES data.

### Program flow
# 0. Define imports
# 1. Define constants
# 2. Define dynamic inputs (e.g. model data, station data)
# Note: The dynamic inputs are currently static, for simplicity, while the algorithm is finalized
# 3. Define initial conditions
# 4. Define data storage variables and other program administrative parameters
# 5. Calculate secondary parameters to support solution algorithm
# 6. Perform solution algorithm
# Note: This is currently an iterative algorithm 

# CRITERIA: Solution criteria for each loop is that each variable doesn't change by > 1% per iteration

### 0. Imports
import time
import glob, sys
import pandas as pd
import numpy as np
from datetime import datetime as dt
from datetime import timedelta as td
import matplotlib.pyplot as plt
import asos_data_reader as asos
import ameriflux_data_processing as af
import mesonet_reader as meso
import sensible_heat_flux as qh
import plotter
import goes_air_temperature.goes_16_air_temperature as g16 # Revised script
import colormap_plot
# import anim

### 1. Constants
rho = 1.225     # Air density at station (kg/m^3) - note: p = rho*R*T
c_p = 1006      # Specific heat capacity of air at station (J/kg-K)
vk = 0.4        # Von Karman constant
g = 9.81        # Gravitational acceleration (m/s^2)
R = 287.05      # Gas constant

# Main function definition

def single_pixel(start_date, end_date, site):
    
    # Import Ameriflux data characteristics for site instrumentation
    # Source: https://ameriflux.lbl.gov/data/measurement-height/
    if site[0:2] == 'US':
        ameriflux_site_data = pd.read_csv('ameriflux_data/BASE_MeasurementHeight_20201218.csv')
        # Try to fetch measurement height (z_r). If not found, use median height of all towers
        if ameriflux_site_data.loc[(ameriflux_site_data['Site_ID'] == site) & (ameriflux_site_data['Variable'] == 'TA_1_1_1')].empty:
            z_r = np.nanmedian(ameriflux_site_data.loc[ameriflux_site_data['Variable'] == 'TA_1_1_1']['Height'])
        else:
            z_r = ameriflux_site_data.loc[(ameriflux_site_data['Site_ID'] == site) & (ameriflux_site_data['Variable'] == 'TA_1_1_1')]['Height']
        print(z_r)
        # Site-specific data
        site_fpath = glob.glob('goes_air_temperature/elevation_data/' + site + '*.csv')[0]
        crd = [float(site_fpath.split('_')[-3]), float(site_fpath.split('_')[-2])]
    else:
        crd = [0, 0]
    
    p_air = 1013.25 # Air pressure, MSLP
    # Define location dictionary with site-specific properties
    # Coordinates and flux tower height above sea level for Bronx, Brooklyn, Manhattan, Queens, Staten Island MesoNET sites
    site_info = {'BRON': [[40.8714, -73.8963], 57.5],
                 'BKLN': [[40.6305, -73.9521], 33.2],
                 'MANH': [[40.7678, -73.9645], 94.8],
                 'QUEE': [[40.7366, -73.8201], 54.6],
                 'STAT': [[40.6021, -74.1504], 33.1]}
                 # 'US-xNQ': [crd, 5] # Ameriflux slot
    
    # Non-iterative quantities    
    bbox = g16.site_info(site_info[site][0])
    goes_grid_info = g16.goes_grid(bbox)    
    nlcd_grid_info = g16.nlcd_grid(site_info[site][0])
    h_0 = g16.roughness_height_grid(site_info[site][0])
    
    # Specify date range of interest
    # Format: YYYYMMDDHHMM
    [start_date, end_date] = [dt.strptime(str(start_date), '%Y%m%d%H%M'),
                                        dt.strptime(str(end_date), '%Y%m%d%H%M')]
    end_date = end_date - td(hours=1)
    date_range = pd.date_range(start_date, end_date, (end_date-start_date).total_seconds()/3600+1)
    
    # Iterate over date range, incrementing by 1 hour
    date_list, T_lst_list, T_air_list, h_0_list = [], [], [], []
    utc_offset = ''
    for t in date_range:
        l, a, h, utc_offset = g16.temperature_data(t, site_info[site][0], site, bbox, goes_grid_info, nlcd_grid_info, h_0)
        date_list.append(t.to_pydatetime())
        T_lst_list.append(l)
        T_air_list.append(a)
        h_0_list.append(h)
    agg_df = pd.DataFrame({'date': date_list, 
                           'T_lst': T_lst_list, 
                           'T_air': T_air_list, 
                           'h_0': h_0_list})
    
    # Access and process data from ASOS FTP for given dates
    asos_data = asos.data_read(start_date, end_date, site_info[site][0])
    agg_df['u_r'] = asos_data['u_r']
    
    # Pull Mesonet data and assign to 'q_h_obs' column in aggregate DataFrame
    if site[0:2] != 'US':
        mesonet_data = meso.csv_reader(start_date, end_date, ("mesonet_data/" + site))
        print(agg_df)
        print(mesonet_data['H'][::2].tolist())
        agg_df['q_h_obs'] = mesonet_data['H'][::2].tolist()
        agg_df['u_star_obs'] = mesonet_data['USTAR'][::2].tolist()
        # agg_df['T_air'] = mesonet_data['T_air'][::2].tolist()
        agg_df['T_air_obs'] = mesonet_data['T_air'][::2].tolist()
    else:
        af_data = af.csv_reader(site, site_info[site][0], start_date, end_date)
        print(af_data)
        agg_df['obs_time'] = af_data['datetime']
        print(agg_df.shape)
        print(af_data.shape)
        agg_df['q_h_obs'] = af_data['H'].tolist()
        agg_df['u_star_obs'] = af_data['u_star'].tolist()
        agg_df['T_air_obs'] = af_data['T_air'].tolist()
    
    C_h_list, u_star_list, L_list, z_0m_list, zeta_list, q_h_list = [], [], [], [], [], []
    for index, row in agg_df.iterrows():
        print(row['date'])
        [C_h, u_star, L, z_0m, zeta, q_h] = qh.q_sens(site_info[site][1], row['h_0'], p_air, row['u_r'], row['T_lst'], row['T_air'])
        print("----------------------------------------------------------------------------------")
        C_h_list.append(C_h)
        u_star_list.append(u_star)
        L_list.append(L)
        z_0m_list.append(z_0m)
        zeta_list.append(zeta)
        q_h_list.append(q_h)
    agg_df['C_h'] = C_h_list
    agg_df['u_star'] = u_star_list
    agg_df['L'] = L_list
    agg_df['z_0m'] = z_0m_list
    agg_df['zeta'] = zeta_list
    agg_df['q_h'] = q_h_list
    del C_h_list, u_star_list, L_list, z_0m_list, zeta_list, q_h_list
    
    # Calculate sensible heat flux error and append to DataFrame
    agg_df['q_h_error'] = 100*((agg_df['q_h'] - agg_df['q_h_obs'])/agg_df['q_h_obs'])
    agg_df['u_star_error'] = 100*((agg_df['u_star'] - agg_df['u_star_obs'])/agg_df['u_star_obs'])
    
    plotter.agg_plot(agg_df, site, savefig=False, t_air_src='hrisko-model')
    
    return agg_df


######################################### COLORMAP WORKING ZONE

### Colormap only valid for New York City metro area
# Non-iterative quantities   
nyc_crd = [40.631762, -73.953678] # Pick Brooklyn as central coordinate for NYC Metro
res = 0.1
bbox = g16.site_info(nyc_crd, res)
goes_grid_info = g16.goes_grid(bbox)   
lats, lons, nc_indx_corners, nc_indx_bounds, lon_plot, lat_plot = goes_grid_info 
nlcd_grid_info = g16.nlcd_grid(nyc_crd, res)
h_0 = g16.roughness_height_grid(nyc_crd, res)

def colormap(start_date, end_date, crd):
    
    # Specify date range of interest
    # Format: YYYYMMDDHHMM
    [start_date, end_date] = [dt.strptime(str(start_date), '%Y%m%d%H%M'),
                                        dt.strptime(str(end_date), '%Y%m%d%H%M')]
    end_date = end_date - td(hours=1)
    date_range = pd.date_range(start_date, end_date, (end_date-start_date).total_seconds()/3600+1)
    
    # Iterate over date range, incrementing by 1 hour    
    p_air = 1013.25 # Air pressure, MSLP
    date_list, T_lst_list, T_air_list, h_0_list, prod_list = [], [], [], [], []
    utc_offset = ''
    for t in date_range:
        l, a, h, utc_offset, prods = g16.temperature_data(t, crd, 'MANH', bbox, goes_grid_info, nlcd_grid_info, h_0)
        date_list.append(t.to_pydatetime())
        T_lst_list.append(l)
        T_air_list.append(a)
        h_0_list.append(h)
        prod_list.append(prods)
    agg_df = pd.DataFrame({'date': date_list, 
                           'T_lst': T_lst_list, 
                           'T_air': T_air_list, 
                           'h_0': h_0_list})
    
    # Access and process data from ASOS FTP for given dates
    asos_data = asos.data_read(start_date, end_date, crd, nyc=True)
    agg_df['u_r'] = asos_data['u_r']
    
    C_h_list, u_star_list, L_list, z_0m_list, zeta_list, q_h_list, u_r_list = [], [], [], [], [], [], []
    for index, row in agg_df.iterrows():
        print(row['date'])
        print('Roughness height: %.2f | Air pressure: %.2f | Wind speed: %.2f | LST: %.2f | AT: %.2f'  % (row['h_0'], p_air, row['u_r'], row['T_lst'], row['T_air']))
        # Note: elevation value is an average value of all sites. Very influential to error.
        [C_h, u_star, L, z_0m, zeta, q_h] = qh.q_sens(54.6, row['h_0'], p_air, row['u_r'], row['T_lst'], row['T_air'])
        print("----------------------------------------------------------------------------------")
        C_h_list.append(C_h)
        u_star_list.append(u_star)
        L_list.append(L)
        z_0m_list.append(z_0m)
        zeta_list.append(zeta)
        q_h_list.append(q_h)
        u_r_list.append(row['u_r'])
    agg_df['C_h'] = C_h_list
    agg_df['u_star'] = u_star_list
    agg_df['L'] = L_list
    agg_df['z_0m'] = z_0m_list
    agg_df['zeta'] = zeta_list
    agg_df['q_h'] = q_h_list
    agg_df['u_r'] = u_r_list
    
    return agg_df, q_h_list, u_r_list, T_air_list, prod_list

### Run algorithm for selected dates
runtime = time.time()

# Initialize plot full of nans
start_date, end_date = [201907281200, 201907281500]
[start_date_dt, end_date_dt] = [dt.strptime(str(start_date), '%Y%m%d%H%M'),
                                        dt.strptime(str(end_date), '%Y%m%d%H%M')]
end_date_dt = end_date_dt - td(hours=1)
date_range = pd.date_range(start_date_dt, end_date_dt, (end_date_dt-start_date_dt).total_seconds()/3600+1)

hrs = int((end_date_dt-start_date_dt).total_seconds()/3600+1)

'''
### Initialize 3D arrays for parameters of interest
# Sensible heat flux
q_h_arr = np.empty([lat_plot.shape[0], lat_plot.shape[1], hrs])
u_arr = np.empty([lat_plot.shape[0], lat_plot.shape[1], hrs])
t_arr = np.empty([lat_plot.shape[0], lat_plot.shape[1], hrs])
prod_list = []

for ii in range(0, lat_plot.shape[0]-1):
    for jj in range(0, lon_plot.shape[1]-1):
        [df, q, u, at, prods] = colormap(start_date, end_date, [lat_plot[ii, jj], lon_plot[ii, jj]])
        q_h_arr[ii, jj] = q
        u_arr[ii, jj] = u
        t_arr[ii, jj] = at
        # Save product data to prod_list only if prod_list is empty
        if not prod_list:
            prod_list = prods
'''
            
### Plot data in 3D array(s)          
# for i in range(0, hrs):
#     colormap_plot.colormap_plot(q_h_arr[:,:,i], nyc_crd, goes_grid_info, prod_list[i], "q_h")

# anim.colormap_plot(q_h_arr, nyc_crd, goes_grid_info, prod_list, "q_h", res)
        
########################################## END COLORMAP WORKING ZONE

# ## Single run
# model_data = single_pixel(201907260000, 201907270000, 'BKLN')

# ## Multiple runs
# station_list = ['BKLN', 'QUEE', 'STAT']

# # Read in GOES-16 Clear Day Log, where dates with continuous clear days are recorded
# rundates = pd.read_csv('goes_air_temperature/aux/goes-16_clear_day_log.csv', header=3, usecols=['Dates of Interest'])
# # Remove intermediate dates or blank rows
# rundates = rundates[rundates['Dates of Interest'].str.len() > 1]
# # Reformat date entries to match desired input formats
# rundates = rundates["Dates of Interest"].tolist()
# rundates = [int(s.replace('/', '') + '0000') for s in rundates]
# rundates = np.reshape(rundates, (int(len(rundates)/2), 2))
# # For each pair of dates in the 'rundates' matrix, run the algorithm for each station
# for station in station_list:
#     # Run 0-22, 23 onwards
#     for i in range(0, 21):
#         model_data = single_pixel(rundates[i][0], rundates[i][1], station)

print('Program runtime: %.4f s' % (time.time()-runtime))