### asos_data_processing branch

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
import pandas as pd
from datetime import datetime as dt
import matplotlib.pyplot as plt
import asos_data_reader as asos
import mesonet_reader as meso
import sensible_heat_flux as qh
import plotter
# import goes_air_temperature.goes_16_air_temperature as g16 # Old script
import goes_air_temperature.goes_16_air_temp as g16 # Revised script

### 1. Constants
rho = 1.225     # Air density at station (kg/m^3) - note: p = rho*R*T
c_p = 1006      # Specific heat capacity of air at station (J/kg-K)
vk = 0.4        # Von Karman constant
g = 9.81        # Gravitational acceleration (m/s^2)
R = 287.05      # Gas constant

# Main function definition
def main(start_date, end_date, site):
    
    p_air = 1013.25 # Air pressure, MSLP
    # Define location dictionary with site-specific properties
    # Height above sea level for Bronx, Brooklyn, Manhattan, Queens, Staten Island MesoNET sites
    site_info = {'BRON': [[11, 26], 57.5, 7.5],
                 'BKLN': [[19, 22], 33.2, 7.5],
                 'MANH': [[15, 22], 94.8, 10],
                 'QUEE': [[15, 28], 54.6, 7.5],
                 'STAT': [[20, 15], 33.1, 5]}
    
    # Specify date range of interest
    # Format: YYYYMMDDHHMM
    [start_date_fmtd, end_date_fmtd] = [dt.strptime(str(start_date), '%Y%m%d%H%M').strftime('%Y-%m-%d %H:%M:%S'),
                                        dt.strptime(str(end_date), '%Y%m%d%H%M').strftime('%Y-%m-%d %H:%M:%S')]
    date_range = pd.date_range(start_date_fmtd, end_date_fmtd, (end_date-start_date)*24/10000+1)
    
    # Iterate over date range, incrementing by 1 hour
    date_list, T_lst_list, T_air_list, h_0_list = [], [], [], []
    for t in date_range:
        t_int = int(t.strftime('%Y%m%d%H%M'))
        l, a, h = g16.temperature_data(t_int, site_info[site][0][0], site_info[site][0][1], site)
        date_list.append(t.to_pydatetime())
        T_lst_list.append(l)
        T_air_list.append(a)
        h_0_list.append(h)
    agg_df = pd.DataFrame({'date': date_list, 
                           'T_lst': T_lst_list, 
                           'T_air': T_air_list, 
                           'h_0': h_0_list})
    # del date_list, T_lst_list, T_air_list, h_0_list
    
    # Access and process data from ASOS FTP for given dates
    asos_data = asos.data_read(start_date, end_date, site)
    agg_df['u_r'] = asos_data['u_r']
    # agg_df['T_air'] = asos_data['T_air']
    agg_df = agg_df[:-1] # Cut off last row to remove hanging timestamp (midnight on next day)
    
    # Pull Mesonet data and assign to 'q_h_obs' column in aggregate DataFrame
    mesonet_data = meso.csv_reader(start_date, end_date, ("mesonet_data/" + site))
    agg_df['q_h_obs'] = mesonet_data['H'][::2].tolist()
    agg_df['u_star_obs'] = mesonet_data['USTAR'][::2].tolist()
    # agg_df['T_air'] = mesonet_data['T_air'][::2].tolist()
    agg_df['T_air_meso'] = mesonet_data['T_air'][::2].tolist()
    
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
    
    plotter.agg_plot(agg_df, site, savefig=True, t_air_src='hrisko-model')
    
    return agg_df

    ## End iterative method

runtime = time.time()
model_data = main(201907260000, 201907290000, 'MANH') 
print('Program runtime: %.4f s' % (time.time()-runtime))