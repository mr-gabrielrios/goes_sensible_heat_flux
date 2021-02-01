'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        Mesonet Reader
Package file path:  ~/vdtn/obs/mesonet/reader.py
Objective:          Scrape relevant data from New York State Mesonet network stations (used for New York City flux validation)
Author:             Gabriel Rios
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

import datetime, numpy as np, os, pandas as pd

##############################################################################################
# END IMPORTS
##############################################################################################
  
##############################################################################################
# Method name:      file_sort
# Method objective: Set key for sorting files chronologically 
#                   (first 8 characters are date in YYYYMMDD format)
# Input(s):         fn [str]
# Outputs(s):       fn[0:8] [str]
##############################################################################################

def file_sort(fn):
    return(fn[0:8])

##############################################################################################
# Method name:      CtoK
# Method objective: Convert temperature in Celsius to Kelvin
# Input(s):         T [float]
# Outputs(s):       T + 273.15 [float]
##############################################################################################

def CtoK(T):
    return T + 273.15

##############################################################################################
# Method name:      csv_reader
# Method objective: Scrape relevant data from CSVs in a given directory for a temporal domain.
# Input(s):         T [float]
# Outputs(s):       T + 273.15 [float]
##############################################################################################

def csv_reader(date_range, data_dir):
    
    # Define CSV columns to extract data from
    # Datetime, sensible heat flux, friction velocity, air temperature (deg C)
    cols = ['datetime', 'H', 'USTAR', 'Tc']
    
    # Initialize empty DataFrame with predefined columns
    data = pd.DataFrame(columns=cols)
    
    # Sort files by date, assuming standard Mesonet filename convention
    # Example: YYYYMMDD_PARAMETER_LOCATION_Parameter_NYSMesonet.csv
    file_list = sorted(os.listdir(data_dir), key=file_sort)
    
    # Specify date range of interest
    # Format: YYYYMMDDHHMM
    date_range = pd.date_range(start=date_range[0], end=date_range[-1], freq='30min')
    # data['datetime'] = date_range
    
    # Iterate through sorted file list and extract data within date range with daily resolution
    for i, file in enumerate(file_list): 
        file_date = datetime.datetime.strptime(os.path.splitext(file)[0][0:8], '%Y%m%d')
        # Reduce datetimes to daily resolution to work on down-filtering
        days = [datetime.datetime.strptime(date_range[0].strftime('%Y%m%d'), '%Y%m%d'),
                datetime.datetime.strptime(date_range[-1].strftime('%Y%m%d'), '%Y%m%d')]
        # Filter files by day - any day within the date range will have a corresponding file
        if days[0] <= file_date <= days[-1]:
            filename = os.path.join(data_dir, file)
            data = data.append(pd.read_csv(filename, usecols=cols))
    
    # Convert date strings to datetime data type
    data['datetime'] = pd.to_datetime(data['datetime'])
    # Filter data entries by full datetime (includes hours and minutes)
    data = data[(data['datetime'] >= date_range[0]) & (data['datetime'] <= date_range[-1])] 
    data = data.reset_index(drop=True)
    
    # Account for missing observation data by inserting nans
    for i, date in enumerate(date_range):
        nanrow = pd.DataFrame([[np.nan] * len(cols)], columns=cols)
        nanrow['datetime'] = date
        if data.loc[i, ['datetime']].item() != date:
            data = pd.concat([data.iloc[:i], nanrow, data.iloc[i:]]).reset_index(drop=True)
        
    # Re-cast numerical strings as floats
    data['H'] = data['H'].astype(float)
    data['USTAR'] = data['USTAR'].astype(float)
    data['Tc'] = data['Tc'].astype(float)
    data['Tc'] = [CtoK(i) for i in data['Tc']]
    # Match parameter names to model parameter names
    data.rename(columns = {'H': 'Q_H',
                           'USTAR': 'u_star',
                           'Tc': 'T_air'}, inplace = True) 
    data = data.iloc[::2].reset_index(drop=True)
    return data
    
if __name__ == "__main__":
    date_range = [datetime.datetime(year=2019, month=7, day=28, hour=5),
                  datetime.datetime(year=2019, month=7, day=29, hour=5)-datetime.timedelta(hours=1)]
    date_range = pd.date_range(start=date_range[0], end=date_range[1], freq='H') 
    data_dir = os.path.join(os.path.dirname(__file__), 'data/BKLN/')
    data = csv_reader(date_range, data_dir)