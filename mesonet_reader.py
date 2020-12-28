### Objective:
# Read data from NYS Mesonet sites

import pandas as pd
import os
from datetime import datetime as dt
from datetime import timedelta
from time_adjust import time_adjust
import numpy as np

### Auxiliary method for file sorting by date in filename
def file_sort(fn):
    return(fn[0:8])

### Auxiliary method for unit conversion
# Convert Celsius to Kelvin
def CtoK(T):
    return T + 273.15

def csv_reader(start_date, end_date, data_dir):

    [start_date, end_date] = [start_date.strftime('%Y%m%d'), end_date.strftime('%Y%m%d')]
    
    data_cols = ['datetime', 'H', 'USTAR', 'Tc'] # Datetime, sensible heat flux, friction velocity, air temperature
    
    # Initialize empty array to hold file data
    data_list = []
    
    # Sort files by date, assuming standard Mesonet filename convention
    # Example: YYYYMMDD_PARAMETER_LOCATION_Parameter_NYSMesonet.csv
    file_list = sorted(os.listdir(data_dir), key=file_sort)
    start_date = int(start_date)
    end_date = int(end_date)
    
    # Specify date range of interest
    # Format: YYYYMMDDHHMM
    # [start_date, end_date] = [201906140000, 201906150000]
    [start_date_fmtd, end_date_fmtd] = [dt.strptime(str(start_date), '%Y%m%d').strftime('%Y-%m-%d %H:%M:%S'),
                                        dt.strptime(str(end_date+1), '%Y%m%d').strftime('%Y-%m-%d %H:%M:%S')]
    date_range = pd.date_range(start_date_fmtd, end_date_fmtd, (end_date+1-start_date)*48+1)
    date_range = [t.to_pydatetime() for t in date_range]
    
    # Iterate through sorted file list and extract data within date range
    for file in file_list:       
        if start_date <= int(os.path.splitext(file)[0][0:8]) <= end_date:
            # print(os.path.splitext(file)[0][0:8])
            filename = data_dir + "/" + file
            print(filename)
            data_list.append(pd.read_csv(filename, usecols=data_cols))
            
    # Concatenate all dataframes in data_list
    data = pd.concat(data_list)
    
    ii = 0 # Parallel iterand for date_range
    idx_list = []
    # Iterate through each date and reformat each entry
    for i in range(0,data.shape[0]):        
        data.iloc[i, 0] = dt.strptime(data.iloc[i, 0], '%Y-%m-%d %H:%M:%S.%f').strftime('%Y-%m-%d %H:%M:%S')
        # If Mesonet data is missing, log the index so that the next loop can fill in missing data with nans
        if str(data.iloc[i, 0]) != str(date_range[ii]):
            print('Missing date: ', data.iloc[i, 0])
            idx_list.append(i)
            ii += 1
        i += 1
        ii += 1
    
    # Fill in missing data with nans for entries with missing timestamps
    for i in idx_list:        
        line = pd.DataFrame({'datetime': date_range[i], 'H': np.nan, 'USTAR': np.nan,'Tc': np.nan}, index=[i])
        data = pd.concat([data.iloc[:(i-1)], line, data.iloc[(i-1):]]).reset_index(drop=True)
    
    # Re-cast datetime strings as datetime objects
    data['datetime'] = pd.to_datetime(data['datetime'])
    # Re-cast numerical strings as floats
    data['H'] = data['H'].astype(float)
    data['USTAR'] = data['USTAR'].astype(float)
    data['Tc'] = data['Tc'].astype(float)
    data['Tc'] = [CtoK(i) for i in data['Tc']]
    data.rename(columns = {'Tc':'T_air'}, inplace = True) 
    
    return data
        
def main(start_date, end_date, data_dir):
    csv_reader(start_date, end_date, data_dir)
    
if __name__ == "main":
    main()