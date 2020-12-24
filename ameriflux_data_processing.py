###################################################################################################
### Ameriflux Data Processing
# Objective:    Extract relevant data for use in sensible heat flux calculation
###################################################################################################

from datetime import datetime as dt
from datetime import timedelta as td
from time_adjust import time_adjust
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

###################################################################################################
### CtoK
# Objective:    Convert Celsius to Kelvin
# Input(s):     Temperature in Celsius (int or float)
# Output(s):    Temperature in Kelvin (float)
###################################################################################################
def CtoK(T):
    return T + 273.15

###################################################################################################
### CSV Reader
# Objective:    Process site-specific raw data and extract relevant parameters into DataFrame
# Input(s):     Site name (str), start date (int or str in YYYYmmDDHHMM format), end date (int or str in YYYYmmDDHHMM format)   
# Output(s):    DataFrame
###################################################################################################
def csv_reader(site, crd, start_date, end_date):
    
    # Adjust datetimes to timezone corresponding to location of interest
    start_date, end_date = [start_date + time_adjust(crd), end_date + time_adjust(crd)]
    
    # Find file corresponding to selected site
    ameriflux_fpath = 'ameriflux_data/*'
    file = glob.glob(ameriflux_fpath + site + '*.csv')[0]
    print(start_date, end_date)
    
    ### DataFrame generation, data filtering, and type conversions
    
    # List of parameters of interest. See Ameriflux README for more information (linked below).
    # https://ftp.fluxdata.org/.ameriflux_downloads/BASE/README_AmeriFlux_BASE.txt
    # Timestamp, sensible heat flux, atmospheric pressure, relative humidity, air temperature, friction velocity, wind speed
    # Use the try/except based on the site-to-site variation in column names for air temperature and wind speed
    try:
        col_list = ['TIMESTAMP_START', 'H', 'PA', 'RH', 'TA_1_1_1', 'USTAR', 'WS_1_1_1']
        data = pd.read_csv(file, header=2, usecols=col_list)
        data = data.rename(columns = {'TIMESTAMP_START': 'datetime', 
                                      'TA_1_1_1': 'T_air', 
                                      'USTAR': 'u_star',
                                      'WS_1_1_1': 'u_r'}) 
    except:
        col_list = ['TIMESTAMP_START', 'H', 'PA', 'RH', 'TA', 'USTAR', 'WS']
        data = pd.read_csv(file, header=2, usecols=col_list)
        data = data.rename(columns = {'TIMESTAMP_START': 'datetime', 
                                      'TA': 'T_air', 
                                      'USTAR': 'u_star',
                                      'WS': 'u_r'}) 
    
    # Remove null sensible heat flux data
    data = data[data > -9999] 
    # Cast times to datetime format
    data['datetime'] = pd.to_datetime(data['datetime'], format='%Y%m%d%H%M') 
    
    # Convert air temperature to Kelvin
    data['T_air'] = CtoK(data['T_air'])
    
    # Date filtering
    # start_date, end_date = [pd.to_datetime(str(start_date), format='%Y%m%d%H%M'), pd.to_datetime(str(end_date), format='%Y%m%d%H%M')]
    data = data[(data['datetime'] > start_date) & (data['datetime'] < end_date)]
    
    date_df = pd.DataFrame()
    date_range = pd.date_range(start_date, end_date, (end_date-start_date).total_seconds()/3600+1)
    
    date_df = date_df.reindex(date_range, fill_value=np.nan)
    date_df = date_df.reset_index().rename(columns={'index': 'datetime'})
    
    data = pd.merge(data, date_df, how='right', on=['datetime'])   
    data = data.sort_values(by=['datetime'])
    data = data.reset_index(drop=True)
    data = data[:-1] # Drop nan row at end due to reindexing
    
    # Plotting for troubleshooting
    # nanmask = np.isfinite(data['H']) # Masking to plot over nans
    # plt.plot(data['datetime'][nanmask], data['H'][nanmask])
    
    return data

def main(site, start_date, end_date):
    csv_reader(site, start_date, end_date)
    
if __name__ == "main":
    main()

# Trial function
# data = csv_reader('US-ONA', [31.19484, -84.46861], dt.strptime(str(201901010000), '%Y%m%d%H%M'), dt.strptime(str(201912310000), '%Y%m%d%H%M'))