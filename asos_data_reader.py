### asos_data_processing branch

### Objective
# The objective of this script is to read relevant parameters from observation data logs from NOAA ASOS stations.

### Inputs: 
# Station name (str)
# Start date (int)
# End date (int)
### Outputs:
# Pandas DataFrame with air temperature (T_air, K), wind speed (u, m/s), and date data

# Status: currently testing single cases to develop generalized algorithm.

##########################################################################################
# Imports
import pandas as pd
import re
import time
from datetime import datetime as dt
from datetime import timedelta
from time_adjust import time_adjust
import os
import numpy as np
import sys

### Auxiliary methods for unit conversion
# Convert Celsius to Kelvin
def CtoK(T):
    return T + 273.15
# Convert knots to meters per second
def KTtoMS(u):
    return u*0.51444

# Boolean to control strings with performance data printing to the console. If True, strings should print.  
str_switch = False

#############################################################################
### ASOS Station Finder
# Objective: Find ASOS station closest to given latitude and longitude
# Input: latitude (float) and longitude (float)
# Output: ASOS station code

# Define crdal ASOS station list filepath. Uncomment the second definition for a path to the online text file
asos_station_fp = os.getcwd() + '/asos_data/asos-stations.txt'
# asos_station_fp = r'https://www.ncdc.noaa.gov/homr/file/asos-stations.txt'

# Define relevant columns
asos_cols = ['CALL', 'NAME', 'LAT', 'LON', 'ELEV', 'UTC']

R = 6378137 # GRS80 semi-major axis of Earth, per GOES-16 PUG-L2, Volume 5, Table 4.2.8
    
## Great circle distance calculator
# Function takes coordinates for point of interest (crd) and ASOS station (asos) and calculates distance
def distance(lat_crd, lon_crd, lat_asos, lon_asos):
    p = np.pi/180
    a = 0.5 - np.cos((lat_asos-lat_crd)*p)/2 + np.cos(lat_crd*p) * np.cos(lat_asos*p) * (1-np.cos((lon_asos-lon_crd)*p))/2
    return 2*R*np.arcsin(np.sqrt(a))

def asos_find(lat_crd, lon_crd):
    t = time.time()
    # Read data into ASOS DataFrame (adf)
    adf = pd.read_fwf(asos_station_fp, usecols=asos_cols).drop([0])
    # Filter stations with null data
    adf = adf[adf['ELEV'].astype(float) > -99999]
    # Read relevant parameters into lists for iterative purposes
    stations, lat_asos, lon_asos = [adf['CALL'].tolist(), adf['LAT'].astype(float).tolist(), adf['LON'].astype(float).tolist()]
    
    dists = 2*np.pi*R # Define arbitrarily large number as an initial condition (rouglhy equal to circumference of Earth)
    station = '' # Initialize empty string for population
    # Iterate over 
    for i in range(1, len(lat_asos)):
        dist = distance(lat_crd, lon_crd, lat_asos[i], lon_asos[i])
        if dist < dists:
            dists = dist
            station = 'K' + stations[i]

    if str_switch:
        print("Closest station: %s, %.2f m away" % (station, dists))
        print('asos_find runtime: %.4f s' % (time.time() - t))
    # asos_find(25.8223, -80.2895) # Function call for troubleshooting purposes. Move out of function if used.
    return station
            
#############################################################################
### Data Read
# Objective: Process data from ASOS .dat files and return DataFrame with selected parameters
# Input: start date (int), end date (int), coordinates (str)
# Output: Pandas DataFrame

def data_read(start_date, end_date, crd):
    
    # Adjust datetimes to timezone corresponding to location of interest
    start_date, end_date = [start_date + time_adjust(crd), end_date + time_adjust(crd)]
    # Generate ints from dates for comparison purposes
    date_str = [dt.strftime(start_date, '%Y%m%d%H%M'),
                dt.strftime(end_date, '%Y%m%d%H%M')]
    
    print(date_str)
    
    # Create file path based on nearest ASOS station
    data_url = 'ftp://ftp.ncdc.noaa.gov/pub/data/asos-fivemin/6401-' + \
    date_str[0][0:4] + '/' + '64010' + asos_find(crd[0], crd[1]) + date_str[0][0:6] + '.dat'
    print(data_url)
    
    # Try connecting to ASOS FTP using data_url.
    try:     
        # Import data to dataframe. Note that column 6 is wind information, column 9 is air temperature
        data = pd.read_table(data_url)
        print('Connected to FTP!')
    except:
        print('Cannot connect to FTP :(')
        sys.exit()
        
    # Backup to re-direct to local files if ASOS FTP is down
    # if crd == "STAT":
    #     data_url = "asos_data/64010KEWR201907.dat"
    # elif crd == "BKLN":
    #     data_url = "asos_data/64010KJFK201907.dat"
    # elif crd == "QUEE" or crd == "BRON" or crd == "MANH":cr
    #     data_url = "asos_data/64010KLGA201907.dat"
        
    
    # Date range where element 0 is start data and 1 is end date
    # Date format: YYYYMMDDHHMM
    date_range = [int(date_str[0]), int(date_str[1])]
    
    # Read data for air temperature, T_air, and wind speed, u
    date_list, station_list, T_list, u_list = [], [], [], []
    
    ## Regex patterns for data mining through the .dat file(s)
    # Air temperature regex: string of 6 characters "(0-9)(0-9)/(0-9)(0-9)" bounded by 2 spaces
    T_pattern = r'\s.?\d\d[+-/].?\d\d\s'
    # Wind speed regex: string of 6 characters "(0-9)(0-9)KT " bounded by 2 numbers and a space
    # Note: This definition ignores gusts
    u_pattern = r"\s\d\d\d\d\d\D"
    # Note: This definition allows the gust becomes the effective wind speed
    # u_pattern = r"\d\d[K][T]\s\d"
    # Sea level pressure: string of 6 characters "(0-9)(0-9)/(0-9)(0-9)" bounded by space and number
    # slp_pattern = r"\d\d[K][T]\s\d"
    
    # Iterate through DataFrame and filter out results
    for row in data.iloc[::1, 0]:
        if date_range[0] <= int(row[13:25]) <= date_range[1]:
            if re.findall(T_pattern, row) and re.findall(u_pattern, row):
                T_air_str = re.findall(T_pattern, row)[0]
                if T_air_str[1] == 'M':
                    T_list.append(-int(T_air_str[2:4]))
                else:
                    T_list.append(int(T_air_str[1:3]))
                u_str = re.findall(u_pattern, row)[0]
                u_list.append(int(u_str[4:6])) # If gusts not accounted for, use this
                # u_list.append(u_str[0:2]) # If gusts accounted for, use this
                station_list.append(row[5:9])
                date_list.append(dt.strptime(row[13:25], '%Y%m%d%H%M'))
                
    # Perform type conversions for further data manipulation
    # [::12]s added to select hourly data (12 entries per hour, 5-minute ASOS data)
    T_list = list(map(int, T_list[::12]))
    T_list = [CtoK(T) for T in T_list]
    u_list = list(map(int, u_list[::12]))
    u_list = [KTtoMS(u) for u in u_list]
    
    # Add data to output DataFrame
    df = pd.DataFrame()
    df['station'] = station_list[::12]
    df['date'] = date_list[::12]
    df['T_air'] = T_list
    df['u_r'] = u_list
    return df

def main(start_date, end_date, crd):
    data_read(start_date, end_date, crd)
    
if __name__ == "__main__":
    main()
