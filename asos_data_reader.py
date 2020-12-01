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
from datetime import datetime 
from datetime import timedelta

# Start script timer
t = time.time()

### Auxiliary methods for unit conversion
# Convert Celsius to Kelvin
def CtoK(T):
    return T + 273.15
# Convert knots to meters per second
def KTtoMS(u):
    return u*0.51444

#############################################################################

def data_read(start_date, end_date, loc):

    # Brooklyn flux data
    # Select data for JFK for June 2019 as a test case
    # TO DO: write algorithm to based FTP data access based on start and end dates specified
    if loc == "STAT":
        data_url = "asos_data/64010KEWR201907.dat"
        # data_url = "ftp://ftp.ncdc.noaa.gov/pub/data/asos-fivemin/6401-2019/64010KEWR201906.dat"
    elif loc == "BKLN":
        data_url = "asos_data/64010KJFK201907.dat"
        # data_url = "ftp://ftp.ncdc.noaa.gov/pub/data/asos-fivemin/6401-2019/64010KJFK201906.dat"
    elif loc == "QUEE" or loc == "BRON" or loc == "MANH":
        # data_url = "ftp://ftp.ncdc.noaa.gov/pub/data/asos-fivemin/6401-2019/64010KLGA201906.dat"
        data_url = "asos_data/64010KLGA201907.dat"
        
    # Import data to dataframe. Note that column 6 is wind information, column 9 is air temperature
    data = pd.read_table(data_url)
    # Date range where element 0 is start data and 1 is end date
    # Date format: YYYYMMDDHHMM
    date_range = [start_date, end_date]
    # Read data for air temperature, T_air, and wind speed, u
    date_list = []
    station_list = []
    T_list = []
    u_list = []
    
    ## Regex patterns for data mining through the .dat file(s)
    # Air temperature regex: string of 6 characters "(0-9)(0-9)/(0-9)(0-9)" bounded by 2 spaces
    T_pattern = r"\s\d\d[^a-z0-9\s]\d\d\s" 
    # Wind speed regex: string of 6 characters "(0-9)(0-9)KT " bounded by 2 numbers and a space
    # Note: This definition ignores gusts
    u_pattern = r"\d[Z]\s\d\d\d\d\d\D"
    # Note: This definition allows the gust becomes the effective wind speed
    # u_pattern = r"\d\d[K][T]\s\d"
    # Sea level pressure: string of 6 characters "(0-9)(0-9)/(0-9)(0-9)" bounded by space and number
    # slp_pattern = r"\d\d[K][T]\s\d"
    
    # Iterate through DataFrame and filter out results
    for row in data.iloc[::1, 0]:
        if date_range[0] <= int(row[13:25]) <= date_range[1]:
            if re.findall(T_pattern, row) and re.findall(u_pattern, row):
                T_air_str = re.findall(T_pattern, row)[0]
                T_list.append(T_air_str[1:3])
                u_str = re.findall(u_pattern, row)[0]
                u_list.append(u_str[6:8]) # If gusts not accounted for, use this
                # u_list.append(u_str[0:2]) # If gusts accounted for, use this
                station_list.append(row[5:9])
                date_list.append(datetime.strptime(row[13:25], '%Y%m%d%H%M'))
    
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
    df['date'] = [i + timedelta(hours=5) for i in df['date']] # Adjust EST to UTC    
    df['T_air'] = T_list
    df['u_r'] = u_list
    return df

def main(start_date, end_date, loc):
    data_read(start_date, end_date, loc)
    
if __name__ == "__main__":
    main()

# Print elapsed time for script function
# print('Time elapsed for iteration: %.3f' % (time.time()-t))