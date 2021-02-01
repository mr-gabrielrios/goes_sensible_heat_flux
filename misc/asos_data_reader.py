'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        ASOS Data Reader
Package file path:  ~/vdtn/asos_data_reader.py
Objective:          Read ASOS data for selected fields from FTP for defined location and time.
Author:             Gabriel Rios
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

import datetime, numpy as np, os, pandas as pd, re, shutil, sys, time, urllib.request

##############################################################################################
# END IMPORTS
##############################################################################################

str_switch = False

##############################################################################################
# Method name:      CtoK
# Method objective: Convert degrees Celsius to degrees Kelvin.
# Input(s):         T [float]
# Outputs(s):       T + 273.15 [float]
##############################################################################################

def CtoK(T):
    return T + 273.15

##############################################################################################
# Method name:      KTtoMS
# Method objective: Convert wind speed from knots to meters per second.
# Input(s):         u [float]
# Outputs(s):       u * 0.51444 [float]
##############################################################################################

def KTtoMS(u):
    return u*0.51444

##############################################################################################
# Method name:      distance
# Method objective: Calculate great circle distance between a point and an ASOS station.
# Input(s):         lat_crd [float], lon_crd [float], lat_asos [float], lon_asos [float]
##############################################################################################
    
def distance(lat_crd, lon_crd, lat_asos, lon_asos):
    # GRS80 semi-major axis of Earth, per GOES-16 PUG-L2, Volume 5, Table 4.2.8
    R = 6378137 
    p = np.pi/180
    a = 0.5 - np.cos((lat_asos-lat_crd)*p)/2 + np.cos(lat_crd*p) * np.cos(lat_asos*p) * (1-np.cos((lon_asos-lon_crd)*p))/2
    return 2*R*np.arcsin(np.sqrt(a))

##############################################################################################
# Method name:      asos_find
# Method objective: Find closest ASOS station to a given coordinate
# Input(s):         lat_crd [float], lon_crd [float]
# Outputs(s):       station [str]
##############################################################################################

def asos_find(lat_crd, lon_crd):
    t = time.time() 
    # Text file listing all ASOS stations with metadata.
    asos_station_fp = r'https://www.ncdc.noaa.gov/homr/file/asos-stations.txt'
    # GRS80 semi-major axis of Earth, per GOES-16 PUG-L2, Volume 5, Table 4.2.8
    R = 6378137
    # Define relevant columns for station location
    asos_cols = ['CALL', 'NAME', 'LAT', 'LON', 'ELEV', 'UTC']
    # Read data into ASOS DataFrame (adf)
    adf = pd.read_fwf(asos_station_fp, usecols=asos_cols).drop([0])
    # Filter stations with null data
    adf = adf[adf['ELEV'].astype(float) > -99999]
    # Read relevant parameters into lists for iterative purposes
    stations, lat_asos, lon_asos = [adf['CALL'].tolist(), adf['LAT'].astype(float).tolist(), adf['LON'].astype(float).tolist()]
    
    # Define arbitrarily large number as an initial condition (rouglhy equal to circumference of Earth)
    dists = 2*np.pi*R 
    # Initialize empty string for population
    station = '' 
    # Iterate over list of statons to find closest station
    for i in range(1, len(lat_asos)):
        dist = distance(lat_crd, lon_crd, lat_asos[i], lon_asos[i])
        if dist < dists:
            dists = dist
            station = 'K' + stations[i]

    if str_switch:
        print("Closest station: %s, %.2f m away" % (station, dists))
        print('asos_find runtime: %.4f s' % (time.time() - t))
    
    # Manual override for work in New York City - KNYC is generally unreliable.
    if station == 'KNYC':
        station = 'KLGA'
        
    return station

##############################################################################################
# Method name:      wfile
# Method objective: Download ASOS data file for month corresponding to given date for the
#                   corresponding location.
# Input(s):         fpath [str], crd [list or tuple]
# Outputs(s):       station [str]
##############################################################################################

def wfile(fpath, crd):
    # Define path where ASOS data will be saved locally
    local_path = os.path.join(os.path.dirname(__file__), 'asos_data/' + '64010' + asos_find(crd[0], crd[1]) + fpath)
    
    print('Reading from {0}...'.format(asos_find(crd[0], crd[1])))
    # If path exists locally, use that. Else, download it and use locally moving forward.
    if os.path.isfile(local_path):
        return local_path
    else:
        # Create ASOS data folder if non-existent
        if not os.path.isdir(os.path.join(os.path.dirname(__file__), local_path.split('/')[-2])):
            # Choose -4 to omit the 3-character file extension
            os.mkdir(os.path.join(os.path.dirname(__file__), local_path.split('/')[-2]))
        url = 'ftp://ftp.ncdc.noaa.gov/pub/data/asos-fivemin/6401-' + fpath[0:4] + '/64010' + asos_find(crd[0], crd[1]) + fpath
        with urllib.request.urlopen(url) as f:
            dat = f.read().decode('utf-8')
            file = open(local_path, 'w')
            file.write(dat)
            file.close()
        return local_path
         
##############################################################################################
# Method name:      data_read
# Method objective: Pull selected data from ASOS file for a given spatial domain and point.
# Input(s):         start_date [datetime], end_date [datetime], crd [list or tuple]
# Outputs(s):       df [Pandas DataFrame]
##############################################################################################

def data_read(date_range, crd, utc_offset):
    
    # Number of ASOS observations an hour
    interval = 12
    freq = ['5min', 'H']
    
    # Adjust datetimes to timezone corresponding to location of interest
    date_range = pd.date_range(start=date_range[0], end=date_range[-1], freq=freq[0]) 
    start_date, end_date = [date_range[0], date_range[-1]]
    
    # Generate strings from dates for comparison purposes
    date_str = [datetime.datetime.strftime(start_date, '%Y%m%d%H%M'),
                datetime.datetime.strftime(end_date, '%Y%m%d%H%M')]
    
    # Initialize DataFrame here to return nan DataFrame in case of failed FTP connection
    df = pd.DataFrame(np.nan, index=date_range, columns=['T_air', 'T_dew', 'u_r', 'p_air'])
    # Set up URL to appropriate data file
    data_url = wfile(date_str[0][0:6] + '.dat', crd)
    
    # Import data to DataFrame 'df'
    try:
        asos_data = pd.read_table(data_url, header=None)
    except:
        print('FTP connection failed. Exiting program...')
        sys.exit()
        
    ## Regex patterns for data mining through the .dat file(s)
    # Air temperature regex: string of 6 characters "(0-9)(0-9)/(0-9)(0-9)" bounded by 2 spaces
    T_pattern = r'\s.?\d\d[+-/].?\d\d\s'
    # Wind speed regex: string of 6 characters "(0-9)(0-9)KT " bounded by 2 numbers and a space
    # Note: This definition ignores gusts
    u_pattern = r"\s\d\d\d\d\d\D"
    # Note: This definition allows the gust becomes the effective wind speed
    # u_pattern = r"\d\d[K][T]\s\d"
    # Air pressure regex: string of 6 characters "SLP(0-9)(0-9)"
    p_pattern = r"[S][L][P]\d\d\d"
    
    # Iterate through all rows in ASOS data file. For dates in file that are within date range, extract data.
    for row in asos_data.iloc[:, 0]:
        if datetime.datetime.strptime(row[13:23], '%Y%m%d%H') in df.index:
            # If temperature pattern is found, extract data.
            if re.findall(T_pattern, row):
                date = datetime.datetime.strptime(row[13:25], '%Y%m%d%H%M')
                
                # Extract air temperature ('M' prefix indicates a negative temperature)
                T_air_str = re.findall(T_pattern, row)[0]                
                if T_air_str[1] == 'M':
                    df.loc[date, 'T_air'] = CtoK(-int(T_air_str[2:4]))
                else:
                    df.loc[date, 'T_air'] = CtoK(int(T_air_str[1:3]))
                    
                # Extract dew point temperature ('M' prefix indicates a negative temperature)
                if T_air_str[-4] == 'M':
                    df.loc[date, 'T_dew'] = CtoK(-int(T_air_str[-3:-1]))
                else:
                    df.loc[date, 'T_dew'] = CtoK(int(T_air_str[-3:-1]))
                    
                # Extract wind speed
                if re.findall(u_pattern, row):
                    u_str = re.findall(u_pattern, row)[0]
                    df.loc[date, 'u_r'] = int(u_str[4:6])
                else:
                    df.loc[date, 'u_r'] = 0
                
                # Extract air pressure
                if re.findall(p_pattern, row):
                    # Convert p_str to pressure in hPa
                    p_temp = 1000 + int(re.findall(p_pattern, row)[0][-3:])/10 
                    df.loc[date, 'p_air'] = int(p_temp)
                else:
                    df.loc[date, 'p_air'] = 1013.25
    
    # Average over all observations to produce hourly, then re-index to set dates to proper indices.
    df = pd.DataFrame(df.values.reshape(-1, interval, df.shape[1]).mean(1), columns=df.columns)
    df['date'] = pd.date_range(start=date_range[0], end=date_range[-1], freq=freq[1])
    df = df.set_index('date')
    
    # Delete ASOS data folder created locally
    shutil.rmtree(os.path.join(os.path.dirname(__file__), data_url.split('/')[-2]))

    return df['u_r'], df['T_dew'], df['p_air']
    
if __name__ == "__main__":
    print("Add troubleshooting options here...")
