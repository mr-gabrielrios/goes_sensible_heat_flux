'''
Script name:        NYS Mesonet File Conversion - mesonet_converter.py
Script objective:   The objective of this script is to generate files to match the format 
                    given by NYS Mesonet for flux data 
'''

### Imports
import os, sys
from datetime import datetime
from datetime import timedelta
import pandas as pd

##########################################################
# File Conversion
# Objective:    Convert files into desired format and generate files accordingly
# Input:        directory name for files of interest (str)
# Output:       none

def file_conversion(fpath):
    
    ## Date range setup
    # Assume that the data date range is in the directory name, with '_' being the separator
    start_date, end_date = [fpath.split('_')[-2][1:],
                            fpath.split('_')[-1][1:]]
    start_date, end_date = [datetime.strptime(start_date, '%Y%m%d'),
                            datetime.strptime(end_date, '%Y%m%d')+timedelta(days=1)]
    
    # Get every day between bounds of the time range
    date_range = pd.date_range(start_date, end_date, (end_date-start_date).total_seconds()/3600/24+1)
    date_range = [t.to_pydatetime() for t in date_range]
    
    ## Station list in NYS Mesonet network, New Yor City metro area only
    station_list = ['BKLN', 'BRON', 'MANH', 'QUEE', 'STAT']
    
    ## Column list for data reading
    # Station ID, datetime, sensible heat flux, friction velocity, air temperature
    data_cols = ['stid', 'datetime', 'H', 'USTAR', 'Tc'] 
    
    ## Filter 1: Create DataFrame with all data filtered by selected columns
    all_data = pd.DataFrame()
    for file in os.listdir(fpath):
        df = pd.read_csv(os.path.join(fpath, file), usecols=data_cols)
        all_data = pd.concat([all_data, df])
    
    # Note:     I know I should've used a dictionary to store different station DataFrames, but this is just \
    #           an auxiliary script I won't use often, if ever again.
    
    ## Filter 2: Create DataFrames with all data filtered by station
    for station in station_list:
        if station == 'BKLN':
            bkln_df = all_data[all_data['stid'] == 'FLUX_BKLN']
            bkln_df['datetime'] = bkln_df['datetime'].astype('str')
        elif station == 'BRON':
            bron_df = all_data[all_data['stid'] == 'FLUX_BRON']
            bron_df['datetime'] = bron_df['datetime'].astype('str')
        elif station == 'MANH':
            manh_df = all_data[all_data['stid'] == 'FLUX_MANH']
            manh_df['datetime'] = manh_df['datetime'].astype('str')
        elif station == 'QUEE':
            quee_df = all_data[all_data['stid'] == 'FLUX_QUEE']
            quee_df['datetime'] = quee_df['datetime'].astype('str')
        elif station == 'STAT':
            stat_df = all_data[all_data['stid'] == 'FLUX_STAT']
            stat_df['datetime'] = stat_df['datetime'].astype('str')
        else:
            print('Station not found! Ending script...')
            sys.exit()
    
    ## Filter 3: For each station's DataFrame, save daily CSVs
    # BKLN
    loc = 'BKLN'
    data_fpath = 'mesonet_data/' + loc
    for day in date_range:
        fname = day.strftime('%Y%m%d') + '_' + 'FLUX_' + loc + '_Flux_NYSMesonet.csv'
        day = day.strftime('%Y-%m-%d')
        date = pd.to_datetime(bkln_df['datetime']).dt.date
        if not bkln_df[date == day].empty:
            temp = bkln_df[date == day]
            # Avoid re-writing copies of the same file
            if os.path.isfile(os.path.join(data_fpath, fname)):
                temp.to_csv(path_or_buf=os.path.join(data_fpath, fname))
    # BRON
    loc = 'BRON'
    data_fpath = 'mesonet_data/' + loc
    for day in date_range:
        fname = day.strftime('%Y%m%d') + '_' + 'FLUX_' + loc + '_Flux_NYSMesonet.csv'
        day = day.strftime('%Y-%m-%d')
        date = pd.to_datetime(bron_df['datetime']).dt.date
        if not bron_df[date == day].empty:
            temp = bron_df[date == day]
            # Avoid re-writing copies of the same file
            if os.path.isfile(os.path.join(data_fpath, fname)):
                temp.to_csv(path_or_buf=os.path.join(data_fpath, fname))
    # MANH
    loc = 'MANH'
    data_fpath = 'mesonet_data/' + loc
    for day in date_range:
        fname = day.strftime('%Y%m%d') + '_' + 'FLUX_' + loc + '_Flux_NYSMesonet.csv'
        day = day.strftime('%Y-%m-%d')
        date = pd.to_datetime(manh_df['datetime']).dt.date
        if not manh_df[date == day].empty:
            temp = manh_df[date == day]
            # Avoid re-writing copies of the same file
            if os.path.isfile(os.path.join(data_fpath, fname)):
                temp.to_csv(path_or_buf=os.path.join(data_fpath, fname))
    # QUEE
    loc = 'QUEE'
    data_fpath = 'mesonet_data/' + loc
    for day in date_range:
        fname = day.strftime('%Y%m%d') + '_' + 'FLUX_' + loc + '_Flux_NYSMesonet.csv'
        day = day.strftime('%Y-%m-%d')
        date = pd.to_datetime(quee_df['datetime']).dt.date
        if not quee_df[date == day].empty:
            temp = quee_df[date == day]
            # Avoid re-writing copies of the same file
            if os.path.isfile(os.path.join(data_fpath, fname)):
                temp.to_csv(path_or_buf=os.path.join(data_fpath, fname))
    # STAT
    loc = 'STAT'
    data_fpath = 'mesonet_data/' + loc
    for day in date_range:
        fname = day.strftime('%Y%m%d') + '_' + 'FLUX_' + loc + '_Flux_NYSMesonet.csv'
        day = day.strftime('%Y-%m-%d')
        date = pd.to_datetime(stat_df['datetime']).dt.date
        if not stat_df[date == day].empty:
            temp = stat_df[date == day]
            # Avoid re-writing copies of the same file
            if os.path.isfile(os.path.join(data_fpath, fname)):
                temp.to_csv(path_or_buf=os.path.join(data_fpath, fname))
    
    return bkln_df, bron_df, manh_df, quee_df, stat_df
                   
fpath = 'mesonet_data/flux_NYC_ALL_s20191001_e20200531'
dates = file_conversion(fpath)