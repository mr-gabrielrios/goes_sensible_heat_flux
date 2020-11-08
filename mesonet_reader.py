### Objective:/
# Read data from NYS Mesonet sites

import pandas as pd
import os
from datetime import datetime

# Inputs for troubleshooting
# start_date = 20190601
# end_date = 20190602
# data_dir = "data"

def file_sort(fn):
    return(fn[0:8])

def csv_reader(start_date, end_date, data_dir):
    data_cols = [2, 23] # Datetime and sensible heat flux
    
    # Initialize empty array to hold file data
    data_list = []
    
    # Sort files by date, assuming standard Mesonet filename convention
    # Example: YYYYMMDD_PARAMETER_LOCAITION_Parameter_NYSMesonet.csv
    file_list = sorted(os.listdir(data_dir), key=file_sort)
    # Iterate through sorted file list and extract data within date range
    for file in file_list:        
        if start_date <= int(os.path.splitext(file)[0][0:8]) <= end_date:
            print(os.path.splitext(file)[0][0:8])
            filename = data_dir + "/" + file
            data_list.append(pd.read_csv(filename, usecols=data_cols))
            
    # Concatenate all dataframes in data_list
    data = pd.concat(data_list)
    
    i = 0 # Initialize iterand for datetime formatting and row assignment
    # Iterate through each date and reformat each entry
    for date in data.iloc[:, 0]:
        data.iloc[i, 0] = datetime.strptime(date, '%Y-%m-%d %H:%M:%S.%f').strftime('%Y-%m-%d %H:%M:%S')
        i += 1
    
    # Re-cast datetime strings as datetime objects
    data['datetime'] = pd.to_datetime(data['datetime'])

    return data
        
def main(start_date, end_date, data_dir):
    csv_reader(start_date, end_date, data_dir)
    
if __name__ == "main":
    main()