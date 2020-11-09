# Use this only in a local application where the sensible heat flux and data processing products 
# are in two different directories
import sys
sys.path.insert(0, '../goes_data_processing')

import goes_data_processing.data_processing as dp
import os
import numpy as np
import pandas as pd
from datetime import datetime as dt
from datetime import timedelta

# Where the .nc files are kept once downloaded using goes_data_processing.access
ncdir = "./goes_data_processing/ncdata/"
# Longitude and latitude for Brooklyn College
lon = -73.9521
lat = 40.6305
bound_sz = .1

def file_sort(fn):
    return(fn[-17:-4])

# def lst_grab():
#     file_list = sorted(os.listdir(ncdir), key=file_sort)
#     file_list = [os.path.join(ncdir, fn) for fn in file_list]
    
#     data_list = dp.goes_data_processor(file_list, lon, lat, bound_sz)
    
#     d = []
#     for entry in data_list:
#         d.append({
#             'datetime': dt.strptime(entry[7][0:-2], '%Y-%m-%d, %H:%M:%S').strftime('%Y-%m-%d %H:%M:%S'),
#             'box_lst': np.average(entry[0])
#         })
#     d = pd.DataFrame(d)
#     print(d)
#     return d

file_list = sorted(os.listdir(ncdir), key=file_sort)
file_list = [os.path.join(ncdir, fn) for fn in file_list]

data_list = dp.goes_data_processor(file_list, lon, lat, bound_sz)

d = []
for entry in data_list:
    d.append({
        'datetime': dt.strptime(entry[7][0:-2], '%Y-%m-%d, %H:%M:%S').strftime('%Y-%m-%d %H:%M:%S'),
        'box_lst': np.average(entry[0])
    })
d = pd.DataFrame(d)
d['datetime'] = pd.to_datetime(d['datetime'])
d['datetime'] = [i - timedelta(hours=5) + timedelta(seconds=39) for i in d['datetime']]
print(d)