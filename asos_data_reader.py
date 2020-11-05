### Objective
# The objective of this script is to read relevant parameters from observation data logs from NOAA ASOS stations.

# Status: currently testing single cases to develop generalized algorithm.

import pandas as pd
import matplotlib.pyplot as plt
import re
import time

# Brooklyn flux data
# Select data for JFK for June 2019 as a test case
data_url = "ftp://ftp.ncdc.noaa.gov/pub/data/asos-fivemin/6401-2019/64010KJFK201906.dat"

# Import data to dataframe. Note that column 6 is wind information, column 9 is air temperature
data = pd.read_table(data_url, sep="\s+", usecols=[1, 6, 9])

# Convert Celsius to Fahrenheit
def CtoF(T):
    return (T*9/5) + 32

# Convert knots to meters per second
def KTtoMS(u):
    return u*0.51444

data.iloc[:,0] = data.iloc[:,0].str[3:15]
data.iloc[:,1] = data.iloc[:,1].str[3:5]
data.iloc[:,2] = data.iloc[:,2].str[0:2]

########## Method 2

# Import data to dataframe. Note that column 6 is wind information, column 9 is air temperature
data = pd.read_table(data_url)

# Read data for air temperature, T_lst, and wind speed, u
T_list = []
u_list = []
# Air temperature regex: string of 6 characters "(0-9)(0-9)/(0-9)(0-9)" bounded by 2 spaces
T_pattern = r"\s\d\d[^a-z0-9\s]\d\d\s" 
# Wind speed regex: string of 6 characters "(0-9)(0-9)/(0-9)(0-9)" bounded by 2 spaces
# Note: if gusts exist, the gust becomes the effective wind speed
u_pattern = r"\d\d[K][T]\s\d"

t = time.time()
for row in data.iloc[:, 0]:
    if re.findall(T_pattern, row):
        T_lst_str = re.findall(T_pattern, row)[0]
        T_list.append(T_lst_str[1:3])
    if re.findall(u_pattern, row):
        u_str = re.findall(u_pattern, row)[0]
        u_list.append(u_str[0:2])
        
T_list = list(map(int, T_list))
T_list = [CtoF(T) for T in T_list]
u_list = list(map(int, u_list))

print('Time elapsed for iteration: %.3f' % (time.time()-t))

# plt.plot(T_list)