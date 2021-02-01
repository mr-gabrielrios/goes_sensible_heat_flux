'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        Observation Functions
Package file path:  ~/vdtn/obs/functions.py
Objective:          Perform various functions for linking model and observational data.
Author:             Gabriel Rios
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

import numpy as np, os, pandas as pd

##############################################################################################
# END IMPORTS
##############################################################################################
  
##############################################################################################
# Method name:      distance
# Method objective: Calculate great circle distance between a point and an observation station.
# Input(s):         lat_crd [float], lon_crd [float], station_lat [float], station_lon [float]
# Output(s):        distance [float]
##############################################################################################
    
def distance(crd_lat, crd_lon, station_lat, station_lon):
    # GRS80 semi-major axis of Earth, per GOES-16 PUG-L2, Volume 5, Table 4.2.8
    R = 6378137 
    p = np.pi/180
    a = 0.5 - np.cos((station_lat-crd_lat)*p)/2 + np.cos(crd_lat*p) * np.cos(station_lat*p) * (1-np.cos((station_lon-crd_lon)*p))/2
    return 2*R*np.arcsin(np.sqrt(a))

##############################################################################################
# Method name:      station_finder
# Method objective: Retrieve path to observation station directory closest to queried point.
# Input(s):         crd [list]
# Output(s):        path [str]
##############################################################################################
  
def station_finder(crd):
    
    # Retrieve latitude and longitude
    lat, lon = crd
    # Initialize empty DataFrame
    stations = pd.DataFrame()
    # Get list of all observation networks, assuming each network has a dedicated directory
    networks = [f.path for f in os.scandir(os.path.dirname(__file__)) 
                if (f.is_dir() and f.path.split('/')[-1] != '__pycache__')]
    # Iterate through list of network paths and grab metadata for network stations from log file
    for network in networks:
        try:
            stations = stations.append(pd.read_csv(os.path.join(network, 'log.csv')))
        except:
            print('No station log here, proceeding to next observation network')
    stations = stations.reset_index(drop=True)
    
    # Search for closest station in stations list
    min_dist = 1e9
    for i in range(0, stations.shape[0]):
        dist = distance(lat, lon, stations.at[i, 'Latitude'], stations.at[i, 'Longitude'])
        if dist < min_dist:
            station = stations.at[i, 'Site']
            min_dist = dist
    
    # Grab file path for the directory belonging to the closest station
    for root, dirs, files in os.walk(os.path.dirname(__file__)):
        for d in dirs:
            if station == d:
                path = os.path.join(root, d)
            
    return path

if __name__ == '__main__':
    fpath = station_finder(40.6305, -73.9521)