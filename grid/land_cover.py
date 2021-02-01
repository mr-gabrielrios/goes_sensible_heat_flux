'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        Land Cover
Package file path:  ~/grid/land_cover.py
Objective:          Generate land cover data distribution over spatial domain defined by 'goes_grid_setup.py'
Author:             Gabriel Rios
Notes:              Currently only accommodates 2016 NLCD data upscaled to 2-km resolution
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

# Library imports
import numpy as np, os, pandas as pd, re, time

##############################################################################################
# END IMPORTS
##############################################################################################

##############################################################################################
# Method name:      land_cover_grid
# Method objective: Provide 2D array representing land cover data corresponding to GOES coordinate grid.
# Input(s):         goes_idx [1D ndarray]
# Outputs(s):       lcd [2D ndarray], water_pixel [2D ndarray]
# Note:             'goes_idx' provided in the following format, with respect to GOES CONUS data
#                   [max latitude index, min latitude index, min longitude index, max longitude index] (per PUG-L2)
#                   Currently only works with 2016 NLCD Land Cover data, upscaled from 30-m to 2-km to match GOES 
##############################################################################################

def land_cover_grid(goes_idx):
    t = time.time()
    
    # Define file path for NLCD distribution CSV file
    nlcd_data_fp = os.path.join(os.path.dirname(__file__), 'nlcd_to_goes_upscaled.csv') 
    
    # Create arrays to help with downselection from complete NLCD dataset based on chosen spatial domain
    nlcd_rows = [int(i) for i in np.linspace(goes_idx[0], goes_idx[1]-1, (goes_idx[1]-goes_idx[0]))]
    nlcd_cols = [int(i) for i in np.linspace(goes_idx[2], goes_idx[3]-1, (goes_idx[3]-goes_idx[2]))]
    
    # Filter raw NLCD dataset using indices corresponding to spatial domain
    nlcd_goes = pd.read_csv(nlcd_data_fp, usecols=nlcd_cols, low_memory=False)
    nlcd_goes = np.array(nlcd_goes[nlcd_goes.index.isin(nlcd_rows)])
    
    # Create 3D array with dimensions latitude (rows), longitude (columns), and land cover class (depth)
    nlcd_upscaled = np.array(np.zeros((np.shape(nlcd_goes)[0],np.shape(nlcd_goes)[1], 20))) # '20' chosen for number of NLCD classes
    for i in range(0,np.shape(nlcd_goes)[0]):
        for j in range(0,np.shape(nlcd_goes)[1]):
            # If data is not available, add nans
            if nlcd_goes[i, j] == 'nan':
                nlcd_ii_jj = np.repeat(np.nan, 20) # '20' chosen for number of NLCD classes
            else:
                nlcd_ii_jj = [float(qq) for qq in re.sub(' +',',',\
                                  str(nlcd_goes[i, j]).replace('[',\
                              '').replace(']','').replace('\n','')).split(',') if qq!='']
            nlcd_upscaled[i, j] = nlcd_ii_jj # Match NLCD data to GOES-16 grid shape
    
    # Record which pixels have total water cover using 2D array with shape matching NLCD grid
    water_pixel = np.zeros(nlcd_goes.shape, dtype=bool)
    for i in range(0, np.shape(nlcd_goes)[0]):
        for j in range(0, np.shape(nlcd_goes)[1]):
            if nlcd_upscaled[i, j][0] > 0.5 or np.isnan(nlcd_upscaled[i, j][0]):
                water_pixel[i, j] = True
    
    # Generic variable define to allow for plug-and-play land cover data in case of eventual downscaling
    lcd = nlcd_upscaled
    
    runtime = time.time() - t
    
    return lcd, water_pixel

if __name__ == '__main__':
    goes_idx = [300, 310, 1809, 1850]
    lcd, ii, jj = land_cover_grid(goes_idx)