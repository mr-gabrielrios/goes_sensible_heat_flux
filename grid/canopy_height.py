'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        Canopy Height
Package file path:  ~/grid/canopy_height.py
Objective:          Generate canopy height distribution over spatial domain defined by 'goes_grid_setup.py'
Author:             Gabriel Rios
Notes:              Currently dependent on NLCD and WRF land cover heights
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

# Library imports
import numpy as np, os, pandas as pd, time

##############################################################################################
# END IMPORTS
##############################################################################################

##############################################################################################
# Method name:      canopy_height_grid
# Method objective: Provide 2D array representing canopy heights corresponding to GOES coordinate grid.
# Input(s):         lcd [3D ndarray]
# Outputs(s):       h_0 [2D ndarray]
# Note:             Canopy height assumptions currently correspond to NLCD land cover data
#                   and WRF assumed heights for urban areas.
##############################################################################################

def canopy_height_grid(lcd):
    
    t = time.time()
    
    # Define file path for NLCD reference heights CSV file
    nlcd_zr_fp = os.path.join(os.path.dirname(__file__), 'nlcd_reference_heights.csv') 
    nlcd_hts = pd.read_csv(nlcd_zr_fp, usecols=['ZR']).to_numpy().ravel()
    
    # Create 2D array matching shape of spatial domain grid, populate with averaged element height data
    h_0 = np.zeros([lcd.shape[0], lcd.shape[1]])
    for i in range(0, h_0.shape[0]):
        for j in range(0, h_0.shape[1]):
            h_0[i, j] = np.dot(lcd[i, j], nlcd_hts)
    
    runtime = time.time() - t
    
    return h_0
    
if __name__ == '__main__':
    print('Add something to test here later.')