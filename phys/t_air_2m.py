'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        Air Temperature Model, 2 meters AGL
Package file path:  ~/phys/t_air_2m.py
Objective:          Calculate air temperature at 2 m AGL for a GOES-16 grid pixel.
Author:             Gabriel Rios
Notes:              Reference Hrisko et al (2020), https://doi.org/10.1016/j.rse.2019.111495
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

# External imports
import csv, glob, numpy as np, os, pandas as pd, time

##############################################################################################
# END IMPORTS
##############################################################################################

def air_temp(domain, pixel, date, goes_idx, pos, lcd, lst):
        
    t = time.time()
    
    # Define custom file name segment with bounding coordinates (NW and SE corners)
    domain_name = str(domain[1]) + 'N_' + str(domain[2]) + 'W_' + str(domain[0]) + 'N_' + str(domain[3]) + 'W'
    
    # Retrieve elevation data file
    fname = glob.glob(os.getcwd() + '/grid/elevation_csvs/' + domain_name + '_elevation.csv')[0]
    goes_elev = np.loadtxt(open(fname, 'r'), delimiter=',')
    
    # Read in Gaussian fit parameters
    col_array, fit_params, gaus_header = [], [], []
    with open(os.path.join(os.path.dirname(__file__), 'gaus_fit_params.csv'), 'r') as fit_params_csv:
        reader = csv.reader(fit_params_csv, delimiter=',')
        for row in reader:
            if gaus_header == []:
                gaus_header = row
                continue
            col_array.append(int(row[0])) # this variable identifies the variable column to multiply
            fit_params.append(row[1:]) # the rest of the fit coefficients
    
    # Collect variables for the Gaussian inputs at each pixel
    input_array = [] 
    if goes_elev[pos[0], pos[1]] > -1:
        # print('Elevation: {0}'.format(goes_elev[pos[0], pos[1]]))
        input_ii_jj = [0, pixel[0], pixel[1], goes_elev[pos[0], pos[1]]]  
    else:   
        # print('No elevation data...')
        input_ii_jj = [0, pixel[0], pixel[1], np.nan]  
    # Roll out land cover list and append to the input parameter list
    input_ii_jj.extend(lcd)
    input_array.append(input_ii_jj)
    
    # This computes the Gaussian variables for the current pixel
    A_0, x_0, sigma, y_0 = [], [], [], []
    for goes_ii in input_array:
        A_ii,x_ii,sigma_ii,y_ii = 0.0,0.0,0.0,0.0
        for aa in range(0,len(col_array)):
            if fit_params[aa][0]!='':
                col_vals = [float(qq) for qq in fit_params[aa][0].replace('[','').replace(']','').split(',')]
                A_ii+=(goes_ii[col_array[aa]]*col_vals[1])+col_vals[2]
            if fit_params[aa][1]!='':
                col_vals = [float(qq) for qq in fit_params[aa][1].replace('[','').replace(']','').split(',')]
                x_ii+=(goes_ii[col_array[aa]]*col_vals[1])+col_vals[2]
            if fit_params[aa][2]!='':
                col_vals = [float(qq) for qq in fit_params[aa][2].replace('[','').replace(']','').split(',')]
                sigma_ii+=(goes_ii[col_array[aa]]*col_vals[1])+col_vals[2]
            if fit_params[aa][3]!='':
                col_vals = [float(qq) for qq in fit_params[aa][3].replace('[','').replace(']','').split(',')]
                y_ii+=(goes_ii[col_array[aa]]*col_vals[1])+col_vals[2]
        A_0.append(A_ii)
        x_0.append(x_ii)
        sigma.append(sigma_ii)
        y_0.append(y_ii)
    
    lst = lst.ravel() # ravel the LST data for processing
    
    # THE GAUSSIAN CALCULATION FUNCTION
    def gaus_func(LST_func, A_0_func, x_0_func, sigma_func, y_0_func, local_t_func):
        return LST_func + y_0_func - (A_0_func*np.exp(-(np.power(np.subtract(local_t_func.hour,x_0_func),2.0))/(2.0*sigma_func**2.0)))
    
    T_air = []
    for val_iter in range(0,len(input_array)):
        T_air_ii = gaus_func(lst[val_iter],A_0[val_iter],x_0[val_iter],sigma[val_iter],y_0[val_iter],date)
        if T_air_ii>330.0 or T_air_ii<180.0:
            T_air.append(np.nan) # if outside reasonable range, nan
        else: 
            T_air.append(T_air_ii)
            
    return T_air[0]

if __name__ == '__main__':
    wdir = os.path.join(os.getcwd(), 'grid/') + 'elevation_csvs'
    print(os.path.dirname(__file__))