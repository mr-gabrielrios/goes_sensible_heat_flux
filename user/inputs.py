'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        Inputs
Package file path:  ~/user/inputs.py
Objective:          Handle user inputs if external package inputs aren't passed
Author:             Gabriel Rios
Notes:              N/A
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

# Library imports
import datetime, pandas, sys
# Internal imports


##############################################################################################
# END IMPORTS
##############################################################################################

##############################################################################################
# Method name:      user_input
# Method objective: Take user input for spatial domain, date range.
# Input(s):         args (bool)
# Outputs(s):       domains [list], date range [datetime]
##############################################################################################

def user_input(domain=False, date_range=False, plot_var=False):
    
    # If spatial domain is not defined, allow user input to do so.
    # Objective:    Support 'plug and play' option with other models by allowing for arguments
    #               to pass into the keywords    
    if domain is False:
        domain = [None] * 4 # Initialize list for geographical domain extent
        user_lats = input('Enter latitudinal extent in degrees N (latitude 1, latitude 2): ') # Gather user input for latitudes
        domain[0], domain[1] = sorted([float(i) for i in user_lats.split(',')])
        if domain[0] <= 24 or domain[1] >= 50: # Catch if requested latitude(s) are outside of CONUS
            print('Product only calculates sensible heat flux for the contiguous United States. Please try running again.')
            sys.exit()    
        user_lons = input('Enter longitudinal extent in degrees E (longitude 1, longitude 2): ') # Gather user input for longitudes
        domain[2], domain[3] = sorted([float(i) for i in user_lons.split(',')])
        if domain[2] <= -125 or domain[3] >= -66:  # Catch if requested longitude(s) are outside of CONUS
            print('Product only calculates sensible heat flux for the contiguous United States. Please try running again.')
            sys.exit()
            
    # If temporal domain is not defined, allow user input to do so.
    # Objective:    Support 'plug and play' option with other models by allowing for arguments
    #               to pass into the keywords    
    if date_range is False:
        date_range = [None] * 2 # Initialize list for time ranges
        freq = 'H' # Initialize default period for 1-hour temporal frequency
        date_range[0] = input('Enter start date for model run (YYYY-mm-dd-HH format): ')
        date_range[1] = input('Enter end date for model run (YYYY-mm-dd-HH format): ')
        date_range = [datetime.datetime.strptime(date, '%Y-%m-%d-%H') for date in date_range]
        date_range = pandas.date_range(start=date_range[0], end=date_range[1], freq=freq) 
        
    # If variable that would like to be plotted isn't defined, let user define it.
    if plot_var is False:
        plot_var = 'Q_H' # Default option for plotted variable
        plot_var = input('Enter name of variable to be plotted (see readme for list of variable names: ')        
        
    return domain, date_range, plot_var