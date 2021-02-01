'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        Variable Names
Package file path:  ~/plot/varnames.py
Objective:          Store metadata for variables used in the calculation of sensible heat flux.
Author:             Gabriel Rios
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

##############################################################################################
# END IMPORTS
##############################################################################################

##############################################################################################
# Method name:      var_metadata
# Method objective: Allow for variable metadata access given a variable name.
# Input(s):         var [str]
# Outputs(s):       metadata [list]
##############################################################################################

def var_metadata(var, plot_var):
    
    if var == 'station_data':
        var = plot_var
    
    T_s = {
        'name': 'T_s',
        'full_name': 'Land surface temperature',
        'units': 'K'}
    
    T_a = {
        'name': 'T_a',
        'full_name': 'Air temperature, 2 m AGL',
        'units': 'K'}
    
    Q_H = {
        'name': 'Q_H',
        'full_name': 'Upward surface sensible heat flux',
        'units': 'W m^{-2}'}
    
    u_r = {
        'name': 'u_r',
        'full_name': 'Wind speed, surface',
        'units': 'm {s^-1}'}
    
    h_0 = {
        'name': 'h_0',
        'full_name': 'Canopy element height',
        'units': 'm'}
    
    z_m = {
        'name': 'z_m',
        'full_name': 'Momentum roughness element height',
        'units': 'm'}
    
    u_star = {
        'name': 'u_star',
        'full_name': 'Friction velocity',
        'units': 'm {s^-1}'}
    
    L = {
        'name': 'L',
        'full_name': 'Obukhov length',
        'units': 'm'}
    
    var_list = [T_s, T_a, Q_H, u_r, h_0, z_m, u_star, L]
    
    for v in var_list:
        if var is v['name']:
            metadata = v
            return metadata  
            

# Use this method for troubleshooting
if __name__ == '__main__':
    print(var_metadata('a_H', 'Q_H'))