### asos_data_processing branch

### Objective
# The objective of this script is to generate the sensible heat flux based on input from GOES data.

### Program flow
# 0. Define imports
# 1. Define constants
# 2. Define dynamic inputs (e.g. model data, station data)
# Note: The dynamic inputs are currently static, for simplicity, while the algorithm is finalized
# 3. Define initial conditions
# 4. Define data storage variables and other program administrative parameters
# 5. Calculate secondary parameters to support solution algorithm
# 6. Perform solution algorithm
# Note: This is currently an iterative algorithm 

# CRITERIA: Solution criteria for each loop is that each variable doesn't change by > 1% per iteration

### 0. Imports
import math
import numpy as np
import time
import pandas as pd
import matplotlib.pyplot as plt
import asos_data_reader as asos
import mesonet_reader as meso
# import sandbox as sdbx

### 1. Constants
rho = 1.225     # Air density at station (kg/m^3) - note: p = rho*R*T
c_p = 1006      # Specific heat capacity of air at station (J/kg-K)
vk = 0.4        # Von Karman constant
g = 9.81        # Gravitational acceleration (m/s^2)
R = 287.05      # Gas constant

### 2. Dynamic inputs
# z_r = 10        # Roughness height
# h_0 = 1         # Elevation
p_air = 1013.25 # Atmospheric pressure at station (hPa)
# u_r = 1.5       # Wind speed at station (m/s)
# T_lst = 295     # Land surface temperature (K)
# T_air = 290     # 2m air temperature at station (K)

def q_sens(z_r, h_0, p_air, u_r, T_lst, T_air):
    print('T_lst: %.3f | T_air: %.3f' % (T_lst, T_air))
    L = 1e15
    u_star = 1
    q_h = 1      
    psi_m = 0      
    psi_h = 0
    C_h = 0.01
    C_d = 1
    
    if u_r == 0:
        u_r = 0.001
        
    def roughness_height(h_i, s_i, N, S_T):
        p = 0
        for i in range(0, N):
            p += h_i*s_i
        return (0.25*p/S_T)
    
    z_0m = roughness_height(h_0, 21600, 200, 4e6) 
    
    C_zil = math.pow(10, -0.4*h_0)
    mu = (1.458e-6)*math.pow(T_air, 3/2)/(T_air + 110.4)
    nu = mu/rho
        
    def theta(T_abs):
        T_pot = T_abs * math.pow(1000/p_air, R/c_p)
        return T_pot
    
    val_dict = {'q_h': [q_h], 
                'C_h': [C_h],
                'C_d': [C_h],
                'z_0m': [z_0m],
                'psi_m': [0],
                'psi_h': [0],
                'u_star': [u_star],
                'L': [L]}
    
    err_val = 1
    err_dict = {'L': [err_val]}
    L_err = val_dict['L'][0]
    
    conv_crit = 0.01
    iter_lim = 10
    
    d_0 = math.exp(0.9793*np.log(h_0)-0.1536)
    z = z_r - d_0
    
    t_start = time.time()
    i = 1
    end_iter = False
    
    while abs(L_err) > conv_crit:
        # print('########## Iteration %d ##########' % i)
    
        # print("Q_h: %.3f | L: %.3f | L_err: %.3f | C_h: %.3f | u*: %.3f " % \
        #   (q_h, L, L_err, C_h, u_star))
            
        zeta = z/L
        # print('Stability parameter: %.4f' % zeta)
        
        if zeta < 1/16:
            x = math.pow(1-16*zeta, 0.25)
        else:
            x = -math.pow(abs(1-16*zeta), 0.25)
        
        if zeta < 0:
            psi_m = 2*np.log((1+x)/2) + np.log((1+math.pow(x,2))/2) - 2*math.atan(x) + math.pi/2
            psi_h = 2*np.log((1+math.pow(x,2))/2)
        else:
            psi_m = -5*zeta
            psi_h = psi_m 
            
        Re_t = z_0m*u_star/nu
        z_t = z_0m*math.exp(-vk*C_zil*math.sqrt(Re_t))    
            
        C_h_n = math.pow(vk, 2)/((np.log(z/z_0m) - psi_m*zeta) * (np.log(z/z_t) - psi_h*zeta))
        val_dict['C_h'].append(C_h_n)
        if i == 1:
            C_h = val_dict['C_h'][-1]
            
        C_d_n = math.pow(vk,2)/math.pow((np.log(z/z_0m)-psi_m*zeta), 2)
        val_dict['C_d'].append(C_d_n)
        if i == 1:
            C_d = val_dict['C_d'][-1]
        
        u_star_n = math.sqrt(C_d*u_r)
        val_dict['u_star'].append(u_star_n)
        
        q_h_n = rho*c_p*C_h*u_r*(theta(T_lst)-theta(T_air))
        val_dict['q_h'].append(q_h_n)
        
        L_n = -rho*c_p*math.pow(u_star,3)*0.5*(theta(T_lst)+theta(T_air))/(vk*g*q_h)
        val_dict['L'].append(L_n)
        
        L_err = abs(L_n-L)/L
        err_dict['L'].append(L_err)
        # print('L(i-1): %.3f; L(i): %.3f; Error: %.6f%%' % (val_dict['L'][-2], val_dict['L'][-1], (L_err*100)))
        # print("--- Iteration end ---")
        # print("Q_h: %.3f | L: %.3f | L_err: %.3f | C_h: %.3f | u*: %.3f " % \
        #   (q_h_n, L_n, L_err, C_h_n, u_star_n))
        
        if abs(L_err) < conv_crit and not end_iter:
            # print('TRIP')
            L_err = 1
            end_iter = True
            
        if i < iter_lim:
            C_h = C_h_n
            C_d = C_d_n
            u_star = u_star_n
            q_h = q_h_n
            L = L_n
            i += 1
        else:
            print('Model did not converge.')
            break
        
    t_elapsed = time.time() - t_start
    
    # print("----------------------------------------------------------------------------------\n")
    # print("Q_h: %.3f | L: %.3f | L_err: %.3f | C_h: %.3f | u*: %.3f " % \
          # (q_h_n, L_n, L_err, C_h_n, u_star_n))
            
    return C_h, C_d, L, z_0m, zeta, q_h

# def main():
    # ### Print statements
    # z_r = 33.2 
    # h_0 = 21.3
    # p_air = 1013.25 
    # u_r = 5.66 # Test case, Upper East Side, 10:45 AM EST
    # T_lst = 286 # Test case, Upper East Side, 10:45 AM EST
    # T_air = 288.7 # Test case, Upper East Side, 10:45 AM EST
    # [C_h, C_d, L, z_0m, zeta, q_h] = q_sens(z_r, h_0, p_air, u_r, T_lst, T_air)
    # print('Sensible heat flux for selected pixel is: %.4f W/m^2' % q_h)
    # print('Stability parameter for the selected pixel is: %.4f' % zeta)
    # print('Obukhov length for the selected pixel is: %.4f m' % L)
    # print('Roughness length for the selected pixel is: %.4f m' % z_0m)
    # print('Heat transfer coefficient for the selected pixel is: %.4f' % C_h)

# Main function definition
def main():
    
    ## Begin 'single data point' method for troubleshooting purposes
    # p_air = 1013.25 
    # T_lst = 300   # Test case
    # # Data point JFK20190601000010306 from 64010KJFK201906.dat
    # u_r = 2.68
    # T_air = 290.102 # Test case, Upper East Side, 10:45 AM EST
    # [C_h, C_d, L, z_0m, zeta, q_h] = q_sens(z_r, h_0, p_air, u_r, T_lst, T_air)
            
    # print('Sensible heat flux for selected pixel is: %.4f W/m^2' % q_h)
    # print('Stability parameter for the selected pixel is: %.4f' % zeta)
    # print('Obukhov length for the selected pixel is: %.4f m' % L)
    # print('Roughness length for the selected pixel is: %.4f m' % z_0m)
    # print('Heat transfer coefficient for the selected pixel is: %.4f' % C_h)
    
    ## End 'single data point' method for troubleshooting purposes
    
    ## Begin iterative method using ASOS data
    
    df = asos.data_read(201906140000, 201906142359)    # Test case, Upper East Side, 10:45 AM EST
    
    # lsts = sdbx.lst_grab()
    u_r = df['u_r']
    T_air = df['T_air']   # Test case, Upper East Side, 10:45 AM EST
        
    # T_lst = 310   # Test case
    z_r = 33.2      # Height above sea level for Brooklyn MesoNET site
    h_0 = 7.5       # Based on Index 2 (Hi-Density Residential) from WRF/URBPARAM.TBL (https://github.com/wrf-model/WRF/blob/master/run/URBPARM.TBL)
    u_list = []
    q_h_list = []
    df['lst'] = d['box_lst']
    df['T_lst'] = [287, 286, 286, 285, 285, 283, 284, 286, 290, 293, 295, 296, 297, np.nan, 297, 299, 297, 295, 293, 291, 288, 286, 286, 285] # Rudimentar
    df['T_air'] = [283, 283, 281, 281, 280, 281, 282, 283, 285, 288, 290, 292, 294, np.nan, 294, 295, 292, 290, 287, 285, 284, 284, 283, 283]
    z_list = []
    
    for index, row in df.iterrows():
        print(row)
        
    for index, row in df.iterrows():
        u_r = row['u_r']
        T_air = row['T_air']
        T_lst = row['T_lst']
        [C_h, C_d, L, z_0m, zeta, q_h] = q_sens(z_r, h_0, p_air, u_r, T_lst, T_air)
        u_list.append(u_r)
        q_h_list.append(q_h)
        z_list.append(zeta)
    
    flux_obs = meso.csv_reader(20190614, 20190615, "data")
        
    fig, ax = plt.subplots(figsize=(8, 6))
    # plt.plot(df['date'], z_list)
    plt.plot(df['date'], q_h_list)
    # plt.plot(df['date'], df['lst'])
    plt.plot(flux_obs['datetime'], flux_obs['H'])
    # ax.set_xticks(ax.get_xticks()[::24])
    plt.xticks(rotation=45)
    plt.title("Sensible heat flux, Brooklyn College")
    plt.ylabel("Q_h [W/m^2]")
    plt.legend(["Model", "Observed"])
    # print(df)

    ## End iterative method
    
if __name__ == "__main__":
    main()