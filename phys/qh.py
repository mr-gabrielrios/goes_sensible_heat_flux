
'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        Sensible Heat Flux Calculation
Package file path:  ~/phys/qh.py
Objective:          Calculate sensible heat flux at a given coordinate and time
Author:             Gabriel Rios
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

import math, numpy as np, sys, time

##############################################################################################
# END IMPORTS
##############################################################################################

# Assumed constants
c_p = 1006      # Specific heat capacity of air at station (J/kg-K)
vk = 0.4        # Von Karman constant
g = 9.81        # Gravitational acceleration (m/s^2)
R = 287.05      # Gas constant
u_0 = 0         # Velocity at ground level

##############################################################################################
# Method name:      theta
# Method objective: Calculate potential and virtual potential temperatures.
# Input(s):         T_air [float], T_dew [float], p_air [float]
# Outputs(s):       theta [float], theta_v [float]
##############################################################################################    

def theta(T_air, T_dew, p_air):
    theta_0 = T_air * math.pow(1000/p_air, R/c_p)
    T_dew = T_dew - 273.15 # Convert from K to C
    theta_v = T_air/(1-0.379*6.11*math.pow(10, 7.5*T_dew/(237.7+T_dew))/p_air)
    # print(theta_0, theta_v)
    return theta_0

##############################################################################################
# Method name:      hfx
# Method objective: Calculate sensible heat flux using an iterative algorithm.
# Input(s):         z_r [float], h_0 [float], p_air [float], u_r [float], T_lst [float], 
#                   T_air [float], T_dew [float]
# Outputs(s):       q_h [float]
##############################################################################################

def hfx(z_r, h_0, p_air, u_r, T_lst, T_air, T_dew):
    # Define initial conditions
    L = 1e8     # Obukhov length (L --> inf as atmosphere stabilizes)
    u_star = 1  # Friction velocity
    q_h = 1     # Sensible heat flux
    psi_m = 1   # Stability parameter, momentum  
    psi_h = 1   # Stability parameter, thermal
    C_h = 1     # Heat transfer coefficient
    
    # Used to prevent a 'division by zero' error
    if h_0 < 0.001:
        h_0 = 0.01
    if u_r == 0:
        u_r = 0.001

    # Calculate air density based on air pressure and air temperature
    rho = p_air*100/(R*T_air)
        
    # Calculate zero-plane displacement height
    # Literature reference: Stanhill, 1969
    z_d = math.exp(0.9793*np.log(h_0)-0.1536)
    z = z_r - z_d

    # Data storage for parameters over iterations
    val_dict = {'q_h': [q_h], 
                'C_h': [C_h],
                'C_d': [C_h],
                'z_0m': [1],
                'psi_m': [0],
                'psi_h': [0],
                'u_star': [u_star],
                'L': [L]}
    err_val = 1
    err_dict = {'q_h': [err_val]}
    q_h_err = val_dict['q_h'][0]
    
    # Iteration control
    conv_crit = 0.01 # Convergence criteria set to 0.1%
    iter_lim = 25
    
    # Sensible heat flux iterative calculation
    i = 1 # Iterand for loop control   
    while abs(q_h_err) > conv_crit:
            
        # Momentum roughness height calculation, in conjunction 
        # with zero-plane displacement calculation outside loop
        # Literature rference: Raupach (1994) as described by Voogt and Grimmond (1999)
        z_0m = h_0*(1-z_d/h_0)*np.exp(-vk*u_r/u_star + 0.193)
        
        # Atmospheric stability parameter
        zeta = z/L    
        if zeta < 0:
            x = math.pow(1-16*zeta, 0.25)
            psi_m = 2*np.log((1+x)/2) + np.log((1+math.pow(x,2))/2) - 2*math.atan(x) + math.pi/2
            psi_h = 2*np.log((1+math.pow(x,2))/2)
        else:
            psi_m = -5*zeta
            psi_h = psi_m 
            
        # Thermal roughness height calculation algorithm
        # Literature reference: Li and Bou-Zeid (2014)
        # Zilitinkevich relationship
        C_zil = math.pow(10, -0.4*h_0)
        # Calculate dynamic viscosity
        mu = (1.458e-6)*math.pow(T_air, 3/2)/(T_air + 110.4)
        # Calculate kinematic viscosity
        nu = mu/rho
        # Calculate Reynolds parameter
        Re_t = z_0m*u_star/nu 
        # Calculate thermal roughness height
        z_t = z_0m*math.exp(-vk*C_zil*math.sqrt(Re_t))          
            
        # Heat transfer coefficient calculation
        C_h_n = math.pow(vk, 2)/((np.log(z/z_0m) - psi_m*zeta) * (np.log(z/z_t) - psi_h*zeta))
        val_dict['C_h'].append(C_h_n)
        if i == 1:
            C_h = val_dict['C_h'][-1]   
        
        # Friction velocity (u*)
        # Literature reference: Kim and Kwong, 2019, Equation 7.
        u_star_n = vk*(u_r-u_0)/(np.log(z/z_0m)-psi_m*(zeta))
        val_dict['u_star'].append(u_star_n)
        
        # Sensible heat flux calculation
        q_h_n = rho*c_p*C_h*u_r*(theta(T_lst, T_dew, p_air)-theta(T_air, T_dew, p_air))
        val_dict['q_h'].append(q_h_n)
        # Sensible heat flux error calculation
        q_h_err = abs(q_h_n-q_h)/q_h
        err_dict['q_h'].append(q_h_err)
        
        # Obukhov length calculation
        L_n = -rho*c_p*math.pow(u_star,3)*0.5*(theta(T_lst, T_dew, p_air)+theta(T_air, T_dew, p_air))/(vk*g*q_h)
        val_dict['L'].append(L_n)  
        
        # Iteration control - stop if iteration number exceeds specified limit
        if i < iter_lim:
            C_h = C_h_n
            u_star = u_star_n
            q_h = q_h_n
            L = L_n
            i += 1
        else:
            print('Model did not converge.')
            break
    
    zeta = z/L
    C_h = math.pow(vk, 2)/((np.log(z/z_0m) - psi_m*zeta) * (np.log(z/z_t) - psi_h*zeta))
    u_star = vk*(u_r-u_0)/(np.log(z/z_0m)-psi_m*(zeta))
    q_h = rho*c_p*C_h*u_r*(theta(T_lst, T_dew, p_air)-theta(T_air, T_dew, p_air))
            
    return q_h

if __name__ == '__main__':
    print('Troubleshooting goes here')