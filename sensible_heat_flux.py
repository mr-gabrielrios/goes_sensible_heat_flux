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

### 1. Constants
rho = 1.225     # Air density at station (kg/m^3) - note: p = rho*R*T
c_p = 1006      # Specific heat capacity of air at station (J/kg-K)
vk = 0.4        # Von Karman constant
g = 9.81        # Gravitational acceleration (m/s^2)
R = 287.05      # Gas constant
u_0 = 0         # Velocity at ground level

### 2. Solution algorithm
def q_sens(z_r, h_0, p_air, u_r, T_lst, T_air):
    
    ### 3. Define initial conditions
    L = 1e8
    u_star = 1
    q_h = 1      
    psi_m = 1     
    psi_h = 1
    C_h = 1
    
    # Calculate air density based on air pressure and air temperature
    rho = p_air*100/(R*T_air)
    
    # Used to prevent a 'division by zero' error
    if u_r == 0:
        u_r = 0.001
     
    # Kondo and Yamazawa (1986) method to determine z_0m
    # def roughness_height(h_i, s_i, N, S_T):
    #     p = 0
    #     for i in range(0, N):
    #         p += h_i*s_i
    #     return (0.25*p/S_T)
    # z_0m = roughness_height(h_0, 21600, 200, 4e6)
    # print('NLCD-sourced roughness height: %.3f | Kondo/Yamazawa derived roughness height: %.3f' % (h_0, z_0m))
    
    # Raupach (1994) method to determine z_d, z_0m
    cd = 7.5 # Per Raupach (1994)
    # IMPORTANT: lambda_f must be picked with rationale
    lambda_f = 1 # Per Raupach (1994). Can range from 1-10.
    z_d = h_0*(1-(1-np.exp(-np.sqrt(cd*lambda_f)))/np.sqrt(cd*lambda_f))
      
    # z_d = math.exp(0.9793*np.log(h_0)-0.1536)
    z = z_r - z_d
    print('z_r: %.3f | z_d: %.3f | h_0: %.3f' % (z_r, z_d, h_0))
    
    def theta(T_abs):
        T_pot = T_abs * math.pow(1000/p_air, R/c_p)
        return T_pot
    
    ### Administrative parameters
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
    
    conv_crit = 0.01 # Convergence criteria set to 0.1%
    iter_lim = 25 # Number chosen based on algorithm by Launiainen (1990)
    
    ### Sensible heat flux iterative calculation
    t_start = time.time()
    i = 1    
    while abs(q_h_err) > conv_crit:
        print('########## Iteration %d ##########' % i)
        print("PRE-ITERATION Q_h: %.3f | q_h_err: %.3f | L: %.3f | C_h: %.3f | u*: %.3f " % \
          (q_h, q_h_err, L, C_h, u_star))
        
        # Momentum roughness height calculation, in conjunction with zero-plane displacement calculation outside loop
        # Per Raupach (1994) as described by Voogt and Grimmond (1999)
        z_0m = h_0*(1-z_d/h_0)*np.exp(-vk*u_r/u_star + 0.193) # Values per Raupach (1994)
        print('z_0m: %.3f m' % z_0m)
        # Atmospheric stability parameter
        zeta = z/L
        print('Stability parameter: %.4f' % zeta)
        
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
            
        # Zilitinkevich relationship, per Li and Zeid (2014)
        C_zil = math.pow(10, -0.4*h_0)
        mu = (1.458e-6)*math.pow(T_air, 3/2)/(T_air + 110.4)
        nu = mu/rho
        Re_t = z_0m*u_star/nu
        z_t = z_0m*math.exp(-vk*C_zil*math.sqrt(Re_t))    
            
        C_h_n = math.pow(vk, 2)/((np.log(z/z_0m) - psi_m*zeta) * (np.log(z/z_t) - psi_h*zeta))
        val_dict['C_h'].append(C_h_n)
        if i == 1:
            C_h = val_dict['C_h'][-1]
        
        # Friction velocity (u*). Reference Kim et al, 2019, Equation 7.
        u_star_n = vk*(u_r-u_0)/(np.log(z/z_0m)-psi_m*(zeta))
        val_dict['u_star'].append(u_star_n)
        
        q_h_n = rho*c_p*C_h*u_r*(theta(T_lst)-theta(T_air))
        val_dict['q_h'].append(q_h_n)
        
        L_n = -rho*c_p*math.pow(u_star,3)*0.5*(theta(T_lst)+theta(T_air))/(vk*g*q_h)
        val_dict['L'].append(L_n)
        
        print("POST-ITERATION Q_h_n: %.3f | L_n: %.3f | C_h_n: %.3f | u*_n: %.3f " % \
          (q_h_n, L_n, C_h_n, u_star_n))
        
        q_h_err = abs(q_h_n-q_h)/q_h
        err_dict['q_h'].append(q_h_err)
            
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
    q_h = rho*c_p*C_h*u_r*(theta(T_lst)-theta(T_air))
    
    print('NLCD-sourced roughness height: %.3f | Raupach derived roughness height: %.3f' % (h_0, z_0m))
    print("Q_h: %.3f | q_h_err: %.3f | L: %.3f | C_h: %.3f | u*: %.3f | psi_h: %.3f " % (q_h, q_h_err, L, C_h, u_star, psi_h))
            
    return C_h, u_star, L, z_0m, zeta, q_h