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
import time
import pandas as pd
import matplotlib.pyplot as plt
import asos_data_reader as asos

### 1. Constants
rho = 1.225     # Air density at station (kg/m^3) - note: p = rho*R*T
c_p = 1006      # Specific heat capacity of air at station (J/kg-K)
vk = 0.4        # Von Karman constant
g = 9.81        # Gravitational acceleration (m/s^2)
R = 287.05      # Gas constant

### 2. Dynamic inputs
# z_r = 10        # Roughness height
# h_0 = 1         # Elevation
# p_air = 1013.25 # Atmospheric pressure at station (hPa)
# u_r = 1.5       # Wind speed at station (m/s)
# T_lst = 295     # Land surface temperature (K)
# T_air = 290     # 2m air temperature at station (K)

def q_sens(z_r, h_0, p_air, u_r, T_lst, T_air):

    ### 3. Initial conditions
    L = 1e15
    u_star = 1
    q_h = 1      
    psi_m = 0      
    psi_h = 0
    C_h = 1
    C_d = 1
    
    # To prevent division-by-zero error for Obukhov length (L) calculation
    if u_r == 0:
        u_r = 0.001
    
    ### 4. Calculate secondary parameters
    
    ## Aerodynamic roughness height as a function of building geometry 
    # Ref. Kondo and Yamazawa (1986) per Stull, pp. 379
    # h_i = building height
    # s_i = building ground surface area
    # N = number of buildings
    # S_T = total surface area of buildings counted
    def roughness_height(h_i, s_i, N, S_T):
        p = 0
        for i in range(0, N):
            p += h_i*s_i
        return (0.25*p/S_T)
    
    # Roughness height for Brooklyn, Flatbush area 
    # with 2km square pixel centered at Brooklyn College
    # Assume that:
    #   average building height is 21.5m (~5 stories)
    #   average block size is 21600 m^2 (80 x 270 m)
    #   number of blocks is 200 (based on street grid within selected pixel)
    #   land area is 4 km^2 (based on GOES-16 2km resolution)
    #   all buildings are of equal height          
    z_0m = roughness_height(h_0, 21600, 200, 4e6) 
        
    ## Potential temperature calculation
    def theta(T_abs):
        T_pot = T_abs * math.pow(1000/p_air, R/c_p)
        return T_pot
    
    ### 5. Program management
    ## Create dictionary to store value history in lists for each variable
    ## Note: values without explicit initial conditions are set at 0
    val_dict = {'q_h': [q_h], 
                'C_h': [C_h],
                'C_d': [C_h],
                'z_0m': [z_0m],
                'psi_m': [0],
                'psi_h': [0],
                'u_star': [u_star],
                'L': [L]}
    
    ## Create dictionary to store error history in lists for each variable
    # Assume initial values for the errors at 100%
    err_val = 1
    err_dict = {'L': [err_val]}
    L_err = val_dict['L'][0]
    
    # Set convergence criteria at 1%
    conv_crit = 0.01
    # Set limit on number of iterations for each loop
    iter_lim = 5
    
    ## Zero-plane displacement height, ref. Kim et al. (2019), Eqn. 3
    d_0 = math.exp(0.9793*math.log(h_0)-0.1536)
    z = z_r - d_0
    
    ### 6. Iterative solution for parameters L
    ## Start timing solution algorithm
    t_start = time.time()
    i = 1 # Iterand for data access and storage purposes
    
    # Iteratively calculate all parameters until L converges
    while L_err > conv_crit:
        print('########## Iteration %d ##########' % i)
        
        ## Stability parameter, zeta
        zeta = z/L
        print('Stability parameter: %.4f' % zeta)
        
        ## Similarity factor
        x = math.pow(1-16*zeta, 0.25)
        
        ## Similarity relation calculations based on atmospheric stability
        if zeta < 0:
            psi_m = 2*math.log((1+x)/2) + math.log((1+math.pow(x,2))/2) - 2*math.atan(x) + math.pi/2
            psi_h = 2*math.log((1+math.pow(x,2))/2)
        else:
            psi_m = -5*zeta
            psi_h = psi_m 
            
        ## Thermal roughness length calculation
        # Zilitinkevich relationship constant, per Chen and Zhang (2009)
        C_zil = math.pow(10, -0.4*h_0)
        # Dynamic viscosity calculation based on the Sutherland Equation
        mu = (1.458e-6)*math.pow(T_air, 3/2)/(T_air + 110.4)
        # Kinematic viscosity of air
        nu = mu/rho
        # Roughness Reynolds number
        Re_t = z_0m*u_star/nu
        # Thermal roughness length
        z_t = z_0m*math.exp(-vk*C_zil*math.sqrt(Re_t))    
            
        ## Heat transfer coefficient for heat flux, or Stanton number
        val_dict['C_h'].append(math.pow(vk,2)/ \
            ((math.log(z/z_0m)-psi_m*zeta)*(math.log(z/z_t)-psi_h*zeta)))
        # On the first iteration, make the initial condition C_h equal to its first calculation
        if i == 1:
            C_h = val_dict['C_h'][-1]
            
        ## Drag coefficient for momentum flux
        val_dict['C_d'].append(math.pow(vk,2)/ math.pow((math.log(z/z_0m)-psi_m*zeta), 2))
        # On the first iteration, make the initial condition C_d equal to its first calculation
        if i == 1:
            C_d = val_dict['C_d'][-1]
        
        ## Friction velocity
        # Ref. Launiainen et al. (1990), Eqn. A7
        val_dict['u_star'].append(math.sqrt(C_d*u_r))
        # print("u*: %.3f" % u_star)
        
        print("rho: %.3f | c_p: %.3f | C_h: %.3f | u_r: %.3f | (theta(T_lst)-theta(T_air)): %.3f" % \
              (rho, c_p, C_h, u_r, theta(T_lst)-theta(T_air)))
        ## Sensible heat flux
        # Ref. Launiainen et al. (1990), Eqn. 5
        val_dict['q_h'].append(rho*c_p*C_h*u_r*(theta(T_lst)-theta(T_air)))
        
        # print(q_h)
        ## Obukhov length (ignoring latent heat flux)
        # Ref. Launiainen et al. (A10), derivation performed by Rios
        val_dict['L'].append(-rho*c_p*math.pow(u_star,3)*0.5*(theta(T_lst)+theta(T_air))/(vk*g*q_h))
        
        # Calculate error between current and previous iteration
        L_err = abs((val_dict['L'][-1]-val_dict['L'][-2])/val_dict['L'][-2])
        err_dict['L'].append(L_err)
        # print('L(i-1): %.3f; L(i): %.3f; Error: %.6f%%' % (val_dict['L'][-2], val_dict['L'][-1], (L_err*100)))
        
        # Control logic to ensure number of iterations is not exceeded
        # Recalculate each parameter
        if i < iter_lim:
            C_h = val_dict['C_h'][-1]
            C_d = val_dict['C_d'][-1]
            u_star = val_dict['u_star'][-1]
            q_h = val_dict['q_h'][-1]
            L = val_dict['L'][-1]
            i += 1
        else:
            print('Model did not converge.')
            break
    t_elapsed = time.time() - t_start
    # print('Time elapsed: %.2f' % t_elapsed)
            
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
    z_r = 33.2      # Height above sea level for Brooklyn MesoNET site
    h_0 = 21.5      # Based on assumed average height of 5 story buildings in Flatbush
    p_air = 1013.25 
    T_lst = 300   # Test case
    # Data point JFK20190601000010306 from 64010KJFK201906.dat
    u_r = 1.78
    T_air = 288.15 # Test case, Upper East Side, 10:45 AM EST
    [C_h, C_d, L, z_0m, zeta, q_h] = q_sens(z_r, h_0, p_air, u_r, T_lst, T_air)
    
    ## End 'single data point' method for troubleshooting purposes
    
    ## Begin iterative method using ASOS data
    
    # df = asos.data_read(201906010000, 201906062359)    # Test case, Upper East Side, 10:45 AM EST
    # u_r = df['u_r']
    # T_air = df['T_air']   # Test case, Upper East Side, 10:45 AM EST
    
    # u_list = []
    # q_h_list = []
    # for index, row in df.iterrows():
    #     u_r = row['u_r']
    #     T_air = row['T_air']
    #     [C_h, C_d, L, z_0m, zeta, q_h] = q_sens(z_r, h_0, p_air, u_r, T_lst, T_air)
    #     u_list.append(u_r)
    #     q_h_list.append(q_h)
        
    # fig, ax = plt.subplots(figsize=(8, 6))
    # plt.plot(df['date'], q_h_list)
    # plt.xticks(rotation=45)
    # print(df)

    ## End iterative method
            
    print('Sensible heat flux for selected pixel is: %.4f W/m^2' % q_h)
    print('Stability parameter for the selected pixel is: %.4f' % zeta)
    print('Obukhov length for the selected pixel is: %.4f m' % L)
    print('Roughness length for the selected pixel is: %.4f m' % z_0m)
    print('Heat transfer coefficient for the selected pixel is: %.4f' % C_h)
    
if __name__ == "__main__":
    main() 