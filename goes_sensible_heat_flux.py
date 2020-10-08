
### Objective
# The objective of this script is to generate the sensible heat flux based on input from GOES data

### Program flow
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
import plot
import time

### 1. Constants
rho = 1.225     # Density of air (kg/m^3), assuming p = 1013.25 hPa and T = 15 degC
c_p = 1006      # Specific heat capacity of air (J/kg-K), assuming p = 1000 hPa and T = 0 degC
vk = 0.4        # von Karman constant
z_r = 2.5         # Reference height (m), typically 2 m, based on air temperature model 
                # per Hrisko et al. (2020)
h_0 = 2.596         # Vegetation canopy height (m), assumption
g = 9.81        # Gravitational acceleration constant (m/s^2)
u_0 = 0.1       # Wind speed at surface (m/s), assumption
R = 287.05      # Gas constant for air (J/kg-K)

### 2. Dynamic inputs
## Station inputs
p_air = 1013.25 # Atmospheric pressure (hPa), assumption - weather station input
u_r = 1         # Wind speed at reference height (m/s), assumption - weather station input
## Model inputs
T_lst = 295     # Land surface temperature (LST) (K), assumption - Hrisko model input
T_air = 290     # Air temperature at reference height (K), assumption - Hrisko model input

### 3. Initial conditions
L = 1e15           # Initial condition for Obukhov length
z_0m = 10        # Initial condition for roughness length
u_star = 1      # Initial condition for friction velocity
C_h = 1         # Initial condition for heat transfer coefficient
q_h = 1         # Initial condition for sensible heat flux

### 4. Data storage and program administration

## Create dictionary to store value history in lists for each variable
## Note: values without explicit initial conditions are set at 0
val_dict = {'q_h': [q_h], 
            'C_h': [C_h],
            'z_0m': [z_0m],
            'psi_m': [0],
            'psi_h': [0],
            'u_star': [u_star],
            'L': [L]}

## Assume initial values for the errors at 100%
err_val = 1
q_h_err, z_0m_err, C_h_err, u_star_err, L_err = [1, 1, 1, 1, 1]
## Create dictionary to store error history in lists for each variable
err_dict = {'q_h': [q_h_err], 
       'C_h': [err_val],
       'z_0m': [err_val],
       'psi_m': [err_val],
       'psi_h': [err_val],
       'u_star': [err_val],
       'L': [err_val]}

# Set convergence criteria at 1%
conv_crit = 0.01
# Set limit on number of iterations for each loop
iter_lim = 5

### 5. Calculate secondary parameters

## Potential temperature calculation
def theta(T_abs):
    T_pot = T_abs * math.pow(1000/p_air, R/c_p)
    return T_pot

## Define potential temperatures
theta_0 = theta(T_lst)
theta_r = theta(T_air)

## Zero-plane displacement height, ref. Kim et al. (2019), Eqn. 3
# h_0 = vegetation height (m)
d_0 = math.exp(0.98*math.pow(math.log(h_0), 2)-0.15)

### 6. Iterative solution for parameters z_0m, q_h, and L
## Start timing solution algorithm
start_time = time.time()

##############################################################################################

## Initialize outer while loop
i = 1 
while (L_err + z_0m_err) > conv_crit:
    print('############################ \nIteration #%d' % i)
           
    sp = (z_r - d_0)/L
    print('Stability parameter: %.3f' % sp)
    
    x = math.pow((1 - 15*sp), 0.25)
    
    if (z_r/L) >= 0:
        psi_m = -5*(z_r - d_0)/L
        psi_h = psi_m
    elif -2 < (z_r/L) < 0:
        psi_m = math.log(((1 + math.pow(x,2))/2) * math.pow((1 + x)/2, 2)) \
            - 2*math.atan(x) + math.pi/2
        psi_h = 2*math.log((1 + math.pow(x, 2))/2)            
    
    # L_n = -rho * c_p * math.pow(u_star, 3) * 0.5 * (theta_0 + theta_r)/(vk * g * q_h)
    L_n = (0.5*(theta_0+theta_r)*math.pow((u_r-u_0),3)/(g*u_r*(theta_0-theta_r))) * \
        ((math.log((z_r-d_0)/z_0m)-psi_h*(sp)) / \
          math.pow((math.log((z_r-d_0)/z_0m)-psi_m*sp), 2)) 
    L_err = abs((L_n - L)/L)
    val_dict['L'].append(L_n)
    err_dict['L'].append(L_err)  
    print('L: %.3f; L error: %.3f' % (L_n, L_err)) 
    
    # Value reset methodology
    # Attempt 1: no automatic reset of variables
    # Attempt 2: reset of variables to predefined, arbitrary IC
    j = 1
    while u_star_err > conv_crit:
        
        u_star_n = (vk*(u_r - u_0))/(math.log((z_r-d_0)/z_0m) - psi_m*sp)
        u_star_err = abs((u_star_n - u_star)/u_star)
        val_dict['u_star'].append(u_star_n)
        err_dict['u_star'].append(u_star_err)    
        print('u_star: %.3f; u_star error: %.3f' % (u_star_n, u_star_err)) 
                
        if j < iter_lim:
            u_star = u_star_n
            j += 1
        else:
            break
    
    z_0m_n = (z_r - d_0) / math.exp((0.4*u_r/u_star) + psi_m*sp)
    # z_0m_n = (z_r - d_0) / \
    # math.exp(0.4*u_r*(math.log((z_r-d_0)/z_0m) - psi_m*sp)/(vk*(u_r-u_0)) + psi_m*sp)
    z_0m_err = abs((z_0m_n - z_0m)/z_0m)
    val_dict['z_0m'].append(z_0m_n)
    err_dict['z_0m'].append(z_0m_err) 
    print('z_0m: %.3f; z_0m error: %.3f' % (z_0m_n, z_0m_err)) 
           
    if i < iter_lim:
        L = L_n
        z_0m = z_0m_n
        i += 1
    else:
        break
    
print('////////////////////')     
    
C_h = math.pow(vk, 2) / ((math.log((z_r-d_0)/z_0m) - psi_m*sp) * \
                        (math.log((z_r-d_0)/z_0m) - psi_h *sp))
           
q_h = rho * c_p * C_h * u_r * (theta_0 - theta_r)
        
##############################################################################################

## Calculate time elapsed for iterative solution
elapsed_time = time.time() - start_time

## Main function definition
def main():
    ### Print statements
    print('Sensible heat flux for selected pixel is: %.2f W/m^2' % q_h)
    print('Obukhov length for the selected pixel is: %.2f m' % L)
    print('Aerodynamic roughness length for the selected pixel is: %.4f' % z_0m)
    print('Friction velocity for the selected pixel is: %.4f' % u_star)
    print('Coefficient for the selected pixel is: %.4f' % C_h)
    print('Stability: %.4f' % sp)
    print('Time elapsed: %.4f s' % elapsed_time)    
    
    ## Plots
    
    # Outer loop parameter plot
    plot.plotter(len(val_dict['L']), ['L (m)', 'L error'], val_dict['L'], err1=err_dict['L'])
    plot.plotter(len(val_dict['z_0m']), ['z_0m (m)', 'z_0m error'], val_dict['z_0m'], err1=err_dict['z_0m'])
    # plot.plotter(len(val_dict['u_star']), ['u* (m/s)', 'u* error'], val_dict['u_star'], err1=err_dict['u_star'])

if __name__ == "__main__":
    main()
    
    

# H = f(C_h)
# --> C_h = f(z_0m, L)
# ----> L = f(u*, H)
# ------> u* = f(z_0m, L)
# ----> z_0m = f(u*, L)
