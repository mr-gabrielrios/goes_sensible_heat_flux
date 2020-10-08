### TRIAL 1 - PROGRAM FLOW SIMILAR TO LITERATURE        

### Objective
# The objective of this script is to generate the sensible heat flux based on input from GOES data

### Program flow
# 1. Define constants
# 2. Define dynamic inputs (e.g. model data, station data)
# Note: The dynamic inputs are currently static, for simplicity, while the algorithm is finalized
# 3. Define initial conditions
# 4. Calculate secondary parameters to support solution algorithm
# 5. Perform solution algorithm
# Note: This is currently an iterative algorithm with a nested loop. Outer loop solves for L, inner loop simultaneously solves u* (friction velocity) and z_0m (aero. roughness length)
# CRITERIA: Solution criteria for each loop is that each variable doesn't change by > 1% per iteration


### 0. Imports
import math
import plot
import time

### 1. Constants
rho = 1.225     # Density of air (kg/m^3), assuming p = 1013.25 hPa and T = 15 degC
c_p = 1006      # Specific heat capacity of air (J/kg-K), assuming p = 1000 hPa and T = 0 degC
k = 0.4         # von Karman constant
z_r = 2.5        # Reference height (m), typically 2 m, based on air temperature model per Hrisko et al. (2020)
                # 3 m chosen to match with Kim et al. (2019), Case GH1
h_0 = 2.60      # Vegetation canopy height (m), assumption
g = 9.81        # Gravitational acceleration constant (m/s^2)
u_0 = .1         # Wind speed at surface (m/s), assumption
R = 287.05      # Gas constant for air (J/kg-K)

### 2. Dynamic inputs
## Station inputs
p_air = 1013.25 # Atmospheric pressure (hPa), assumption - weather station input
u_r = 0.3         # Wind speed at reference height (m/s), assumption - weather station input
## Model inputs
T_lst = 303     # Land surface temperature (LST) (K), assumption - Hrisko model input
T_air = 290    # Air temperature at reference height (K), assumption - Hrisko model input

### 3. Initial conditions
L = 100         # Initial condition for Obukhov length
u_star = 1      # Initial condition for friction velocity
z_0m = .1       # Initial condition for aerodynamic roughness length
q_h = 100       # Initial condition for sensible heat flux
x = 0.1         # Initial condition for parameter that's a function of the Richardson number
psi_m = 0       # Initial condition for similarity functions for stability
psi_h = 0
C_h = 1         # Initial condition for heat transfer coefficient

### 4. Calculate secondary parameters

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

## Create dictionary to store value history in lists for each variable
## Note: values without explicit initial conditions are set at 0
val_dict = {'q_h': [q_h], 
            'C_h': [0],
            'z_0m': [z_0m],
            'psi_m': [0],
            'psi_h': [0],
            'u_star': [u_star],
            'L': [L]}

## Assume initial values for the errors at 100%
err_val = 1
q_h_err, z_0m_err, u_star_err, C_h_err, L_err = [1, 1, 1, 1, 1]
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
iter_limit = 100

### 5. Iterative solution for parameters z_0m, q_h, and L
## Start timing solution algorithm
start_time = time.time()
# Initialize iterand at 1 to allow for a list item with index 0 to be called
i = 1 
## Initialize outer while loop
while L_err > conv_crit:
    
    print('#####################################')
    print('Iteration %d' % i)
    print('Stability parameter: %.4f' % ((z_r-d_0)/L))
    
    ######
    print('z_0m: %.4f | z_0m_err: %.4f' % (z_0m, z_0m_err))
    z_0m_n = (z_r - d_0) / math.exp(0.4*u_r/u_star + psi_m*((z_r-d_0)/L))
    # z_0m_n = (z_r - d_0) / \
    #     math.exp(0.4*u_r*(math.log((z_r-d_0)/z_0m) - psi_m*((z_r-d_0)/L))/(k*(u_r-u_0)) + psi_m*((z_r-d_0)/L))          
    val_dict['z_0m'].append(z_0m_n)
    z_0m_err = abs((z_0m_n - z_0m)/z_0m)
    err_dict['z_0m'].append(z_0m_err)
    ######
    
    j = 1
    while q_h_err > conv_crit:
        print('----------------------------')
        print('Sub-iteration %d' % j)
        
        ######
        print('u_star %.4f | u_star_err: %.4f' % (u_star, u_star_err))
        u_star_n = (k*(u_r - u_0))/(math.log((z_r-d_0)/z_0m) - psi_m*((z_r - d_0)/L))  
        u_star_err = abs((u_star_n - u_star)/u_star)  
        val_dict['u_star'].append(u_star)
        err_dict['u_star'].append(u_star_err)
        ######
        
        ######
        print('C_h %.4f | C_h_err: %.4f' % (C_h, C_h_err))
        C_h_n = math.pow(k, 2) / ((math.log((z_r-d_0)/z_0m) - psi_m * (z_r-d_0)/L) * \
                                    (math.log((z_r-d_0)/z_0m) - psi_h * (z_r-d_0)/L))
        C_h_err = abs((C_h_n - C_h)/C_h)
        val_dict['C_h'].append(C_h_n)
        err_dict['C_h'].append(C_h_err)
        ######
        
        ######
        print('q_h %.4f | q_h_err: %.4f' % (q_h, q_h_err))
        q_h_n = rho * c_p * C_h * u_r * (theta_0 - theta_r)
        q_h_err = abs((q_h_n - q_h)/q_h)
        val_dict['C_h'].append(q_h_n)
        err_dict['C_h'].append(q_h_err)
        ######
        
        # Cold break to prevent run-on solution
        if j < iter_limit:
            u_star = u_star_n
            C_h = C_h_n
            q_h = q_h_n
            j += 1
        else:
            break 
    
    ######        
    x = math.pow((1 - 15*(z_r - d_0)/L), 0.25)
    
    if (z_r/L) >= 0:
        psi_m = -5*(z_r - d_0)/L
        psi_h = psi_m
    elif (z_r/L) < 0:
        psi_m = math.log(((1 + math.pow(x,2))/2) * math.pow((1 + x)/2, 2)) \
            - 2*math.atan(x) + math.pi/2
        psi_h = 2*math.log((1 + math.pow(x, 2))/2)
    ######
        
    ######    
    print('L: %.4f | L_err: %.4f' % (L, L_err))
    L_n = -rho * c_p * math.pow(u_star, 3) * 0.5*(theta_0 + theta_r) / (k * g * q_h)
    # L_n = (0.5*(theta_0+theta_r)*math.pow((u_r-u_0),3)/(g*u_r*(theta_0-theta_r))) * \
    #     ((math.log((z_r-d_0)/z_0m)-psi_h*((z_r-d_0)/L)) / \
    #      math.pow((math.log((z_r-d_0)/z_0m)-psi_m*((z_r-d_0)/L)), 2)) 
    L_err = abs((L_n - L)/L)
    val_dict['L'].append(L_n)   
    err_dict['L'].append(L_err) 
    ######
    
    # Cold break to prevent run-on solution
    if i < iter_limit:
        z_0m = z_0m_n
        L = L_n
        i += 1
        print('\n')
    else:
        break   
   
## Calculate time elapsed for iterative solution
elapsed_time = time.time() - start_time

## Main function definition
def main():
    ### Print statements
    print('Sensible heat flux for selected pixel is: %.2f W/m^2' % q_h)
    print('Obukhov length for the selected pixel is: %.2f m' % L)
    print('Roughness length for the selected pixel is: %.2f m' % z_0m)
    print('Time elapsed: %.4f s' % elapsed_time)   
    
    # plot.plotter(len(val_dict['q_h']), ['Q_h (W/m^2)', 'Q_h error'], val_dict['q_h'], err1=err_dict['q_h'])    
    # plot.plotter(len(val_dict['L']), ['L (m)', 'L error'], val_dict['L'], err1=err_dict['L'])    
    # plot.plotter(len(val_dict['z_0m']), ['z_0m (m)', 'z_0m error'], val_dict['z_0m'], err1=err_dict['z_0m'])

if __name__ == "__main__":
    main()  
