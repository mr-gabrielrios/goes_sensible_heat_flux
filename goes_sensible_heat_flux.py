
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
vk = 0.4         # von Karman constant
z_r = 2.3         # Reference height (m), typically 2 m, based on air temperature model per Hrisko et al. (2020)
                # 3 m chosen to match with Kim et al. (2019), Case GH1
h_0 = 2.6      # Vegetation canopy height (m), assumption
g = 9.81        # Gravitational acceleration constant (m/s^2)
u_0 = 0.1         # Wind speed at surface (m/s), assumption
R = 287.05      # Gas constant for air (J/kg-K)

### 2. Dynamic inputs
## Station inputs
p_air = 1013.25 # Atmospheric pressure (hPa), assumption - weather station input
u_r = 2         # Wind speed at reference height (m/s), assumption - weather station input
## Model inputs
T_lst = 300     # Land surface temperature (LST) (K), assumption - Hrisko model input
T_air = 290     # Air temperature at reference height (K), assumption - Hrisko model input

### 3. Initial conditions
L = 1        # Initial condition for Obukhov length
u_star = 1      # Initial condition for friction velocity
z_0m = 1       # Initial condition for aerodynamic roughness length
q_h = 1e3       # Initial condition for sensible heat flux
x = 1         # Initial condition for parameter that's a function of the Richardson number
psi_m = 0       # Initial condition for similarity functions for stability
psi_h = 0
C_h = 1         # Initial condition for the heat transfer coefficient

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
            'C_h': [C_h],
            'z_0m': [z_0m],
            'psi_m': [0],
            'psi_h': [0],
            'u_star': [u_star],
            'L': [L]}

## Assume initial values for the errors at 100%
err_val = 1
q_h_err, z_0m_err, u_star_err, L_err = [1, 1, 1, 1]
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
iter_lim = 10
conv_bool = False

### 5. Iterative solution for parameters z_0m, q_h, and L
## Start timing solution algorithm
start_time = time.time()
# Initialize iterand at 1 to allow for a list item with index 0 to be called
i = 1 
## Initialize outer while loop
while q_h_err > conv_crit:
    print('#####################################')
    print('Iteration #%d' % i)
    L_err = 1
    j = 1
    while (L_err) > conv_crit:        
        
        print('j-iteration #%d' % j)
        # Print list of initial sub-iteration values (keep commented except for troubleshooting)
        print('Stability parameter: %.4f' % ((z_r-d_0)/L))
        print('L: %.4f ' % (L))
                   
        ## Obukhov length calculation, ref. Kim et al. (2019), Eqn. 8
        # L_n = -rho * c_p * math.pow(u_star, 3) * 0.5*(theta_0 + theta_r) / (vk * g * q_h)
        ## Obukhov length calculation, modified version of Kim et al. (2019), Eqn. 8 by Rios
        ## Note: this modification removes dependence on u_star to allow for recursive calculation
        L_n = (0.5*(theta_0+theta_r)*math.pow((u_r-u_0),3)/(g*u_r*(theta_0-theta_r))) * \
            ((math.log((z_r-d_0)/z_0m)-psi_h*((z_r-d_0)/L)) / \
             math.pow((math.log((z_r-d_0)/z_0m)-psi_m*((z_r-d_0)/L)), 2))        
        
        ## Append sub-iteration values 
        val_dict['L'].append(L_n)     
        # Print list of final sub-iteration values (keep commented except for troubleshooting)
        # print('L_n: %.4f' % (L_n))   
        
        # Calculate and print errors
        # L_err = abs((val_dict['L'][j] - val_dict['L'][j-1])/val_dict['L'][j-1])
        L_err = abs((L_n - L)/L)
        err_dict['L'].append(L_err)
        
        print('L error: %.4f' % L_err)       
                 
        # Limit nested loop to 20 iterations
        if j < iter_lim:
            # Reset parameters to new values
            L = L_n       
            j += 1
        else:
            print("LIMIT REACHED, NO CONVERGENCE")
            break   
    
    # Print lists for the inner loop values
    # print('### Values (z_0m - 1st line, L - 2nd line) ###')
    # print(*val_dict['z_0m'], sep=', ')
    # print(*val_dict['L'], sep=', ')
    # print('### Errors (z_0m - 1st line, L - 2nd line) ###')
    # print(*err_dict['z_0m'], sep=', ')
    # print(*err_dict['L'], sep=', ')
    # Print final z_0m value
    # print('z_0m = %.4f | L = %.4f' % (z_0m, L))    
    # break # This serves to isolate the inner loop for troubleshooting    
    
    ## Constant calculation as a function of Ri, ref. Kim et al. (2019), Eqn. 6
    x = math.pow((1 - 15*(z_r - d_0)/L), 0.25)
    ## Alternate calculation of constant as f(Ri), ref. Pearlmutter et al. (2004), Eqn. 4
    # T = (T_lst-T_air)/2 # Assume temperature is mean of ref. ht. temperature and LST
    # Ri = g*((T_lst-T_air)/z_r)/(T*math.pow((u_r-u_0)/z_r, 2))
    # x = math.pow((1 - 16*Ri), 0.25)
    
    ## Stability similarity functions, ref. Kim et al. (2019), Eqn. 5, 6 and Garratt, Ch 3.3
    if (z_r/L) >= 0:
        psi_m = -5*(z_r - d_0)/L
        psi_h = psi_m
    elif (z_r/L) < 0:
        psi_m = math.log(((1 + math.pow(x,2))/2) * math.pow((1 + x)/2, 2)) \
            - 2*math.atan(x) + math.pi/2
        psi_h = 2*math.log((1 + math.pow(x, 2))/2) # Potential source of inaccuracy/error
    
    #############################################
    print('--------------------------------')
    z_0m_err = 1
    k = 1
    while z_0m_err > conv_crit:
        
        print('k-iteration #%d' % k)
        print('z_0m: %.4f' % z_0m)
        # print('z_0m: %.4f ' % (z_0m))  
        
        ## Aerodynamic roughness length, ref. Kim et al. (2019), Eqn. 4 and 
        ## Brummer et al. (2002), Eqn. 1
        # u_star = friction velocity (m/s)
        z_0m_n = (z_r - d_0) / math.exp((0.4*u_r/u_star) + psi_m*((z_r-d_0)/L))  
        ## Aerodynamic roughness length, modified version of Kim et al. (2019) and Brummer et al. (2002) by Rios
        ## Note: this modification removes dependence on u_star to allow for recursive calculation
        
        # z_0m_n = (z_r - d_0) / \
            # math.exp(0.4*u_r*(math.log((z_r-d_0)/z_0m) - psi_m*((z_r-d_0)/L))/(vk*(u_r-u_0)) + psi_m*((z_r-d_0)/L))          
                
        ## Append sub-iteration values
        val_dict['z_0m'].append(z_0m_n)
        
        # Print list of final sub-iteration values (keep commented except for troubleshooting)
        # print('z_0m_n: %.4f' % (z_0m_n)) 
        
        # Calculate and print errors
        # z_0m_err = abs((val_dict['z_0m'][i] - val_dict['z_0m'][i-1])/val_dict['z_0m'][i-1])
        z_0m_err = abs((z_0m_n - z_0m)/z_0m)
        err_dict['z_0m'].append(z_0m_err) 
        
        print('z_0m error: %.4f' % z_0m_err)
                
        # Limit nested loop to 20 iterations
        if k < iter_lim:
            # Reset parameters to new values
            z_0m = z_0m_n
            k += 1
        else:
            print("LIMIT REACHED, NO CONVERGENCE")
            break 
    
    #####################################################################
    
    print('-----------------------------------')
    
    ### u_star loop
    m = 1
    u_star_err = 1
    while u_star_err > conv_crit:
        
        print('m-iteration #%d' % m)
        print('u_star: %.4f' % u_star)
        
        ## Friction velocity calculation, ref. Kim et al. (2019), Eqn. 7
        # u_0 =  wind speed at surface (m/s), assumption
        u_star_n = (vk*(u_r - u_0))/(math.log((z_r-d_0)/z_0m) - psi_m*((z_r - d_0)/L)) 
        u_star_err = abs((u_star_n - u_star)/u_star)
        val_dict['u_star'].append(u_star_n) 
        err_dict['u_star'].append(u_star_err)  
                            
        print('u_star error: %.4f' % u_star_err)   
        
        # Limit nested loop to 20 iterations
        if m < iter_lim:
            u_star = u_star_n
            m += 1
        else:
            print("LIMIT REACHED, NO CONVERGENCE")
            break   
    
    print('----------------------------------------------------')

    ### C_h loop
    n = 1
    C_h_err = 1
    while C_h_err > conv_crit:
        
        print('n-iteration #%d' % n)
        print('C_h: %.4f' % C_h)

        ## Exchange coefficient calculation, ref. Kim et al. (2019), Eqn. 2
        # d_0 = zero-plane displacement height (m)
        # z_0m = aerodynamic roughness length (m)
        # psi_m = atmospheric stability parameter for momentum
        # psi_h = atmospheric stability parameter for sensible heat flux
        # L = Obukhov length (m)
        C_h = math.pow(vk, 2) / ((math.log((z_r-d_0)/z_0m) - psi_m * (z_r-d_0)/L) * \
                                (math.log((z_r-d_0)/z_0m) - psi_h * (z_r-d_0)/L))
        val_dict['C_h'].append(C_h)  
        C_h_err = abs((val_dict['C_h'][i] - val_dict['C_h'][i-1])/val_dict['C_h'][i-1])      
        err_dict['C_h'].append(q_h_err)
        
        print('C_h error: %.4f' % C_h_err)  
        
        # Limit nested loop to 20 iterations
        if n < iter_lim:
            C_h = u_star_n
            n += 1
        else:
            print("LIMIT REACHED, NO CONVERGENCE")
            break   
        
    ## Sensible heat flux calculation, ref. Kim et al. (2019), Eqn. 1
    # q_h = sensible heat flux (W/m^2)
    # C_h = sensible heat exchange coefficient
    # u_r = wind speed at reference height (m/s)
    # theta_0 = potential temperature at surface (K)
    # theta_r = potential temperature at reference height (K) 
    q_h = rho * c_p * C_h * u_r * (theta_0 - theta_r)
    val_dict['q_h'].append(q_h)     
    q_h_err = abs((val_dict['q_h'][i] - val_dict['q_h'][i-1])/val_dict['q_h'][i-1])   
    err_dict['q_h'].append(q_h_err)
    
    print('q_h error: %.4f' % q_h_err)    
    
    print('\n')
    # Cold break to prevent run-on solution
    if i < iter_lim:
        i += 1
    else:
        i += 1           
        print("LIMIT REACHED, NO CONVERGENCE")
        break
    
## Calculate time elapsed for iterative solution
elapsed_time = time.time() - start_time

## Main function definition
def main():
    ### Print statements
    print('Sensible heat flux for selected pixel is: %.2f W/m^2' % q_h)
    print('Obukhov length for the selected pixel is: %.2f m' % L)
    print('Time elapsed: %.4f s' % elapsed_time)    
    
    ## Plots
    
    # Outer loop parameter plot
    plot.plotter(i, ['Q_h (W/m^2)', 'Q_h error'], val_dict['q_h'], err1=err_dict['q_h'])
    plot.plotter(len(val_dict['u_star']), ['u* (m/s)', 'u* error'], val_dict['u_star'], err1=err_dict['u_star'])
    plot.plotter(len(val_dict['L']), ['L (m)', 'L error'], val_dict['L'], err1=err_dict['L'])
    plot.plotter(len(val_dict['z_0m']), ['z_0m (m)', 'z_0m error'], val_dict['z_0m'], err1=err_dict['z_0m'])
    plot.plotter(len(val_dict['C_h']), ['C_h', 'C_h error'], val_dict['C_h'], err1=err_dict['C_h'])

if __name__ == "__main__":
    main()