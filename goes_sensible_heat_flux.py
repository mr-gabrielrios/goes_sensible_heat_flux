### Objective
# The objective of this script is to generate the sensible heat flux based on input from GOES data

### Imports
import math
import plot
import time

### Constants
rho = 1.225     # Density of air (kg/m^3), assuming p = 1013.25 hPa and T = 15 degC
c_p = 1006      # Specific heat capacity of air (J/kg-K), assuming p = 1000 hPa and T = 0 degC
k = 0.4         # von Karman constant
z_r = 2         # Reference height (m), typically 2 m, based on air temperature model per Hrisko et al. (2020)
h_0 = 1         # Vegetation canopy height (m), assumption
g = 9.81        # Gravitational acceleration constant (m/s^2)
u_0 = 0         # Wind speed at surface (m/s), assumption
R = 287.05      # Gas constant for air (J/kg-K)

### Station inputs
p_air = 1013.25 # Atmospheric pressure (hPa), assumption - weather station input
u_r = 0.2         # Wind speed at reference height (m/s), assumption - weather station input

### Model inputs
T_lst = 340     # Land surface temperature (LST) (K), assumption - Hrisko model input
T_air = 290     # Air temperature at reference height (K), assumption - Hrisko model input

### Initial conditions
L = 1e6         # Initial condition for Obukhov length
u_star = u_0    # Initial condition for friction velocity
z_0m = .1       # Initial condition for aerodynamic roughness length
q_h = 10        # Initial condition for sensible heat flux

### Potential temperature calculation
def theta(T_abs):
    T_pot = T_abs * math.pow(1000/p_air, R/c_p)
    return T_pot

### Define potential temperatures
theta_0 = theta(T_lst)
theta_r = theta(T_air)

### Zero-plane displacement height, ref. Kim et al. (2019), Eqn. 3
# h_0 = vegetation height (m)
d_0 = math.exp(0.98*math.pow(math.log(h_0), 2)-0.15)

### Constant calculation as a function of Ri, ref. Kim et al. (2019), Eqn. 6
x = math.pow((1 - 15*(z_r - d_0)/L), 0.25)
### Alternate calculation of constant as f(Ri), ref. Pearlmutter et al. (2004), Eqn. 4
# T = (T_lst-T_air)/2 # Assume temperature is mean of ref. ht. temperature and LST
# Ri = g*((T_lst-T_air)/z_r)/(T*math.pow((u_0-u_r)/z_r, 2))
# x = math.pow((1 - 16*Ri), 0.25)

# Assume initial values for the errors
L_err = 1
z_0m_err = 1
i = 0

# Store iteration history of iteratively-solved variables
L_err_list = [L_err]
L_list = [L]
z_0m_err_list = [z_0m_err]
u_star_list = [u_star]
z_0m_list = [z_0m]
iter_list = [i]

# Set convergence criteria at 1%
conv_crit = 0.01
print('L error: %.2f' % L_err)


start_time = time.time()
### Iterative solution for parameters z_0m, q_h, and L
while L_err > conv_crit:
    
    print('Iteration #%d' % i)
    
    ### Stability similarity functions, ref. Kim et al. (2019), Eqn. 5, 6 and Garratt, Ch 3.3
    if (z_r/L) >= 0:
        psi_m = -5*(z_r - d_0)/L
        psi_h = psi_m
    elif (z_r/L) < 0:
        psi_m = math.log(((1 + math.pow(x,2))/2) * math.pow((1 + x)/2, 2)) \
            - 2*math.atan(x) + math.pi/2
        psi_h = 2*math.log((1 + math.pow(x, 2))/2) # Potential source of inaccuracy/error
    
    ### Write nested loop to simult. solve u_star and z_0m
    j = 0 # Nested loop iterand    
    while z_0m_err > conv_crit:
         ### Friction velocity calculation, ref. Kim et al. (2019), Eqn. 7
        # u_0 =  wind speed at surface (m/s), assumption
        u_star = (k*(u_r - u_0))/(math.log((z_r-d_0)/z_0m) - psi_m*((z_r - d_0)/L))
        
        ### Aerodynamic roughness length, ref. Kim et al. (2019), Eqn. 4
        # u_star = friction velocity (m/s)
        z_0m = z_r / math.exp((0.4*u_r/u_star) + psi_m*((z_r-d_0)/L)) # Should I forget the summation?
        
        z_0m_list.append(z_0m)
        u_star_list.append(u_star)
        
        z_0m_err = abs((z_0m_list[j] - z_0m_list[j-1])/z_0m_list[j-1])        
        
        z_0m = z_0m*z_0m_err
        
        # Limit nested loop to 50 iterations
        if j < 50:
            j += 1
        else:
            break
        
        print("u*: %.4f; z_0m: %.4f; z_0m_err: %.4f" % (u_star, z_0m, z_0m_err))
        
    ### Exchange coefficient calculation, ref. Kim et al. (2019), Eqn. 2
    # d_0 = zero-plane displacement height (m)
    # z_0m = aerodynamic roughness length (m)
    # psi_m = atmospheric stability parameter for momentum
    # psi_h = atmospheric stability parameter for sensible heat flux
    # L = Obukhov length (m)
    C_h = math.pow(k, 2) / ((math.log((z_r-d_0)/z_0m) - psi_m * (z_r-d_0)/L) * \
                            (math.log((z_r-d_0)/z_0m) - psi_h * (z_r-d_0)/L))
        
    ### Sensible heat flux calculation, ref. Kim et al. (2019), Eqn. 1
    # q_h = sensible heat flux (W/m^2)
    # C_h = sensible heat exchange coefficient
    # u_r = wind speed at reference height (m/s)
    # theta_0 = potential temperature at surface (K)
    # theta_r = potential temperature at reference height (K) 
    q_h = rho * c_p * C_h * u_r * (theta_0 - theta_r)
    
    ### Obukhov length calculation, ref. Kim et al. (2019), Eqn. 8
    L_new = -rho * c_p * math.pow(u_star, 3) * 0.5*(theta_0 + theta_r) / (k * g * q_h)
    
    ### Calculate error between last L value and the current one
    L_err = abs((L_new-L)/L)
    # Append to list
    L_err_list.append(L_err)
    L_list.append(L)
    
    ### Print iteration results
    print('L: %.4f  |  L_new: %.4f m  |  Error: %.4f  |  Previous: %.4f' % \
          (L, L_new, L_err_list[i], L_err_list[i-1]))
    
    ### Obukhov length re-calculation
    # Initial condition
    L = L_new
    
    iter_list.append(i)
    
    # Cold break at i = 51 to prevent run-on solution
    if i > 100:
        break
    else:
        i += 1    
    
### Calculate time elapsed for iterative solution
# CONUS Target: 1 min, 
elapsed_time = time.time() - start_time

### Main function defintion
def main():
    ### Print statements
    print('Sensible heat flux for selected pixel is: %.2f W/m^2' % q_h)
    print('Obukhov length for the selected pixel is: %.2f m' % L)
    print('Time elapsed: %.4f s' % elapsed_time)    
    # plot.plotter(iter_list, L_list, L_err_list)

if __name__ == "__main__":
    main()