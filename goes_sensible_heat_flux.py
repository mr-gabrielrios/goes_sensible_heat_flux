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
z_r = 2         # reference height (m), typically 2 m, based on air temperature model per Hrisko et al. (2020)
h_0 = 1         # vegetation canopy height (m), assumption
g = 9.81        # gravitational acceleration constant (m/s^2)
u_0 = 0         # wind speed at surface (m/s), assumption
R = 287.05      # Gas constant for air (J/kg-K)

### Station inputs
p_air = 1013.25 # Atmospheric pressure (hPa), assumption - weather station input
u_r = 1         # wind speed at reference height (m/s), assumption - weather station input

### Model inputs
T_lst = 350     # land surface temperature (LST) (K), assumption - Hrisko model input
T_air = 285     # air temperature at reference height (K), assumption - Hrisko model input

### Initial conditions
L = 1000       # Initial condition for Obukhov length
z_0m = .1        # Initial condition for aerodynamic roughness length
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

err_list = []
L_list = []
err = 2 # Assume an initial value for the error

conv_crit = 0.01
print('Error: %.2f' % err)

iter_list = []
i = 0

start_time = time.time()
### Iterative solution for parameters z_0m, q_h, and L
while err > conv_crit:
    
    print('Iteration #%d' % i)
    print('z_r: %.2f  |  d_0: %.2f  |  L: %.2f' % (z_r, d_0, L))
    
    ### Constant calculation as a function of Ri, ref. Kim et al. (2019), Eqn. 6
    x = math.pow((1 - 15*(z_r - d_0)/L), 0.25)
    ### Alternate calculation of constant as f(Ri), ref. Pearlmutter et al. (2004), Eqn. 4
    # T = (T_lst-T_air)/2 # Assume temperature is mean of ref. ht. temperature and LST
    # Ri = g*((T_lst-T_air)/z_r)/(T*(u_0-u_r)/z_r)
    # x = math.pow((1 - 16*Ri), 0.25)
    
    ### Stability similarity functions, ref. Kim et al. (2019), Eqn. 5, 6 and Garratt, Ch 3.3
    if (z_r/L) >= 0:
        psi_m = -5*(z_r - d_0)/L
        psi_h = psi_m
    elif (z_r/L) < 0:
        psi_m = math.log(((1 + math.pow(x,2))/2) * math.pow((1 + x)/2, 2)) \
            - 2*math.atan(x) + math.pi/2
        psi_h = 2*math.log((1 + math.pow(x, 2))/2)
    
    ### Friction velocity calculation, ref. Kim et al. (2019), Eqn. 7
    # u_0 =  wind speed at surface (m/s), assumption
    u_star = (k*(u_r - u_0))/(math.log((z_r-d_0)/z_0m) - psi_m*((z_r - d_0)/L))
    
    ### Aerodynamic roughness length, ref. Kim et al. (2019), Eqn. 4
    # u_star = friction velocity (m/s)
    z_0m = z_r / math.exp((0.4*u_r/u_star) + psi_m*((z_r-d_0)/L))
    
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
    err = abs((L_new-L)/L)
    # Append to list
    err_list.append(err)
    L_list.append(L)
    
    print('L: %.4f  |  L_new: %.4f m  |  Error: %.4f' % (L, L_new, err))
       
    if i == 0:
        L = L_new
    else:
        if (err_list[i] < err_list[i-1]):
            L = L + L_new*math.pow(err, 2)
            print("subtracting")
        else:
            L = L - L_new*math.pow(err, 2)
            print("adding")
    
    if i > 50:
        break
    
    iter_list.append(i)
    i += 1    
    
elapsed_time = time.time() - start_time

### Main function defintion
def main():
    ### Print statements
    print('Sensible heat flux for selected pixel is: %.2f W/m^2' % q_h)
    print('Obukhov length for the selected pixel is: %.2f m' % L)
    print('Time elapsed: %.4f s' % elapsed_time)    
    plot.plotter(iter_list, L_list, err_list)

if __name__ == "__main__":
    main()
