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
import plot
import time

### 1. Constants
rho = 1.225
c_p = 1006   
vk = 0.4
z_r = 10
h_0 = 1
g = 9.81
R = 287.05 

### 2. Dynamic inputs
p_air = 1013.25 
u_r = 1
T_lst = 295
T_air = 290

### 3. Initial conditions
L = 1e15
u_star = 1
z_0m = 2 # Per Stull, Fig 9.6
z_t = z_0m/7
q_h = 1      
psi_m = 0      
psi_h = 0
C_h = 1   
C_d = 1   

### 4. Calculate secondary parameters

## Potential temperature calculation
def theta(T_abs):
    T_pot = T_abs * math.pow(1000/p_air, R/c_p)
    return T_pot

### 5. Program management
## Create dictionary to store value history in lists for each variable
## Note: values without explicit initial conditions are set at 0
val_dict = {'q_h': [q_h], 
            'C_h': [C_h],
            'z_0m': [z_0m],
            'psi_m': [0],
            'psi_h': [0],
            'u_star': [u_star],
            'L': [L]}

## Create dictionary to store error history in lists for each variable
# Assume initial values for the errors at 100%
err_val = 1
err_dict = {'q_h': [err_val], 
       'C_h': [err_val],
       'z_0m': [err_val],
       'psi_m': [err_val],
       'psi_h': [err_val],
       'u_star': [err_val],
       'L': [err_val]}
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
start_time = time.time()
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
        
    ## Heat transfer coefficient for heat flux, or Stanton number
    C_h_n = math.pow(vk,2)/ \
        ((math.log(z/z_0m)-psi_m*zeta)*(math.log(z/z_t)-psi_h*zeta))
    
    ## Drag coefficient for momentum flux
    C_d_n = math.pow(vk,2)/ math.pow((math.log(z/z_0m)-psi_m*zeta), 2)
    
    ## Friction velocity
    # Ref. Launiainen et al. (1990), Eqn. A7
    u_star_n = math.sqrt(C_d*u_r)
    
    ## Sensible heat flux
    # Ref. Launiainen et al. (1990), Eqn. 5
    q_h_n = rho*c_p*C_h*u_r*(theta(T_lst)-theta(T_air))
    
    ## Obukhov length (ignoring latent heat flux)
    # Ref. Launiainen et al. (A10), derivation performed by Rios
    L_n = -rho*c_p*math.pow(u_star,3)*0.5*(theta(T_lst)+theta(T_air))/(vk*g*q_h)
    
    # Calculate error between current and previous iteration
    L_err = abs((L_n-L)/L)
    print('L: %.3f; L_n: %.3f; L_err: %.3f' % (L, L_n, L_err))
    
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

## Main function definition
def main():
    ### Print statements
    print('Sensible heat flux for selected pixel is: %.4f W/m^2' % q_h)
    print('Obukhov length for the selected pixel is: %.4f m' % L)
    print('Roughness length for the selected pixel is: %.4f m' % z_0m)
    print('Heat transfer coefficient for the selected pixel is: %.4f' % C_h)

if __name__ == "__main__":
    main() 