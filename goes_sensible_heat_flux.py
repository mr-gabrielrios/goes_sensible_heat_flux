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
rho = 1.225
c_p = 1006   
k = 0.4
z_r = 2
h_0 = 1
g = 9.81
u_0 = .1
R = 287.05 

### 2. Dynamic inputs
p_air = 1013.25 
u_r = 2.5
T_lst = 295
T_air = 290

### 3. Initial conditions
L = 1e15
u_star = 1
z_0m = 1.1
q_h = 1      
psi_m = 0      
psi_h = 0
C_h = 1   

### 4. Calculate secondary parameters

## Potential temperature calculation
def theta(T_abs):
    T_pot = T_abs * math.pow(1000/p_air, R/c_p)
    return T_pot

## Mean calculation for series of non-list items
def mean(*args):
    sigma = 0
    for arg in args:
        sigma += arg
    return (sigma/len(args))

## Define potential temperatures
theta_0 = theta(T_lst)
theta_r = theta(T_air)

## Zero-plane displacement height, ref. Kim et al. (2019), Eqn. 3
d_0 = math.exp(0.9793*math.log(h_0)-0.1536)

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
q_h_err, z_0m_err, u_star_err, C_h_err, L_err = [1, 1, 1, 1, 1]
## Create dictionary to store error history in lists for each variable
err_dict = {'q_h': [q_h_err], 
       'C_h': [err_val],
       'z_0m': [err_val],
       'psi_m': [err_val],
       'psi_h': [err_val],
       'u_star': [err_val],
       'L': [err_val]}

u_star_err_max = [1, 1, 1]

# Set convergence criteria at 1%
conv_crit = 0.01
ilc = False
iter_limit = 100

start_time = time.time()
i = 1 

print(mean(1, 2, 3))

#################################################################

while L_err > conv_crit:
    zeta = (z_r-d_0)/L
    zeta_0 = (z_r-d_0)/z_0m
    print('Stability parameter: %.4f' % zeta)
    
    x = math.pow(1-16*zeta, 0.25)
    if zeta < 0:
        psi_m = 2*math.log((1+x)/2) + math.log((1+math.pow(x,2))/2) - 2*math.atan(x) + math.pi/2
        psi_h = 2*math.log((1+math.pow(x,2))/2)
    else:
        psi_m = -5*zeta
        psi_h = psi_m    

    if ilc:
        z_0m_err, u_star_err = [0, 0]
    else:
        z_0m_err, u_star_err = [1, 1]
        z_0m, u_star = [val_dict['z_0m'][0], val_dict['u_star'][0]]        
    
    j = 1
    while z_0m_err + u_star_err > conv_crit:
        print ('\t Sub-iteration %d parameters: z_0m = %.3f; u_star = %.3f' % (j, z_0m, u_star))
        ######
        z_0m_n = (z_r - d_0) / math.exp(0.4*u_r/u_star + psi_m*zeta)
        z_0m_err = abs((z_0m_n - z_0m)/z_0m)     
        val_dict['z_0m'].append(z_0m_n)
        err_dict['z_0m'].append(z_0m_err) 
        ######
        
        ######
        print('\t u_star %.4f | u_star_err: %.4f' % (u_star, u_star_err))
        u_star_n = k*(u_r - u_0)/(math.log(zeta_0) - psi_m*zeta)  
        u_star_err = abs((u_star_n - u_star)/u_star)  
        val_dict['u_star'].append(u_star)
        err_dict['u_star'].append(u_star_err)
        ######
        
        if j == 1:
            u_star_err_max.append(u_star_err)
        
        print('\t ------')
        
        if j < iter_limit:
            z_0m = z_0m_n
            u_star = u_star_n
            j += 1
        else:
            break
    
    print('z_0m: %.4f | z_0m_err: %.4f' % (z_0m, z_0m_err))
    print('u_star %.4f | u_star_err: %.4f' % (u_star, u_star_err)) 
    
    if (u_star_err_max[-1] < max(u_star_err_max)) and (u_star_err_max[-2] < u_star_err_max[-1]):
        ilc = True
    
    if not(ilc):
        print('Number of sub-iterations to convergence: %d' % j)
    ######
    print('C_h %.4f | C_h_err: %.4f' % (C_h, C_h_err))
    C_h_n = math.pow(k, 2) / ((math.log(zeta_0) - psi_m*zeta) * (math.log(zeta_0) - psi_h*zeta))
    C_h_err = abs((C_h_n - C_h)/C_h)
    val_dict['C_h'].append(C_h_n)
    err_dict['C_h'].append(C_h_err)
    ######    
    
    ######
    print('q_h %.4f | q_h_err: %.4f' % (q_h, q_h_err))
    q_h_n = rho * c_p * C_h * u_r * (theta_0 - theta_r)
    q_h_err = abs((q_h_n - q_h)/q_h)
    val_dict['q_h'].append(q_h_n)
    err_dict['q_h'].append(q_h_err)
    ######
    
    ######    
    print('L: %.4f | L_err: %.4f' % (L, L_err))
    L_n = -rho * c_p * math.pow(u_star, 3) * 0.5*(theta_0 + theta_r) / (k * g * q_h)
    L_err = abs((L_n - L)/L)
    val_dict['L'].append(L_n)   
    err_dict['L'].append(L_err) 
    ###### 
    
    if i < iter_limit:
        C_h = C_h_n
        q_h = q_h_n
        L = L_n
        i += 1
    else:
        break
    
    print('###### END ITERATION %d' % i)

#################################################################

## Calculate time elapsed for iterative solution
elapsed_time = time.time() - start_time

## Main function definition
def main():
    ### Print statements
    print('Sensible heat flux for selected pixel is: %.4f W/m^2' % q_h)
    print('Obukhov length for the selected pixel is: %.4f m' % L)
    print('Roughness length for the selected pixel is: %.4f m' % z_0m)
    plot.plotter(len(val_dict['L']), ['L (m)', 'L error'], val_dict['L'], err_data=err_dict['L'], length_plot=True) 
    # plot.plotter(len(val_dict['z_0m']), ['z_0m (m)', 'z_0m error'], val_dict['z_0m'], err_data=err_dict['z_0m'])
    # plot.plotter(len(val_dict['u_star']), ['u_star (m/s)', 'u_star error'], val_dict['u_star'], err_data=err_dict['u_star'])  
    plot.plotter(len(val_dict['C_h']), ['C_h', 'C_h error'], val_dict['C_h'], err_data=err_dict['C_h']) 
    plot.plotter(len(val_dict['q_h']), ['q_h (W/m^2', 'q_h error'], val_dict['q_h'], err_data=err_dict['q_h']) 
    # plot.plotter(len(u_star_err_max), ['u_star_error'], u_star_err_max) 

if __name__ == "__main__":
    main()  

# print('\t Stability parameter: %.4f' % ((z_r-d_0)/L))
                    
# x = math.pow((1 - 15*(z_r - d_0)/L), 0.25)

# if (z_r/L) >= 0:
#     psi_m = -5*(z_r - d_0)/L
#     psi_h = psi_m
# elif (z_r/L) < 0:
#     psi_m = math.log(((1 + math.pow(x,2))/2) * math.pow((1 + x)/2, 2)) \
#         - 2*math.atan(x) + math.pi/2
#     psi_h = 2*math.log((1 + math.pow(x, 2))/2)
# ######
    
# ######    
# print('\t L: %.4f | L_err: %.4f' % (L, L_err))
# L_n = -rho * c_p * math.pow(u_star, 3) * 0.5*(theta_0 + theta_r) / (k * g * q_h)
# # L_n = -(0.5*(theta_0+theta_r)*math.pow((u_r-u_0),3)/(g*u_r*(theta_0-theta_r))) * ((math.log((z_r-d_0)/z_0m)-psi_h*((z_r-d_0)/L)) / math.pow((math.log((z_r-d_0)/z_0m)-psi_m*((z_r-d_0)/L)), 2)) 
# L_err = abs((L_n - L)/L)
# val_dict['L'].append(L_n)   
# err_dict['L'].append(L_err) 
# ######  

# ######
# print('\t z_0m: %.4f | z_0m_err: %.4f' % (z_0m, z_0m_err))
# z_0m_n = (z_r - d_0) / math.exp(0.4*u_r/u_star + psi_m*((z_r-d_0)/L))
# # z_0m_n = (z_r - d_0) / math.exp(0.4*u_r*(math.log((z_r-d_0)/z_0m) - psi_m*((z_r-d_0)/L))/(k*(u_r-u_0)) + psi_m*((z_r-d_0)/L))     
# z_0m_err = abs((z_0m_n - z_0m)/z_0m)     
# val_dict['z_0m'].append(z_0m_n)
# err_dict['z_0m'].append(z_0m_err)
# ######

# ######
# print('u_star %.4f | u_star_err: %.4f' % (u_star, u_star_err))
# u_star_n = (k*(u_r - u_0))/(math.log((z_r-d_0)/z_0m) - psi_m*((z_r - d_0)/L))  
# u_star_err = abs((u_star_n - u_star)/u_star)  
# val_dict['u_star'].append(u_star)
# err_dict['u_star'].append(u_star_err)
# ######

# ######
# print('C_h %.4f | C_h_err: %.4f' % (C_h, C_h_err))
# C_h_n = math.pow(k, 2) / ((math.log((z_r-d_0)/z_0m) - psi_m * (z_r-d_0)/L) * \
#                             (math.log((z_r-d_0)/z_0m) - psi_h * (z_r-d_0)/L))
# C_h_err = abs((C_h_n - C_h)/C_h)
# val_dict['C_h'].append(C_h_n)
# err_dict['C_h'].append(C_h_err)
# ######

# ######
# print('q_h %.4f | q_h_err: %.4f' % (q_h, q_h_err))
# q_h_n = rho * c_p * C_h * u_r * (theta_0 - theta_r)
# q_h_err = abs((q_h_n - q_h)/q_h)
# val_dict['C_h'].append(q_h_n)
# err_dict['C_h'].append(q_h_err)
# ######