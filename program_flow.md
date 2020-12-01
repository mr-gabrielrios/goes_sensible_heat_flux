# Sensible Heat Flux Product - Program Flow

The objective of this document is to catalogue the program flow of the sensible heat flux package.

## Output(s)
- Sensible heat flux, (Q<sub>H</sub>)
    - Sensible heat flux for each GOES pixel. Intent is to build upon [Josh Hrisko's air temperature model](https://github.com/makerportal).

## Input(s)
- Land surface,skin temperature (LST), (T<sub>LST</sub>)
- Air temperature (T<sub>air</sub>)
- Wind speed, (u<sub>r</sub>)
- Heat transfer coefficient, (C<sub>H</sub>)
- Specific heat capacity of air, (c<sub>p</sub>)
- Density of air, &rho;

## Technical Algorithm
The outline below details the calculation structure for this package. The lowest level of each branch is either an input (by the user or from an external source) or a variable that is being solved for iteratively. 
 
``` 
q_H (iterative)
├── T_lst
│   ├── GOES-16 L2+ LSTC
├── T_air
│   ├── GOES-16 L2+ LSTC
│   ├── Numerical model
├── C_H
│   ├── z_0m
│   │   ├── NLCD
│   ├── z_0t
│   │   ├── z_0m (see above for inputs)
│   │   ├── u_star
│   │   │   ├── C_D
│   │   │   │   ├── z
│   │   │   │   │   ├── z_r
│   │   │   │   │   ├── d_0
│   │   │   │   │   │   ├── h_0
│   │   │   │   ├── z_0m
│   │   │   │   │   ├── NLCD
│   │   │   │   ├── zeta
│   │   │   │   │   ├── z
│   │   │   │   │   │   ├── z_r
│   │   │   │   │   │   ├── d_0
│   │   │   │   │   ├── L (iterative)
│   │   │   │   ├── psi_m
│   │   │   │   │   ├── x
│   │   │   │   │   │   ├── zeta (see above for inputs)
│   │   │   ├── u_r
│   │   │   │   ├── ASOS
│   │   ├── mu
│   │   ├── rho
│   ├── L (iterative)
│   │   ├── rho
│   │   ├── c_p
│   │   ├── u_star (see above for inputs)
│   │   ├── T_lst (see above for inputs)
│   │   ├── T_air (see above for inputs)
│   │   ├── q_H (iterative)
├── c_p
├── rho
``` 
- **T**<sub>**LST**</sub> is taken from the GOES-16 Level 2 Land Surface Temperature (LST) product for the continental United States (CONUS).
- **T**<sub>**air**</sub> is taken from an air temperature model based on the GOES-16 LST product.
- **C**<sub>**H**</sub> is calculated as a function of multiple parameters, as listed above.
    - **z**<sub>**0m**</sub> (momentum roughness height) is a function of the National Land Cover Database (NLCD).
    - **z**<sub>**0t**</sub> (thermal roughness height) is a function of the National Land Cover Database (NLCD).
    - **u_star** is a scaling parameter for wind speed.
    - **C**<sub>**D**</sub> is the momentum drag coefficient.
- **zeta** (&zeta;) is the stability parameter, whose value dictates the convective state of the boundary layer. Where `zeta < 0`, the boundary layer is typically unstable, and where `zeta > 0`, it is typically stable.
    - **z** is the height from the reference point, **z**<sub>**r**</sub>, to the zero-plane displacement height, **d**<sub>**0**</sub>, which itself is a function of the roughness element height, **h**<sub>**0**</sub>.
    - **L** is the Obukhov length.
    
### Sensible Heat Flux Algorithm
The sensible heat flux calculation iteratively solves for `L` to calculate `q_H`. 

