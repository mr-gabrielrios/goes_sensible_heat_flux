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
The outline below details the calculation structure for this package.
``` 
{
Q_H
├── T_lst
├── T_air
├── C_h
│   ├── z_0m
│   │   ├── NLCD
│   ├── z_0t
│   │   ├── z_0m, u_star, mu, rho
│   ├── L
│   │   ├── rho
│   │   ├── c_p
│   │   ├── u_star
│   │   ├── T_lst
│   │   ├── T_air
│   │   ├── q_h
├── c_p
├── rho
} 
``` 