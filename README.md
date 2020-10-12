# GOES Sensible Heat Flux Algorithm

The objective of this package is to provide an algorithm for calculating the sensible heat flux using GOES satellite data as a base input, processed through [Josh Hrisko's](https://github.com/makerportal) models for land surface temperature (LST) and air temperature at a height of 2m from the elevation of the reference location.

## 'kim_method' Branch
This branch holds a rudimentary algorithm based on the set of equations proposed by Kim et al. (2019) [see publication here](https://doi.org/10.3390/atmos10070363). The failure with this algorithm is the lack of a proper aerodynamic roughness length (z_0m) calculation, leaving it to be iteratively solved. Although the publication provides a range of roughness length values for different test sites, the origin of these values was never defined. 

This branch is obsolete and is being stored here for future reference.