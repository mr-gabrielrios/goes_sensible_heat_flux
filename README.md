# GOES Sensible Heat Flux Algorithm

The objective of this script is to provide an algorithm for calculating the sensible heat flux using GOES satellite data as a base input, processed through [Josh Hrisko's](https://github.com/makerportal) models for land surface temperature (LST) and air temperature at a height of 2m from the elevation of the reference location.

## Technical notes
The algorithm in this script is loosely based on Launiainen et al. (1995) ([link to publication](https://www.sciencedirect.com/science/article/pii/026698389090021W)). Note that Kim et al. (2019) ([link to publication](https://doi.org/10.3390/atmos10070363)) used a simpler implementation of this algorithm. Convergence has been shown for both unstable and stable atmospheric conditions (`zeta < 0, zeta > 0, zeta = 0`).

### Known Issues
* No connection to dynamic inputs (e.g. MesoNET heat flux stations for validation, LST/air temperatures from Josh Hrisko's models, geographic data)
* Sensitivity to certain initial conditions (C_h, C_d) to be explored
* Implementation not user-friendly (at all)
* No bounds on parameter values, allowing for this script to go haywire

Feel free to point anything out to me that you feel can be improved upon! I'm always looking to improve my programming form, adopt best practices, or try new algorithms. Thanks!