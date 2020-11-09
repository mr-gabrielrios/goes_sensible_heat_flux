# GOES Sensible Heat Flux Algorithm

The objective of this script is to provide an algorithm for calculating the sensible heat flux (Q_h) using GOES satellite data as a base input, processed through [Josh Hrisko's](https://github.com/makerportal) product for land surface temperature (LST) and air temperature at a height of 2m from the elevation of the reference location.

## Technical notes
The algorithm in this script is loosely based on Launiainen et al. (1995) ([link to publication](https://www.sciencedirect.com/science/article/pii/026698389090021W)). Note that Kim et al. (2019) ([link to publication](https://doi.org/10.3390/atmos10070363)) used a simpler implementation of this algorithm. Convergence has been shown for both unstable and stable atmospheric conditions (`zeta < 0, zeta > 0, zeta = 0`).

### Next Steps
* Input from Hrisko product for LST and air temperature
* Robust algorithm for canopy height as a function of surrounding roughness elements
* Validation at key urban test points ([NYS Mesonet](http://www.nysmesonet.org/) stations)
  - Manhattan (CUNY Hunter College)
  - Bronx (CUNY Lehman College)
  - Brooklyn (CUNY Brooklyn College)
  - Queens (CUNY Queens College)
  - Staten Island (CUNY Staten Island)

### Known Issues
* No connection to dynamic inputs (e.g. MesoNET heat flux stations for validation, LST/air temperatures from Josh Hrisko's models, geographic data)
* Rudimentary momentum roughness length (`z_0m`) estimation based on buildings in Manhattan
* Implementation not user-friendly (at all)
* No bounds on parameter values, allowing for this script to go haywire

### Dependencies
* `math`
* `time`

Feel free to point anything out to me that you feel can be improved upon! I'm always looking to improve my programming form, adopt best practices, or try new algorithms. Thanks!