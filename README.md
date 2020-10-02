# GOES Sensible Heat Flux

The objective of this package is to provide an algorithm for calculating the sensible heat flux using GOES satellite data as a base input, processed through [Josh Hrisko's](https://github.com/makerportal) models for land surface temperature (LST) and air temperature at a height of 2m from the elevation of the reference location.

### Known Issues
* No connection to dynamic inputs (e.g. MesoNET heat flux stations for validation, LST/air temperatures from Josh Hrisko's models, geographic data)
* Convergence difficulties for specific parameters
* Implementation not user-friendly (at all)
* No bounds on parameter values, allowing for this script to go haywire

Feel free to point anything out to me that you feel can be improved upon! I'm always looking to improve my programming form, adopt best practices, or try new algorithms. Thanks!