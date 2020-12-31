# GOES Sensible Heat Flux Algorithm

The objective of this script is to provide an algorithm for calculating the sensible heat flux (Q_h) using GOES satellite data as a base input, processed through [Josh Hrisko's](https://github.com/makerportal) product for land surface temperature (LST) and air temperature at a height of 2m from the elevation of the reference location.

## Sample Outputs

- New York, 07/28/2019 @ 12:00 to 23:59 UTC
![Ah damn it, this didn't work, huh?'](https://github.com/mr-gabrielrios/goes_sensible_heat_flux/blob/main/plots/q_h20190728NYC.gif)

## Technical notes
The algorithm in this script is loosely based on Launiainen et al. (1995) ([link to publication](https://www.sciencedirect.com/science/article/pii/026698389090021W)). Note that Kim et al. (2019) ([link to publication](https://doi.org/10.3390/atmos10070363)) used a simpler implementation of this algorithm. Convergence has been shown for both unstable and stable atmospheric conditions (`&zeta; < 0, &zeta; > 0, &zeta; = 0`).

Zero-plane displacement `z_d` and momentum roughness height `z_0m` are calculated using the methods described in Raupach (1994) and carried out in Grimmond and Oke (1999).

Please see `program_flow.md` for a low-level description of the algorithm.

### Next Steps
* Validation at key urban test points ([NYS Mesonet](http://www.nysmesonet.org/) stations)
  - Manhattan (CUNY Hunter College)
  - Bronx (CUNY Lehman College)
  - Brooklyn (CUNY Brooklyn College)
  - Queens (CUNY Queens College)
  - Staten Island (CUNY Staten Island)
* Validation at locations outside of New York City at locations with varying land cover types, per NLCD 2016
* Spatial & temporal visualizations for sensible heat flux at user-defined locations
* Creation of error analysis suite with respect to observation networks (ASOS, NYS Mesonet, Ameriflux)

### Known Issues
* GitHub files are omitting auxiliary files including critical inputs such as elevation data, upscaled NLCD data, and supplemental observation data from the NOAA ASOS and New York State Mesonet network sites in New York City. This was done in the interest of preserving storage space.
* Implementation not user-friendly (at all)
* No bounds on parameter values, allowing for this script to go haywire

### Dependencies
* `math`
* `numpy`
* `time`
* `os`
* `datetime`
* `google.cloud`
* `matplotlib`
* `cartopy`
* `netCDF4`
* `rasterio`

Feel free to point anything out to me that you feel can be improved upon! I'm always looking to improve my programming form, adopt best practices, or try new algorithms. Thanks!