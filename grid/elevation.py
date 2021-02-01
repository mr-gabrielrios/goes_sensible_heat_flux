'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        Elevation
Package file path:  ~/grid/elevation.py
Objective:          Generate gridded elevation data for selected spatial domain.
Author:             Gabriel Rios
Note:               Credit to Josh Hrisko (@maker-portal) for developing large parts of this code.
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

import numpy as np, os, pyproj, elevation, rasterio, time

##############################################################################################
# END IMPORTS
##############################################################################################

def elevation_generator(lats, lons, domain):
    
    t = time.time()
    
    # Define custom file name segment with bounding coordinates (NW and SE corners)
    domain_name = str(domain[1]) + 'N_' + str(domain[2]) + 'W_' + str(domain[0]) + 'N_' + str(domain[3]) + 'W'
        
    # Define names to dependent directories
    # (one for rasterio .tifs, the other for corresponding CSV data)
    tif_dir, csv_dir = 'elevation_tifs/', 'elevation_csvs/'
    tif_path = os.path.join(os.path.dirname(__file__), tif_dir)
    csv_path = os.path.join(os.path.dirname(__file__), csv_dir)    
            
    # Create the file name for the output .csv file
    csv_fname =  domain_name + '_elevation.csv'
    
    # Create new elevation .tif and .csv if the file doesn't already exist for this domain
    if not os.path.isfile(os.path.join(csv_path, csv_fname)):
                
        # If directories don't exist, create them
        if not os.path.isdir(tif_path):
            os.mkdir(tif_path)
        if not os.path.isdir(csv_path):
            os.mkdir(csv_path)
        
        # Customized domain list to match elevation package convention
        elev_domain = [domain[2], domain[0], domain[3], domain[1]] 
        # Retrieve elevation data from SRTM and return a .tif of the spatial domain.
        # If SRTM3 (90m resolution) doesn't work, try SRTM1 (30m resolution)
        try:
            print('\t Using SRTM3...')
            elevation.clip(bounds=elev_domain, output=(tif_path + domain_name + '.tif'), product='SRTM3')
        except:
            print('\t Using SRTM1, will take longer than SRTM3...')
            elevation.clip(bounds=elev_domain, output=(tif_path + domain_name + '.tif'), product='SRTM1')
        elevation.clean()
        
        # Raster image processing
        dem_raster = rasterio.open(tif_path + domain_name + '.tif')     
        src_crs = dem_raster.crs # Get projection info
        utm = pyproj.Proj(src_crs) # Pass CRS of image from rasterio
        lonlat = pyproj.Proj(init='epsg:4326') 
        src_height, src_width = dem_raster.shape
        T0 = rasterio.transform.from_bounds(domain[2], domain[0], domain[3], domain[1], src_width, src_height)
        x_ll, y_ll = T0*((pyproj.transform(lonlat, utm, domain[0], domain[1])))
        x_ul, y_ul = T0*((pyproj.transform(lonlat, utm, domain[0], domain[3])))
        x_lr, y_lr = T0*((pyproj.transform(lonlat, utm, domain[2], domain[1])))
        x_ur, y_ur = T0*((pyproj.transform(lonlat, utm, domain[2], domain[3])))
        cols, rows = np.meshgrid(np.arange(src_width), np.arange(src_height))
        eastings, northings = T0*(cols, rows)
        cols, rows = [], []
        # Transformation from raw projection to WGS84 lat/lon values
        xvals, yvals = (eastings, northings)
        
        # Retrieve actual elevation data and corresponding coordinates from DEM
        elev = dem_raster.read(1)
        lats_elev, lons_elev = (pyproj.transform(utm, lonlat, xvals, yvals))
        
        # Ravel each dataset for processing
        elev = elev.ravel()
        lats_elev = lats_elev.ravel()
        lons_elev = lons_elev.ravel()
        
        # Find the elevation value closes to each GOES-16 pixel
        goes_elev = [elev[np.argmin(np.abs(np.subtract(ii, lons_elev)) + np.abs(np.subtract(jj, lats_elev)))] 
                     for ii,jj in zip(lons.ravel(), lats.ravel())]
        # Reshape data to match GOES-16 grid
        goes_elev = np.reshape(goes_elev, np.shape(lats))
        # Remove the .tif file
        os.remove(tif_path + domain_name + '.tif') 
        # Save the .csv to the corresponding directory
        np.savetxt(os.path.join(csv_path, csv_fname), goes_elev, delimiter=",")
    
        return goes_elev
    
    runtime = time.time() - t
