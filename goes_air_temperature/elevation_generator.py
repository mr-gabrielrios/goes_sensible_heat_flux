### Import statements
import glob,os,pyproj,elevation,rasterio,time,sys
import numpy as np
from netCDF4 import Dataset
import pandas as pd

### Override to delete all elevation data, if necessary. CAUTION: Will delete all .csv files in directory.

delete_files = False
if delete_files:
    files = glob.glob('elevation_data/*')
    for f in files:
        os.remove(f)
if delete_files:
    sys.exit() # Forced end to script to ensure

def elevation_generator(site_lon, site_lat, res):
    bbox = np.array([site_lon - res, site_lat - res, site_lon + res, site_lat + res])
    ########################################################################
    ### Build latitude and longitude grid for GOES-16
    ########################################################################
    
    t = time.time()
    
    goes_grid_file = os.getcwd() + '/aux/GOESR_ABI_CONUS_East.nc'
    grid_dataset = Dataset(goes_grid_file)
    lats,lons = grid_dataset.variables['Latitude'][:].data,grid_dataset.variables['Longitude'][:].data
    
    
    # city bounds indices based on lat/lon bbox
    nc_indx_corners = np.array([np.unravel_index(np.argmin(np.abs(np.subtract(lons,
                                                                    bbox[0]))+np.abs(np.subtract(lats,bbox[1]))),np.shape(lats)),
                     np.unravel_index(np.argmin(np.abs(np.subtract(lons,
                                                                   bbox[0]))+np.abs(np.subtract(lats,bbox[3]))),np.shape(lats)),
                     np.unravel_index(np.argmin(np.abs(np.subtract(lons,
                                                                   bbox[2]))+np.abs(np.subtract(lats,bbox[1]))),np.shape(lats)),
                     np.unravel_index(np.argmin(np.abs(np.subtract(lons,
                                                                   bbox[2]))+np.abs(np.subtract(lats,bbox[3]))),np.shape(lats))])
    # clip indices for lat/lon of city and GOES-16 data
    nc_indx_bounds = [np.min([nc_indx_corners[0][0],nc_indx_corners[1][0],nc_indx_corners[2][0],nc_indx_corners[3][0]]),
                     np.max([nc_indx_corners[0][0],nc_indx_corners[1][0],nc_indx_corners[2][0],nc_indx_corners[3][0]]),
                     np.min([nc_indx_corners[0][1],nc_indx_corners[1][1],nc_indx_corners[2][1],nc_indx_corners[3][1]]),
                     np.max([nc_indx_corners[0][1],nc_indx_corners[1][1],nc_indx_corners[2][1],nc_indx_corners[3][1]])]
    
    # Actual GOES-16 coordinate grid
    lon_plot = lons[nc_indx_bounds[0]:nc_indx_bounds[1],nc_indx_bounds[2]:nc_indx_bounds[3]]
    lat_plot = lats[nc_indx_bounds[0]:nc_indx_bounds[1],nc_indx_bounds[2]:nc_indx_bounds[3]]
    
    print('--------------------------------------------------')
    print("Lat/lon data pull and clip elapsed time: %.2f s" % (time.time() - t))
    print('--------------------------------------------------')
    
    
    ########################################################################
    # Get elevation data
    # From STM DEM
    # From elevation api STRM3 = 90m resolution elevation DEM
    ########################################################################
    
    t = time.time()
    
    elevation.clip(bounds=bbox, output=(os.getcwd().replace(' ','\ '))+'/elevation_tifs/'+site_name.lower().replace(' ','_')+'.tif',product='SRTM1')
    elevation.clean()
    dem_raster = rasterio.open('./elevation_tifs/'+site_name.lower().replace(' ','_')+'.tif') 
    
    src_crs = dem_raster.crs # get projection info
    utm = pyproj.Proj(src_crs) # Pass CRS of image from rasterio
    lonlat = pyproj.Proj(init='epsg:4326') 
    src_shape = src_height, src_width = dem_raster.shape
    T0 = rasterio.transform.from_bounds(bbox[0], bbox[1], bbox[2], bbox[3], src_width, src_height)
    x_ll,y_ll = T0*((pyproj.transform(lonlat,utm,bbox[0],bbox[1])))
    x_ul,y_ul = T0*((pyproj.transform(lonlat,utm,bbox[0],bbox[3])))
    x_lr,y_lr = T0*((pyproj.transform(lonlat,utm,bbox[2],bbox[1])))
    x_ur,y_ur = T0*((pyproj.transform(lonlat,utm,bbox[2],bbox[3])))
    x_span = [x_ll,x_ul,x_lr,x_ur]
    y_span = [y_ll,y_ul,y_lr,y_ur]
    cols, rows = np.meshgrid(np.arange(src_width),np.arange(src_height))
    eastings,northings = T0*(cols,rows)
    cols,rows = [],[]
    # transformation from raw projection to WGS84 lat/lon values
    xvals,yvals = (eastings,northings)
    
    elev = dem_raster.read(1) # actual elevation data from DEM
    lats_elev,lons_elev = (pyproj.transform(utm,lonlat,xvals,yvals)) # actual lats/lons for elevation data
    
    elev = elev.ravel() # ravel the elevation for processing
    lats_elev = lats_elev.ravel() # ravel for processing
    lons_elev = lons_elev.ravel() # ravel for processing
    
    # this is somewhat of a computationally expensive process, but it finds the elevation nearest to each GOES-16 pixel
    goes_elev = [elev[np.argmin(np.abs(np.subtract(ii,lons_elev))+np.abs(np.subtract(jj,\
                                      lats_elev)))] for ii,jj in zip(lon_plot.ravel(),lat_plot.ravel())]
    goes_elev = np.reshape(goes_elev,np.shape(lat_plot)) # reshape to match GOES-16 grid
    
    os.remove('./elevation_tifs/'+site_name.lower().replace(' ','_')+'.tif') # remove the .tif file
    
    csv_fname = site_name + '_' + str(site_lat) + '_' + str(site_lon) + '_elevation.csv'
    np.savetxt(os.getcwd() + '/elevation_data/' + csv_fname, goes_elev, delimiter=",")
    
    print('--------------------------------------------------')
    print("Elevation data access time: %.2f s" % (time.time() - t))
    print('--------------------------------------------------')

########################################################################
### Define location(s) of interest
########################################################################

site_fpath = '../ameriflux_data/' + 'GOES-16 CONUS Ameriflux Network Locations - Scatter Map Data.csv'
data = pd.read_csv(site_fpath)

extent = [-125, -66, 20, 50] # CONUS limiting coordinates

bbox = np.zeros(shape=(data.shape[0], 4))
res = 0.25
    
for i, site_id in enumerate(data['Site Id'].tolist()):
    site_name = site_id
    site_lat, site_lon = [data['Latitude (degrees)'].tolist()[i], data['Longitude (degrees)'].tolist()[i]]
    
    if site_name == 'US-ONA' and (extent[0] <= site_lon <= extent[1] or extent[2] <= site_lat <= extent[3]):
        elevation_generator(site_lon, site_lat, res)
    
