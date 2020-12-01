### Import statements
import os,datetime,timezonefinder,pytz,pyproj,csv,re,elevation,rasterio,requests,time
from google.cloud import storage
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.img_tiles as tiles
import matplotlib.colors as colors
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import numpy as np
from netCDF4 import Dataset
import warnings
import pandas as pd

# Suppress warnings
warnings.filterwarnings("ignore")

##################################################################################################################
### Begin program
##################################################################################################################

cwd = os.getcwd() + '/goes_air_temperature/'
aux_dir = '/aux/'
# Boolean to control strings with performance data printing to the console. If True, strings should print.  
str_switch = False
# Boolean to control plotting of temperature product. If True, plots should print.  
plot_switch = False

######################################################################################################
### Google Cloud Storage Credentials
# Objective: Connect to Google Cloud Storage

t = time.time()
os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "GOOGLE_AUTH_CREDS_CUERG.json" # local google auth credentials
client = storage.Client() # start the storage client
bucket_16 = client.get_bucket('gcp-public-data-goes-16') # call the GOES-16 storage bucket
#
# The Google Storage repository is formatted as follows:
# Example file:
# https://console.cloud.google.com/storage/browser/_details/gcp-public-data-goes-16/ABI-L2-LSTC/2020/
# 111/21/OR_ABI-L2-LSTC-M6_G16_s20201112101130_e20201112103503_c20201112104424.nc
# directory conventions: 
# data repository/product-type/year/year-day/hour/filename.nc
#

if str_switch:
    print("Google Cloud Storage data access elapsed time: %.2f s" % (time.time() - t))

# Bounding box definition for New York - use only for debugging
bbox = [-74.50564249616093, 40.24587322565815, -73.45027336215357, 41.165290318055426]

## ################## END Google Cloud Storage Credentials #########################

######################################################################################################
### Basemapper Plot (to be replaced by plotter.py)
# Objective: Grab LST data from GOES-16 ABI L2+ on Google Cloud

def basemapper():
    
    t = time.time()

    layout_map = tiles.QuadtreeTiles()
    fig = plt.figure(figsize=(12,9))
    # Create a GeoAxes in the tile's projection.
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    # Limit the extent of the map to a small longitude/latitude range.
    ax.set_extent([bbox[0],bbox[2],bbox[1],bbox[3]], crs=ccrs.PlateCarree())
    # Add the Stamen data at zoom level 8.
    ax.add_image(layout_map, 10)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=1, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    
    if str_switch:
        print("Basemapper plotter elapsed time: %.2f s" % (time.time() - t))

    return fig, ax

## ################## END Basemapper #########################

######################################################################################################
### GOES-16 Grid Build
t = time.time()

goes_grid_file = cwd + aux_dir + 'GOESR_ABI_CONUS_East.nc'
grid_dataset = Dataset(goes_grid_file)
lats,lons = grid_dataset.variables['Latitude'][:].data,grid_dataset.variables['Longitude'][:].data

# Boundary grid indices based on lat/lon bounding box
nc_indx_corners = np.array([np.unravel_index(np.argmin(np.abs(np.subtract(lons, bbox[0])) + np.abs(np.subtract(lats,bbox[1]))), np.shape(lats)), \
                            np.unravel_index(np.argmin(np.abs(np.subtract(lons,bbox[0]))+np.abs(np.subtract(lats,bbox[3]))),np.shape(lats)), \
                            np.unravel_index(np.argmin(np.abs(np.subtract(lons,bbox[2]))+np.abs(np.subtract(lats,bbox[1]))),np.shape(lats)), \
                            np.unravel_index(np.argmin(np.abs(np.subtract(lons,bbox[2]))+np.abs(np.subtract(lats,bbox[3]))),np.shape(lats))])
# clip indices for lat/lon of city and GOES-16 data
nc_indx_bounds = [np.min([nc_indx_corners[0][0],nc_indx_corners[1][0],nc_indx_corners[2][0],nc_indx_corners[3][0]]),
                 np.max([nc_indx_corners[0][0],nc_indx_corners[1][0],nc_indx_corners[2][0],nc_indx_corners[3][0]]),
                 np.min([nc_indx_corners[0][1],nc_indx_corners[1][1],nc_indx_corners[2][1],nc_indx_corners[3][1]]),
                 np.max([nc_indx_corners[0][1],nc_indx_corners[1][1],nc_indx_corners[2][1],nc_indx_corners[3][1]])]

########## These are the real GOES-16 lat/lon coords that correspond with the city boundaries ############
lon_plot = lons[nc_indx_bounds[0]:nc_indx_bounds[1],nc_indx_bounds[2]:nc_indx_bounds[3]]
lat_plot = lats[nc_indx_bounds[0]:nc_indx_bounds[1],nc_indx_bounds[2]:nc_indx_bounds[3]]

if str_switch:
    print("Lat/lon data pull and clip elapsed time: %.2f s" % (time.time() - t))      
        
## ################## END GOES-16 Grid Build #########################

######################################################################################################
### NLCD Grid Setup
# Objective: Read through upscaled NLCD CSV, filter out elements corresponding to bound box, and average
# Averaging methodology:    pull WRF reference heights, match to corresponding NLCD class, 
#                           get dot product of each pixel, then get average

t = time.time()

nlcd_goes = []

nlcd_data_fp = cwd + aux_dir + 'nlcd_to_goes_upscaled.csv' # Define file path for NLCD distribution CSV file
nlcd_zr_fp = cwd + aux_dir + 'nlcd_reference_heights.csv' # Define file path for NLCD reference heights CSV file

nlcd_rows = np.linspace(nc_indx_bounds[0], nc_indx_bounds[1]-1, (nc_indx_bounds[1]-nc_indx_bounds[0]))
nlcd_rows = [int(i) for i in nlcd_rows]
nlcd_cols = np.linspace(nc_indx_bounds[2], nc_indx_bounds[3]-1, (nc_indx_bounds[3]-nc_indx_bounds[2]))
nlcd_cols = [int(i) for i in nlcd_cols]

nlcd_goes = pd.read_csv(nlcd_data_fp, usecols=nlcd_cols)
nlcd_goes = nlcd_goes[nlcd_goes.index.isin(nlcd_rows)]

nlcd_goes = np.array(nlcd_goes)

nlcd_goes_clip = nlcd_goes[nc_indx_bounds[0]:nc_indx_bounds[1],
                            nc_indx_bounds[2]:nc_indx_bounds[3]] 
nlcd_upscaled = np.array(np.zeros((np.shape(nlcd_goes)[0],np.shape(nlcd_goes)[1],20)))
for nlcd_ii in range(0,np.shape(nlcd_goes)[0]):
    for nlcd_jj in range(0,np.shape(nlcd_goes)[1]):
        if nlcd_goes[nlcd_ii][nlcd_jj]=='nan':
            nlcd_ii_jj = np.repeat(np.nan,20)
        else:
            nlcd_ii_jj = [float(qq) for qq in re.sub(' +',',',\
                              str(nlcd_goes[nlcd_ii][nlcd_jj]).replace('[',\
                          '').replace(']','').replace('\n','')).split(',') if qq!='']

        # if len(nlcd_ii_jj)!=20:
        #     print('ISSUE WITH NLCD')
        nlcd_upscaled[nlcd_ii][nlcd_jj] = nlcd_ii_jj # match NLCD data to GOES-16 grid shape

# Filter out water locations
water_locs_ii,water_locs_jj = [],[]
for water_ii in range(0,np.shape(nlcd_goes)[0]):
    for water_jj in range(0,np.shape(nlcd_goes)[1]):
        if nlcd_upscaled[water_ii][water_jj][0]>0.5 or np.isnan(nlcd_upscaled[water_ii][water_jj][0]):
            water_locs_ii.append(water_ii)
            water_locs_jj.append(water_jj)

# Use DataFrame to read in the NLCD WRF reference heights
nlcd_hts = pd.read_csv(nlcd_zr_fp, usecols=['ZR']).to_numpy().ravel()
h_0 = np.zeros([nlcd_upscaled.shape[0], nlcd_upscaled.shape[1]])
for i in range(0, h_0.shape[0]):
    for j in range(0, h_0.shape[1]):
        h_0[i, j] = np.dot(nlcd_upscaled[i, j], nlcd_hts)
        
# Plot to show roughness heights for selected plot area, based on upscaled NLCD values
fig, ax = plt.subplots()
im = ax.pcolormesh(nlcd_cols, nlcd_rows, h_0)
im_cbar = fig.colorbar(im, ax = ax)
im_cbar.ax.get_yaxis().labelpad = 15
im_cbar.ax.set_ylabel('h_0 [m]', rotation=270)
plt.gca().invert_yaxis()
plt.gca().set_aspect('equal', adjustable='box')
plt.title('Roof Level based on Upscaled NLCD Values')
plt.show()

if str_switch:
    print("NLCD upscaled value access/plot time: %.2f s" % (time.time() - t))
    
## ################## END NLCD Grid Setup #########################

######################################################################################################
### Elevation Grid Setup
# Objective: Read through upscaled elevation CSV
# Progress notes: create upscaled grid for CONUS. Only NYC grid exists to support short-term efforts

t = time.time()

# IMPORTANT:    Only 2 degree square centered around Manhattan is saved to elevation.csv. For non-New York
#               locales, elevation.csv must be generated for them            
dem_raster = rasterio.open(cwd + 'elevation_tifs/elevation.tif')

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

if str_switch:
    print("Elevation data access time: %.2f s" % (time.time() - t))
    
######################################################################################################
### GOES-16 LST Data Access and Plotting

def temperature_data(query_date, lat_pixel, lon_pixel, site):
    query_date = datetime.datetime.strptime(str(query_date), '%Y%m%d%H%M').strftime('%Y%j%H%M%S%f')
    
    # Start program timer
    start_time = time.time()
    ######################################################################################################
    ### GOES-16 LST Data Access and Plotting
    # Objective: Grab LST data from GOES-16 ABI L2+ on Google Cloud
    t = time.time()
    
    goes_folder_local = './goes_data/' # where the temporary goes file will be stored
    
    tf = timezonefinder.TimezoneFinder() # we need UTC and local time - this will get us both
    # Get the tz-database-style time zone name (e.g. 'America/Vancouver') or None
    timezone_str = tf.certain_timezone_at(lat=bbox[1]+(bbox[3]-bbox[1])/2.0, lng=bbox[0]+(bbox[2]-bbox[0])/2.0) # get time zone for city center
    
    if timezone_str is None:
        print("Could not determine the time zone")
    else:
        # Display the current time in that time zone
        timezone = pytz.timezone(timezone_str)
        dt = datetime.datetime.strptime(query_date, '%Y%j%H%M%S%f')
        print ("The queried time in %s is %s" % (timezone_str, dt + timezone.utcoffset(dt)))
        
    prod = 'LST' # Level L-2 Product identifier
    product_ID = 'ABI-L2-'+prod+'C/' # ABI, L2 product, CONUS
    
    local_t = dt + timezone.utcoffset(dt) # local time with UTC offset (need this for Air temp algorithm)
    t_search = datetime.datetime.strptime(query_date,'%Y%j%H%M%S%f') # UTC time for searching for GOES-16 data files - this can be changed to anything
    
    goes_16_date_vec = []
    goes_16_file_vec = []
    while goes_16_file_vec==[]:
        yr = '{:0004d}'.format(t_search.year) 
        yday = '{:003d}'.format(t_search.timetuple().tm_yday)
        thour = '{:02d}'.format(t_search.hour)
    
        prefix_16 = product_ID+yr+'/'+yday+'/'+thour # this is the relative folder where the GOES-16 data will be
        blobs_16 = bucket_16.list_blobs(prefix=prefix_16) # get all files in prefix
    
        for kk_16,blob_16 in enumerate(blobs_16): # loop through files
    
            start_time_16 = blob_16.name.split('_')[3][1:]
            start_datetime_16 = datetime.datetime.strptime(start_time_16,'%Y%j%H%M%S%f')
    
            if len(blob_16.name.split('_'))<=5:
                end_time_16 = blob_16.name.split('_')[4][1:-3]
                end_datetime_16 = datetime.datetime.strptime(end_time_16,'%Y%j%H%M%S%f')
            else:
                end_time_16 = blob_16.name.split('_')[4][1:]
                end_datetime_16 = datetime.datetime.strptime(end_time_16,'%Y%j%H%M%S%f')
    
            goes_16_date_vec.append(end_datetime_16)
            goes_16_file_vec.append(blob_16)
        
        if goes_16_file_vec==[]:
            t_search = t_search - datetime.timedelta(1/24) # if there are no files at that time, look an hour behind
        
    curr_filename = (goes_16_file_vec[0].name).split('/')[-1] # the the last file in search
    if os.path.isfile(goes_16_file_vec[0].name+curr_filename):
        pass
    else:
        # download the file to the local GOES-16 repository
        (goes_16_file_vec[0]).download_to_filename(cwd+goes_folder_local+curr_filename)
    
    try:
        prod_set = Dataset(cwd+goes_folder_local+curr_filename) # load LST data
        os.remove(cwd+goes_folder_local+curr_filename)
    except:
        # if there was an error, remove the faulty file and redownload
        os.remove(cwd+goes_folder_local+curr_filename) 
        (goes_16_file_vec[0]).download_to_filename(cwd+goes_folder_local+curr_filename)
        prod_set = Dataset(cwd+goes_folder_local+curr_filename)
        os.remove(cwd+goes_folder_local+curr_filename)
    
    prod_names = []
    for prod_name_ii in prod_set.variables:
        prod_names.append(prod_name_ii) # get all variable names in GOES-16 file
    
    if str_switch:
        print("GOES-16 LST data access time: %.2f s" % (time.time() - t))
    
    ## ################## END GOES-16 LST Data Access and Plotting #########################
    
    ######################################################################################################
    ### GOES-16 LST Data Plotting
    # Objective: Plot LST data from GOES-16 ABI L2+ on Google Cloud
    
    t = time.time()
    
    ######################## COMMENTED CODE
    lon_adjust = 0.5*np.diff(lons[nc_indx_bounds[0]:nc_indx_bounds[1],
                              nc_indx_bounds[2]-1:nc_indx_bounds[3]],axis=1) # this shifts lons just for pcolormesh (lower-left corner plot shift)
    lat_adjust = 0*np.diff(lats[nc_indx_bounds[0]:nc_indx_bounds[1]+1,
                              nc_indx_bounds[2]:nc_indx_bounds[3]],axis=0) # this shifts lats just for pcolormesh (lower-left corner plot shift)
    ######################## END COMMENTED CODE
    
    # Get the actual LST data:
    prod_vals = ((prod_set.variables[prod_names[0]])[nc_indx_bounds[0]:nc_indx_bounds[1],nc_indx_bounds[2]:nc_indx_bounds[3]]).data
    # get the data quality flag data
    dqf = ((prod_set.variables[prod_names[1]])[nc_indx_bounds[0]:nc_indx_bounds[1],nc_indx_bounds[2]:nc_indx_bounds[3]]).data
    
    # filtering LST data:
    prod_vals[dqf!=0] = np.nan # any non-clear/valid pixels are nans now
    prod_vals[water_locs_ii,water_locs_jj] = np.nan # any water pixels not captured by the GOES-16 algorithm are nulled as well
    
    # LST Plot
    if plot_switch:
        fig,ax = basemapper() # plot the basemap for visualization
        norm = colors.BoundaryNorm(boundaries=np.linspace(280, 310, 10), ncolors=256)
        p1 = ax.pcolormesh(lon_plot-lon_adjust,lat_plot-lat_adjust,prod_vals,zorder=999,alpha=0.6,transform=ccrs.PlateCarree(), norm=norm) # plotting the product data
        
        cb = fig.colorbar(p1) # colorbar
        cb.set_label(r'%s' % (prod_names[0])) # LST label
        plt.title('{0} from {1} {2} - {3} UTC'.format(prod_set.variables[prod_names[0]].long_name,
                                                  prod_set.time_coverage_start.split('T')[0],
                                                  prod_set.time_coverage_start.split('T')[1].replace('Z',''),
                                              prod_set.time_coverage_end.split('T')[1].replace('Z',''))) # title with date/tie
        
        sat_time_end = datetime.datetime.strptime(prod_set.time_coverage_end,'%Y-%m-%dT%H:%M:%S.%fZ') # end of satellite scan
    
    if str_switch:
        print("GOES-16 LST data plot time: %.2f s" % (time.time() - t))
    
    ## ################## END GOES-16 LST Data Plotting #########################
    
    ######################################################################################################
    ### Gaussian Fit Coefficient Definition
    # Objective: Calculate Gaussian fit coefficients based on generated data
    
    t = time.time()
    
    col_array,fit_params,gaus_header = [],[],[]
    with open(cwd + aux_dir + 'gaus_fit_params.csv','r') as fit_params_csv:
        reader = csv.reader(fit_params_csv,delimiter=',')
        for row in reader:
            if gaus_header==[]:
                gaus_header = row
                continue
            col_array.append(int(row[0])) # this variable identifies the variable column to multiply
            fit_params.append(row[1:]) # the rest of the fit coefficients
    
    # here, we loop through and collect variables for the gaussian inputs
    input_array = []               
    for ii in range(0,np.shape(lon_plot)[0]):
        for jj in range(0,np.shape(lon_plot)[1]):
            if goes_elev[ii,jj]>-1.0:
                input_ii_jj = [0.0,lat_plot[ii,jj],lon_plot[ii,jj],goes_elev[ii,jj]]
            else:
                input_ii_jj = [0.0,lat_plot[ii,jj],lon_plot[ii,jj],np.nan]
            input_ii_jj.extend(nlcd_upscaled[ii,jj])
        
            input_array.append(input_ii_jj) # this contains the 24 input variables (null, lat, lon, elev, 20 NLCD)
            
    # This computes the Gaussian variables for each pixel in the city
    A_0,x_0,sigma,y_0 = [],[],[],[]
    for goes_ii in input_array:
        A_ii,x_ii,sigma_ii,y_ii = 0.0,0.0,0.0,0.0
        for aa in range(0,len(col_array)):
            if fit_params[aa][0]!='':
                col_vals = [float(qq) for qq in fit_params[aa][0].replace('[','').replace(']','').split(',')]
                A_ii+=(goes_ii[col_array[aa]]*col_vals[1])+col_vals[2]
            if fit_params[aa][1]!='':
                col_vals = [float(qq) for qq in fit_params[aa][1].replace('[','').replace(']','').split(',')]
                x_ii+=(goes_ii[col_array[aa]]*col_vals[1])+col_vals[2]
            if fit_params[aa][2]!='':
                col_vals = [float(qq) for qq in fit_params[aa][2].replace('[','').replace(']','').split(',')]
                sigma_ii+=(goes_ii[col_array[aa]]*col_vals[1])+col_vals[2]
            if fit_params[aa][3]!='':
                col_vals = [float(qq) for qq in fit_params[aa][3].replace('[','').replace(']','').split(',')]
                y_ii+=(goes_ii[col_array[aa]]*col_vals[1])+col_vals[2]
        A_0.append(A_ii)
        x_0.append(x_ii)
        sigma.append(sigma_ii)
        y_0.append(y_ii)
        
    if str_switch:
        print("Gaussian fit coefficient generation runtime: %.2f s" % (time.time() - t))
    
    ## ################## END Gaussian Fit Coefficient Definition #########################
    
    ######################################################################################################
    ### Air Temperature Calculation
    # Objective: Calculate air temperature
    
    t = time.time()
    
    lst = prod_vals.ravel() # ravel the LST data for processing
    
    # THE GAUSSIAN CALCULATION FUNCTION
    def gaus_func(LST_func,A_0_func,x_0_func,sigma_func,y_0_func,local_t_func):
        return LST_func + y_0_func - (A_0_func*np.exp(-(np.power(np.subtract(local_t_func.hour,x_0_func),2.0))/(2.0*sigma_func**2.0)))
    
    
    T_air = []
    for val_iter in range(0,len(input_array)):
        T_air_ii = gaus_func(lst[val_iter],A_0[val_iter],x_0[val_iter],sigma[val_iter],y_0[val_iter],local_t)
        if T_air_ii>330.0 or T_air_ii<180.0:
            T_air.append(np.nan) # if outside reasonable range, nan
        else: 
            T_air.append(T_air_ii)
    
    T_air = np.reshape(T_air,np.shape(lon_plot))
    pixel_res = 3
    
    area_lst = np.nanmean(prod_vals[(lat_pixel-pixel_res):(lat_pixel+pixel_res+1), (lon_pixel-pixel_res):(lon_pixel+pixel_res+1)])
    area_T_air = np.nanmean(T_air[(lat_pixel-pixel_res):(lat_pixel+pixel_res+1), (lon_pixel-pixel_res):(lon_pixel+pixel_res+1)])
        
    # Air Temperature Plot
    if plot_switch:
        fig2,ax2 = basemapper() # plot the air temp basemap
        norm = colors.BoundaryNorm(boundaries=np.linspace(270, 320, 10), ncolors=256)
        
        p2 = ax2.pcolormesh(lon_plot-lon_adjust,lat_plot-lat_adjust,T_air,zorder=999,alpha=0.6,transform=ccrs.PlateCarree(), norm=norm) # plotting the air temp data
        
        cb = fig2.colorbar(p2)
        cb.set_label('$T_{air}$ [K]')
        plt.title('Air Temperature from {1} {2} - {3} UTC'.format(prod_set.variables[prod_names[0]].long_name,
                                                  prod_set.time_coverage_start.split('T')[0],
                                                  prod_set.time_coverage_start.split('T')[1].replace('Z',''),
                                              prod_set.time_coverage_end.split('T')[1].replace('Z','')))
        
        plt.show()
    
    if str_switch:
        print("Air temperature calculation runtime: %.2f s" % (time.time() - t))
    
    print('--------------------------------------------------')    
    print('Total GOES_16_LST_to_Tair runtime: %.3f s \n' % (time.time() - start_time))
    
    return area_lst, area_T_air, h_0[lat_pixel, lon_pixel]