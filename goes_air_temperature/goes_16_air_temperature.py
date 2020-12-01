import os,datetime,timezonefinder,pytz,csv,time
from google.cloud import storage
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as colors
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import warnings
import grid_processes

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
### Basemapper Plot (to be replaced by plotter.py)
# Objective: Grab LST data from GOES-16 ABI L2+ on Google Cloud

def basemapper():
    t = time.time()
    
    fig = plt.figure(figsize=(12,9))
    # Create a GeoAxes in the tile's projection.
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    # Limit the extent of the map to a small longitude/latitude range.
    ax.set_extent([bbox[0],bbox[2],bbox[1],bbox[3]], crs=ccrs.PlateCarree())
    # Add the Stamen data at zoom level 8.
    # ax.add_image(layout_map, 10)
    ax.coastlines(resolution='10m')

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    
    if str_switch:
        print('--------------------------------------------------')
        print("Basemapper plotter elapsed time: %.2f s" % (time.time() - t))
    return fig,ax

## ################## END Basemapper #########################

######################################################################################################
### Google Cloud Storage Credentials
# Objective: Connect to Google Cloud Storage

t = time.time()
os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "GOOGLE_AUTH_CREDS_CUERG.json" # local google auth credentials
client = storage.Client() # start the storage client
bucket_16 = client.get_bucket('gcp-public-data-goes-16') # call the GOES-16 storage bucket
# #
# # The Google Storage repository is formatted as follows:
# # Example file:
# # https://console.cloud.google.com/storage/browser/_details/gcp-public-data-goes-16/ABI-L2-LSTC/2020/
# # 111/21/OR_ABI-L2-LSTC-M6_G16_s20201112101130_e20201112103503_c20201112104424.nc
# # directory conventions: 
# # data repository/product-type/year/year-day/hour/filename.nc

if str_switch:
    print('--------------------------------------------------')
    print("Google Cloud Storage data access elapsed time: %.2f s" % (time.time() - t))

## ################## END Google Cloud Storage Credentials #########################

t = time.time()

# GOES-16 constants to be removed once structural integration is performed
H = 42164160
r_pol = 6356752.31414
r_eq = 6378137
lambda_0 = -1.308996939 
# 0.3 degree bound size selected to cover all 4 Mesonet sites
[central_lon, central_lat, bound_sz] = [-73.9521, 40.6305, 0.3]
[x0, x1, y1, y0] = grid_processes.grid_grab(central_lon, central_lat, bound_sz, H, r_pol, r_eq, lambda_0)
bbox = [central_lon - bound_sz,
        central_lat - bound_sz,
        central_lon + bound_sz,
        central_lat + bound_sz]

if str_switch:
    print('Coordinate box:')
    print(bbox)
    print('Corresponding indices:')
    print(x0, x1, y1, y0)

######################################################################################################
### GOES-16 Grid Build

t = time.time()
goes_grid_file = cwd + aux_dir + r'GOESR_ABI_CONUS_East.nc'
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

# These are the real GOES-16 lat/lon coords that correspond with the city boundaries
lon_plot = lons[nc_indx_bounds[0]:nc_indx_bounds[1],nc_indx_bounds[2]:nc_indx_bounds[3]]
lat_plot = lats[nc_indx_bounds[0]:nc_indx_bounds[1],nc_indx_bounds[2]:nc_indx_bounds[3]]

if str_switch:
    print(nc_indx_bounds)
    print('--------------------------------------------------')
    print("GOES-16 grid build elapsed time: %.2f s" % (time.time() - t))

######################################################################################################
### NLCD Grid Setup
# Objective: Read through upscaled NLCD CSV, filter out elements corresponding to bound box, and average
# Averaging methodology:    pull WRF reference heights, match to corresponding NLCD class, 
#                           get dot product of each pixel, then get average

t = time.time()

nlcd_data_fp = cwd + aux_dir + 'nlcd_to_goes_upscaled.csv' # Define file path for NLCD distribution CSV file
nlcd_zr_fp = cwd + aux_dir + 'nlcd_reference_heights.csv' # Define file path for NLCD reference heights CSV file

# Reset all variables
nlcd_rows, nlcd_cols, nlcd_goes, nlcd_hts, nlcds = [], [], [], [], []

# Create continuous vectors for rows and columns for indexing the CSV
nlcd_rows = np.linspace(nc_indx_bounds[0], nc_indx_bounds[1], (nc_indx_bounds[1]-nc_indx_bounds[0]+1))
nlcd_rows = [int(i) for i in nlcd_rows]
nlcd_cols = np.linspace(nc_indx_bounds[2], nc_indx_bounds[3], (nc_indx_bounds[3]-nc_indx_bounds[2]+1))
nlcd_cols = [int(i) for i in nlcd_cols]
# Use DataFrame to read in and filter the CSV to the user-selected box
nlcd_goes = pd.read_csv(nlcd_data_fp, usecols=nlcd_cols)
nlcd_goes = nlcd_goes[nlcd_goes.index.isin(nlcd_rows)]
nlcd_goes = nlcd_goes.to_numpy() # Revert to NumPy array for numerical operations on each element

# Use DataFrame to read in the NLCD WRF reference heights
nlcd_hts = pd.read_csv(nlcd_zr_fp, usecols=['ZR']).to_numpy().ravel()

nlcds = np.array(np.zeros((nlcd_goes.shape[0],nlcd_goes.shape[1],20))) # Initialize empty 2D array of same size and shape as bound box

temp = []
for i in range(0, nlcd_goes.shape[0]):
    # print('Row: %d' % i)
    for j in range(0, nlcd_goes.shape[1]):
        # Catch all the strings of lists
        if type(nlcd_goes[i, j]) is str:
            # print('Row %d, col %d' % (i, j))
            temp = nlcd_goes[i,j].split()
            temp = [k.replace('[','').replace(']','') for k in temp]
            temp = [float(k) for k in temp if k]
            # temp = np.dot(nlcd_hts, temp)/len(temp)
            # if np.isnan(temp):
            #     temp = 0
            nlcds[i, j] = temp
        # Everything else is nan, as verified by trialed conditionals to look for other data types
        else:
            nlcds[i, j] = 0

if str_switch:
    print("Grid size: ", nlcds.shape)    
    print('--------------------------------------------------')
    print("NLCD upscaled value access/plot time: %.2f s" % (time.time() - t))

h_0 = np.zeros([nlcds.shape[0], nlcds.shape[1]])
for i in range(0, h_0.shape[0]):
    for j in range(0, h_0.shape[1]):
        h_0[i, j] = np.dot(nlcds[i, j], nlcd_hts)

# Plot to show roughness heights for selected plot area, based on upscaled NLCD values
if not plot_switch:
    fig, ax = plt.subplots()
    im = ax.pcolormesh(nlcd_cols, nlcd_rows, h_0)
    im_cbar = fig.colorbar(im, ax = ax)
    im_cbar.ax.get_yaxis().labelpad = 15
    im_cbar.ax.set_ylabel('h_0 [m]', rotation=270)
    plt.gca().invert_yaxis()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title('Roof Level based on Upscaled NLCD Values')
    plt.show()

## ################## END NLCD Grid Setup #########################

######################################################################################################
### Elevation Grid Setup
# Objective: Read through upscaled elevation CSV
# Progress notes: create upscaled grid for CONUS. Only NYC grid exists to support short-term efforts

t = time.time()

# IMPORTANT:    Only 2 degree square centered around Manhattan is saved to elevation.csv. For non-New York
#               locales, elevation.csv must be generated for them              
goes_elev = np.empty(lon_plot.shape)
goes_elev_fp = cwd + aux_dir + 'elevation_nyc.csv'
with open(goes_elev_fp, 'r') as file:
    rdr = csv.reader(file, delimiter=',')    
    rdr_data = []    
    for row in rdr:
        rdr_data.append(row)  
    # Clip as necessary based on grid bounds
    rdr_data = rdr_data[nc_indx_bounds[0]:nc_indx_bounds[1]][nc_indx_bounds[2]:nc_indx_bounds[3]]
    # Assign elevation valeus to matrix
    rc = 0
    for row in rdr_data:
        goes_elev[rc] = row        
        rc += 1

if str_switch:
    print('--------------------------------------------------')
    print("Elevation data access time: %.2f s" % (time.time() - t))

## ################## END Elevation Grid Setup #########################

def temperature_data(query_date, lat_pixel, lon_pixel, site):
    query_date = datetime.datetime.strptime(str(query_date), '%Y%m%d%H%M').strftime('%Y%j%H%M%S%f')
    # print('Querying %s at %s' % (site, datetime.datetime.strptime(query_date, '%Y%j%H%M%S%f').strftime('%Y-%m-%d %H:%M:%S')))
    
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
        dt = datetime.datetime.utcnow()
        
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
        print('--------------------------------------------------')
        print("GOES-16 LST data access time: %.2f s" % (time.time() - t))
    
    ## ################## END GOES-16 LST Data Access and Plotting #########################
    
    ######################################################################################################
    ### GOES-16 LST Data Plotting
    # Objective: Plot LST data from GOES-16 ABI L2+ on Google Cloud
    
    t = time.time()
    
    lon_adjust = 0.5*np.diff(lons[nc_indx_bounds[0]:nc_indx_bounds[1],
                              nc_indx_bounds[2]-1:nc_indx_bounds[3]],axis=1) # this shifts lons just for pcolormesh (lower-left corner plot shift)
    lat_adjust = 0.5*np.diff(lats[nc_indx_bounds[0]:nc_indx_bounds[1]+1,
                              nc_indx_bounds[2]:nc_indx_bounds[3]],axis=0) # this shifts lats just for pcolormesh (lower-left corner plot shift)
    
    # Get the actual LST data:
    prod_vals = ((prod_set.variables[prod_names[0]])[nc_indx_bounds[0]:nc_indx_bounds[1],nc_indx_bounds[2]:nc_indx_bounds[3]]).data
    # get the data quality flag data
    dqf = ((prod_set.variables[prod_names[1]])[nc_indx_bounds[0]:nc_indx_bounds[1],nc_indx_bounds[2]:nc_indx_bounds[3]]).data
    
    # filtering LST data:
    prod_vals[dqf!=0] = np.nan # any non-clear/valid pixels are nans now
    
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
                                              prod_set.time_coverage_end.split('T')[1].replace('Z',''))) # title with date/time
    
    if str_switch:
        print('--------------------------------------------------')
        print("GOES-16 LST data plot time: %.2f s" % (time.time() - t))
    
    ## ################## END GOES-16 LST Data Plotting #########################
    
    ######################################################################################################
    ### Gaussian Fit Coefficient Definition
    # Objective: Calculate Gaussian fit coefficients based on generated data
    
    t = time.time()
    
    col_array,fit_params,gaus_header = [],[],[]
    with open(cwd+aux_dir+'gaus_fit_params.csv','r') as fit_params_csv:
        reader = csv.reader(fit_params_csv,delimiter=',')
        for row in reader:
            if gaus_header==[]:
                gaus_header = row
                continue
            col_array.append(int(row[0])) # this variable identifies the variable column to multiply
            fit_params.append(row[1:]) # the rest of the fit coefficients
    
    input_array = []    
    for ii in range(0,lon_plot.shape[0]):
        for jj in range(0,lon_plot.shape[1]):
            if goes_elev[ii,jj]>-1.0:
                input_ii_jj = [0.0,lat_plot[ii,jj],lon_plot[ii,jj],goes_elev[ii,jj]]
            else:
                input_ii_jj = [0.0,lat_plot[ii,jj],lon_plot[ii,jj],np.nan]
            input_ii_jj.extend(nlcds[ii,jj])
        
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
        print('--------------------------------------------------')
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
    
    # for i in range(T_air[lat_pixel, lon_pixel-pixel_res])
    
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
        print('--------------------------------------------------')
        print("Air temperature calculation runtime: %.2f s" % (time.time() - t))
        
    return area_lst, area_T_air, h_0[lat_pixel, lon_pixel]

# Uncomment this line for in-file troubleshooting, to run indepedently of sensible heat product
# temperature_data(201906231100, 11, 10, 'BKLN')