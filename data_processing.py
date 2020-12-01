### Script objective: output N-dimensional NumPy array with desired data
# Inputs (and data types):
    # ABI product level (int)
    # ABI dataset (str)
    # Start time (long)
    # End time (long)
    # Working directory (string)
# Outputs (and data types):
    # N-dimensional NumPy array (ndarray)
  
##############################################################################
    
import os
import os.path
from netCDF4 import Dataset
import numpy as np
import time
import grid_processes as gp

### Data gathering algorithm
### Function objective: list all .nc files that satisfy input conditions
# Inputs (and data types)
    # ABI product level (int)
    # ABI dataset (str)
    # Start time (long)
    # End time (long)
    # Root directory with .nc files (str)
# Outputs (and data types):
    # List of complying files (list)      
# Note: goes-data-storage.py must be run before this

def gather_data(gcs_times, ncdir):
    
    t = time.time()
    
    date_length = 14 # Length of GOES filename time markers
    date_markers = ["_s", "_e"] # Markers for GOES filename start and end dates
    
    file_list = [] # Create container for list of files within specified timespan
    
    # Iterate over all .nc files in root directory and add those within 
    # the specified timespan to file list
    for root, dirs, files in os.walk(ncdir):
        for fn in [f for f in files if f.endswith(".nc")]:
            # Get indices where start and end markers, _s and _e, end
            start_idx = fn.find(date_markers[0]) + len(date_markers[0])
            end_idx = fn.find(date_markers[1]) + len(date_markers[1])
            # Grab substrings with dates and convert them to integers for comparison
            fst = int(fn[start_idx:start_idx + date_length]) # File start time
            fet = int(fn[end_idx:end_idx + date_length]) # File end time
            # If the specified start time is before the iterand file start time
            # and greater than the iterand file end time, grab it
            if int(gcs_times[0]) <= fst and int(gcs_times[1]) >= fet:
                file_list.append(os.path.join(root, fn))
        
    def file_sort(fn):
        return(fn[-17:-4])    
    
    file_list = sorted(file_list, key=file_sort)
    
    print('-----------------------')
    print("gather_data() runtime: %.3f" % (time.time()-t))
    print('-----------------------')
    
    return file_list

##############################################################################

### GOES satellite projection and data processing algorithm
### Function objective: project GOES satellite data onto grid of lat/lon coordinates, provide output data
# Inputs (and data types)
    # Root directory with .nc files (str)
    # File iterand (int)
# Outputs (and data types):
    # Grid of lon/lat degree values (np.meshgrid)
    # Data array (np.ndarray)
    # Dataset abbreviated name (str)
    # Dataset full name (str)
    # Data units (str)
    # Data time (datetime) 
# Note: Intended to be iterated over for list of .nc files

def goes_img_nav(nc_file, *args):
    
    start_time = time.time()
    
    t = time.time()    
    # Access .nc file
    working_dir = os.getcwd() # Save current working directory
    g16_data_file = nc_file  # To be modified when user specifies day and hour
    print("Access data time: %.3f" % (time.time()-t))

    t = time.time()    
    # Designate dataset
    g16nc = Dataset(g16_data_file, 'r')  # Create dataset from .nc file
    dataset_name = [i for i in g16nc.variables][0]  # Gets name of relevant dataset
    # dataset_name = "LSTC"  # Gets name of relevant dataset
    print("Designate data time: %.3f" % (time.time()-t))

    t = time.time()    
    # .nc file metadata
    dataset_long_name = g16nc.variables[dataset_name].long_name
    # dataset_long_name = "ABI L2+ Land Surface Temperature"
    data_units = g16nc.variables[dataset_name].units
    data_time = ((g16nc.time_coverage_end).replace('T', ', ')).replace('Z', '')
    print("Metadata time: %.3f" % (time.time()-t))

    t = time.time()    
    # GOES-R projection info
    proj_info = g16nc.variables['goes_imager_projection']
    lon_origin = proj_info.longitude_of_projection_origin
    H = proj_info.perspective_point_height + proj_info.semi_major_axis
    r_eq = proj_info.semi_major_axis
    r_pol = proj_info.semi_minor_axis    
    lambda_0 = lon_origin * np.pi / 180  # Longitude of origin converted to radians
    print("Projection info time: %.3f" % (time.time()-t))

    t = time.time()    
    # GOES-R grid info
    
    if args:
        [central_lon, central_lat, bound_sz] = args[:]
        [x0, x1, y1, y0] = gp.grid_grab(central_lon, central_lat, bound_sz, H, r_pol, r_eq, lambda_0)
        lat_rad_1d = g16nc.variables['x'][x0:x1]
        lon_rad_1d = g16nc.variables['y'][y0:y1]
        data = g16nc.variables[dataset_name][y0:y1, x0:x1]
    else:
        lat_rad_1d = g16nc.variables['x'][:]
        lon_rad_1d = g16nc.variables['y'][:]
        data = g16nc.variables[dataset_name][:]
    print("Grid info read time: %.3f" % (time.time()-t))

    # Close .nc file
    g16nc.close()
    g16nc = None
    os.chdir(working_dir)  # Revert to script directory

    t = time.time()    
    # Create meshgrid
    lat_rad, lon_rad = np.meshgrid(lat_rad_1d, lon_rad_1d)  # x and y (reference PUG, Section 4.2.8.1)   
    print("Mesh grid time: %.3f" % (time.time()-t))

    t = time.time()    
    # Latitude/longitude projection calculation from satellite radian angle vectors (reference PUG, Section 4.2.8.1)
    a = np.power(np.sin(lat_rad), 2) + np.power(np.cos(lat_rad), 2) * (
                np.power(np.cos(lon_rad), 2) + np.power(np.sin(lon_rad), 2) * np.power(r_eq, 2) / np.power(r_pol, 2))
    b = (-2) * H * np.cos(lat_rad) * np.cos(lon_rad)
    c = np.power(H, 2) - np.power(r_eq, 2)
    r_s = ((-b) - np.sqrt(np.power(b, 2) - 4 * a * c)) / (2 * a)  # distance from satellite to surface

    s_x = r_s * np.cos(lat_rad) * np.cos(lon_rad)
    s_y = (-1 * r_s) * np.sin(lat_rad)
    s_z = r_s * np.cos(lat_rad) * np.sin(lon_rad)    
    print("Projection time: %.3f" % (time.time()-t))

    ########################################################################################

    t = time.time()    
    # Transform radian latitude and longitude values to degrees
    lat_deg = (180 / np.pi) * np.arctan(
        (np.power(r_eq, 2) / np.power(r_pol, 2)) * s_z / (np.sqrt((np.power((H - s_x), 2)) + np.power(s_y, 2))))
    lon_deg = (180 / np.pi) * (lambda_0 - np.arctan(s_y / (H - s_x)))
    print("Transformation time: %.3f" % (time.time()-t))
    
    t = time.time()    
    # Create mask to filter out data where values are default
    data_mask = np.where(data.data != 65535)
    # Generate ndarray with relevant data for future processing/stats
    output_data = np.where(data.data != 65535, data, float('NaN'))
    
    # Create bound box for plotting purposes
    bound_box = [np.min(lon_deg),
                  np.max(lon_deg),
                  np.min(lat_deg),
                  np.max(lat_deg)]
    
    print("Masking time: %.3f" % (time.time()-t))
    print('-----------------------')
    print("goes_img_nav() runtime: %.3f" % (time.time()-start_time))
    print('-----------------------')
    return output_data, lat_deg, lon_deg, data, dataset_name, dataset_long_name, data_units, data_time, bound_box

##############################################################################

### GOES satellite projection and data processing algorithm
### Function objective: project GOES satellite data onto grid of lat/lon coordinates, provide output data
# Inputs (and data types)
# Outputs (and data types):
# Note: Intended to be iterated over for list of .nc files

def goes_data_processor(file_list, *args):
    data_list = []
    
    for file in file_list:
        print(file)
        output_data, lat_deg, lon_deg, data, dataset_name, dataset_long_name, \
            data_units, data_time, bound_box = goes_img_nav(file, *args)
        data_list.append([output_data, lat_deg, lon_deg, data, dataset_name, dataset_long_name, \
            data_units, data_time, bound_box])        
    return data_list

##############################################################################

### 
###
# Inputs (and data types)
# Outputs (and data types):
    
def data_stats(data_list):    
    i = 0
    for file in data_list: 
        print(np.nanmin(file[0]), np.nanmax(file[0]))
        if i == 0:
            min_val = np.nanmin(file[0])
            max_val = np.nanmax(file[0])
        else:
            if np.nanmin(file[0]) < min_val:
                min_val = np.nanmin(file[0])
            elif np.nanmax(file[0]) > max_val:
                max_val = np.nanmax(file[0])     
        i += 1
            
    print(min_val, max_val)
    

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #
### OLD SCRIPT

# ### Script objective: output N-dimensional NumPy array with desired data
# # Inputs (and data types):
#     # ABI product level (int)
#     # ABI dataset (str)
#     # Start time (long)
#     # End time (long)
#     # Working directory (string)
# # Outputs (and data types):
#     # N-dimensional NumPy array (ndarray)
  
# ##############################################################################
    
# import os
# import os.path
# from netCDF4 import Dataset
# import numpy as np
# import time

# ### Data gathering algorithm
# ### Function objective: list all .nc files that satisfy input conditions
# # Inputs (and data types)
#     # ABI product level (int)
#     # ABI dataset (str)
#     # Start time (long)
#     # End time (long)
#     # Root directory with .nc files (str)
# # Outputs (and data types):
#     # List of complying files (list)      
# # Note: goes-data-storage.py must be run before this

# def gather_data(gcs_times, ncdir):
    
#     t = time.time()
    
#     date_length = 14 # Length of GOES filename time markers
#     date_markers = ["_s", "_e"] # Markers for GOES filename start and end dates
    
#     file_list = [] # Create container for list of files within specified timespan
    
#     # Iterate over all .nc files in root directory and add those within 
#     # the specified timespan to file list
#     for root, dirs, files in os.walk(ncdir):
#         for fn in [f for f in files if f.endswith(".nc")]:
#             # Get indices where start and end markers, _s and _e, end
#             start_idx = fn.find(date_markers[0]) + len(date_markers[0])
#             end_idx = fn.find(date_markers[1]) + len(date_markers[1])
#             # Grab substrings with dates and convert them to integers for comparison
#             fst = int(fn[start_idx:start_idx + date_length]) # File start time
#             fet = int(fn[end_idx:end_idx + date_length]) # File end time
#             # If the specified start time is before the iterand file start time
#             # and greater than the iterand file end time, grab it
#             if int(gcs_times[0]) <= fst and int(gcs_times[1]) >= fet:
#                 file_list.append(os.path.join(root, fn))
        
#     def file_sort(fn):
#         return(fn[-17:-4])    
    
#     file_list = sorted(file_list, key=file_sort)
    
#     print('-----------------------')
#     print("gather_data() runtime: %.3f" % (time.time()-t))
#     print('-----------------------')
    
#     return file_list

# ##############################################################################

# ### GOES satellite projection and data processing algorithm
# ### Function objective: project GOES satellite data onto grid of lat/lon coordinates, provide output data
# # Inputs (and data types)
#     # Root directory with .nc files (str)
#     # File iterand (int)
# # Outputs (and data types):
#     # Grid of lon/lat degree values (np.meshgrid)
#     # Data array (np.ndarray)
#     # Dataset abbreviated name (str)
#     # Dataset full name (str)
#     # Data units (str)
#     # Data time (datetime) 
# # Note: Intended to be iterated over for list of .nc files

# def goes_img_nav(nc_file, *args):
    
#     start_time = time.time()
    
#     t = time.time()    
#     # Access .nc file
#     working_dir = os.getcwd() # Save current working directory
#     g16_data_file = nc_file  # To be modified when user specifies day and hour
#     print("Access data time: %.3f" % (time.time()-t))

#     t = time.time()    
#     # Designate dataset
#     g16nc = Dataset(g16_data_file, 'r')  # Create dataset from .nc file
#     dataset_name = [i for i in g16nc.variables][0]  # Gets name of relevant dataset
#     # dataset_name = "LSTC"  # Gets name of relevant dataset
#     print("Designate data time: %.3f" % (time.time()-t))

#     t = time.time()    
#     # .nc file metadata
#     dataset_long_name = g16nc.variables[dataset_name].long_name
#     # dataset_long_name = "ABI L2+ Land Surface Temperature"
#     data_units = g16nc.variables[dataset_name].units
#     data_time = ((g16nc.time_coverage_end).replace('T', ', ')).replace('Z', '')
#     print("Metadata time: %.3f" % (time.time()-t))

#     t = time.time()    
#     # GOES-R projection info
#     proj_info = g16nc.variables['goes_imager_projection']
#     lon_origin = proj_info.longitude_of_projection_origin
#     H = proj_info.perspective_point_height + proj_info.semi_major_axis
#     r_eq = proj_info.semi_major_axis
#     r_pol = proj_info.semi_minor_axis
#     print("Projection info time: %.3f" % (time.time()-t))

#     t = time.time()    
#     # GOES-R grid info
    
#     lat_rad_1d = g16nc.variables['x'][:]
#     lon_rad_1d = g16nc.variables['y'][:]
#     data = g16nc.variables[dataset_name][:]
#     # print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
#     # print(g16nc.variables['x'][1809:1889])
#     # print(len(g16nc.variables['x'][1809:1889]))
#     # print(g16nc.variables['y'][282:354])
#     # print(len(g16nc.variables['y'][282:354]))
#     # print(data.shape)
#     # print(data[282:354, 1809:1889])
#     # print(np.amax(data[282:354, 1809:1889]))
#     # print(data[282:354, 1809:1889].shape)
#     # print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
#     print("Grid info read time: %.3f" % (time.time()-t))

#     # Close .nc file
#     g16nc.close()
#     g16nc = None
#     os.chdir(working_dir)  # Revert to script directory

#     t = time.time()    
#     # Create meshgrid
#     lat_rad, lon_rad = np.meshgrid(lat_rad_1d, lon_rad_1d)  # x and y (reference PUG, Section 4.2.8.1)   
#     print('Lat rad min and max:')    
#     print(max(lat_rad_1d), min(lat_rad_1d)) 
#     print('Lon rad min and max:') 
#     print(max(lon_rad_1d), min(lon_rad_1d))
#     print("Mesh grid time: %.3f" % (time.time()-t))

#     t = time.time()    
#     # Latitude/longitude projection calculation from satellite radian angle vectors (reference PUG, Section 4.2.8.1)
#     lambda_0 = lon_origin * np.pi / 180  # Longitude of origin converted to radians
#     a = np.power(np.sin(lat_rad), 2) + np.power(np.cos(lat_rad), 2) * (
#                 np.power(np.cos(lon_rad), 2) + np.power(np.sin(lon_rad), 2) * np.power(r_eq, 2) / np.power(r_pol, 2))
#     b = (-2) * H * np.cos(lat_rad) * np.cos(lon_rad)
#     c = np.power(H, 2) - np.power(r_eq, 2)
#     r_s = ((-b) - np.sqrt(np.power(b, 2) - 4 * a * c)) / (2 * a)  # distance from satellite to surface

#     s_x = r_s * np.cos(lat_rad) * np.cos(lon_rad)
#     s_y = (-1 * r_s) * np.sin(lat_rad)
#     s_z = r_s * np.cos(lat_rad) * np.sin(lon_rad)    
#     print("Projection time: %.3f" % (time.time()-t))

#     ########################################################################################

#     t = time.time()    
#     # Transform radian latitude and longitude values to degrees
#     lat_deg = (180 / np.pi) * np.arctan(
#         (np.power(r_eq, 2) / np.power(r_pol, 2)) * s_z / (np.sqrt((np.power((H - s_x), 2)) + np.power(s_y, 2))))
#     lon_deg = (180 / np.pi) * (lambda_0 - np.arctan(s_y / (H - s_x)))
#     print('Latitude extrema')
#     print(np.min(lat_deg), np.max(lat_deg))
#     print('Longitude extrema')
#     print(np.min(lon_deg), np.max(lon_deg))
#     print("Transformation time: %.3f" % (time.time()-t))
    
#     t = time.time()    
#     # Create mask to filter out data where values are default
#     data_mask = np.where(data.data != 65535)
#     # Generate ndarray with relevant data for future processing/stats
#     output_data = np.where(data.data != 65535, data, float('NaN'))
    
#     # Create bound box for plotting purposes
#     bound_box = [np.min(lon_deg[data_mask]),
#                   np.max(lon_deg[data_mask]),
#                   np.min(lat_deg[data_mask]),
#                   np.max(lat_deg[data_mask])]
    
#     # If statement filters out data beyond user-defined window
#     if args:
#         [central_lon, central_lat, bound_sz] = args[:]
#         # Create window from user input
#         bound_box = [central_lon - bound_sz,
#                       central_lon + bound_sz,
#                       central_lat - bound_sz,
#                       central_lat + bound_sz]        
#         # Generate lons for points within window boundaries
#         lon_deg_idx = np.where(np.logical_and(lon_deg>=central_lon-bound_sz,lon_deg<=central_lon+bound_sz), lon_deg, np.nan)
#         # Generate lats for points within window boundaries
#         lat_deg_idx = np.where(np.logical_and(lat_deg>=central_lat-bound_sz,lat_deg<=central_lat+bound_sz), lat_deg, np.nan)
#         # Generate ndarray with relevant data for future processing/stats within
#         # window boundaries
#         output_data_idx = np.where(np.logical_and(~np.isnan(lon_deg_idx), ~np.isnan(lat_deg_idx)), output_data, np.nan)
#         # Re-filter data, arbitrary upper limit set at 1000
#         output_data = np.where(~np.isnan(output_data_idx), output_data, np.nan)
#         output_data = output_data[~np.isnan(output_data)]
         
#         print("Masking time: %.3f" % (time.time()-t))
#         print('-----------------------')
#         print("goes_img_nav() runtime: %.3f" % (time.time()-start_time))
#         print('-----------------------')
        
#         return output_data, lat_deg, lon_deg, data, dataset_name, dataset_long_name, data_units, data_time, bound_box
    
#     # Else, plot CONUS or Full-Disk image
#     else:
#         print("Masking time: %.3f" % (time.time()-t))
#         print('-----------------------')
#         print("goes_img_nav() runtime: %.3f" % (time.time()-start_time))
#         print('-----------------------')
#         return output_data, lat_deg, lon_deg, data, dataset_name, dataset_long_name, data_units, data_time, bound_box

# ##############################################################################

# ### GOES satellite projection and data processing algorithm
# ### Function objective: project GOES satellite data onto grid of lat/lon coordinates, provide output data
# # Inputs (and data types)
# # Outputs (and data types):
# # Note: Intended to be iterated over for list of .nc files

# def goes_data_processor(file_list, *args):
#     data_list = []
    
#     for file in file_list:
#         print(file)
#         output_data, lat_deg, lon_deg, data, dataset_name, dataset_long_name, \
#             data_units, data_time, bound_box = goes_img_nav(file, *args)
#         data_list.append([output_data, lat_deg, lon_deg, data, dataset_name, dataset_long_name, \
#             data_units, data_time, bound_box])        
#     return data_list

# ##############################################################################

# ### 
# ###
# # Inputs (and data types)
# # Outputs (and data types):
    
# def data_stats(data_list):    
#     i = 0
#     for file in data_list: 
#         print(np.nanmin(file[0]), np.nanmax(file[0]))
#         if i == 0:
#             min_val = np.nanmin(file[0])
#             max_val = np.nanmax(file[0])
#         else:
#             if np.nanmin(file[0]) < min_val:
#                 min_val = np.nanmin(file[0])
#             elif np.nanmax(file[0]) > max_val:
#                 max_val = np.nanmax(file[0])     
#         i += 1
            
#     print(min_val, max_val)