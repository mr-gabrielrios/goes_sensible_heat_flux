'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        Land Surface Temperature
Package file path:  ~/phys/t_lst.py
Objective:          Read LST from GOES-16 Level 2 product for spatial pixel.
Author:             Gabriel Rios
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

# External imports
from netCDF4 import Dataset
import datetime, numpy as np, os, pytz, time, timezonefinder

##############################################################################################
# END IMPORTS
##############################################################################################

##############################################################################################
# Method name:      data_download
# Method objective: Download GOES-16 data from Google Cloud for selected time.
# Input(s):         domain center [list], date [datetime], bucket [Google Client object]
# Outputs(s):       goes_products [netCDF file object], goes_varnames [list]
##############################################################################################

def data_download(domain_center, date, bucket):
    # Create directory where temporary GOES netCDF file will be stored
    dir_name = 'goes_data'
    # If directory doesn't exist, make it
    if not os.path.isdir(os.path.join(os.path.dirname(__file__), dir_name)):
        os.mkdir(os.path.join(os.path.dirname(__file__), dir_name))
    goes_dir = os.path.join(os.path.dirname(__file__), dir_name)
    
    prod = 'LST' # Level L-2 Product identifier
    product_ID = 'ABI-L2-' + prod + 'C/' # ABI, L2 product, CONUS
    
    goes_16_date_vec, goes_16_file_vec = [], []
    
    print(date.strftime('%Y-%m-%d %H:%M:%S'))
    
    while goes_16_file_vec == []:
        yr = '{:0004d}'.format(date.year) 
        yday = '{:003d}'.format(date.timetuple().tm_yday)
        thour = '{:02d}'.format(date.hour)
    
        # this is the relative folder where the GOES-16 data will be
        prefix_16 = product_ID + yr + '/' + yday + '/' + thour 
        blobs_16 = bucket.list_blobs(prefix=prefix_16) # get all files in prefix
    
        for kk_16, blob_16 in enumerate(blobs_16): # loop through files
    
            if len(blob_16.name.split('_')) <= 5:
                end_time_16 = blob_16.name.split('_')[4][1:-3]
                end_datetime_16 = datetime.datetime.strptime(end_time_16,'%Y%j%H%M%S%f')
            else:
                end_time_16 = blob_16.name.split('_')[4][1:]
                end_datetime_16 = datetime.datetime.strptime(end_time_16,'%Y%j%H%M%S%f')
    
            goes_16_date_vec.append(end_datetime_16)
            goes_16_file_vec.append(blob_16)
        
        if goes_16_file_vec == []:
            date = date - datetime.timedelta(hours=1) # if there are no files at that time, look an hour behind
        
    curr_filename = (goes_16_file_vec[0].name).split('/')[-1] # the the last file in search
    if os.path.isfile(goes_16_file_vec[0].name+curr_filename):
        pass
    else:
        # download the file to the local GOES-16 repository
        (goes_16_file_vec[0]).download_to_filename(goes_dir+curr_filename)
    
    try:
        # Load LST data
        goes_product = Dataset(goes_dir+curr_filename) 
        # Remove LST data from local directory
        os.remove(goes_dir+curr_filename)
    except:
        # If there was an error, remove the faulty file and re-download
        os.remove(goes_dir+curr_filename) 
        (goes_16_file_vec[0]).download_to_filename(goes_dir+curr_filename)
        # Load LST data
        goes_product = Dataset(goes_dir+curr_filename)
        # Remove LST data from local directory
        os.remove(goes_dir+curr_filename)
    
    goes_varnames = []
    for i in goes_product.variables:
        # Get all variable names in GOES-16 file
        goes_varnames.append(i)

    return goes_product, goes_varnames

##############################################################################################
# Method name:      lst
# Method objective: Get land surface temperature from GOES-16 data for select pixel.
# Input(s):         goes_products [netCDF file object], goes_varnames [list], date [datetime], goes_idx [list]
# Outputs(s):       T_s [float]
##############################################################################################

def lst(goes_product, goes_varnames, date, goes_idx):
    t = time.time()
    
    # Get the actual LST data
    T_s = ((goes_product.variables[goes_varnames[0]])[goes_idx[0]:goes_idx[1], goes_idx[2]:goes_idx[3]]).data
    # Retrieve the data quality flag
    dqf = ((goes_product.variables[goes_varnames[1]])[goes_idx[0]:goes_idx[1], goes_idx[2]:goes_idx[3]]).data
    
    # Identify any water or non-clear/valid pixels and set to nan
    # T_s[water_locs_ii,water_locs_jj] = np.nan
    T_s[dqf != 0] = np.nan
    
    runtime = time.time() - t
    
    return T_s
    
if __name__ == '__main__':
    # Test to see if LST algorithm is working properly
    def gcp_bucket():
        from google.cloud import storage
        
        # Reference local Google Credentials
        os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "./GOOGLE_AUTH_CREDS_CUERG.json"
        # Start the storage client
        client = storage.Client() 
        # Call the GOES-16 storage bucket
        bucket = client.get_bucket('gcp-public-data-goes-16') 
    
        return bucket
    
    bucket = gcp_bucket()
    pixel = [40.6305, -73.9521]
    date = datetime.datetime(year=2019, month=7, day=28, hour=15)
    goes_idx = [307, 344, 1829, 1871]
    T_s = lst(pixel, date, goes_idx, bucket)