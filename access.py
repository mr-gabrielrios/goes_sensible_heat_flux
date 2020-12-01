import os
from datetime import timedelta
from datetime import datetime
import time

# Objective: specify GOES-16 data type (L1 or L2) and date to prompt download 
#            from Google Cloud Storage
# Input: user prompts for desired data product and date
# Output: strings for Google Cloud Storage and local directories for desired data
# Note: GOES only uploads L2 products hourly

def create_path_dirs(**kwargs):
            
    t = time.time()
    
    # Grab date from user input - input_bool controls while loop
    input_bool = [False, False]
    # Generate empty list to hold dates
    input_dates = 2*[None]
    
    ### Pre-determined input (useful for functions)
    if kwargs:
        for key, value in kwargs.items():
            if key == "goes_level":
                input_product_level = 'L' + str(value)  
            elif key == "goes_product":
                input_dataset = str(value)
            elif key == "start_date":
                input_dates[0] = str(value)
            elif key == "end_date":
                input_dates[1] = str(value)               
    
    ### User input (useful for singular use) 
    else: 
        # Grab product level desired by user
        input_product_level = input('Enter a GOES-16 ABI product level: ')
        input_product_level = 'L' + input_product_level
    
        # Grab data type desired by user
        input_dataset = input('Enter a GOES-16 ABI dataset: ')
        
        # Get start date from user input
        while input_bool[0] == False:  # Ensures only good data is input
            date_str = input('Enter a start date in YYYY-MM-DD-HH format:\n HH ranges from 00 to 23 ')
            if len(date_str) == 13 and date_str[4] == '-' and date_str[7] == '-' and date_str[10] == '-':
                input_bool[0] = True
            else:
                print('\n Whoa there skipper, read the instructions and try again!')
                continue
        input_dates[0] = date_str.replace('-', '')
            
        # Get end date from user input
        while input_bool[1] == False:  # Ensures only good data is input
            date_str = input('Enter an end date in YYYY-MM-DD-HH format:\n HH ranges from 00 to 23 ')
            if len(date_str) == 13 and date_str[4] == '-' and date_str[7] == '-' and date_str[10] == '-':
                input_bool[1] = True
            else:
                print('\n Whoa there skipper, read the instructions and try again!')
                continue
        input_dates[1] = date_str.replace('-', '')
    
    # This is the root directory for all locally-stored GOES-16 data
    g16dir = r'ncdata/' + input_dataset + '/'
    isdir = os.path.isdir(g16dir) # Boolean to check if ncdata folder already exists
    if not os.path.exists(g16dir): # If it doesn't exist, create it as a subfolder
        os.makedirs(os.path.join(os.getcwd(), g16dir))
    
    # Generate web path string
    gcs_path = "ABI-" + input_product_level + "-" + input_dataset
    # Initialize list of GCS bucket paths and times
    gcs_times = [None] * len(input_dates)

    # Generate 3-digit day of year from input_date
    i = 0 # Iterand for for-loop
    for date in input_dates:
        doy = str(datetime(int(date[0:4]), int(date[4:6]),
                           int(date[6:8])).timetuple().tm_yday)
        # Prepend leading zeros to day of the year
        doy = (3 - len(doy)) * "0" + doy
        # Append zeroes for minutes, seconds, and millisecond
        gcs_times[i] = date[0:4] + doy + date[8:10] + "00000"
        i += 1
          
    print(gcs_times)

    print('-----------------------')
    print("create_path_dirs() runtime: %.3f" % (time.time()-t))
    print('-----------------------')

    return gcs_path, gcs_times, g16dir

# Objective: specify GOES-16 data type (L1 or L2) and date to prompt download 
#            Google Cloud Storage
# Input:     user prompts for desired data product and date
# Output:    .nc file downloaded to specified local directory
# Note: GOES only uploads L2 products hourly

def gcs_data_access(gcs_path, gcs_times, local_dir):
    t = time.time()
    
    date_length = 14 # Length of GOES filename time markers
    date_markers = ["_s", "_e"] # Markers for GOES filename start and end dates
    
    # Loop to find minimum and maximum timestamps for downloaded files
    file_times = []
    i = 0
    for file in os.listdir(local_dir):        
        start_idx = file.find(date_markers[0]) + len(date_markers[0])
        end_idx = file.find(date_markers[1]) + len(date_markers[1])
        fst = int(file[start_idx:start_idx + date_length]) # File start time
        fet = int(file[end_idx:end_idx + date_length]) # File end time        
        if i == 0:
            file_times.append(fst)
            file_times.append(fet)
        else:
            if fst < file_times[0]:
                file_times[0] = fst
            if fet > file_times[1]:
                file_times[1] = fet
        i += 1
    
    gcs_times_dl = [gcs_times[0], gcs_times[1]] # Download timestamp bounds
    dl = True # When True, download files
    if file_times:
        # Logic to handle which files to download if one of the dates is beyond the already-downloaded file timestamps
        dl = False # When True, download files
        if int(gcs_times[0]) >= file_times[0] and int(gcs_times[1]) <= file_times[1]:
            print('No download necessary, moving to data processing...')
        elif int(gcs_times[0]) <= file_times[0] and int(gcs_times[1]) <= file_times[1]: # Earlier date requested
            gcs_times_dl[1] = str(file_times[0])
            dl = True
        elif int(gcs_times[0]) >= file_times[0] and int(gcs_times[1]) >= file_times[1]: # Later date requested
            gcs_times_dl[0] = str(file_times[1])
            dl = True
        elif int(gcs_times[0]) <= file_times[0] and int(gcs_times[1]) >= file_times[1]: # Both earlier and later date requested
            gcs_times_dl = gcs_times
            dl = True
        else:
            dl = True
    
    print('Download times: ', gcs_times_dl)
    dt = datetime.strptime(gcs_times_dl[1][0:-1], '%Y%j%H%M%S') - datetime.strptime(gcs_times_dl[0][0:-1], '%Y%j%H%M%S')
    dt = timedelta.total_seconds(dt)/3600
    print(dt)
    
    if dl and dt > 1:  
        # 1. Connect to Google Cloud Storage and authenticate for data access
        from google.cloud import storage
        # Reference credentials on local network
        os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = \
            r"GOOGLE_AUTH_CREDS_CUERG.json"
        client = storage.Client()
    
        # 2. Generate paths for Google Cloud Repo (source) and local directory (destination)
        # NOTE: Please only use for L2 products at the moment. Functionality for L1 products to be added.
        # Get Google Cloud path to read from and local path to write to
        # Connect to GOES-16/-17 data bucket
        bucket = client.get_bucket('gcp-public-data-goes-16')
        # Get .nc files that match with user-designate data
        blobs_16 = bucket.list_blobs(prefix=gcs_path)
        
        for blob in blobs_16:
            # Get indices where start and end markers, _s and _e, end
            start_idx = blob.name.find(date_markers[0]) + len(date_markers[0])
            end_idx = blob.name.find(date_markers[1]) + len(date_markers[1])
            # Grab substrings with dates and convert them to integers for comparison
            fst = int(blob.name[start_idx:start_idx + date_length]) # File start time
            fet = int(blob.name[end_idx:end_idx + date_length]) # File end time
            # If the specified start time is before the iterand file start time
            # and greater than the iterand file end time, grab it
            if int(gcs_times_dl[0]) <= fst and int(gcs_times_dl[1]) >= fet:
                blob_fn = local_dir + r"\\" + blob.name.split('/')[-1]
                # Download the .nc file
                blob.download_to_filename(blob_fn)
                print("File " + blob_fn + " downloaded!")
    
    print('-----------------------')
    print("gcs_data_access() runtime: %.3f" % (time.time()-t))
    print('-----------------------')