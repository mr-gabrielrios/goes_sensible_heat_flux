'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        Time Adjust
Package file path:  ~/misc/time_adjust.py
Objective:          Returns the timedelta between the pixel coordinates and UTC
Author:             Gabriel Rios
Notes:              N/A
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

import datetime, pytz, timezonefinder

##############################################################################################
# END IMPORTS
##############################################################################################

##############################################################################################
# Method name:      time_adjust
# Method objective: Take user input for spatial domain, date range.
# Input(s):         lat [float], lon [float]
# Outputs(s):       utc_offset [timedelta]
##############################################################################################

def time_adjuster(lat, lon):
    # Get timezone string for coordinates 
    # Example: crd = [40.17759, -112.45244] returns 'America/Denver'
    loc = timezonefinder.TimezoneFinder().timezone_at(lat=lat, lng=lon)
    # Dummy date to get time offset
    date = datetime.datetime(year=1970, month=1, day=1, hour=0, minute=0)
    tz = pytz.timezone(loc)
    # Note: Arizona doesn't observe DST, so Arizona locations may experience timeshifts of 1 hour
    utc_offset = tz.utcoffset(date, is_dst=True)
    
    return utc_offset

if __name__ == '__main__':
    # Troubleshooting: generate random coordinate within CONUS and print UTC offsetv
    import random
    lat = random.uniform(24, 50)
    lon = random.uniform(-125, -66)
    print('({0:.3f} N, {1:.3f} W)'.format(lat, lon))
    print(time_adjuster(lat, lon))

