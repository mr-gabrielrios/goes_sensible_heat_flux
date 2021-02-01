'''
Package:            GOES-16 Sensible Heat Flux Numerical Model
Script name:        GOES Grid Setup
Package file path:  ~/grid/goes_grid_setup.py
Objective:          Create grid to use as spatial basis for the package
Author:             Gabriel Rios
Notes:              Credit to Josh Hrisko @maker-portal for initial development of this code
'''

##############################################################################################
# BEGIN IMPORTS
##############################################################################################

import numpy as np, os, time
from scipy import spatial
from netCDF4 import Dataset

##############################################################################################
# END IMPORTS
##############################################################################################

##############################################################################################
# Method name:      goes_grid
# Method objective: Provide 2D arrays representing coordinate grids for selected spatial domain.
# Input(s):         domain [list] - format: (lat min, lat max, lon min, lon max)
# Outputs(s):       lats [2D ndarray], lons [2D ndarray], goes_idx [1D ndarray]
##############################################################################################

def goes_grid(domain):
    t = time.time()
    
    # Retrieve GOES-16 coordinate data from dedicated netCDF file
    goes_grid_file = os.path.join(os.path.dirname(__file__), 'GOESR_ABI_CONUS_East.nc')
    grid_dataset = Dataset(goes_grid_file)
    goes_lats, goes_lons = grid_dataset.variables['Latitude'][:].data, grid_dataset.variables['Longitude'][:].data
    
    # Boundary grid indices based on defined spatial domain
    nc_indx_corners = np.array([np.unravel_index(np.argmin(np.abs(np.subtract(goes_lons, domain[2])) + np.abs(np.subtract(goes_lats,domain[0]))), np.shape(goes_lats)), \
                                np.unravel_index(np.argmin(np.abs(np.subtract(goes_lons,domain[2]))+np.abs(np.subtract(goes_lats,domain[1]))),np.shape(goes_lats)), \
                                np.unravel_index(np.argmin(np.abs(np.subtract(goes_lons,domain[3]))+np.abs(np.subtract(goes_lats,domain[0]))),np.shape(goes_lats)), \
                                np.unravel_index(np.argmin(np.abs(np.subtract(goes_lons,domain[3]))+np.abs(np.subtract(goes_lats,domain[1]))),np.shape(goes_lats))])
        
    # Create index clip for coordinate and data clipping
    # Format: [max latitude index, min latitude index, min longitude index, max longitude index] (per PUG-L2)
    goes_idx = [np.min([nc_indx_corners[0][0],nc_indx_corners[1][0],nc_indx_corners[2][0],nc_indx_corners[3][0]]),
                np.max([nc_indx_corners[0][0],nc_indx_corners[1][0],nc_indx_corners[2][0],nc_indx_corners[3][0]]),
                np.min([nc_indx_corners[0][1],nc_indx_corners[1][1],nc_indx_corners[2][1],nc_indx_corners[3][1]]),
                np.max([nc_indx_corners[0][1],nc_indx_corners[1][1],nc_indx_corners[2][1],nc_indx_corners[3][1]])]
    
    # Clip coordinates with regards to defined spatial domain
    lons = goes_lons[goes_idx[0]:goes_idx[1], goes_idx[2]:goes_idx[3]]
    lats = goes_lats[goes_idx[0]:goes_idx[1], goes_idx[2]:goes_idx[3]]
    
    runtime = time.time() - t    
    
    return lats, lons, goes_idx

##############################################################################################
# Method name:      crd_idx
# Method objective: Return data grid index given a latitude and longitude
# Input(s):         crd [list or tuple], lats [2D ndarray], lons [2D ndarray]
# Outputs(s):       i [int], j [int]
##############################################################################################

def crd_idx(crd, lats, lons):
    _, index = spatial.KDTree(np.c_[lats.ravel(), lons.ravel()]).query(crd)
    i = index // lats.shape[1] - 1 # Get row number
    j = index % lons.shape[1] - 1 # Get column number
    return i, j

##############################################################################################
# Method name:      grid_alignment
# Method objective: Translate data from a WRF grid to align with the GOES grid.
# Input(s):         goes_grid [list], wrf_grid [list]
# Outputs(s):       goes_data [2D ndarray]
##############################################################################################

def grid_alignment(goes_grid, wrf_grid):
    
    goes_lons, goes_lats = goes_grid
    wrf_lons, wrf_lats, wrf_data = wrf_grid
    goes_data = np.empty(np.shape(goes_lons))
            
    tree = spatial.KDTree(np.c_[wrf_lons, wrf_lats])
    for i in range(0, np.shape(goes_data)[0]):
        for j in range(0, np.shape(goes_data)[1]):
            goes_lon, goes_lat = [goes_lons[i, j], goes_lats[i, j]]
            dist, ii = tree.query([(goes_lon, goes_lat)])
            goes_data[i, j] = wrf_data.ravel()[ii]
            
    return goes_data

# Use this method for troubleshooting
if __name__ == '__main__': 
    # Produce map with inset focused on spatial domain as script test
    import cartopy.crs as ccrs
    import cartopy.io.img_tiles as cimgt
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    
    # Define plot data
    domain = [40.50, 40.90, -74.15, -73.75]
    lats, lons, goes_idx = goes_grid(domain)
    _lats = [np.amin(lats), np.amax(lats)]
    _lons = [np.amin(lons), np.amax(lons)]
    # Coordinate-to-index test
    i, j = crd_idx([40.849035, -73.875716], lats, lons)
    
    # Define projections
    proj = ccrs.PlateCarree()
    # proj = ccrs.Orthographic(central_latitude=sum(_lats)/2, central_longitude=sum(_lons)/2)
    fig, ax = plt.subplots(subplot_kw={'projection': proj}, dpi=300)
    # Miscellaneous plot formatting
    zoom = 0.1
    ax.set_extent([_lons[0] - zoom, _lons[1] + zoom, _lats[0] - zoom, _lats[1] + zoom])
    ax.coastlines('10m')
    # Add the satellite background data at zoom level 12. Left uncommented to save troubleshooting time.
    # sat_img = cimgt.GoogleTiles(style='satellite')
    # ax.add_image(sat_img, 12)
    # Draw user-defined spatial domain for reference
    box = patches.Rectangle(xy=(_lons[0], _lats[0]), 
                           width=(_lons[1] - _lons[0]), 
                           height=(_lats[1] - _lats[0]), 
                           lw=1, 
                           edgecolor='r', 
                           facecolor='none',
                           transform=ccrs.PlateCarree())
    ax.add_patch(box)
    # Create filler data. Note that pcolormesh data dimensions must 1 smaller along each axis than the X and Y arrays.
    grid_data = np.zeros([lats.shape[0]-1, lats.shape[1]-1])
    # Highlight cell selected by coordinate-to-index test
    grid_data[i, j] = 1
    # Plot
    ax.pcolormesh(lons, lats, grid_data, cmap='Greys', linewidth=0.5, alpha=0.2, edgecolor='k', transform=ccrs.PlateCarree())
    # Draw and format gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.1, linestyle='-')
    gl.top_labels = False
    gl.right_labels = False
    