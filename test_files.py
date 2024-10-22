import xarray as xr
import matplotlib.pyplot as plt
import glob
import numpy as np

# Find all NetCDF files in the current directory
nc_files = glob.glob('*.nc')

# Open the first file
ds = xr.open_dataset(nc_files[3])
print(ds)

year = 5
## min and max of the connectivity matrix for the selected year
print(np.nanmin(ds.connectivity.isel(year=year)), np.nanmax(ds.connectivity.isel(year=year)))

# Plotting the connectivity matrix for the first year for 2 simulations
ds.isel(year=year).connectivity.plot.imshow(col='scenario', col_wrap=2)
plt.show()
