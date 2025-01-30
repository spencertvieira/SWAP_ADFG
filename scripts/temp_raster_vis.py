import xarray as xr
import rioxarray as rxr
import os

# Define file paths
input_file = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/raw/monthly/2m_temperature_2010_monthly.nc"
output_dir = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/dump"
output_file = os.path.join(output_dir, "2m_temperature_2010_month1.tif")

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Open dataset
ds = xr.open_dataset(input_file)

# Select the first month's data (valid_time index 0)
t2m = ds["t2m"].isel(valid_time=0)

# Convert to a GeoTIFF-compatible format
t2m = t2m.rio.write_crs("EPSG:4326")  # Assuming lat/lon in degrees, adjust if necessary

# Save as GeoTIFF
t2m.rio.to_raster(output_file)

print(f"GeoTIFF saved to: {output_file}")
