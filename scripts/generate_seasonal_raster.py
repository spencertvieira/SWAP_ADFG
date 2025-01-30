import xarray as xr
import rioxarray
from rasterio.enums import Resampling
import pandas as pd
import os

# Define seasons
SEASONS = {
    "DJF": [12, 1, 2],
    "MAM": [3, 4, 5],
    "JJA": [6, 7, 8],
    "SON": [9, 10, 11],
    "Annual": list(range(1, 13)),
    "MJJAS": [5, 6, 7, 8, 9],
    "ONDJFMAM": [10, 11, 12, 1, 2, 3, 4, 5],
}

# Mapping of human-readable variable names to NetCDF variable names
VARIABLE_MAP = {
    "2m_temperature": "t2m",
    "total_precipitation": "tp",
    "snowfall": "sf",
    "evaporation": "e",
    "10m_u_component_of_wind": "u10",
    "10m_v_component_of_wind": "v10",
    "sea_surface_temperature": "sst",
}

# Variables that require unit conversion
CONVERT_TO_MM = ["total_precipitation", "snowfall", "evaporation"]
CONVERT_TO_C = ["2m_temperature", "sea_surface_temperature"]

def generate_seasonal_raster(netcdf_file, variable, season, output_dir):
    """ Generates a seasonal raster from a NetCDF file. """

    # Validate inputs
    if not os.path.exists(netcdf_file):
        raise FileNotFoundError(f"NetCDF file not found: {netcdf_file}")

    if variable not in VARIABLE_MAP:
        raise ValueError(f"Variable '{variable}' not recognized.")

    if season not in SEASONS:
        raise ValueError(f"Season '{season}' is not valid.")

    netcdf_variable_name = VARIABLE_MAP[variable]
    os.makedirs(output_dir, exist_ok=True)

    print(f"Processing file: {netcdf_file}")

    # Load the dataset
    ds = xr.open_dataset(netcdf_file)
    data = ds[netcdf_variable_name]

    print("Original Data Summary:")
    print(data)

    print("Value range before processing:", data.min().values, data.max().values)

    # Detect the correct time dimension
    time_dim = None
    for possible_time in ["valid_time", "time"]:
        if possible_time in data.coords:
            time_dim = possible_time
            break

    if time_dim is None:
        raise KeyError("No valid time coordinate found in the dataset.")

    print(f"Using time dimension: {time_dim}")

    # Ensure `valid_time` is treated as a datetime index
    if time_dim == "valid_time":
        data = data.assign_coords(valid_time=pd.to_datetime(data["valid_time"].values))

    # Debugging: Check available time values
    print(f"Checking available {time_dim} values: {data[time_dim].values}")
    print(f"Checking available {time_dim} months: {data[time_dim].dt.month.values}")

    # Select the months corresponding to the season
    months = SEASONS[season]

    if hasattr(data[time_dim], "dt"):
        seasonal_subset = data.sel({time_dim: data[time_dim].dt.month.isin(months)})
    else:
        print(f"Warning: {time_dim} does not support datetime accessor `dt`. Skipping processing.")
        return

    # Compute seasonal mean
    seasonal_stat = seasonal_subset.mean(dim=time_dim, skipna=True)

    # Convert units if necessary
    if variable in CONVERT_TO_C:
        seasonal_stat -= 273.15  # Convert Kelvin to Celsius

    if variable in CONVERT_TO_MM:
        seasonal_stat *= 1000 * 31  # Convert meters to mm and adjust for days

    print("Value range after unit conversion:", seasonal_stat.min().values, seasonal_stat.max().values)

    # Assign CRS if missing
    if seasonal_stat.rio.crs is None:
        print("Assigning CRS EPSG:4326 to NetCDF data.")
        seasonal_stat = seasonal_stat.rio.write_crs("EPSG:4326")

    # Reproject to EPSG:3338 using nearest-neighbor resampling
    seasonal_stat_3338 = seasonal_stat.rio.reproject(
        "EPSG:3338",
        resolution=(1000, 1000),
        resampling=Resampling.nearest
    )

    print("Value range after reprojection:", seasonal_stat_3338.min().values, seasonal_stat_3338.max().values)

    # Output raster filename
    year = netcdf_file.split("_")[-2]  # Extract year from filename
    output_raster = os.path.join(output_dir, f"{variable}_{year}_{season}.tif")

    # Save raster
    print(f"Saving raster to {output_raster}")
    seasonal_stat_3338.rio.to_raster(
        output_raster,
        nodata=-9999,  # Explicitly set nodata value
        dtype="float32"
    )
    print(f"Raster saved: {output_raster}")

    # Debugging output raster
    from rasterio import open as rio_open
    with rio_open(output_raster) as src:
        print("Output raster info:")
        print("Shape:", src.shape)
        print("CRS:", src.crs)
        print("Bounds:", src.bounds)
        print("Value range:", src.read(1).min(), src.read(1).max())


# Hardcoded parameters for a single run (for troubleshooting)
NETCDF_FILE = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/raw/monthly/2m_temperature_2010_monthly.nc"
VARIABLE = "2m_temperature"
SEASON = "DJF"
OUTPUT_DIR = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/dump"

generate_seasonal_raster(NETCDF_FILE, VARIABLE, SEASON, OUTPUT_DIR)
