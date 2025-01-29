import xarray as xr
import rioxarray
from rasterio.enums import Resampling
import os

# Mapping of variable names to NetCDF variable names
VARIABLE_MAP = {
    "2m_temperature": "t2m",
    "total_precipitation": "tp",
}

# Variables to be converted from meters to millimeters
CONVERT_TO_MM = ["total_precipitation"]

def generate_monthly_rasters(netcdf_dir, variable, output_dir, year_range):
    netcdf_variable_name = VARIABLE_MAP.get(variable, variable)
    os.makedirs(output_dir, exist_ok=True)

    years = range(year_range[0], year_range[1] + 1)

    for year in years:
        netcdf_file = os.path.join(netcdf_dir, f"{variable}_{year}.nc")
        if not os.path.exists(netcdf_file):
            print(f"NetCDF file for year {year} not found: {netcdf_file}")
            continue

        print(f"Processing file: {netcdf_file}")
        data = xr.open_dataset(netcdf_file, chunks={"latitude": 100, "longitude": 100})[netcdf_variable_name]

        print("Original Data:")
        print(data)
        print("Value range before processing:", data.min().values, data.max().values)

        # Convert temperature from Kelvin to Celsius
        if variable == "2m_temperature":
            data -= 273.15

        # Convert cumulative variables
        if variable in CONVERT_TO_MM:
            if "m" in data.attrs.get("units", "").lower():
                data *= 1000

        print("Value range after unit conversion:", data.min().values, data.max().values)

        # Handle nodata values
        nodata_value = data.encoding.get("_FillValue", float("nan"))
        data = data.where(data != nodata_value)

        # Set CRS for the NetCDF (assumes EPSG:4326 for ERA5)
        data = data.rio.write_crs("EPSG:4326")

        # Reproject to EPSG:3338 using nearest-neighbor resampling
        data_3338 = data.rio.reproject(
            "EPSG:3338",
            resolution=(1000, 1000),
            resampling=Resampling.nearest
        )

        print("Value range after reprojection:", data_3338.min().values, data_3338.max().values)

        for month in range(1, 13):
            monthly_data = data_3338.sel(valid_time=data_3338["valid_time.month"] == month)

            print(f"Monthly Data for {year}-{month:02d}:")
            print("Shape:", monthly_data.shape)
            print("Value range:", monthly_data.min().values, monthly_data.max().values)

            if monthly_data.size == 0 or monthly_data.isnull().all():
                print(f"No valid data for {year}-{month:02d}")
                continue

            output_raster = os.path.join(output_dir, f"{variable}_{year}_{month:02d}.tif")

            # Write the raster with explicit nodata and data type
            print(f"Saving raster for {year}-{month:02d} to {output_raster}")
            monthly_data.rio.to_raster(
                output_raster,
                nodata=nodata_value or -9999,  # Use source nodata or default to -9999
                dtype="float32"
            )
            print(f"Raster saved: {output_raster}")

            # Debug the written raster
            from rasterio import open as rio_open
            with rio_open(output_raster) as src:
                print("Output raster info:")
                print("Shape:", src.shape)
                print("CRS:", src.crs)
                print("Bounds:", src.bounds)
                print("Value range:", src.read(1).min(), src.read(1).max())


# Hardcoded parameters
NETCDF_DIR = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/raw"
VARIABLE = "total_precipitation"
OUTPUT_DIR = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/processed"
YEAR_RANGE = (2000, 2001)

generate_monthly_rasters(NETCDF_DIR, VARIABLE, OUTPUT_DIR, YEAR_RANGE)
