import cdsapi
import os
import argparse
import geopandas as gpd
import logging
import xarray as xr

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Variables that should be summed instead of averaged
SUM_VARIABLES = ["total_precipitation", "snowfall", "evaporation"]

def get_bounding_boxes(shapefile_path, source_crs="EPSG:3338", target_crs="EPSG:4326"):
    """
    Extracts the bounding box from a shapefile and reprojects it to WGS84 if needed.
    """
    try:
        gdf = gpd.read_file(shapefile_path)
        logging.info(f"Successfully read input shapefile: {shapefile_path}")

        if gdf.crs is None:
            raise ValueError(f"Shapefile {shapefile_path} does not have a CRS. Assign one before proceeding.")

        # Log the input CRS
        logging.info(f"Input shapefile CRS: {gdf.crs}")

        # Reproject to WGS84 if needed
        if gdf.crs.to_string() != target_crs:
            logging.info(f"Reprojecting bounding box from {gdf.crs} to {target_crs}")
            gdf = gdf.to_crs(target_crs)

        # Extract bounding box in lat/lon format (minx, miny, maxx, maxy)
        minx, miny, maxx, maxy = gdf.total_bounds
        north, west, south, east = maxy, minx, miny, maxx  # Reorder for ERA5 API

        logging.info(f"Bounding box extracted: North={north}, West={west}, South={south}, East={east}")
        return [north, west, south, east]

    except Exception as e:
        logging.error(f"Error in get_bounding_boxes: {e}")
        raise

def download_era5_hourly(variable, start_year, end_year, area, hourly_dir, daily_dir, custom_params=None):
    """
    Download ERA5 hourly reanalysis data, store as yearly NetCDF files, and compute daily averages.
    """
    try:
        c = cdsapi.Client()

        for year in range(start_year, end_year + 1):
            hourly_file = f"{hourly_dir}/{variable}_{year}_hourly.nc"
            daily_file = f"{daily_dir}/{variable}_{year}_daily.nc"

            if os.path.exists(hourly_file) and os.path.exists(daily_file):
                logging.info(f"Hourly and daily files already exist for {variable} {year}. Skipping.")
                continue

            if not os.path.exists(hourly_file):
                logging.info(f"Downloading hourly ERA5 data for {variable}, {year}...")

                # ERA5 request parameters
                params = {
                    'product_type': 'reanalysis',
                    'variable': variable,
                    'year': str(year),
                    'month': [f"{i:02d}" for i in range(1, 13)],  # All months
                    'day': [f"{i:02d}" for i in range(1, 32)],  # All days
                    'time': [f"{i:02d}:00" for i in range(0, 24)],  # Hourly data
                    'format': 'netcdf',
                    'area': area,  # [North, West, South, East]
                }

                if custom_params:
                    params.update(custom_params)

                c.retrieve('reanalysis-era5-single-levels', params, hourly_file)
                logging.info(f"Saved hourly file: {hourly_file}")

            # Convert hourly data to daily averages
            if not os.path.exists(daily_file):
                convert_to_daily(hourly_file, daily_file, variable)

    except Exception as e:
        logging.error(f"Error in download_era5_hourly: {e}")
        raise

def convert_to_daily(hourly_file, daily_file, variable):
    """
    Converts hourly ERA5 NetCDF data to daily averages or sums (for precip, snowfall, evap).
    """
    try:
        logging.info(f"Converting {hourly_file} to daily data...")

        # Load dataset
        ds = xr.open_dataset(hourly_file)

        # Convert hourly time series to daily data
        if variable in SUM_VARIABLES:
            ds_daily = ds.resample(valid_time="1D").sum()  # Sum for precipitation, snowfall, evaporation
            logging.info(f"Applying SUM for {variable}")
        else:
            ds_daily = ds.resample(valid_time="1D").mean()  # Mean for other variables
            logging.info(f"Applying MEAN for {variable}")

        # Save as daily NetCDF file
        ds_daily.to_netcdf(daily_file, mode="w", format="NETCDF4")

        logging.info(f"Saved daily file: {daily_file}")

        # Close datasets
        ds.close()
        ds_daily.close()

    except Exception as e:
        logging.error(f"Error converting hourly to daily: {e}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download ERA5 hourly reanalysis, compute daily averages, and store in NetCDF.")
    parser.add_argument("--variable", required=True, help="ERA5 variable to download.")
    parser.add_argument("--start_year", type=int, required=True, help="Start year of data to download.")
    parser.add_argument("--end_year", type=int, required=True, help="End year of data to download.")
    parser.add_argument("--shapefile", required=True, help="Path to the shapefile for bounding box extraction.")
    args = parser.parse_args()

    try:
        # Get the correctly reprojected bounding box
        bounding_box_wgs84 = get_bounding_boxes(args.shapefile)
        logging.info(f"Bounding Box for ERA5 API (WGS84): {bounding_box_wgs84}")

        # Ensure output directories exist
        hourly_dir = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/raw/hourly"
        daily_dir = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/raw/daily"
        os.makedirs(hourly_dir, exist_ok=True)
        os.makedirs(daily_dir, exist_ok=True)

        # Download ERA5 hourly data and compute daily averages
        download_era5_hourly(
            variable=args.variable,
            start_year=args.start_year,
            end_year=args.end_year,
            area=bounding_box_wgs84,
            hourly_dir=hourly_dir,
            daily_dir=daily_dir
        )

    except Exception as e:
        logging.error(f"Script failed: {e}")
