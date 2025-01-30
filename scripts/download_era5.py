import cdsapi
import os
import argparse
import geopandas as gpd
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def get_bounding_boxes(shapefile_path, source_crs="EPSG:3338", target_crs="EPSG:4326"):
    """
    Extracts the bounding box from a shapefile and reprojects it to WGS84 if needed.
    """
    try:
        # Load the shapefile
        gdf = gpd.read_file(shapefile_path)
        logging.info(f"Successfully read input shapefile: {shapefile_path}")

        # Ensure the input CRS matches expected
        if gdf.crs is None:
            raise ValueError(f"Shapefile {shapefile_path} does not have a CRS. Assign a CRS before proceeding.")

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

def download_era5_monthly(variable, start_year, end_year, area, output_dir, custom_params=None):
    """
    Download ERA5 monthly averaged data for a given variable and year range.
    """
    try:
        c = cdsapi.Client()

        for year in range(start_year, end_year + 1):
            file_name = f"{output_dir}/{variable}_{year}_monthly.nc"
            if os.path.exists(file_name):
                logging.info(f"File already exists: {file_name}. Skipping.")
                continue

            # ERA5 request parameters
            params = {
                'product_type': 'monthly_averaged_reanalysis',
                'variable': variable,
                'year': str(year),
                'month': [f"{i:02d}" for i in range(1, 13)],  # All months
                'time': '00:00',
                'format': 'netcdf',
                'area': area,  # [North, West, South, East]
            }

            if custom_params:
                params.update(custom_params)

            logging.info(f"Downloading {variable} for {year} (monthly averages)...")
            c.retrieve('reanalysis-era5-single-levels-monthly-means', params, file_name)
            logging.info(f"Saved to {file_name}")

    except Exception as e:
        logging.error(f"Error in download_era5_monthly: {e}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download ERA5 monthly averages for a specific variable and year range.")
    parser.add_argument("--variable", required=True, help="ERA5 variable to download.")
    parser.add_argument("--start_year", type=int, required=True, help="Start year of data to download.")
    parser.add_argument("--end_year", type=int, required=True, help="End year of data to download.")
    parser.add_argument("--shapefile", required=True, help="Path to the shapefile for bounding box extraction.")
    args = parser.parse_args()

    try:
        # Get the correctly reprojected bounding box
        bounding_box_wgs84 = get_bounding_boxes(args.shapefile)
        logging.info(f"Bounding Box for ERA5 API (WGS84): {bounding_box_wgs84}")

        # Ensure output directory exists
        output_dir = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/raw/monthly"
        os.makedirs(output_dir, exist_ok=True)

        # Download the data for the given year range
        download_era5_monthly(
            variable=args.variable,
            start_year=args.start_year,
            end_year=args.end_year,
            area=bounding_box_wgs84,
            output_dir=output_dir,
        )

    except Exception as e:
        logging.error(f"Script failed: {e}")
