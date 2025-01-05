import cdsapi
import os
import argparse
from bounding_box import get_bounding_boxes

def download_era5(variable, years, area, output_dir, custom_params=None):
    """
    Download ERA5 reanalysis data for a given variable and year range.

    Args:
        variable (str): ERA5 variable name.
        years (list): List of years to download data for.
        area (list): Bounding box in [North, West, South, East] format.
        output_dir (str): Directory to save downloaded files.
        custom_params (dict, optional): Custom parameters for the CDS API.
    """
    c = cdsapi.Client()

    for year in years:
        file_name = f"{output_dir}/{variable}_{year}.nc"
        if os.path.exists(file_name):
            print(f"File already exists: {file_name}. Skipping download.")
            continue

        params = {
            'product_type': 'monthly_averaged_reanalysis',
            'variable': variable,
            'year': year,
            'month': [f"{i:02d}" for i in range(1, 13)],
            'time': '00:00',
            'format': 'netcdf',
            'area': area,  # [North, West, South, East]
        }

        if custom_params:
            params.update(custom_params)

        print(f"Downloading {variable} for {year}...")
        c.retrieve('reanalysis-era5-single-levels-monthly-means', params, file_name)
        print(f"Saved to {file_name}")


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Download ERA5 data for a specific variable and year range.")
    parser.add_argument("--variable", required=True, help="Climate variable to download.")
    parser.add_argument("--start_year", type=int, required=True, help="Start year of the data range.")
    parser.add_argument("--end_year", type=int, required=True, help="End year of the data range.")
    args = parser.parse_args()

    # Define parameters
    shapefile_path = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/resources/shapefiles/BiogeographicRegions_AKSWAP_ABRedited/BiogeographicRegions_AKSWAP_ABRedited.shp"

    # Get the bounding boxes
    bounding_box_wgs84, _ = get_bounding_boxes(shapefile_path)
    print("Bounding Box for ERA5 API (WGS84):", bounding_box_wgs84)

    # Output directory
    output_dir = "../data/raw/"
    os.makedirs(output_dir, exist_ok=True)

    # Download data
    download_era5(
        variable=args.variable,
        years=range(args.start_year, args.end_year + 1),
        area=bounding_box_wgs84,
        output_dir=output_dir,
    )

