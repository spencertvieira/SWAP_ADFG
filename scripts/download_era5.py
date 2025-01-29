import cdsapi
import os
import argparse
import geopandas as gpd

def get_bounding_boxes(shapefile_path, target_crs="EPSG:4326"):
    """
    Extracts the bounding box from a shapefile and reprojects it to WGS84 if needed.
    """
    gdf = gpd.read_file(shapefile_path)

    if gdf.crs.to_string() != target_crs:
        gdf = gdf.to_crs(target_crs)

    bounds = gdf.total_bounds  # [minx, miny, maxx, maxy]
    north, west, south, east = bounds[3], bounds[0], bounds[1], bounds[2]
    
    return [north, west, south, east]

def download_era5_monthly(variable, years, area, output_dir, custom_params=None):
    """
    Download ERA5 monthly averaged data for a given variable and year range.
    """
    c = cdsapi.Client()

    for year in years:
        file_name = f"{output_dir}/{variable}_{year}_monthly.nc"
        if os.path.exists(file_name):
            print(f"File already exists: {file_name}. Skipping.")
            continue
        
        #monthly
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

        print(f"Downloading {variable} for {year} (monthly averages)...")
        c.retrieve('reanalysis-era5-single-levels-monthly-means', params, file_name)
        print(f"Saved to {file_name}")

        #daily
        # params = {
        #     'product_type': 'reanalysis',
        #     'variable': variable,  # e.g., 'total_precipitation'
        #     'year': year,
        #     'month': [f"{i:02d}" for i in range(1, 13)],  # All months
        #     'day': [f"{i:02d}" for i in range(1, 32)],    # All days
        #     'time': '00:00',  # One time step per day for daily totals
        #     'format': 'netcdf',
        #     'area': area,  # [North, West, South, East]
        # }
        # if custom_params:
        #     params.update(custom_params)

        # print(f"Downloading {variable} for {year}...")
        # c.retrieve('reanalysis-era5-single-levels', params, file_name)
        # print(f"Saved to {file_name}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download ERA5 monthly averages for a specific variable and year range.")
    parser.add_argument("--variable", required=True, help="ERA5 variable to download.")
    parser.add_argument("--start_year", type=int, required=True, help="Start year of data.")
    parser.add_argument("--end_year", type=int, required=True, help="End year of data.")
    parser.add_argument("--shapefile", required=True, help="Path to the shapefile for bounding box extraction.")
    args = parser.parse_args()

    bounding_box_wgs84 = get_bounding_boxes(args.shapefile)
    print("Bounding Box for ERA5 API (WGS84):", bounding_box_wgs84)

    output_dir = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/raw/monthly/monthly_era5_AK"
    os.makedirs(output_dir, exist_ok=True)

    download_era5_monthly(
        variable=args.variable,
        years=range(args.start_year, args.end_year + 1),
        area=bounding_box_wgs84,
        output_dir=output_dir,
    )
