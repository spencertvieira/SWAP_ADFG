import xarray as xr
import geopandas as gpd
import rioxarray
from rasterstats import zonal_stats
import pandas as pd
import os
import argparse

# Define seasons
SEASONS = {
    "DJF": [12, 1, 2],
    "MAM": [3, 4, 5],
    "JJA": [6, 7, 8],
    "SON": [9, 10, 11],
    "Annual": list(range(1, 13)),  # All months
    "Oct-May": [10, 11, 12, 1, 2, 3, 4, 5],  # Snowfall
    "May-Sep": [5, 6, 7, 8, 9],  # Evaporation
}

# Mapping of human-readable variable names to NetCDF variable names
VARIABLE_MAP = {
    "2m_temperature": "t2m",
    "total_precipitation": "tp",
    "snowfall": "sf",
    "evaporation": "e",
    "10m_u_component_of_wind": "u10",
    "10m_v_component_of_wind": "v10",
}

def calculate_seasonal_means(netcdf_files, shapefile_path, start_year, end_year, variable_type, output_csv):
    """
    Calculate seasonal and annual statistics for polygons in a shapefile.

    Args:
        netcdf_files (list): List of NetCDF file paths.
        shapefile_path (str): Path to the shapefile in EPSG:3338.
        start_year (int): Start year of the period.
        end_year (int): End year of the period.
        variable_type (str): Type of variable ("mean" or "sum").
        output_csv (str): Path to save the output CSV.
    """
    # Load the shapefile
    gdf = gpd.read_file(shapefile_path)

    # Ensure the shapefile is in EPSG:3338
    if gdf.crs.to_epsg() != 3338:
        gdf = gdf.to_crs(epsg=3338)

    # Initialize a results storage
    all_stats = []
    seasonal_data = {season: None for season in SEASONS.keys()}
    units = None  # Initialize units variable

    for netcdf_file in netcdf_files:
        # Extract variable name and year
        filename = os.path.basename(netcdf_file)
        variable_name, year_with_ext = filename.rsplit("_", 1)
        year = int(year_with_ext.split(".")[0])  # Extract the year

        # Map the human-readable variable name to the NetCDF variable name
        netcdf_variable_name = VARIABLE_MAP.get(variable_name, variable_name)

        if start_year <= year <= end_year:
            print(f"Processing file: {netcdf_file} for year {year}")

            # Load NetCDF and extract the variable
            ds = xr.open_dataset(netcdf_file)

            if netcdf_variable_name not in ds.variables:
                raise KeyError(f"Variable '{netcdf_variable_name}' not found in {netcdf_file}. Available variables: {list(ds.variables.keys())}")

            data = ds[netcdf_variable_name]

            # Extract units from metadata
            if units is None:
                units = data.attrs.get("units", "unknown")

            # Set CRS for the NetCDF (assumes EPSG:4326 from ERA5)
            data = data.rio.write_crs("EPSG:4326")

            # Reproject to EPSG:3338
            data_3338 = data.rio.reproject("EPSG:3338")

            # Group data by season and aggregate
            for season, months in SEASONS.items():
                seasonal_subset = data_3338.sel(valid_time=data_3338["valid_time.month"].isin(months))
                if variable_type == "mean":
                    seasonal_mean = seasonal_subset.mean(dim="valid_time")
                elif variable_type == "sum":
                    seasonal_mean = seasonal_subset.sum(dim="valid_time")
                
                # Accumulate data for multi-year averaging
                if seasonal_data[season] is None:
                    seasonal_data[season] = seasonal_mean
                else:
                    seasonal_data[season] += seasonal_mean

    # Average over the years in the period
    for season in SEASONS.keys():
        seasonal_data[season] /= (end_year - start_year + 1)

        # Save the seasonal raster to a temporary GeoTIFF
        temp_raster = f"temp_raster_{season}.tif"
        seasonal_data[season].rio.to_raster(temp_raster)

        # Perform zonal statistics
        stats = zonal_stats(
            gdf,
            temp_raster,
            stats=["mean"],  # Weighted average
            all_touched=True,
            geojson_out=True,
        )

        # Add time, period, and units information
        for feature, region_stats in zip(gdf.itertuples(), stats):
            all_stats.append({
                "Region": feature.Region,  # Assuming region name is in the "Region" column
                "Period": f"{start_year}-{end_year}",
                "Season": season,
                "Mean": region_stats["properties"]["mean"],
                "Units": units
            })

        # Cleanup temporary files
        if os.path.exists(temp_raster):
            os.remove(temp_raster)

    # Convert results to a DataFrame
    results_df = pd.DataFrame(all_stats)

    # Save to CSV
    results_df.to_csv(output_csv, index=False)
    print(f"Seasonal and annual statistics saved to: {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate seasonal and annual statistics.")
    parser.add_argument("--variable", required=True, help="Variable to process (e.g., 2m_temperature).")
    parser.add_argument("--start_year", type=int, required=True, help="Start year of the period.")
    parser.add_argument("--end_year", type=int, required=True, help="End year of the period.")
    parser.add_argument("--shapefile", required=True, help="Path to the shapefile.")
    parser.add_argument("--output_csv", required=True, help="Path to save the output CSV.")
    parser.add_argument("--variable_type", choices=["mean", "sum"], required=True, help="Type of variable: mean or sum.")
    args = parser.parse_args()

    # Input NetCDF files
    netcdf_files = [f"/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/raw/{args.variable}_{year}.nc" for year in range(args.start_year, args.end_year + 1)]

    # Calculate seasonal means
    calculate_seasonal_means(
        netcdf_files,
        args.shapefile,
        args.start_year,
        args.end_year,
        args.variable_type,
        args.output_csv
    )

