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
    "MJJAS": [5, 6, 7, 8, 9],  # Custom season for evaporation
    "ONDJFMAM": [10, 11, 12, 1, 2, 3, 4, 5],  # Custom season for snowfall
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

# Variables to be converted from meters to millimeters
CONVERT_TO_MM = ["total_precipitation", "snowfall", "evaporation"]

def calculate_seasonal_means_yearly(netcdf_dir, variable, shapefile_path, season, output_csv):
    """
    Calculate seasonal and annual statistics for a single variable and append to a CSV.

    Args:
        netcdf_dir (str): Directory containing NetCDF files.
        variable (str): Variable to process.
        shapefile_path (str): Path to the shapefile in EPSG:3338.
        season (str): Season to process (e.g., "DJF", "MAM").
        output_csv (str): Path to save or append the output CSV.
    """
    # Load the shapefile
    gdf = gpd.read_file(shapefile_path)

    # Ensure the shapefile is in EPSG:3338
    if gdf.crs.to_epsg() != 3338:
        gdf = gdf.to_crs(epsg=3338)

    netcdf_variable_name = VARIABLE_MAP.get(variable, variable)
    units = None  # Initialize units variable
    all_stats = []

    # Collect NetCDF files for this variable
    netcdf_files = [
        os.path.join(netcdf_dir, f"{variable}_{year}.nc")
        for year in range(1950, 2024)
    ]

    for netcdf_file in netcdf_files:
        year = int(os.path.basename(netcdf_file).rsplit("_", 1)[-1].split(".")[0])
        print(f"Processing file: {netcdf_file} for year {year}")

        # Load NetCDF and extract the variable
        ds = xr.open_dataset(netcdf_file)

        if netcdf_variable_name not in ds.variables:
            raise KeyError(
                f"Variable '{netcdf_variable_name}' not found in {netcdf_file}. Available variables: {list(ds.variables.keys())}"
            )

        data = ds[netcdf_variable_name]

        # Extract units from metadata
        if units is None:
            units = data.attrs.get("units", "unknown").replace(" ", "_")  # Clean up unit name

        # Convert temperature from Kelvin to Celsius if processing temperature
        if variable == "2m_temperature":
            data -= 273.15
            units = "C"

        # Convert cumulative variables to daily totals and scale to monthly cumulative
        if variable in CONVERT_TO_MM:
            if "m" in data.attrs.get("units", "").lower():
                data *= 1000  # Convert meters to millimeters
                units = "mm"
            elif "mm" in data.attrs.get("units", "").lower():
                units = "mm"
            else:
                raise ValueError(f"Unexpected units for {variable}: {data.attrs.get('units', 'unknown')}")

            # Multiply by 31 to approximate monthly cumulative value
            data *= 31

        # Set CRS for the NetCDF (assumes EPSG:4326 from ERA5)
        data = data.rio.write_crs("EPSG:4326")

        # Reproject to EPSG:3338
        data_3338 = data.rio.reproject("EPSG:3338")

        # Select the season's months
        months = SEASONS[season]
        seasonal_subset = data_3338.sel(valid_time=data_3338["valid_time.month"].isin(months))

        # Calculate statistics
        seasonal_stat = (
            seasonal_subset.mean(dim="valid_time")
            if variable == "2m_temperature"
            else seasonal_subset.sum(dim="valid_time")
        )

        # Perform zonal statistics
        stats = zonal_stats(
            gdf,
            seasonal_stat.values,  # Use the numpy array representation of the xarray DataArray
            affine=seasonal_stat.rio.transform(),  # Get the affine transform from rioxarray
            stats=["mean", "min", "max"],  # Include desired statistics
            all_touched=True,
            geojson_out=True,
            nodata=seasonal_stat.encoding.get("_FillValue", None)  # Use _FillValue if it exists, otherwise None
        )

        # Add time, season, and variable information
        for feature, region_stats in zip(gdf.itertuples(), stats):
            all_stats.append({
                "Region": feature.Region,
                "Year": year,
                "Season": season,
                f"Mean_{variable}_{units}": region_stats["properties"]["mean"],
                f"Min_{variable}_{units}": region_stats["properties"]["min"],
                f"Max_{variable}_{units}": region_stats["properties"]["max"],
            })

    # Convert results to a DataFrame
    results_df = pd.DataFrame(all_stats)

    # Append to or create the output CSV file
    if os.path.exists(output_csv):
        existing_df = pd.read_csv(output_csv)
        results_df = pd.concat([existing_df, results_df], ignore_index=True)

    results_df.to_csv(output_csv, index=False)
    print(f"Results for {variable}, {season} saved to: {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate yearly seasonal and annual statistics.")
    parser.add_argument("--netcdf_dir", required=True, help="Directory containing NetCDF files.")
    parser.add_argument("--variable", required=True, help="Variable to process (e.g., 2m_temperature).")
    parser.add_argument("--shapefile", required=True, help="Path to the shapefile.")
    parser.add_argument("--season", required=True, choices=SEASONS.keys(), help="Season to process (e.g., DJF).")
    parser.add_argument("--output_csv", required=True, help="Path to save or append the output CSV.")
    args = parser.parse_args()

    calculate_seasonal_means_yearly(
        args.netcdf_dir,
        args.variable,
        args.shapefile,
        args.season,
        args.output_csv,
    )
