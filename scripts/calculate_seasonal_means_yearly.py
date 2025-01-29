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
    "Annual": list(range(1, 13)),
    "MJJAS": [5, 6, 7, 8, 9],
    "ONDJFMAM": [10, 11, 12, 1, 2, 3, 4, 5],
}

# Valid variable-season combinations
SEASON_VARIABLE_COMBINATIONS = {
    "Annual": ["2m_temperature", "total_precipitation", "10m_u_component_of_wind", "10m_v_component_of_wind", "sea_surface_temperature"],
    "DJF": ["2m_temperature", "total_precipitation", "10m_u_component_of_wind", "10m_v_component_of_wind", "sea_surface_temperature"],
    "MAM": ["2m_temperature", "total_precipitation", "10m_u_component_of_wind", "10m_v_component_of_wind", "sea_surface_temperature"],
    "JJA": ["2m_temperature", "total_precipitation", "10m_u_component_of_wind", "10m_v_component_of_wind", "sea_surface_temperature"],
    "SON": ["2m_temperature", "total_precipitation", "10m_u_component_of_wind", "10m_v_component_of_wind", "sea_surface_temperature"],
    "MJJAS": ["evaporation"],
    "ONDJFMAM": ["snowfall"],
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

def ensure_crs_alignment(data, gdf, target_crs="EPSG:3338"):
    """Ensures CRS alignment between a NetCDF dataset and a GeoDataFrame."""
    if data.rio.crs is None:
        print("Assigning CRS EPSG:4326 to NetCDF data.")
        data = data.rio.write_crs("EPSG:4326")
    
    if data.rio.crs != target_crs:
        print(f"Reprojecting NetCDF data to {target_crs}.")
        data = data.rio.reproject(target_crs)

    if gdf.crs is None:
        raise ValueError("Shapefile does not have a CRS. Please assign one.")
    if gdf.crs.to_epsg() != int(target_crs.split(":")[1]):
        print(f"Reprojecting shapefile to {target_crs}.")
        gdf = gdf.to_crs(target_crs)

    return data, gdf

def calculate_seasonal_means_yearly(netcdf_dir, variable, shapefile_path, season, output_csv):
    """Calculates seasonal and annual statistics for a given variable and stores results in CSV."""
    if variable not in SEASON_VARIABLE_COMBINATIONS.get(season, []):
        print(f"Skipping: Variable '{variable}' is not valid for season '{season}'.")
        return

    gdf = gpd.read_file(shapefile_path)
    netcdf_variable_name = VARIABLE_MAP.get(variable, variable)
    units = None
    all_stats = []

    netcdf_files = sorted([
        f for f in os.listdir(netcdf_dir)
        if f.startswith(f"{variable}_") and f.endswith("_monthly.nc")
    ])

    for netcdf_file in netcdf_files:
        try:
            year = int(netcdf_file.split("_")[-2])
        except ValueError:
            print(f"Skipping invalid filename format: {netcdf_file}")
            continue

        print(f"Processing file: {netcdf_file} for year {year}")
        ds = xr.open_dataset(os.path.join(netcdf_dir, netcdf_file))
        data = ds[netcdf_variable_name]

        # Ensure CRS alignment
        data, gdf = ensure_crs_alignment(data, gdf)

        # Convert missing values to NaN
        data = xr.where(data != data.attrs.get("_FillValue", -9999), data, float("nan"))

        if units is None:
            units = data.attrs.get("units", "unknown").replace(" ", "_")

        # Convert units
        if variable in ["evaporation", "snowfall", "total_precipitation"]:
            data *= 1000 * 31  # Convert meters to mm and adjust for days
            units = "mm"
        elif variable in ["2m_temperature", "sea_surface_temperature"]:
            data -= 273.15  # Convert Kelvin to Celsius
            units = "C"

        # Detect the correct time dimension
        time_dim = None
        for possible_time in ["valid_time", "time"]:
            if possible_time in data.coords:
                time_dim = possible_time
                break

        if time_dim is None:
            raise KeyError("No valid time coordinate found in the dataset.")

        # Ensure valid_time is treated as a datetime index
        if time_dim == "valid_time":
            data["valid_time"] = pd.to_datetime(data["valid_time"].values)

        # Select months correctly
        months = SEASONS[season]
        seasonal_subset = data.sel({time_dim: data[time_dim].dt.month.isin(months)})
        seasonal_stat = seasonal_subset.mean(dim=time_dim, skipna=True)

        # Perform zonal statistics
        stats = zonal_stats(
            gdf, seasonal_stat.values, affine=seasonal_stat.rio.transform(),
            stats=["mean", "min", "max"], all_touched=True, geojson_out=True,
            nodata=float("nan")
        )

        # Append statistics for each region
        for feature, region_stats in zip(gdf.itertuples(), stats):
            all_stats.append({
                "Region": feature.Region, "Year": year, "Season": season,
                f"Mean_{variable}": region_stats["properties"]["mean"],
                f"Min_{variable}": region_stats["properties"]["min"],
                f"Max_{variable}": region_stats["properties"]["max"],
            })

    # Save results to CSV
    pd.DataFrame(all_stats).to_csv(output_csv, index=False)
    print(f"Results saved to: {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate seasonal statistics from NetCDF data.")
    parser.add_argument("--netcdf_dir", required=True, help="Directory containing NetCDF files.")
    parser.add_argument("--variable", required=True, help="Variable to process.")
    parser.add_argument("--shapefile", required=True, help="Path to shapefile.")
    parser.add_argument("--season", required=True, choices=SEASONS.keys(), help="Season to process.")
    parser.add_argument("--output_csv", required=True, help="CSV output file.")
    args = parser.parse_args()

    calculate_seasonal_means_yearly(args.netcdf_dir, args.variable, args.shapefile, args.season, args.output_csv)
