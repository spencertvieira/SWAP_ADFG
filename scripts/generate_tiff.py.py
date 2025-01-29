import xarray as xr
import rioxarray
import argparse

def ensure_crs(data, target_crs="EPSG:3338"):
    """
    Ensures the CRS of the dataset is aligned with the target CRS.
    If CRS is missing, assigns EPSG:4326 and reprojects to the target CRS.

    Args:
        data (xr.DataArray): The input data array.
        target_crs (str): The target CRS to ensure.

    Returns:
        xr.DataArray: The reprojected dataset.
    """
    if data.rio.crs is None:
        print("Assigning CRS EPSG:4326 to NetCDF data.")
        data = data.rio.write_crs("EPSG:4326")
    
    if data.rio.crs != target_crs:
        print(f"Reprojecting NetCDF data to {target_crs}.")
        data = data.rio.reproject(target_crs)
    
    return data

def generate_tiff(netcdf_file, variable, year, month, output_tiff):
    """
    Generate a single-band GeoTIFF for a specific variable, month, and year.

    Args:
        netcdf_file (str): Path to the NetCDF file.
        variable (str): Variable name to extract (e.g., "tp").
        year (int): The year to extract data for.
        month (int): The month to extract data for.
        output_tiff (str): Path to save the output GeoTIFF.
    """
    print(f"Opening NetCDF file: {netcdf_file}")
    ds = xr.open_dataset(netcdf_file)
    data = ds[variable]

    print(f"Filtering data for year {year} and month {month}")
    data_filtered = data.sel(time=data["time"].dt.year == year)
    data_filtered = data_filtered.sel(time=data_filtered["time"].dt.month == month)

    if data_filtered.time.size == 0:
        raise ValueError(f"No data found for year {year} and month {month} in {netcdf_file}.")

    # Calculate the average for the specified month across time (if needed)
    print("Calculating the average for the selected month.")
    data_avg = data_filtered.mean(dim="time", skipna=True)

    # Ensure CRS and reproject to EPSG:3338
    print("Ensuring CRS alignment.")
    data_avg = ensure_crs(data_avg)

    # Save to GeoTIFF
    print(f"Saving data to GeoTIFF: {output_tiff}")
    data_avg.rio.to_raster(output_tiff)

    print(f"GeoTIFF saved at {output_tiff}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a single-band GeoTIFF for a specific variable, month, and year.")
    parser.add_argument("--netcdf_file", required=True, help="Path to the NetCDF file.")
    parser.add_argument("--variable", required=True, help="Variable to extract (e.g., 'tp').")
    parser.add_argument("--year", required=True, type=int, help="Year to filter data.")
    parser.add_argument("--month", required=True, type=int, help="Month to filter data (1-12).")
    parser.add_argument("--output_tiff", required=True, help="Path to save the output GeoTIFF.")
    args = parser.parse_args()

    generate_tiff(
        args.netcdf_file,
        args.variable,
        args.year,
        args.month,
        args.output_tiff,
    )
