import xarray as xr
import rioxarray
import argparse
import os

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

def detect_time_coordinate(ds):
    """
    Detect the time coordinate in the dataset.
    """
    possible_time_coords = ["valid_time", "time", "datetime"]
    for coord in possible_time_coords:
        if coord in ds.coords:
            return coord
    raise KeyError("No valid time coordinate found in the dataset.")

def generate_tiff(netcdf_file, variable, year, month, output_tiff):
    """
    Generate a single-band TIFF file for one month of data using the time coordinate.
    """
    print(f"Opening NetCDF file: {netcdf_file}")
    ds = xr.open_dataset(netcdf_file)
    data = ds[variable]

    # Detect the time coordinate
    time_coord = detect_time_coordinate(ds)
    print(f"Using time coordinate: {time_coord}")

    print(f"Filtering data for year {year} and month {month}")
    data_filtered = data.sel({time_coord: data[time_coord].dt.year == year})
    data_filtered = data_filtered.sel({time_coord: data_filtered[time_coord].dt.month == month})

    if data_filtered[time_coord].size == 0:
        raise ValueError(f"No data found for year {year} and month {month}.")

    # Ensure CRS
    data_filtered = ensure_crs(data_filtered)

    # Ensure the output directory exists
    output_dir = os.path.dirname(output_tiff)
    if not os.path.exists(output_dir):
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir)

    # Save the filtered data to a TIFF
    print(f"Saving data to {output_tiff}")
    data_filtered.squeeze().rio.to_raster(output_tiff, compress="LZW")
    print(f"TIFF saved at {output_tiff}")


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
