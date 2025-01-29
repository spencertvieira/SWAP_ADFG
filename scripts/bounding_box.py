import geopandas as gpd
import os
from shapely.geometry import box

def save_bounding_box_as_shapefile(shapefile_path, output_dir):
    """
    Save the bounding box of a shapefile as a new shapefile.

    Args:
        shapefile_path (str): Path to the input shapefile.
        output_dir (str): Directory where the bounding box shapefile will be saved.
    """
    # Check if the shapefile exists
    if not os.path.exists(shapefile_path):
        raise FileNotFoundError(f"Shapefile not found: {shapefile_path}")

    # Load the shapefile
    gdf = gpd.read_file(shapefile_path)

    # Calculate the bounding box
    minx, miny, maxx, maxy = gdf.total_bounds  # [min_x, min_y, max_x, max_y]
    bbox_geom = box(minx, miny, maxx, maxy)  # Create a shapely box geometry

    # Create a GeoDataFrame for the bounding box
    bbox_gdf = gpd.GeoDataFrame({"geometry": [bbox_geom]}, crs=gdf.crs)

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Save the bounding box as a shapefile
    bbox_shapefile_path = os.path.join(output_dir, "bounding_box.shp")
    bbox_gdf.to_file(bbox_shapefile_path)

    print(f"Bounding box shapefile saved to: {bbox_shapefile_path}")


if __name__ == "__main__":
    # Input shapefile and output directory
    shapefile_path = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/resources/BiogeographicRegions_AKSWAP_ABRedited_edited_clipped_again_again_again_again_again_again.shp"
    output_dir = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/resources"

    # Save the bounding box as a shapefile
    try:
        save_bounding_box_as_shapefile(shapefile_path, output_dir)
    except FileNotFoundError as e:
        print(e)
    except Exception as e:
        print(f"An error occurred: {e}")