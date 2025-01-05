import geopandas as gpd
import os

def get_bounding_boxes(shapefile_path):
    """
    Generate bounding boxes in WGS84 and EPSG:3338 for a given shapefile.

    Args:
        shapefile_path (str): Path to the shapefile.

    Returns:
        tuple: Bounding boxes in [North, West, South, East] format:
            - WGS84 (EPSG:4326)
            - EPSG:3338 (meters)
    """
    # Check if the shapefile exists
    if not os.path.exists(shapefile_path):
        raise FileNotFoundError(f"Shapefile not found: {shapefile_path}")

    # Load the shapefile
    gdf = gpd.read_file(shapefile_path)

    # Reproject to WGS84 (EPSG:4326)
    gdf_wgs84 = gdf.to_crs(epsg=4326)
    minx_wgs84, miny_wgs84, maxx_wgs84, maxy_wgs84 = gdf_wgs84.total_bounds
    bbox_wgs84 = [maxy_wgs84, minx_wgs84, miny_wgs84, maxx_wgs84]  # [North, West, South, East]

    # Reproject to EPSG:3338 (Alaska Albers Equal Area)
    gdf_epsg3338 = gdf.to_crs(epsg=3338)
    minx_epsg3338, miny_epsg3338, maxx_epsg3338, maxy_epsg3338 = gdf_epsg3338.total_bounds
    bbox_epsg3338 = [maxy_epsg3338, minx_epsg3338, miny_epsg3338, maxx_epsg3338]  # [North, West, South, East]

    return bbox_wgs84, bbox_epsg3338


if __name__ == "__main__":
    # Ensure the path is correct (absolute or relative)
    shapefile_path = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/resources/shapefiles/BiogeographicRegions_AKSWAP_ABRedited/BiogeographicRegions_AKSWAP_ABRedited.shp"

    # Get the bounding boxes
    try:
        bbox_wgs84, bbox_epsg3338 = get_bounding_boxes(shapefile_path)
        print("Bounding Box for ERA5 API (WGS84):", bbox_wgs84)
        print("Bounding Box (EPSG:3338):", bbox_epsg3338)
    except FileNotFoundError as e:
        print(e)
    except Exception as e:
        print(f"An error occurred: {e}")

