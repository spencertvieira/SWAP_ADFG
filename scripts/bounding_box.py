import geopandas as gpd
from shapely.geometry import box
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Define input and output file paths
INPUT_SHP = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/resources/BiogeographicRegions_AKSWAP_ABRedited_edited_clipped_again_again_again_again_again_again.shp"
OUTPUT_SHP = "/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/resources/bounding_box.shp"

def save_bounding_box(input_shapefile, output_shapefile, expected_crs="EPSG:3338"):
    """
    Extracts and saves the bounding box of a shapefile as a new shapefile.
    Keeps the bounding box in the original CRS (EPSG:3338).

    Args:
        input_shapefile (str): Path to the input shapefile.
        output_shapefile (str): Path to save the output bounding box shapefile.
        expected_crs (str): Expected CRS of the input shapefile (default: EPSG:3338).
    """
    try:
        # Read the input shapefile
        gdf = gpd.read_file(input_shapefile)
        logging.info(f"Successfully read input shapefile: {input_shapefile}")

        # Check if the input shapefile has a CRS
        if gdf.crs is None:
            raise ValueError("Shapefile does not have a CRS. Please assign one before proceeding.")

        # Log the CRS of the input shapefile
        logging.info(f"Input shapefile CRS: {gdf.crs}")

        # Ensure input is in expected CRS
        if gdf.crs.to_string() != expected_crs:
            raise ValueError(f"Unexpected CRS {gdf.crs}. Expected {expected_crs}. Check input file.")

        # Extract bounding box
        minx, miny, maxx, maxy = gdf.total_bounds
        logging.info(f"Bounding box extracted: minx={minx}, miny={miny}, maxx={maxx}, maxy={maxy}")

        # Check if the bounding box coordinates are reasonable for Alaska in EPSG:3338
        if not (-2500000 <= minx <= 1600000 and 300000 <= miny <= 2800000):
            logging.warning("Bounding box coordinates are outside the expected range for Alaska in EPSG:3338.")
            logging.warning("Expected approximate range for Alaska in EPSG:3338:")
            logging.warning("minx: -2,500,000 to 1,600,000")
            logging.warning("miny: 300,000 to 2,800,000")

        # Create bounding box geometry
        bbox_geom = gpd.GeoDataFrame(geometry=[box(minx, miny, maxx, maxy)], crs=gdf.crs)

        # Save the bounding box in the original CRS
        bbox_geom.to_file(output_shapefile)
        logging.info(f"Bounding box shapefile saved to: {output_shapefile}")

        # Provide instructions for visualizing the bounding box
        logging.info("\nTo visualize the bounding box:")
        logging.info(f"1. Open {output_shapefile} in QGIS or another GIS tool.")
        logging.info("2. Set the project CRS to EPSG:3338 (Project > Properties > CRS).")
        logging.info("3. Overlay the bounding box on a base map of Alaska to verify alignment.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    save_bounding_box(INPUT_SHP, OUTPUT_SHP)