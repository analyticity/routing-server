import geopandas as gpd
import pandas as pd


def load_traffic_data(path: str) -> gpd.GeoDataFrame:
    """
    Load traffic data from a provided GeoJSON file, perform preprocessing, convert CRS and return a GeoDataFrame.

    Args:
        path (str): Path to the GeoJSON file containing traffic data.

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame containing the processed traffic data.
    """
    # Load the traffic data from a GeoJSON file
    traffic_data = gpd.read_file(path)

    # Check if the data is in EPSG:4326
    if traffic_data.crs is None:
        raise ValueError("No CRS found in the file")
    # Convert coordinates to CRS EPSG:32633
    traffic_data.to_crs("EPSG:32633", inplace=True)

    # Convert pubMillis to datetime
    traffic_data["started"] = pd.to_datetime(traffic_data["pubMillis"], unit="ms")
    traffic_data["lastupdated"] = pd.to_datetime(traffic_data["lastupdated"], unit="ms")
    traffic_data["date"] = traffic_data["started"].dt.date

    # If finished, calculate duration
    traffic_data["finished"] = traffic_data["finished"].astype(bool)
    traffic_data.loc[traffic_data["finished"], "duration"] = (
        traffic_data["lastupdated"] - traffic_data["started"]
    ).dt.total_seconds()
    # If not, default to 1 hour
    traffic_data.loc[~traffic_data["finished"], "duration"] = 3600
    # If longer than 4 hours, set to 4 hours
    traffic_data.loc[traffic_data["duration"] > 14400, "duration"] = 14400

    columns_to_keep = [
        "id",
        "country",
        "city",
        "street",
        "date",
        "started",
        "lastupdated",
        "duration",
        "delay",
        "speed",
        "geometry",
    ]

    for col in traffic_data.columns:
        if col not in columns_to_keep:
            traffic_data.drop(columns=col, inplace=True)

    # Save traffic data to a GeoJSON file
    traffic_data.to_file("data/processed_traffic.geojson", driver="GeoJSON")

    return traffic_data
