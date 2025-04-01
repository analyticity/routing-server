import geopandas as gpd
import networkx as nx
from geopy.distance import geodesic
from shapely.geometry import LineString, MultiLineString, Point

from osm import osm_data_for_area


def create_graph_from_base(routing_base: gpd.GeoDataFrame) -> nx.MultiDiGraph:
    """
    Create a directed graph from the routing base GeoDataFrame.

    Args:
        routing_base (gpd.GeoDataFrame): GeoDataFrame containing the routing base data.

    Returns:
        nx.MultiDiGraph: Directed graph representing the routing base.
    """
    graph = nx.MultiDiGraph()

    # Define speed limits for different road types in m/s
    speeds = {
        "motorway": 27.78,  # 100 km/h
        "motorway_link": 22.22,  # 80 km/h
        "trunk": 22.22,  # 80 km/h
        "trunk_link": 16.67,  # 60 km/h
        "primary": 16.67,  # 60 km/h
        "primary_link": 13.89,  # 50 km/h
        "secondary": 13.89,  # 50 km/h
        "secondary_link": 11.11,  # 40 km/h
        "tertiary": 11.11,  # 40 km/h
        "tertiary_link": 8.33,  # 30 km/h
        "residential": 8.33,  # 30 km/h
        "living_street": 5.56,  # 20 km/h
        "unclassified": 8.33,  # 30 km/h
        "service": 5.56,  # 20 km/h
    }

    for _, row in routing_base.iterrows():
        geometry = row.geometry
        tags = row.get("tags", {})

        road_type = tags.get("highway", "unclassified")  # Default to unclassified
        speed = speeds.get(road_type, 8.33)  # Default to 30 km/h

        if isinstance(geometry, (LineString, MultiLineString)):
            coords = list(geometry.coords)
            for i in range(len(coords) - 1):
                start = coords[i]
                end = coords[i + 1]

                distance = geodesic(
                    (start[1], start[0]), (end[1], end[0])
                ).meters  # geodesic takes (lat, lon)
                weight = (
                    distance / speed if speed > 0 else distance
                )  # Avoid division by zero

                # Add edges to the graph
                graph.add_edge(Point(start), Point(end), weight=weight)
                # Add reverse edge if not one-way
                if "oneway" not in tags or tags["oneway"] == "no":
                    graph.add_edge(Point(end), Point(start), weight=weight)

    return graph


def get_routing_base(area: str) -> gpd.GeoDataFrame:
    """
    Retrieve current map data for a given area from OpenStreetMap and convert it to a GeoDataFrame.

    Args:
        area (str): Name of the area to query, e.g., 'Brno', 'South Moravian Region', 'Czech Republic'.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame containing the routing base data.
    """
    osm_data = osm_data_for_area(area)
    routing_base = gpd.GeoDataFrame.from_features(osm_data["features"])
    return routing_base
