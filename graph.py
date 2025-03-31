import geopandas as gpd
import networkx as nx
from geopy.distance import geodesic
from shapely.geometry import LineString, MultiLineString

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

    for _, row in routing_base.iterrows():
        geometry = row.geometry
        tags = row.get("tags", {})
        if isinstance(geometry, (LineString, MultiLineString)):
            coords = list(geometry.coords)
            for i in range(len(coords) - 1):
                start = coords[i]
                end = coords[i + 1]
                graph.add_edge(start, end, weight=geodesic(start, end).meters)
                if "oneway" not in tags or tags["oneway"] == "no":
                    graph.add_edge(end, start, weight=geodesic(end, start).meters)

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
