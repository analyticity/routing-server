# graph.py

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
        id = row.get("id", "")
        name = row.get("name", "")
        highway = row.get(
            "highway", "unclassified"
        )  # Type of road, default to "unclassified"
        oneway = row.get("oneway", "no")  # Default to "no"
        if oneway != "yes":
            oneway = "no"
        speed = speeds.get(highway, 8.33)  # Default to 30 km/h

        if isinstance(geometry, (LineString, MultiLineString)):
            coords = list(geometry.coords)
            for i in range(len(coords) - 1):
                start = coords[i]
                end = coords[i + 1]

                length = geodesic(
                    (start[1], start[0]), (end[1], end[0])
                ).meters  # geodesic takes (lat, lon)
                traversal_time = (
                    length / speed if speed > 0 else length
                )  # Avoid division by zero

                # Add edges to the graph
                graph.add_edge(
                    Point(start),
                    Point(end),
                    geometry=LineString([Point(start), Point(end)]),
                    id=id,
                    name=name,
                    highway=highway,
                    length=length,
                    traversal_time=traversal_time,
                    oneway=oneway,
                )
                # Add reverse edge if not one-way
                if oneway == "no":
                    graph.add_edge(
                        Point(end),
                        Point(start),
                        geometry=LineString([Point(end), Point(start)]),
                        id=id,
                        name=name,
                        highway=highway,
                        length=length,
                        traversal_time=traversal_time,
                        oneway=oneway,
                    )

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
    routing_base = gpd.GeoDataFrame.from_features(osm_data["features"], crs="OGC:CRS84")
    tags_to_extract = ["name", "highway", "oneway"]
    for tag in tags_to_extract:
        routing_base[tag] = routing_base["tags"].apply(lambda x: x.get(tag, ""))
    
    # Debug for visualization on geojson.io
    routing_base["stroke"] = routing_base.apply(
        lambda x: "#{:06x}".format(hash(x["geometry"]) % 0xFFFFFF), axis=1
    )
    routing_base["stroke-width"] = 10
    
    routing_base = routing_base.drop(columns=["tags"])

    # Save the routing base to a GeoJSON file
    routing_base.to_file("data/sample_routing_base.geojson", driver="GeoJSON")

    return routing_base
