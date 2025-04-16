# graph.py

from typing import Tuple

import geopandas as gpd
import networkx as nx
from shapely.geometry import LineString, MultiLineString, Point

from osm import osm_data_for_area


def find_projected_point(graph: nx.MultiDiGraph, point: Point) -> Tuple[Point, dict]:
    """
    Finds the closest edge to a point and projects the point onto that edge's linestring.

    Args:
        graph (nx.MultiDiGraph): The graph with edges containing geometry data.
        point (Point): The point to project.

    Returns:
        Tuple[Point, dict]: The projected point and edge data {'u': u, 'v': v, 'edge': edge},
        or (None, None) if no valid edges are found.
    """
    if not graph.edges():
        return None, None

    closest_edge = None
    min_distance = float("inf")

    for u, v, edge in graph.edges(data=True):
        geometry = edge.get("geometry")
        if isinstance(geometry, LineString):
            distance = point.distance(geometry)
            if distance < min_distance:
                min_distance = distance
                closest_edge = {"u": u, "v": v, "edge": edge}

    if closest_edge is None:
        return None, None

    geometry = closest_edge["edge"]["geometry"]
    projected_point = geometry.interpolate(geometry.project(point))
    return projected_point, closest_edge


def split_edge_at_point(
    graph: nx.MultiDiGraph, projected_point: Point, edge_data: dict
) -> Point:
    """
    Splits an edge at the projected point, replacing it with two new edges in the graph.
    Maintains original edge properties, recalculating length and traversal time.

    Args:
        graph (nx.MultiDiGraph): The graph to modify.
        projected_point (Point): The point where the edge should be split.
        edge_data (dict): Edge data containing 'u', 'v', and 'edge' with properties.

    Returns:
        Point: The projected point (new node) added to the graph.
    """
    p = projected_point
    u, v = edge_data["u"], edge_data["v"]
    edge = edge_data["edge"]
    
    # Check if the edge is already split
    if graph.has_node(p):
        return p

    # Extract edge properties
    edge_id = edge.get("id")
    name = edge.get("name", "")
    highway = edge.get("highway", "unclassified")
    oneway = edge.get("oneway", "no")
    speed = edge.get("speed", 8.33)  # Default 30 km/h
    is_blocked = edge.get("is_blocked", False)

    # Create new segments
    line1 = LineString([u, p])
    line2 = LineString([p, v])
    length1 = u.distance(p)
    length2 = p.distance(v)
    traversal_time1 = length1 / speed if speed > 0 else length1
    traversal_time2 = length2 / speed if speed > 0 else length2

    # Update graph
    graph.add_node(p)
    graph.add_edge(
        u,
        p,
        geometry=line1,
        id=edge_id,
        name=name,
        highway=highway,
        length=length1,
        traversal_time=traversal_time1,
        oneway=oneway,
        is_blocked=is_blocked,
    )
    graph.add_edge(
        p,
        v,
        geometry=line2,
        id=edge_id,
        name=name,
        highway=highway,
        length=length2,
        traversal_time=traversal_time2,
        oneway=oneway,
        is_blocked=is_blocked,
    )

    # Remove original edge
    graph.remove_edge(u, v)

    # Handle bidirectional edges
    if oneway == "no":
        graph.add_edge(
            p,
            u,
            geometry=line1.reverse(),
            id=edge_id,
            name=name,
            highway=highway,
            length=length1,
            traversal_time=traversal_time1,
            oneway=oneway,
            is_blocked=is_blocked,
        )
        graph.add_edge(
            v,
            p,
            geometry=line2.reverse(),
            id=edge_id,
            name=name,
            highway=highway,
            length=length2,
            traversal_time=traversal_time2,
            oneway=oneway,
            is_blocked=is_blocked,
        )
        if graph.has_edge(v, u):
            graph.remove_edge(v, u)

    return p


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
        is_blocked = row.get("is_blocked", False)

        if isinstance(geometry, (LineString, MultiLineString)):
            coords = list(geometry.coords)
            for i in range(len(coords) - 1):
                start_coord = coords[i]
                end_coord = coords[i + 1]

                start_node = Point(start_coord)
                end_node = Point(end_coord)
                
                length = start_node.distance(end_node)
                
                traversal_time = (
                    length / speed if speed > 0 else length / speeds.get("unclassified")
                ) 

                # Add edges to the graph
                graph.add_edge(
                    start_node,
                    end_node,
                    geometry=LineString([start_node, end_node]),
                    id=id,
                    name=name,
                    highway=highway,
                    length=length,
                    traversal_time=traversal_time,
                    oneway=oneway,
                    is_blocked=is_blocked,
                )
                # Add reverse edge if not one-way
                if oneway == "no":
                    graph.add_edge(
                        end_node,
                        start_node,
                        geometry=LineString([end_node, start_node]),
                        id=id,
                        name=name,
                        highway=highway,
                        length=length,
                        traversal_time=traversal_time,
                        oneway=oneway,
                        is_blocked=is_blocked,
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
    # Convert to WGS84
    routing_base = routing_base.to_crs("EPSG:32633")
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
    # routing_base.to_file("data/sample_routing_base.geojson", driver="GeoJSON")

    return routing_base
