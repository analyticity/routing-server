# graph.py

import math
import random
from copy import deepcopy
from datetime import date, datetime
from typing import Tuple

import geopandas as gpd
import networkx as nx
from shapely.geometry import LineString, MultiLineString, Point

from osm import osm_data_for_area
from traffic import update_graph_with_traffic


def create_graph_cache():
    """
    Creates a cache for traffic graphs to avoid recomputing them multiple times.

    Returns:
        dict: A dictionary to cache traffic graphs, where keys are tuples of start and end dates.
    """
    traffic_graph_cache: dict[tuple[datetime, datetime], nx.MultiDiGraph] = {}
    return traffic_graph_cache

def get_graph_with_traffic_cached(
    base_graph: nx.MultiDiGraph,
    cache: dict[tuple[date, date], nx.MultiDiGraph],
    edge_jam_overlaps: dict,
    start_date: date,
    end_date: date,
) -> nx.MultiDiGraph:
    """
    Retrieves a cached graph with traffic data for the specified date range.

    Args:
        base_graph (nx.MultiDiGraph): The base graph to be modified.
        cache (dict[tuple[date, date], nx.MultiDiGraph]): Cache for storing traffic graphs.
        edge_jam_overlaps (dict): Dictionary of edge-jam overlaps.
        start_date (date): Start date of the traffic data.
        end_date (date): End date of the traffic data.

    Returns:
        nx.MultiDiGraph: The graph with traffic data applied for the specified date range.
    """
    key = (start_date, end_date)
    if key in cache:
        return cache[key]

    # Get edge-jam overlaps for the specified date range
    available_dates = sorted(edge_jam_overlaps.keys())
    date_range_dates = [d for d in available_dates if start_date <= d <= end_date]
    date_range_overlaps = []
    if date_range_dates:
        for date in date_range_dates:
            date_range_overlaps.extend(edge_jam_overlaps[date])
    num_days = len(date_range_dates)
    print(f"Number of days with traffic in range: {num_days}")
    if num_days > 0:
        graph_copy = deepcopy(base_graph)
        graph_copy = update_graph_with_traffic(graph_copy, date_range_overlaps, num_days)
        cache[key] = graph_copy
        return graph_copy
    else:
        cache[key] = base_graph
        return base_graph


def preprocess_alt(
    graph: nx.MultiDiGraph, num_landmarks: int = 8, weight: str = "traversal_time"
) -> list:
    """
    Performs ALT preprocessing by selecting landmarks and calculating traversal times, saving them in the graph.

    Args:
        graph: The nx.MultiDiGraph to preprocess.
        num_landmarks: How many landmarks to select. Defaults to 8, an optimal number for a medium-sized city.
        weight: The edge attribute representing cost. Defaults to 'traversal_time'.

    Returns:
        A list of the selected landmark nodes.
    """
    all_nodes = list(graph.nodes())
    if not all_nodes:
        print("Graph has no nodes, skipping ALT preprocessing")
        return []
    if len(all_nodes) <= num_landmarks:
        print(
            f"Graph has fewer nodes ({len(all_nodes)}) than requested landmarks ({num_landmarks}). Using all nodes as landmarks"
        )
        landmark_nodes = list(all_nodes)
    else:
        # Determine bounding box of the graph
        min_x = min(n.x for n in all_nodes)
        max_x = max(n.x for n in all_nodes)
        min_y = min(n.y for n in all_nodes)
        max_y = max(n.y for n in all_nodes)

        width = max_x - min_x
        height = max_y - min_y

        # Determine grid dimensions aiming for roughly square cells
        num_cols = max(
            1,
            (
                math.ceil(math.sqrt(num_landmarks * width / height))
                if height > 0
                else num_landmarks
            ),
        )
        num_rows = max(1, math.ceil(num_landmarks / num_cols))
        # Ensure there is enough cells, adjust if calculation is too low due to aspect ratio
        while num_rows * num_cols < num_landmarks:
            num_rows += 1

        cell_width = width / num_cols
        cell_height = height / num_rows

        # Assign nodes to regions
        nodes_in_region = {}  # (row, col): list of nodes
        for node in all_nodes:
            # Calculate column index, clamping to [0, num_cols-1]
            col = (
                min(num_cols - 1, math.floor((node.x - min_x) / cell_width))
                if cell_width > 0
                else 0
            )
            # Calculate row index, clamping to [0, num_rows-1]
            row = (
                min(num_rows - 1, math.floor((node.y - min_y) / cell_height))
                if cell_height > 0
                else 0
            )
            # Ensure non-negative indices
            col = max(0, col)
            row = max(0, row)

            region_key = (row, col)
            if region_key not in nodes_in_region:
                nodes_in_region[region_key] = []
            nodes_in_region[region_key].append(node)

        # Select landmarks
        selected_landmarks_set = set()
        # Get list of regions that contain nodes
        available_regions = [
            (r, c) for (r, c), nodes in nodes_in_region.items() if nodes
        ]
        random.shuffle(available_regions)  # Shuffle order for region processing

        # Pick one landmark from distinct available regions (up to num_landmarks)
        for r, c in available_regions:
            if len(selected_landmarks_set) >= num_landmarks:
                break
            # Select a random node from this region's list
            candidates = nodes_in_region[(r, c)]
            chosen_node = random.choice(candidates)
            selected_landmarks_set.add(chosen_node)

        print(f"Selected {len(selected_landmarks_set)} landmarks")

        # If more landmarks needed, pick random nodes from the remaining pool
        if len(selected_landmarks_set) < num_landmarks:
            # Create a pool of all nodes not yet selected
            remaining_nodes_pool = [
                node for node in all_nodes if node not in selected_landmarks_set
            ]
            random.shuffle(remaining_nodes_pool)

            needed = num_landmarks - len(selected_landmarks_set)
            # Add needed number of landmarks from the shuffled remaining pool
            selected_landmarks_set.update(remaining_nodes_pool[:needed])

        landmark_nodes = list(selected_landmarks_set)

    # Initialize storage on nodes
    for node in all_nodes:
        graph.nodes[node]["landmark_traversal_time"] = {}

    # Calculate traversal times for each landmark
    # Pre-create a reversed graph for efficiency
    reversed_graph = graph.reverse(copy=False)

    for i, landmark in enumerate(landmark_nodes):
        print(
            f"\tCalculating traversal times for landmark {i + 1}/{len(landmark_nodes)} = {100 * (i + 1) / len(landmark_nodes):.2f}%"
        )

        # Traversal times from landmark
        traversal_time_from_L = nx.single_source_dijkstra_path_length(
            graph, landmark, weight=weight
        )
        for node, traversal_time in traversal_time_from_L.items():
            if landmark not in graph.nodes[node]["landmark_traversal_time"]:
                graph.nodes[node]["landmark_traversal_time"][
                    landmark
                ] = {}  # Ensure dict exists
            graph.nodes[node]["landmark_traversal_time"][landmark][
                "to"
            ] = traversal_time

        # Traversal times to landmark using reversed graph
        traversal_time_to_L = nx.single_source_dijkstra_path_length(
            reversed_graph, landmark, weight=weight
        )
        for node, traversal_time in traversal_time_to_L.items():
            if landmark not in graph.nodes[node]["landmark_traversal_time"]:
                graph.nodes[node]["landmark_traversal_time"][
                    landmark
                ] = {}  # Ensure dict exists
            graph.nodes[node]["landmark_traversal_time"][landmark][
                "from"
            ] = traversal_time

    return landmark_nodes


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
                    traversal_time_with_traffic=traversal_time,
                    oneway=oneway,
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
                        traversal_time_with_traffic=traversal_time,
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
    if not osm_data:
        raise ValueError(f"No OSM data found for area: {area}")
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
    routing_base["stroke-width"] = 5

    routing_base = routing_base.drop(columns=["tags"])

    # Save the routing base to a GeoJSON file
    # routing_base.to_file("data/sample_routing_base.geojson", driver="GeoJSON")

    return routing_base
