# routing.py

from datetime import datetime
from heapq import heappop, heappush
from typing import Literal

import networkx as nx
from shapely import LineString, Point

from graph import find_projected_point, split_edge_at_point


def straightline_heuristic(a: Point, b: Point) -> float:
    """
    Calculate the minimum possible traversal time between two points.
    The heuristic must never overestimate the actual traversal time.

    Args:
        a (Point): Start point.
        b (Point): End point.

    Returns:
        float: Estimated traversal time in seconds.
    """
    max_speed = 27.78  # 100 km/h
    meters = a.distance(b)  # Distance in meters
    return meters / max_speed


def alt_heuristic(a: Point, b: Point, graph: nx.MultiDiGraph, landmarks: list) -> float:
    """
    ALT heuristic calculation.

    Args:
        a (Point): Start point.
        b (Point): End point.
        graph: The graph containing precomputed landmark traversal times.
        landmarks: List of landmark nodes used in preprocessing.

    Returns:
        float: Estimated traversal time in seconds.
    """
    # Retrieve precomputed traversal times for a and b
    a_traversal_times = graph.nodes[a].get("landmark_traversal_time", {})
    b_traversal_times = graph.nodes[b].get("landmark_traversal_time", {})

    max_heuristic_value = 0.0

    for landmark in landmarks:
        # Get traversal times, defaulting to infinity if landmark or direction is missing
        a_traversal_time_data = a_traversal_times.get(landmark, {})
        b_traversal_time_data = b_traversal_times.get(landmark, {})

        traversal_time_L_b = b_traversal_time_data.get("to", float("inf"))
        traversal_time_L_a = a_traversal_time_data.get("to", float("inf"))
        traversal_time_a_L = a_traversal_time_data.get("from", float("inf"))
        traversal_time_b_L = b_traversal_time_data.get("from", float("inf"))

        # Calculate potential heuristic values using triangle inequality
        if traversal_time_L_b != float("inf") and traversal_time_L_a != float("inf"):
            h1 = traversal_time_L_b - traversal_time_L_a
            max_heuristic_value = max(max_heuristic_value, h1)

        if traversal_time_a_L != float("inf") and traversal_time_b_L != float("inf"):
            h2 = traversal_time_a_L - traversal_time_b_L
            max_heuristic_value = max(max_heuristic_value, h2)

    # Ensure heuristic is non-negative
    return max(0.0, max_heuristic_value)


def astar_route(
    graph: nx.MultiDiGraph, source_node: Point, destination_node: Point
) -> LineString:
    """
    Find the shortest path between two nodes in a graph using A* algorithm.

    Args:
        graph (nx.MultiDiGraph): Graph to search in.
        source_node (Point): Source node.
        destination_node (Point): Destination node.

    Returns:
        LineString: LineString representing the path from source to destination,
                    or empty LineString if no path is found.
    """
    weight_property = "traversal_time"  # Graph's edge attribute for weight

    open_set = []  # (f_score, counter, node)
    came_from = {}
    g_score = {source_node: 0}
    f_score = {node: float("inf") for node in graph.nodes()}
    f_score[source_node] = straightline_heuristic(source_node, destination_node)
    counter = 0
    heappush(open_set, (f_score[source_node], counter, source_node))

    while open_set:
        f, _, current = heappop(open_set)

        if f > f_score[current]:
            continue

        if current == destination_node:
            path = []
            temp = current
            while temp in came_from:
                path.append(temp)
                temp = came_from[temp]
            path.append(source_node)
            path.reverse()
            return LineString(path)

        for neighbor in graph.neighbors(current):
            weight = graph[current][neighbor][0].get(weight_property, 1)
            tentative_g = g_score.get(current, float("inf")) + weight

            if tentative_g < g_score.get(neighbor, float("inf")):
                came_from[neighbor] = current
                g_score[neighbor] = tentative_g
                f = tentative_g + straightline_heuristic(neighbor, destination_node)
                f_score[neighbor] = f
                counter += 1
                heappush(open_set, (f, counter, neighbor))

    return LineString([])


def alt_route(
    graph: nx.MultiDiGraph, landmarks: list, source_node: Point, destination_node: Point
) -> LineString:
    """
    Find the shortest path between two nodes using A* with ALT heuristic.

    Args:
        graph: Graph to search in with traffic and landmark heuristic data.
        landmarks: List of landmark nodes used for the heuristic.
        source_node: Source node.
        destination_node: Destination node.

    Returns:
        LineString: LineString representing the path from source to destination,
                    or empty LineString if no path is found.
    """
    weight_property = "traversal_time"  # Graph's edge attribute for weight

    open_set = []  # (f_score, counter, node)
    came_from = {}
    g_score = {source_node: 0}
    f_score = {node: float("inf") for node in graph.nodes()}
    f_score[source_node] = alt_heuristic(
        source_node, destination_node, graph, landmarks
    )
    counter = 0
    heappush(open_set, (f_score[source_node], counter, source_node))

    while open_set:
        f, _, current = heappop(open_set)
        if f > f_score[current]:
            continue

        if current == destination_node:
            path = []
            temp = current
            while temp in came_from:
                path.append(temp)
                temp = came_from[temp]
            path.append(source_node)
            path.reverse()
            return LineString(path)

        for neighbor in graph.neighbors(current):
            weight = graph[current][neighbor][0].get(weight_property, 1)
            tentative_g = g_score.get(current, float("inf")) + weight

            if tentative_g < g_score.get(neighbor, float("inf")):
                came_from[neighbor] = current
                g_score[neighbor] = tentative_g
                f = tentative_g + alt_heuristic(
                    neighbor, destination_node, graph, landmarks
                )
                f_score[neighbor] = f
                counter += 1
                heappush(open_set, (f, counter, neighbor))

    return LineString([])


def find_nearest_node(graph: nx.MultiDiGraph, point: Point) -> Point:
    """
    Find the nearest node in a graph to a given point.

    Args:
        graph (nx.MultiDiGraph): Graph to search in.
        point (Point): Point to find the nearest node to.

    Raises:
        ValueError: If the graph has no nodes.

    Returns:
        Point: The nearest node in the graph to the given point.
    """
    return (
        min(graph.nodes(), key=lambda n: point.distance(n)) if graph.nodes() else None
    )


def find_route(
    graph: nx.MultiDiGraph,
    source_coord: Point,
    destination_coord: Point,
    algorithm: Literal["astar", "alt"] = "astar",
    landmarks: list = None,
) -> tuple[LineString, int, int, list]:
    """
    Find the shortest path between two coordinates in a graph.

    Args:
        graph (nx.MultiDiGraph): Graph to search in.
        source_coord (Point): Source coordinate.
        destination_coord (Point): Destination coordinate.
        algorithm (Literal["astar", "alt"], optional): Algorithm to use for routing, A* or ALT.
        landmarks (list, optional): List of landmark nodes for ALT heuristic.

    Returns:
        tuple: (LineString, total_length, total_traversal_time, street_segments)
            - LineString: LineString representing the path from source to destination
            - total_length: Total length of the path in meters
            - total_traversal_time: Total traversal time in seconds
            - street_segments: List containing street segments with coordinates and names and traffic status
    """
    if algorithm == "astar":
        # Use projected point for source and destination, splitting the edges
        source_node, source_edge = find_projected_point(graph, source_coord)
        destination_node, destination_edge = find_projected_point(
            graph, destination_coord
        )

        if source_node and destination_node:
            source_node = split_edge_at_point(graph, source_node, source_edge)
            destination_node = split_edge_at_point(
                graph, destination_node, destination_edge
            )

            try:
                path = astar_route(graph, source_node, destination_node)
            except nx.NetworkXNoPath or nx.NodeNotFound:
                return []
        else:
            return []

    elif algorithm == "alt":
        # If landmarks are not provided, raise an error
        if landmarks is None:
            raise ValueError("Landmarks must be provided for ALT algorithm.")

        # Use nearest node for source and destination without splitting the edges
        source_node = find_nearest_node(graph, source_coord)
        destination_node = find_nearest_node(graph, destination_coord)

        if source_node and destination_node:
            try:
                path = alt_route(graph, landmarks, source_node, destination_node)
            except nx.NetworkXNoPath or nx.NodeNotFound:
                return []
        else:
            return []
    else:
        raise ValueError("Invalid algorithm specified. Use 'astar' or 'alt'.")

    path = LineString(list(path.coords))
    street_segments = []

    # Calculate traversal time of the path
    total_length = 0.0
    total_traversal_time = 0.0
    
    for i in range(len(path.coords) - 1):
        start = Point(path.coords[i])
        end = Point(path.coords[i + 1])
        edge_data = graph.get_edge_data(start, end, key=0)
        if edge_data:
            segment_length = edge_data.get("length", 0.0)
            segment_time = edge_data.get("traversal_time", 0.0)
            street_name = edge_data.get("name", "")
            
            # Add segment to dictionary
            street_segments.append({
                "street_name": street_name,
                "path": [path.coords[i], path.coords[i + 1]],
                "severity": "2" if edge_data.get("is_penalized_by_traffic", False) else "0",
            })
            
            total_length += segment_length
            total_traversal_time += segment_time

    return path, total_length, total_traversal_time, street_segments
