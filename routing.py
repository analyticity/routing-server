# routing.py

from datetime import datetime
from heapq import heappop, heappush

import networkx as nx
from geopy.distance import geodesic
from shapely import LineString, Point


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
    meters = geodesic((a.y, a.x), (b.y, b.x)).meters  # geodesic takes (lat, lon)
    return meters / max_speed


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
        LineString: LineString representing the path from source to destination.
    Raises:
        nx.NetworkXNoPath: If no path exists between source and destination.
        nx.NodeNotFound: If source or destination node is not found in the graph.
    """
    weight_property = "traversal_time"  # Graph's edge attribute for weight

    open_set = []  # (f_score, counter, node)
    heappush(open_set, (0, 0, source_node))
    came_from = {}
    g_score = {source_node: 0}
    f_score = {source_node: straightline_heuristic(source_node, destination_node)}
    counter = 0

    while open_set:
        f, _, current = heappop(open_set)
        if current == destination_node:
            path = []
            while current in came_from:
                path.append(current)
                current = came_from[current]
            path.append(source_node)
            path.reverse()
            return LineString(path)

        for neighbor in graph.neighbors(current):
            weight = graph[current][neighbor][0].get(weight_property, 1)
            tentative_g = g_score[current] + weight
            if neighbor not in g_score or tentative_g < g_score[neighbor]:
                came_from[neighbor] = current
                g_score[neighbor] = tentative_g
                f = tentative_g + straightline_heuristic(neighbor, destination_node)
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
    start_date: datetime.date,
    end_date: datetime.date,
) -> LineString:
    """
    Find the shortest path between two coordinates in a graph.

    Args:
        graph (nx.MultiDiGraph): Graph to search in.
        source_coord (Point): Source coordinate.
        destination_coord (Point): Destination coordinate.
        start_date (datetime.date): Start date for traffic data.
        end_date (datetime.date): End date for traffic data.

    Returns:
        list: List of nodes in the path from source to destination.
    """
    source_node = find_nearest_node(graph, source_coord)
    destination_node = find_nearest_node(graph, destination_coord)

    try:
        path = astar_route(graph, source_node, destination_node)
    except nx.NetworkXNoPath:
        return []
    except nx.NodeNotFound:
        return []

    # Add source and destination coordinates to the path to complete
    path = LineString([source_coord] + list(path.coords) + [destination_coord])

    return path
