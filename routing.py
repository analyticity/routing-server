import networkx as nx
from shapely import Point


def find_nearest_node(graph: nx.MultiDiGraph, point: Point) -> str:
    """
    Find the nearest node in a graph to a given point.

    Args:
        graph (nx.MultiDiGraph): Graph to search in.
        point (Point): Point to find the nearest node to.

    Returns:
        str: Coordinate of the nearest node.
    """
    node = None
    min_distance = float("inf")
    for n in graph.nodes():
        distance = point.distance(Point(n))
        if distance < min_distance:
            min_distance = distance
            node = n
    return node


def find_route(
    graph: nx.MultiDiGraph,
    source_coord: Point,
    destination_coord: Point,
    # start_date: datetime,
    # end_date: datetime,
) -> list:
    """
    Find the shortest path between two coordinates in a graph.

    Args:
        graph (nx.MultiDiGraph): Graph to search in.
        source_coord (Point): Source coordinate.
        destination_coord (Point): Destination coordinate.

    Returns:
        list: List of nodes in the path from source to destination.
    """

    source_node = find_nearest_node(graph, source_coord)
    destination_node = find_nearest_node(graph, destination_coord)

    try:
        path = nx.shortest_path(graph, source_node, destination_node, weight="weight")
    except nx.NetworkXNoPath:
        return []
    except nx.NodeNotFound:
        return []

    path_coords = []

    path_coords.append(source_node)
    path_coords.extend(path)
    path_coords.append(destination_node)

    return path_coords
