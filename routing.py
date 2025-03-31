import networkx as nx
from shapely import LineString, Point


def find_nearest_node(graph: nx.MultiDiGraph, point: Point) -> str:
    """
    Find the nearest node in a graph to a given point.

    Args:
        graph (nx.MultiDiGraph): Graph to search in.
        point (Point): Point to find the nearest node to.

    Returns:
        str: Coordinate of the nearest node.
    """
    # Todo: Create a new node on the graph which will be the closest to the point going perpendicular to the street
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
) -> LineString:
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
        path = [Point(coord) for coord in path]
    except nx.NetworkXNoPath:
        return []
    except nx.NodeNotFound:
        return []

    path.insert(0, source_coord)
    path.append(destination_coord)

    return LineString(path)
