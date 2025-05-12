# main.py
from copy import deepcopy
from datetime import datetime

import geopandas as gpd
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from shapely import Point

from graph import create_graph_from_base, get_routing_base, preprocess_alt
from models import RoutingCoordRequestBody
from routing import find_route
from traffic import get_edge_jam_overlaps, load_traffic_data, update_graph_with_traffic

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

AREA = "Brno"

base = get_routing_base(AREA)
unmodified_graph = create_graph_from_base(base)

# Create a template graph with traffic penalties
template_graph = deepcopy(unmodified_graph)
# Store original traversal times for quick restoration
original_data = {}
for u, v, key, data in template_graph.edges(data=True, keys=True):
    original_data[(u, v, key)] = data.copy()

landmarks = preprocess_alt(graph=unmodified_graph)
traffic = load_traffic_data("data/partial_processed_jams.geojson")

# Precalculate all edge-jam overlaps for the entire dataset
all_edge_jam_overlaps = get_edge_jam_overlaps(
    graph=unmodified_graph,
    traffic=traffic,
)

# Store overlaps by date for quick access
edge_jam_overlaps_by_date = {}
for (u, v, key), jam_rows in all_edge_jam_overlaps:
    for date in jam_rows["date"].unique():
        if date not in edge_jam_overlaps_by_date:
            edge_jam_overlaps_by_date[date] = []
        # Filter jam rows for this date
        date_jam_rows = jam_rows[jam_rows["date"] == date]
        if not date_jam_rows.empty:
            edge_jam_overlaps_by_date[date].append(((u, v, key), date_jam_rows))


@app.post("/find_route_by_coord")
async def find_route_by_coord(body: RoutingCoordRequestBody):
    """Finds the route between two coordinates in the given timeframe."""
    source_coord = Point((body.src_coord))
    destination_coord = Point((body.dst_coord))
    start_date = datetime.strptime(body.from_time, "%Y-%m-%d").date()
    end_date = datetime.strptime(body.to_time, "%Y-%m-%d").date()
    use_traffic = bool(body.use_traffic)

    print(f"Source: {source_coord}")
    print(f"Destination: {destination_coord}")
    print(f"From: {start_date}")
    print(f"To: {end_date}")
    print(f"Use traffic: {use_traffic}")

    # Convert coordinates from CRS EPSG:4326 to CRS EPSG:32633
    source = gpd.GeoSeries([source_coord], crs="EPSG:4326").to_crs("EPSG:32633").iloc[0]
    destination = (
        gpd.GeoSeries([destination_coord], crs="EPSG:4326").to_crs("EPSG:32633").iloc[0]
    )

    if use_traffic:
        # Get all dates in range that have traffic data
        available_dates = sorted(edge_jam_overlaps_by_date.keys())
        date_range_dates = [d for d in available_dates if start_date <= d <= end_date]

        if date_range_dates:
            # Get all overlaps for the date range
            date_range_overlaps = []
            for date in date_range_dates:
                date_range_overlaps.extend(edge_jam_overlaps_by_date[date])

            # Restore original data
            for (u, v, key), data in original_data.items():
                for attr, value in data.items():
                    template_graph.edges[u, v, key][attr] = value

            # Update graph with collected overlaps
            date_range = len(date_range_dates)
            graph = update_graph_with_traffic(
                graph=template_graph,
                edge_jam_overlaps=date_range_overlaps,
                date_range=date_range,
            )
        else:
            graph = template_graph
    else:
        graph = unmodified_graph

    route, length, time, streets = find_route(
        graph=graph,
        source_coord=source,
        destination_coord=destination,
        algorithm="alt",
        landmarks=landmarks,
    )

    if route:
        # Convert route back to WGS84 for output
        route = gpd.GeoSeries([route], crs="EPSG:32633").to_crs("EPSG:4326")[0]
        path = [[coord[1], coord[0]] for coord in route.coords]
        streets_coord = {"street_name": "", "path": path, "color": "blue"}
        return {
            "streets_coord": [streets_coord],
            "route": [],
            "src_street": "",
            "dst_street": "",
        }
    else:
        # No route found
        return {
            "streets_coord": [],
            "route": [],
            "src_street": "",
            "dst_street": "",
        }
