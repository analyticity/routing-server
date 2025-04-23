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
landmarks = preprocess_alt(graph=unmodified_graph)
traffic = load_traffic_data("data/processed_jams.geojson")

days = sorted(traffic["date"].unique())[-3:]
traffic_subset = traffic[(traffic["date"] >= days[0]) & (traffic["date"] <= days[-1])]

graph = deepcopy(unmodified_graph)
edge_jam_overlaps = get_edge_jam_overlaps(
    graph=graph,
    traffic=traffic_subset,
)
new_graph = update_graph_with_traffic(
    graph=graph,
    edge_jam_overlaps=edge_jam_overlaps,
    date_range=len(traffic_subset["date"].unique()),
)


@app.post("/find_route_by_coord")
async def find_route_by_coord(body: RoutingCoordRequestBody):
    """Finds the route between two coordinates in the given timeframe."""
    source_coord = Point((body.src_coord))
    destination_coord = Point((body.dst_coord))
    start_date = datetime.strptime(body.from_time, "%Y-%m-%d").date()
    end_date = datetime.strptime(body.to_time, "%Y-%m-%d").date()

    print(f"Source: {source_coord}")
    print(f"Destination: {destination_coord}")
    print(f"From: {start_date}")
    print(f"To: {end_date}")

    # Convert coordinates from CRS EPSG:4326 to CRS EPSG:32633
    source = gpd.GeoSeries([source_coord], crs="EPSG:4326").to_crs("EPSG:32633").iloc[0]
    destination = (
        gpd.GeoSeries([destination_coord], crs="EPSG:4326").to_crs("EPSG:32633").iloc[0]
    )

    route = find_route(graph=graph, source_coord=source, destination_coord=destination, algorithm="alt", landmarks=landmarks)

    if route:
        # Convert route back to WGS84 for output
        route = gpd.GeoSeries([route], crs="EPSG:32633").to_crs("EPSG:4326")[0]
        return {
            "streets_coord": [],
            "route": list(route.coords),
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
