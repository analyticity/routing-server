# main.py
import os
from copy import deepcopy
from datetime import datetime, timedelta

import geopandas as gpd
from dotenv import load_dotenv
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from shapely import Point

from graph import (
    create_graph_cache,
    create_graph_from_base,
    get_graph_with_traffic_cached,
    get_routing_base,
    preprocess_alt,
)
from models import RoutingCoordRequestBody
from routing import find_route, prepare_route_response
from traffic import get_edge_jam_overlaps, load_jam_data_from_db, preprocess_jams

load_dotenv()

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

AREA = "Brno"  # Area to be used for routing
db_config = {
    "host": os.getenv("DB_HOST"),
    "port": int(os.getenv("DB_PORT")),
    "user": os.getenv("DB_USER"),
    "password": os.getenv("DB_PASSWORD"),
    "dbname": os.getenv("DB_NAME"),
}

base = get_routing_base(AREA, local=True)
unmodified_graph = create_graph_from_base(base)
graph = deepcopy(unmodified_graph)
landmarks = preprocess_alt(graph=unmodified_graph)
jams = load_jam_data_from_db(db_config)
traffic = preprocess_jams(jams)
edge_jam_overlaps = get_edge_jam_overlaps(
    graph=unmodified_graph,
    traffic=traffic,
)
traffic_graph_cache = create_graph_cache()
# Cache the last week of traffic data
get_graph_with_traffic_cached(
        base_graph=unmodified_graph,
        cache=traffic_graph_cache,
        edge_jam_overlaps=edge_jam_overlaps,
        start_date=datetime.now().date() - timedelta(days=7),
        end_date=datetime.now().date(),
)

@app.post("/find_route_by_coord")
async def find_route_by_coord(body: RoutingCoordRequestBody):
    """Finds the route between two coordinates in the given timeframe."""
    algorithm = "alt"  # Default to ALT
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

    current_graph = get_graph_with_traffic_cached(
        base_graph=unmodified_graph,
        cache=traffic_graph_cache,
        edge_jam_overlaps=edge_jam_overlaps,
        start_date=start_date,
        end_date=end_date,
    )     

    route, length, time_with_traffic, time_without_traffic, streets = find_route(
        graph=current_graph,
        source_coord=source,
        destination_coord=destination,
        algorithm=algorithm,
        landmarks=landmarks,
        use_traffic=use_traffic,
    )

    if route:
        return prepare_route_response(
            route, length, time_with_traffic, time_without_traffic, streets
        )
    else:
        # No route found
        return {
            "streets_coord": [],
            "route": [],
            "src_street": "",
            "dst_street": "",
            "length": 0.0,
            "time_with_traffic": 0.0,
            "time_without_traffic": 0.0,
        }
