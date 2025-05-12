# main.py
from copy import deepcopy
from datetime import datetime

import geopandas as gpd
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from shapely import LineString, Point

from graph import create_graph_from_base, get_routing_base, preprocess_alt, create_graph_cache, get_graph_with_traffic_cached 
from models import RoutingCoordRequestBody
from routing import find_route
from traffic import get_edge_jam_overlaps, load_traffic_data

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
graph = deepcopy(unmodified_graph)
landmarks = preprocess_alt(graph=unmodified_graph)
traffic = load_traffic_data("data/partial_processed_jams.geojson")
edge_jam_overlaps = get_edge_jam_overlaps(
    graph=unmodified_graph,
    traffic=traffic,
)
traffic_graph_cache = create_graph_cache()


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

    graph = get_graph_with_traffic_cached(
        base_graph=unmodified_graph,
        cache=traffic_graph_cache,
        edge_jam_overlaps=edge_jam_overlaps,
        start_date=start_date,
        end_date=end_date,
    )     

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

        # Convert each street segment to WGS84
        for segment in streets:
            linestring = LineString(segment["path"])
            linestring = gpd.GeoSeries([linestring], crs="EPSG:32633").to_crs("EPSG:4326")[0]
            segment["path"] = [[coord[1], coord[0]] for coord in linestring.coords]
            if segment["severity"] == "0":
                segment["color"] = "green"
            elif segment["severity"] == "1":
                segment["color"] = "orange"
            elif segment["severity"] == "2":
                segment["color"] = "red"
            else:
                segment["color"] = "green"
        
        # Find the first street and last street names
        src_street, dst_street = "", ""
        for segment in streets:
            if segment["street_name"] != "":
                src_street = segment["street_name"]
                break
        for segment in reversed(streets):
            if segment["street_name"] != "":
                dst_street = segment["street_name"]
                break
        
        return {
            "streets_coord": streets,
            "route": path,
            "src_street": src_street,
            "dst_street": dst_street,
            "length": length,
            "time": time,
        }
    else:
        # No route found
        return {
            "streets_coord": [],
            "route": [],
            "src_street": "",
            "dst_street": "",
        }
