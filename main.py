"""
Filename: main.py
Author: Matyáš Strelec
Email: xstrel03@vut.cz
Date: 2025-05-20
Project: Routing Server
"""

from copy import deepcopy
from datetime import datetime, timedelta
import asyncio
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
from contextlib import asynccontextmanager

# Globals
AREA = "Brno"

unmodified_graph = None
graph = None
landmarks = None
traffic_graph_cache = None
edge_jam_overlaps = None
traffic = None


@asynccontextmanager
async def lifespan(app: FastAPI):
    global \
        unmodified_graph, \
        graph, \
        traffic_graph_cache, \
        edge_jam_overlaps, \
        date_edge_jam_overlaps, \
        landmarks

    load_dotenv()

    base = get_routing_base(AREA, local=True)
    unmodified_graph = create_graph_from_base(base)
    graph = deepcopy(unmodified_graph)
    landmarks = preprocess_alt(graph=unmodified_graph, num_landmarks=1)
    jams = load_jam_data_from_db()
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

    asyncio.create_task(schedule_daily_traffic_refresh())
    asyncio.create_task(cache_date_pairs())
    yield


app = FastAPI(lifespan=lifespan)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


async def cache_single_pair(start_date, end_date):
    """
    Caches a single pair of dates for traffic data.
    """
    get_graph_with_traffic_cached(
        base_graph=unmodified_graph,
        cache=traffic_graph_cache,
        edge_jam_overlaps=edge_jam_overlaps,
        start_date=start_date,
        end_date=end_date,
    )


async def cache_date_pairs():
    """
    Caches date pairs for the last 7 days.
    """
    today = datetime.now().date()
    tasks = []

    # Last 7 days
    for offset_end in range(7):
        end_date = today - timedelta(days=offset_end)
        for offset_start in range(offset_end + 1, 7):
            start_date = today - timedelta(days=offset_start)
            tasks.append(cache_single_pair(start_date, end_date))

    # Run up to N concurrently
    semaphore = asyncio.Semaphore(100)

    async def sem_task(task):
        async with semaphore:
            await task

    await asyncio.gather(*(sem_task(t) for t in tasks))
    print("Cached date pairs for traffic data.")


async def schedule_daily_traffic_refresh():
    """
    Schedules a daily refresh of traffic data at 1 AM.
    """
    while True:
        now = datetime.now()
        # 1 AM
        next_run = now.replace(hour=1, minute=0, second=0, microsecond=0)
        if next_run <= now:
            next_run += timedelta(days=1)
        await asyncio.sleep((next_run - now).total_seconds())
        await refresh_traffic()


async def refresh_traffic():
    """
    Refreshes the traffic data and updates the graph cache.
    """
    global edge_jam_overlaps
    jams = load_jam_data_from_db()
    traffic = preprocess_jams(jams)
    edge_jam_overlaps = get_edge_jam_overlaps(
        graph=unmodified_graph,
        traffic=traffic,
    )
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
