from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from shapely import Point

from graph import create_graph_from_base, get_routing_base
from models import RoutingCoordRequestBody
from routing import find_route

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


base = get_routing_base("Brno")
graph = create_graph_from_base(base)


@app.post("/find_route_by_coord")
async def find_route_by_coord(body: RoutingCoordRequestBody):
    source_coord = Point((body.src_coord))
    destination_coord = Point((body.dst_coord))
    # start_date = body.from_time
    # end_date = body.to_time

    print(f"Source: {source_coord}")
    print(f"Destination: {destination_coord}")
    # print(f"From: {start_date}")
    # print(f"To: {end_date}")

    route = find_route(
        graph,
        source_coord,
        destination_coord,
        # start_date,
        # end_date
    )

    return {
        "streets_coord": [],
        "route": list(route.coords),
        "src_street": "",
        "dst_street": "",
    }
