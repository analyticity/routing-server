# osm.py

import osm2geojson
import requests


def osm_data_for_area(area: str) -> dict:
    """
    Queries OpenStreetMap (OSM) Overpass API for 'highway' data in a given area.

    Args:
        area (str): Name of the area to query, e.g., 'Brno', 'South Moravian Region', 'Czech Republic'.

    Returns:
        dict: Response data from the Overpass API in JSON format.
    """
    url = "https://overpass-api.de/api/interpreter"
    fallback_url = "https://maps.mail.ru/osm/tools/overpass/api/interpreter"

    headers = {"Content-Type": "application/x-www-form-urlencoded"}

    # All roads in the area for car use which are publicly accessible
    query = f"""
    [out:json][timeout:60];
    area[name="{area}"]->.searchArea;
    (
        way["highway"~"motorway|motorway_link|trunk|trunk_link|primary|primary_link|secondary|secondary_link|tertiary|tertiary_link|unclassified|residential|living_street|service"]["access"!~"private|no"]["motor_vehicle"!~"no"]["vehicle"!~"private"](area.searchArea);
    );
    out body;
    >;
    out skel qt;
    """

    # Make the request and ensure successful response
    try:
        response = requests.post(url, data=query, headers=headers)
        response.raise_for_status()
    except requests.exceptions.RequestException:
        # If the first request fails, try the fallback URL
        response = requests.post(fallback_url, data=query, headers=headers)
        response.raise_for_status()
    response_json = response.json()

    # Convert the OSM format to a GeoJSON
    return osm2geojson.json2geojson(response_json)
