"""Build data/coast.json (land polygons + coastline) for the FAL Wing map from OpenStreetMap.
OSM data (c) OpenStreetMap contributors, ODbL. Run by .github/workflows/coast.yml."""
import json, sys, time, urllib.parse, urllib.request
from shapely.geometry import LineString, box, Point, mapping
from shapely.ops import unary_union, polygonize, linemerge

S, W, N, E = 50.00, -5.25, 50.28, -4.85          # data box (lat/lon)
SEA_SEEDS = [(-4.95, 50.10), (-5.00, 50.03)]     # (lon, lat) points known to be open sea
TOL = 0.00006                                     # simplify tolerance, degrees (~5 m)

q = f'[out:json][timeout:120];way["natural"="coastline"]({S},{W},{N},{E});out geom;'
eps = ["https://overpass-api.de/api/interpreter",
       "https://overpass.kumi.systems/api/interpreter",
       "https://maps.mail.ru/osm/tools/overpass/api/interpreter"]
data = None
for attempt in range(3):
    for ep in eps:
        try:
            req = urllib.request.Request(ep, data=urllib.parse.urlencode({"data": q}).encode(),
                                         headers={"User-Agent": "FAL Wing coastline builder (github.com/bobmorek/Falwing)"})
            with urllib.request.urlopen(req, timeout=180) as r:
                data = json.load(r)
            if data.get("elements"): break
        except Exception as e:
            print("fail", ep, e, file=sys.stderr)
    if data and data.get("elements"): break
    time.sleep(20)
if not data or not data.get("elements"):
    sys.exit("no coastline data")

frame = box(W, S, E, N)
lines = [LineString([(g["lon"], g["lat"]) for g in el["geometry"]])
         for el in data["elements"] if el.get("geometry") and len(el["geometry"]) > 1]
coast = linemerge(unary_union(lines)).intersection(frame)
noded = unary_union([coast, frame.exterior])
polys = list(polygonize(noded))
seeds = [Point(p) for p in SEA_SEEDS]
land = [p for p in polys if not any(p.contains(s) for s in seeds)]
print(f"ways={len(lines)} polygons={len(polys)} land={len(land)}", file=sys.stderr)

def rnd(o):
    if isinstance(o, (list, tuple)):
        return [rnd(x) for x in o] if not (len(o) == 2 and all(isinstance(v, float) for v in o)) else [round(o[0], 5), round(o[1], 5)]
    return o

land_g = unary_union(land).simplify(TOL, preserve_topology=True)
coast_g = coast.simplify(TOL, preserve_topology=True)
out = {"bbox": [W, S, E, N],
       "attribution": "Coastline (c) OpenStreetMap contributors, ODbL",
       "land": rnd(mapping(land_g)),
       "coast": rnd(mapping(coast_g))}
json.dump(out, open("data/coast.json", "w"), separators=(",", ":"))
print("bytes", len(json.dumps(out, separators=(",", ":"))), file=sys.stderr)
