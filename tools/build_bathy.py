"""Build data/bathy.json for FAL Wing: seabed depth grid from the EMODnet Bathymetry DTM (WCS).
Data: EMODnet Bathymetry Consortium, used under CC BY 4.0. Run by .github/workflows/bathy.yml."""
import base64, io, json, sys, urllib.request
import numpy as np, tifffile

W, S, E, N = -5.30, 49.93, -4.85, 50.28
RES = 1/960                                  # 1/16 arc-minute, the DTM's native grid
nx, ny = round((E-W)/RES), round((N-S)/RES)
urls = [
  ("https://ows.emodnet-bathymetry.eu/wcs?service=WCS&version=1.0.0&request=GetCoverage"
   f"&coverage=emodnet:mean&crs=EPSG:4326&bbox={W},{S},{E},{N}&width={nx}&height={ny}&format=image/tiff"),
  ("https://ows.emodnet-bathymetry.eu/wcs?service=WCS&version=2.0.1&request=GetCoverage"
   f"&coverageId=emodnet__mean&subset=Lat({S},{N})&subset=Long({W},{E})&format=image/tiff"),
]
arr = None
for u in urls:
    try:
        with urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "FAL Wing (github.com/bobmorek/Falwing)"}), timeout=180) as r:
            raw = r.read()
        arr = tifffile.imread(io.BytesIO(raw)).astype("float32")
        print("ok", u[:90], arr.shape, float(np.nanmin(arr)), float(np.nanmax(arr)), file=sys.stderr)
        break
    except Exception as e:
        print("fail", u[:90], e, file=sys.stderr)
if arr is None: sys.exit("no bathymetry")
if arr.ndim == 3: arr = arr[..., 0] if arr.shape[-1] < 5 else arr[0]
arr[~np.isfinite(arr)] = 0
arr[arr < -1e4] = 0
depth = np.clip(-arr, 0, 127.5)               # metres below mean sea level; land -> 0
q = np.round(depth * 2).astype("uint8")       # 0.5 m steps
h, w = q.shape
out = {"bbox": [W, S, E, N], "nx": int(w), "ny": int(h), "scale": 0.5, "rowsFrom": "north",
       "attribution": "Bathymetry: EMODnet Bathymetry Consortium, CC BY 4.0",
       "depth_b64": base64.b64encode(q.tobytes()).decode()}
json.dump(out, open("data/bathy.json", "w"), separators=(",", ":"))
print("grid", w, h, "max depth", float(depth.max()), file=sys.stderr)
