"""Ground roughness z0 from ESA WorldCover 2021 (10 m) for an OpenFOAM ABL case.

Writes, inside the case:
  constant/boundaryData/<patch>/points   sample points on a regular grid (z = 0)
  constant/boundaryData/<patch>/0/z0     z0 per point          -> per-face list (tools/z0_faces.py)
  constant/boundaryData/<patch>/0/albedo, 0/sensibleFraction   -> ground materials for solar_illumination.py
  geo/z0_map.png, geo/z0_map_zoom.png    land cover + z0 + building roofs (alignment check)
  geo/inlet_z0_sectors.csv               log-mean z0 of the upstream fetch per inflow direction (for z0In)

Mesh coordinates must be a Transverse Mercator centred on (lat0, lon0), x = east, y = north, metres, z up,
ground at z = 0.
See README.md.

Usage:
  python make_z0.py --case <case> --lat0 <lat> --lon0 <lon> --radius <m> [options]
"""
import argparse
import math
import os
import urllib.request
import numpy as np
import rasterio
from pyproj import Transformer
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, LogNorm
from matplotlib.collections import PolyCollection

# WorldCover v200 class -> (name, z0 [m], colour); z0 after Wieringa (1992) / Blocken (2015) Table 1
CLASSES = {
    10: ("tree cover", 1.0, "#1c6b2c"),
    20: ("shrubland", 0.25, "#9a8c3a"),
    30: ("grassland", 0.03, "#b9d77a"),
    40: ("cropland", 0.10, "#e8d36a"),
    50: ("built-up", 1.5, "#c0392b"),
    60: ("bare", 0.005, "#bdb6a8"),
    70: ("snow/ice", 0.0002, "#ffffff"),
    80: ("water", 0.0002, "#2a78d6"),
    90: ("herb. wetland", 0.03, "#5fb3a6"),
    95: ("mangroves", 0.5, "#0f8f6f"),
    100: ("moss/lichen", 0.03, "#d9d0a8"),
}
NODATA_Z0 = CLASSES[80][1]          # WorldCover has no data over open sea -> water

# Surface for the Solar tool: class -> (albedo, sensible fraction of absorbed shortwave)
# Typical values (Oke et al. 2017, Urban Climates); written next to z0 as albedo / sensibleFraction
SURFACE = {
    10: (0.18, 0.4),   # trees: shaded, transpiring
    20: (0.20, 0.6),
    30: (0.23, 0.5),   # grass
    40: (0.20, 0.5),
    50: (0.15, 0.9),   # built-up seen from above: roofs + roads + some green
    60: (0.20, 0.7),   # bare soil
    70: (0.60, 0.3),
    80: (0.07, 0.15),  # water: absorbs, mostly evaporation / storage
    90: (0.15, 0.3),
    95: (0.15, 0.3),
    100: (0.20, 0.5),
}
# Built-up pixels between explicit buildings keep the built-up surface values. A separate asphalt value
# (0.12, 1.0) here gave a sharp 50 m halo of different albedo around the explicit buildings.
SURFACE_EXPLICIT = SURFACE[50]
URL = "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/{}"

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--case", default=".", help="OpenFOAM case directory")
ap.add_argument("--lat0", type=float, required=True, help="projection origin latitude")
ap.add_argument("--lon0", type=float, required=True, help="projection origin longitude")
ap.add_argument("--stl", default="geo/Buildings.stl", help="explicit buildings, relative to case")
ap.add_argument("--patch", default="Terrain", help="ground patch name")
ap.add_argument("--radius", type=float, required=True,
                help="half-size of the sampled area, m: at least the farthest ground point from the origin")
ap.add_argument("--dx", type=float, default=20.0, help="sample spacing, m")
ap.add_argument("--explicit-dist", type=float, default=50.0,
                help="built-up pixels closer than this to an explicit building get --z0-explicit")
ap.add_argument("--z0-explicit", type=float, default=0.1, help="z0 between explicit buildings, m")
ap.add_argument("--fetch", type=float, default=5000.0, help="upstream fetch outside the domain for inlet z0, m")
ap.add_argument("--sectors", type=int, default=32, help="number of inflow sectors for inlet z0")
ap.add_argument("--tile-dir", default=os.path.expanduser("~/worldcover"), help="WorldCover tile cache")
a = ap.parse_args()
os.chdir(a.case)

to_ll = Transformer.from_crs(f"+proj=tmerc +lat_0={a.lat0} +lon_0={a.lon0} +k=1 +x_0=0 +y_0=0 +datum=WGS84",
                             "EPSG:4326", always_xy=True)


# ---- WorldCover tiles (3 x 3 deg, named by lower-left corner), downloaded once into --tile-dir
def tile_name(lat, lon):
    la, lo = 3 * math.floor(lat / 3), 3 * math.floor(lon / 3)
    return (f"ESA_WorldCover_10m_2021_v200_{'N' if la >= 0 else 'S'}{abs(la):02d}"
            f"{'E' if lo >= 0 else 'W'}{abs(lo):03d}_Map.tif")


_open = {}


def tile(name):
    if name not in _open:
        path = os.path.join(a.tile_dir, name)
        if not os.path.exists(path):
            os.makedirs(a.tile_dir, exist_ok=True)
            print(f"downloading {name} ...", flush=True)
            try:
                urllib.request.urlretrieve(URL.format(name), path + ".part")
            except Exception as e:        # tile over open ocean does not exist
                print(f"  not available ({e}); treated as water")
                _open[name] = None
                return None
            os.replace(path + ".part", path)
        _open[name] = rasterio.open(path)
    return _open[name]


def landcover(x, y):
    """WorldCover class at mesh points (x, y); 0 = no data."""
    lon, lat = to_ll.transform(np.asarray(x, float), np.asarray(y, float))
    out = np.zeros(lon.shape, np.uint8)
    keys = np.char.add((3 * np.floor(lat / 3)).astype(int).astype(str),
                       np.char.add("_", (3 * np.floor(lon / 3)).astype(int).astype(str)))
    for k in np.unique(keys):
        m = keys == k
        src = tile(tile_name(lat[m][0], lon[m][0]))
        if src is None:
            continue
        rows, cols = rasterio.transform.rowcol(src.transform, lon[m], lat[m])
        rows, cols = np.asarray(rows), np.asarray(cols)
        r0, r1, c0, c1 = rows.min(), rows.max() + 1, cols.min(), cols.max() + 1
        block = src.read(1, window=((r0, r1), (c0, c1)))
        out[m] = block[rows - r0, cols - c0]
    return out


def z0_of(cls):
    z = np.full(cls.shape, NODATA_Z0)
    for k, (_, z0, _) in CLASSES.items():
        z[cls == k] = z0
    return z


# ---- explicit buildings (binary or ASCII STL)
def read_stl(path):
    with open(path, "rb") as f:
        head = f.read(84)
        n = int(np.frombuffer(head[80:84], "<u4")[0])
        if os.path.getsize(path) == 84 + 50 * n:
            dt = np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
            return np.frombuffer(f.read(), dt, count=n)["v"].astype(float)
    txt = open(path).read().split()
    idx = [i for i, t in enumerate(txt) if t == "vertex"]
    v = np.array([[float(txt[i + 1]), float(txt[i + 2]), float(txt[i + 3])] for i in idx])
    return v.reshape(-1, 3, 3)


tris = read_stl(a.stl)
bxy = np.unique(np.round(tris[:, :, :2].reshape(-1, 2)), axis=0)
tree = cKDTree(bxy)
nz = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
roofs = tris[np.abs(nz[:, 2]) > 0.5 * np.linalg.norm(nz, axis=1) + 1e-12][:, :, :2]

# ---- sample grid on the ground
g = np.arange(-a.radius, a.radius + a.dx / 2, a.dx)
X, Y = np.meshgrid(g, g)
inside = np.hypot(X, Y) <= a.radius
x, y = X[inside], Y[inside]
cls = landcover(x, y)
z0 = z0_of(cls)
near_bldg = tree.query(np.c_[x, y], distance_upper_bound=a.explicit_dist)[0] < np.inf
z0[(cls == 50) & near_bldg] = a.z0_explicit

# ---- write boundaryData
out = os.path.join("constant", "boundaryData", a.patch)
os.makedirs(os.path.join(out, "0"), exist_ok=True)


def foam_header(cls_name, obj):
    return ("FoamFile\n{\n    version 2.0;\n    format ascii;\n"
            f"    class {cls_name};\n    object {obj};\n}}\n\n")


with open(os.path.join(out, "points"), "w") as f:
    f.write(foam_header("vectorField", "points") + f"{len(x)}\n(\n")
    np.savetxt(f, np.c_[x, y, np.zeros_like(x)], fmt="(%.1f %.1f %.1f)")
    f.write(")\n")
alb = np.full(cls.shape, SURFACE[80][0]); sens = np.full(cls.shape, SURFACE[80][1])
for k, (al, sf) in SURFACE.items():
    alb[cls == k], sens[cls == k] = al, sf
alb[(cls == 50) & near_bldg], sens[(cls == 50) & near_bldg] = SURFACE_EXPLICIT
for name, vals in (("z0", z0), ("albedo", alb), ("sensibleFraction", sens)):
    with open(os.path.join(out, "0", name), "w") as f:
        f.write(foam_header("scalarField", name) + f"{len(vals)}\n(\n")
        np.savetxt(f, vals, fmt="%g")
        f.write(")\n")

# ---- inlet z0 per sector: log-mean over the fetch just outside the domain
rows = []
for s in range(a.sectors):
    d = s * 360 / a.sectors                                   # wind FROM direction, deg clockwise from N
    th = np.radians(d + np.linspace(-0.5, 0.5, 9) * 360 / a.sectors)
    rr = np.arange(a.radius, a.radius + a.fetch, 50.0)
    T, Rr = np.meshgrid(th, rr)
    zz = z0_of(landcover((Rr * np.sin(T)).ravel(), (Rr * np.cos(T)).ravel()))
    rows.append((d, float(np.exp(np.log(zz).mean()))))
os.makedirs("geo", exist_ok=True)
with open(os.path.join("geo", "inlet_z0_sectors.csv"), "w") as f:
    f.write("dir_from_deg,z0_logmean_m\n")
    f.writelines(f"{d:g},{z:.4f}\n" for d, z in rows)

# ---- report
print(f"origin lat/lon {a.lat0:.6f}, {a.lon0:.6f}; {len(tris)} building triangles ({a.stl})")
print(f"{len(x)} sample points, dx = {a.dx:g} m, radius = {a.radius:.0f} m -> {out}")
for k, (name, zv, _) in CLASSES.items():
    n = int((cls == k).sum())
    if n:
        print(f"  {name:14s} z0 = {zv:<7g} {100 * n / len(cls):5.1f} %")
if (cls == 0).any():
    print(f"  no data (water) z0 = {NODATA_Z0:<7g} {100 * (cls == 0).mean():5.1f} %")
print(f"  built-up within {a.explicit_dist:g} m of explicit buildings -> {a.z0_explicit} m: "
      f"{100 * ((cls == 50) & near_bldg).mean():.1f} %")
print("inlet z0 (log-mean of fetch) per sector:", ", ".join(f"{d:g}:{z:.2f}" for d, z in rows))

# ---- figures
fig, axs = plt.subplots(1, 2, figsize=(14, 7), facecolor="#fcfcfb")
keys = sorted(CLASSES)
cmap = ListedColormap([CLASSES[k][2] for k in keys])
norm = BoundaryNorm([k - 0.5 for k in keys] + [keys[-1] + 0.5], cmap.N)
C = np.full(X.shape, np.nan); C[inside] = np.where(cls == 0, 80, cls)
Z = np.full(X.shape, np.nan); Z[inside] = z0
ext = (g[0], g[-1], g[0], g[-1])
axs[0].imshow(C, origin="lower", extent=ext, cmap=cmap, norm=norm, interpolation="nearest")
axs[0].set_title("ESA WorldCover 2021 + explicit buildings")
im = axs[1].imshow(Z, origin="lower", extent=ext, cmap="viridis", norm=LogNorm(2e-4, 2), interpolation="nearest")
axs[1].set_title(f"z0 on patch {a.patch}, m")
fig.colorbar(im, ax=axs[1], shrink=0.8)
for ax in axs:
    ax.add_collection(PolyCollection(roofs, facecolor="#0b0b0b", edgecolor="none", linewidth=0))
    ax.set_xlabel("x (east), m"); ax.set_ylabel("y (north), m"); ax.set_aspect("equal")
handles = [plt.Line2D([], [], ls="", marker="s", color=CLASSES[k][2], label=CLASSES[k][0])
           for k in keys if (cls == k).any()]
axs[0].legend(handles=handles, loc="lower left", fontsize=7, framealpha=0.9)
fig.tight_layout()
fig.savefig(os.path.join("geo", "z0_map.png"), dpi=130)

lo, hi = tris[:, :, :2].reshape(-1, 2).min(axis=0) - 200, tris[:, :, :2].reshape(-1, 2).max(axis=0) + 200
for ax in axs:
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1])
fig.set_size_inches(16, 7.5)
fig.tight_layout()
fig.savefig(os.path.join("geo", "z0_map_zoom.png"), dpi=150)
print("written geo/z0_map.png, geo/z0_map_zoom.png, geo/inlet_z0_sectors.csv")
