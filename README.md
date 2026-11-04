<p align="center"><img src="docs/images/logo.png" width="480" alt="z0-worldcover-openfoam logo"></p>

# z0-worldcover-openfoam

![z0 per ground face from ESA WorldCover, TU Delft example](docs/images/z0_ground_3d.png)

Ground roughness **z0** for OpenFOAM atmospheric boundary layer (ABL) cases, taken from the
**ESA WorldCover 2021** land-cover map (10 m). Each ground face gets the z0 of the land under it:
trees, grass, water and built-up areas each have their own value, and buildings that are already
in the mesh are not counted twice. The same run writes ground **albedo** and **sensible fraction**
per point, which the solar ray-casting tool
[solar-raycast-openfoam](https://github.com/Jupihawx/solar-raycast-openfoam) reads.

- OpenFOAM (openfoam.com): tested with v2412 and v2606.
- Python ≥ 3.10.
- Status: the tools run and are stable on meshes from 0.5 M to 34 M cells. The z0 values are table
  values per land-cover class, not a calibrated or validated roughness map; you can change them in the
  table at the top of `make_z0.py` (see [How z0 is assigned](#how-z0-is-assigned)).

| File | What it does |
|---|---|
| `make_z0.py` | Downloads WorldCover tiles (once), samples them on a grid in mesh coordinates, writes z0, albedo, sensibleFraction and alignment figures |
| `tools/z0_faces.py` | Turns the z0 map into a **per-face list in polyMesh face order** for the wall functions |
| `tools/decompose_z0.py` | After `decomposePar`: writes each processor its own slice of the z0 list (works around a v2412 decomposePar issue) |
| `tools/check_wall_lists.py` | Puts the z0 list (and what the solver used) on the ground faces as a `.vtp` for ParaView |
| `example/` | TU Delft campus (OSM buildings, STL included) and a neutral `simpleFoam` case, 855 k cells |

---

## Conventions

| | |
|---|---|
| Mesh axes | x = East, y = North, z = up, metres, **ground at z = 0** |
| Projection | Transverse Mercator centred on (`lat0`, `lon0`), `+k=1`, WGS84: the geometry's origin is at (`lat0`, `lon0`), x/y in metres from it |
| Wind direction `dirFrom` | meteorological: the direction the wind blows **from**, degrees clockwise from North. `flowDir = (−sin dirFrom, −cos dirFrom, 0)` |
| `Zref` | 100 m (ERA5 100 m wind) |
| `z0In` | the **only** inlet roughness: one value per site, or per direction sector (see [Inlet roughness](#inlet-roughness-z0in)) |
| Patches | `Terrain` (ground, from blockMesh) and `Buildings` (snappyHexMesh) are the walls the tools work on. The example adds `inlet`, `outlet`, `sides`, `top` (box aligned with the wind) |

If your geometry was moved, rotated or scaled after export, the lat/lon no longer lines up: undo the
transform, or pass the lat/lon of the new origin.

## Install

```bash
python3 -m venv ~/venvs/z0
~/venvs/z0/bin/python -m pip install -r requirements.txt
```

WorldCover tiles (35–100 MB each) are downloaded on the first run and cached in `~/worldcover`.

## Inputs

- `<case>/geo/Buildings.stl`: the buildings that are in the mesh, in the same coordinates as the mesh.
- `--radius <m>` (command line): half-size of the area to cover, at least the farthest ground point from the origin.
- `--lat0`, `--lon0` (command line): the lat/lon of the mesh origin (0, 0). The example's `Allrun.sh`
  reads them from `geo/site.json`.

## Outputs

| Output | Use |
|---|---|
| `constant/boundaryData/Terrain/points` | Sample points, regular grid (`--dx`, default 20 m), z = 0 |
| `constant/boundaryData/Terrain/0/z0` | z0 per point |
| `constant/boundaryData/Terrain/0/albedo`, `0/sensibleFraction` | Ground materials for the solar tool (`--material-map`) |
| `constant/wallLists/Terrain.z0` | Per-face z0 list from `tools/z0_faces.py` |
| `geo/z0_map.png`, `geo/z0_map_zoom.png` | Land cover and z0 with the building roofs on top: **alignment check** |
| `geo/inlet_z0_sectors.csv` | Log-mean z0 of the land **upstream** of the domain per inflow direction (for `z0In`) |
| printed summary | Share of each land-cover class, inlet z0 per sector |

## Quickstart

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc        # your OpenFOAM install
cd <case>                                               # blockMesh (and snappyHexMesh) already run
python3 <repo>/make_z0.py --lat0 <lat0> --lon0 <lon0> --radius <m>   # -> constant/boundaryData/Terrain, geo/*.png
python3 <repo>/tools/z0_faces.py . --patch Terrain      # -> constant/wallLists/Terrain.z0
```

Check `geo/z0_map_zoom.png` (the building roofs must sit on red built-up pixels; parks, canals and
lakes must line up):

![Alignment check: land cover and z0 with the building roofs on top](docs/images/z0_map.png)

Then set the wall functions:

```
// 0/nut
Terrain
{
    type            atmNutkWallFunction;
    z0              #include "$FOAM_CASE/constant/wallLists/Terrain.z0";
    boundNut        true;
    value           uniform 0;
}
// 0/epsilon
Terrain
{
    type            atmEpsilonWallFunction;
    z0              #include "$FOAM_CASE/constant/wallLists/Terrain.z0";     // the same list
    value           uniform $epsilonInlet;
}
```

In parallel:

```bash
decomposePar -force
python3 <repo>/tools/decompose_z0.py . Terrain      # REQUIRED in v2412, see below
mpirun -np 4 simpleFoam -parallel
```

### The per-face list and the decomposePar step

The list belongs to one mesh: **rerun `z0_faces.py` after every
remesh** (it reads `constant/polyMesh`, or the latest `<time>/polyMesh` if snappyHexMesh wrote one;
`--time` to choose).

In OpenFOAM v2412, `decomposePar` does **not** split the `z0` list of `atmNutkWallFunction` /
`atmEpsilonWallFunction`: every processor receives the start of the full list, so the ground gets
scrambled roughness (a woven pattern in `nut`, `alphat`, `epsilon`). Ordinary lists (e.g. `value`,
`refValue`) are split correctly. `tools/decompose_z0.py` writes each processor its own slice using
`processor*/constant/polyMesh/faceProcAddressing` (ascii and binary fields). `reconstructPar` has the
same problem in reverse: the reconstructed `nut` z0 is scrambled, while the flow fields are fine.
Not checked against the OpenFOAM source; re-test on newer versions.

## How z0 is assigned

Default values below; you can change them in the `CLASSES` (z0) and `SURFACE` (albedo, sensible
fraction) tables at the top of `make_z0.py`, then rerun it (no remeshing needed).

| WorldCover class | z0 [m] (Wieringa 1992) | albedo | sensible fraction |
|---|---|---|---|
| 10 tree cover | 1.0 | 0.18 | 0.4 |
| 20 shrubland | 0.25 | 0.20 | 0.6 |
| 30 grassland | 0.03 | 0.23 | 0.5 |
| 40 cropland | 0.10 | 0.20 | 0.5 |
| 50 built-up | 1.5 | 0.15 | 0.9 |
| 60 bare | 0.005 | 0.20 | 0.7 |
| 70 snow/ice | 0.0002 | 0.60 | 0.3 |
| 80 water, and no data (open sea) | 0.0002 | 0.07 | 0.15 |
| 90 herbaceous wetland | 0.03 | 0.15 | 0.3 |
| 95 mangroves | 0.5 | 0.15 | 0.3 |
| 100 moss/lichen | 0.03 | 0.20 | 0.5 |

Albedo and sensible fraction (share of the absorbed shortwave that heats the air) are typical values
after Oke et al. (2017), *Urban Climates*.

**Explicit buildings are not counted twice.** Built-up pixels within `--explicit-dist` (50 m) of a
building in `Buildings.stl` get `--z0-explicit` (0.1 m), because the mesh already resolves those
buildings' drag. Their albedo and sensible fraction stay at the built-up values. Trees, grass and water between the buildings
keep their own z0: they are not in the mesh.

## Inlet roughness `z0In`

The ground z0 above covers the mesh. The wind arriving at the inlet has crossed the land **outside**
the mesh, so the inlet profile (U, k, epsilon) needs its own roughness: `z0In` in
`0/include/ABLConditions`.

`make_z0.py` estimates it for every wind direction and writes `geo/inlet_z0_sectors.csv`:

| `dir_from_deg` | `z0_logmean_m` |
|---|---|
| wind direction (from), every 11.25° | typical z0 of the land upwind, m |

For each direction it samples the land cover in a wedge from the domain edge to 5 km further upwind
(`--fetch`) and averages the z0 of those samples as `exp(mean(ln z0))`. The log average is used because
the wind profile depends on `ln z0`; a plain average would be dominated by a few rough pixels.

To use it: find the row closest to your wind direction and copy its value into `z0In`. Example, wind
from the south-west (the CSV row `225,0.08` gives 0.08 m; the example case uses 0.1):

```
// 0/include/ABLConditions
dirFrom   225;
z0In      0.08;
```

## Options of `make_z0.py`

| Option | Default | Meaning |
|---|---|---|
| `--case` | `.` | Case directory |
| `--lat0`, `--lon0` | required | Projection origin |
| `--stl` | `geo/Buildings.stl` | Explicit buildings, relative to the case |
| `--patch` | `Terrain` | Ground patch name: `boundaryData/<patch>` |
| `--radius` | required | Half-size of the sampled area, m: at least the farthest ground point from the origin |
| `--dx` | 20 | Sample spacing, m (WorldCover is 10 m: finer gains nothing) |
| `--explicit-dist` | 50 | Built-up pixels closer than this to a building count as explicit, m |
| `--z0-explicit` | 0.1 | z0 between explicit buildings, m |
| `--fetch` | 5000 | Upstream strip for the inlet z0, m |
| `--sectors` | 32 | Inflow sectors in the CSV |
| `--tile-dir` | `~/worldcover` | Tile cache |

## RAM and runtime

Example (TU Delft, 855 k cells, 127 k ground faces, 8 ranks): `make_z0.py` with `--dx 10` over a 2.4 km radius (~180 k points) plus `z0_faces.py` and `decompose_z0.py` took ~6 s together; the first run adds ~1 min for the tile download. Full example timings: [`example/README.md`](example/README.md).

`z0_faces.py` reads the full `points` and `faces` files into memory: plan RAM of a few times their size on disk.

## Example: TU Delft campus

`example/` is a real site (~1.3 km of TU Delft campus and housing from OpenStreetMap, 855 k cells) that
runs end to end in about 6.5 minutes on a desktop: see [`example/README.md`](example/README.md). The same case is
the example of the solar repository.

![Wind speed 10 m above ground, neutral simpleFoam run, wind from the south-west](docs/images/wind_10m.png)

## Limitations

- WorldCover 2021 is classified from satellite data at 10 m: class edges can be off, and 2021 may predate
  recent construction.
- z0 per class is a table value.
- `atmNutkWallFunction` uses `log((y + z0)/z0)`. With large z0 (≥ 1 m), keep the first cell reasonably
  tall there, or the wall function sits deep inside the roughness layer.
- No displacement height (`d = 0`).

## Data and attribution

- **ESA WorldCover 2021 v200**: © ESA WorldCover project 2021 / Contains modified Copernicus Sentinel
  data (2021) processed by ESA WorldCover consortium. Licence CC BY 4.0. Zanaga et al. (2022),
  doi:[10.5281/zenodo.7254221](https://doi.org/10.5281/zenodo.7254221). Tiles are downloaded at run time,
  not stored in this repository.
- **OpenStreetMap** (example geometry): © OpenStreetMap contributors, Open Database Licence (ODbL 1.0),
  <https://www.openstreetmap.org/copyright>. `example/case/geo/Buildings.stl` is derived from it and
  stays under the ODbL; if you share or adapt it, keep this attribution.

## Licence and citation

Code: MIT, see [LICENSE](LICENSE). To cite, see [CITATION.cff](CITATION.cff).
