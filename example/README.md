# Example: TU Delft campus (neutral ABL)

A real site, so that WorldCover gives meaningful land cover: the TU Delft campus and the housing to its
east, ~1.3 km × 1.3 km of buildings, origin 52.000825 N, 4.377 E. The land cover mixes built-up areas,
trees, lawns and canals, so z0 varies from 0.0002 m to 1.5 m.

| | |
|---|---|
| Geometry | `geo/Buildings.stl`, shipped: OSM buildings (extent lon 4.36745–4.38655, lat 51.99439–52.00726), cleaned (small sheds and broken triangles removed). 406 solids, tallest 90 m (EWI) |
| Domain | box aligned with the wind: buildings within 1000 m of the origin, plus 300 m upstream and to the sides, 900 m downstream, 300 m high. **Smaller than the COST 732 guideline** (5 H / 15 H / 6 H with H = 90 m) to keep the example short |
| Mesh | blockMesh 20 m cells, rotated to `dirFrom` with `transformPoints`. snappyHexMesh: the whole ground at 10 m (first cell 10 m tall everywhere, no refinement boxes), 5 m within 10 m of buildings, 2.5 m on building surfaces |
| Inflow | 5 m/s at 100 m from 225° (south-west), illustrative, not a measured day. `z0In` 0.1 m from `geo/inlet_z0_sectors.csv` (225° sector: 0.11 m) |
| Ground | z0 per face from WorldCover (`make_z0.py` + `tools/z0_faces.py`) |
| Solver | `simpleFoam`, k-epsilon, rough atm wall functions on ground and buildings |

## Step by step

Linux or WSL2 (Ubuntu). Commands start in the repository root.

**1. OpenFOAM (openfoam.com; tested with v2412 and v2606).** Install it from <https://develop.openfoam.com/Development/openfoam/-/wikis/precompiled>
(Ubuntu: `apt` package `openfoam2412-default`), then load it in every new shell:

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
```

**2. Python environment:**

```bash
python3 -m venv ~/venvs/z0
~/venvs/z0/bin/python -m pip install -r requirements.txt
```

**3. Run the example** (internet needed on the first run: one 55 MB WorldCover tile):

```bash
cd example/case
PYTHON=~/venvs/z0/bin/python NPROCS=4 bash Allrun.sh
```

`NPROCS` = MPI ranks (1 = serial). `Allrun.sh` does, with a log per step in `logs/`:

| # | Step | Writes |
|---|---|---|
| 1 | read the origin from `geo/site.json` | |
| 2 | `blockMesh`, `transformPoints` (box turned to `dirFrom`), `snappyHexMesh`, `checkMesh` | `constant/polyMesh` |
| 3 | `make_z0.py`: WorldCover → z0, albedo, sensible fraction | `constant/boundaryData/Terrain`, `geo/z0_map*.png`, `geo/inlet_z0_sectors.csv` |
| 4 | `tools/z0_faces.py`: z0 per ground face | `constant/wallLists/Terrain.z0` |
| 5 | `decomposePar`, `tools/decompose_z0.py` (v2412 z0 workaround) | `processor*/` |
| 6 | `simpleFoam` (300 iterations), `reconstructPar` | `300/`, `postProcessing/probes/` |

**4. Look at the results**

```bash
touch case.foam          # open in ParaView: latest time, patches Terrain and Buildings
```

See [Check](#check) below for what to look for. `bash Allclean.sh` removes everything generated
(keeps `geo/Buildings.stl` and `geo/site.json`).

## Runtime (measured)

16-core WSL2 machine, 8 MPI ranks, OpenFOAM v2606, WorldCover tile already cached:

| Step | Time |
|---|---|
| WorldCover tile (first run only) | ~1 min, one 55 MB tile, cached in `~/worldcover` |
| blockMesh + rotate + snappyHexMesh (serial) | ~2.5 min, 855 k cells; `checkMesh`: 8 faces with skewness up to 4.8 (limit 4), otherwise OK |
| `make_z0.py` + `z0_faces.py` + `decompose_z0.py` | ~6 s (127 k ground faces) |
| simpleFoam | 0.79 s per iteration. In a 500-iteration run the probes were steady to 4 digits from ~200, so `endTime` is 300 (~4 min) |
| **Total** | **~6.5 min** (401 s) |

## Check

- `geo/z0_map_zoom.png`: roofs (black) on built-up pixels (red); lawns and canals where they are on a map.
- `logs/make_z0.log`: land-cover shares and the inlet z0 per sector.
- `logs/decompose_z0.log`: `wrote z0 for N of N faces`.
- ParaView (`touch case.foam`): `nut` on `Terrain` follows the land cover with no stripes along processor
  boundaries.
- `postProcessing/probes/0/U`: the probes (upstream at 10 m, above the site at 100 m) settle.

## Another wind direction

Change `dirFrom` (and `z0In` if the fetch changes, see `geo/inlet_z0_sectors.csv`) in
`0/include/ABLConditions`, delete `constant/polyMesh` and rerun: the box is rotated to the new direction,
and the z0 list must be regenerated for the new mesh (`Allrun.sh` does both).

## Your own site

1. A building STL in metres, z up, ground at z = 0, x = east and y = north from a known origin
   (Transverse Mercator centred on that lat/lon). Remove buildings smaller than the mesh resolution;
   they only make slivers.
2. Put the origin in `geo/site.json`, set `Rb` (radius of the building area) and the margins in
   `system/blockMeshDict`, and a `locationInMesh` above the buildings in `system/snappyHexMeshDict`.

## Licence of the geometry

`geo/Buildings.stl` is derived from OpenStreetMap: © OpenStreetMap contributors, available under the
Open Database Licence (ODbL 1.0), <https://www.openstreetmap.org/copyright>. It is a derived database
under the ODbL, separate from the MIT-licensed code; if you share or adapt it, keep this attribution and
the ODbL terms.
