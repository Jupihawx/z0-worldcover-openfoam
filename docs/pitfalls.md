# Pitfalls

Things that went wrong while using these tools on OpenFOAM v2412 (34 M cells, 9.3 M ground faces,
32 cores) and on the small example, and what fixed them.

## z0 list scrambled after decomposePar (v2412)

**Symptom:** a woven or striped pattern on the ground in `nut`, `alphat`, `epsilon` (and, in heated
runs, wall temperature and wall heat flux), following the processor boundaries.

**Cause:** `decomposePar` copies the `z0` list of `atmNutkWallFunction` / `atmEpsilonWallFunction` to
every processor without splitting it. Each processor reads the first N values of the global list.
Ordinary patch lists (`value`, `refValue`, `valueFraction`) are split correctly. Found with
`tools/check_wall_lists.py`: the input list was smooth, the `nut` z0 the solver used was scrambled.
Not checked against the OpenFOAM source; worth a report at develop.openfoam.com and a re-test on
newer versions.

**Fix:** run `tools/decompose_z0.py <case> Terrain` right after `decomposePar`. Confirmed on a
2000-iteration heated run: the pattern in `alphat` is gone. `reconstructPar` scrambles the reconstructed
`nut` z0 in the same way, so inspect processor files or the downstream fields, not reconstructed z0.

## Smooth wall functions on buildings

`nutkWallFunction` (smooth) gave `nut = 0` at low y+ on the building walls, which killed turbulent
transport there (and wall heat transfer in heated runs). Use `atmNutkWallFunction` /
`atmEpsilonWallFunction` with a small facade z0 (0.05 m) on `Buildings` too.

## Wall results depend on the first-cell height

Wall shear, `alphat` and wall heat flux go through the wall functions
(`nut ∝ y / ln((y + z0)/z0)`, `epsilon ∝ 1/y`), so they change where the ground cell height changes:
the edge of a refinement box is visible in `alphat`, `epsilon` and wall temperature. Interpret wall
results only where the ground cells have one height (inside the dense box), or refine the whole ground
patch to one level, or add snappy layers.

## Slivers against the ground

STL details smaller than the building cell size (notches, narrow courtyards and corridors) produce
sliver cells against the ground. Options: clean the footprints (2D close/open), `gapLevel`, or finer
refinement levels. The example drops footprints below 15 m² and simplifies outlines to 0.3 m.

## Large z0 and low first cells

`atmNutkWallFunction` uses `log((y + z0)/z0)`. With z0 ≥ 1 m (trees, built-up) and a first cell centre
of 1–2 m, the wall function sits inside the roughness layer. Keep the first cell reasonably tall where
z0 is large; between explicit buildings `make_z0.py` already lowers z0 to 0.1 m.

## Geometry not where WorldCover says it is

If the roofs in `geo/z0_map_zoom.png` do not sit on the red built-up pixels, the lat/lon is wrong or the
geometry was moved, rotated or scaled after export. A shift of a few metres is within WorldCover's
accuracy; tens of metres is not.

## Inlet z0 from local water

The sector CSV averages a 5 km strip just outside the domain. A reservoir or canal in that strip gives a
low z0 that does not represent the tens of kilometres that shape the inlet profile up to 100 m.
Keep one `z0In` per site unless the fetch really changes (open sea on one side).

## Windows line endings

Shell scripts edited on Windows fail in WSL/Linux with `$'\r': command not found`.
`sed -i 's/\r$//' *.sh`. The repository's `.gitattributes` keeps them LF.
