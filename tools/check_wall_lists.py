#!/usr/bin/env python3
"""Put the per-face lists the solver actually uses on the wall faces, for ParaView.

    python3 tools/check_wall_lists.py <case> <time> [--patches Terrain Buildings]
                                       (reconstructed time, or 0 before running)

Reads from <time>/{T,nut,epsilon,alphat} whatever exists of: refValue, valueFraction, value (T),
z0 and value (nut, epsilon), value (alphat), plus the input lists constant/wallLists/<patch>.z0 and
solarIllumination/T_solar/<patch>.{Tsolair,z0}, and writes <case>/check_<time>.vtp with one cell
array per list. Open it in ParaView: a list that shows a woven / scrambled pattern while the input
list (file_*) is smooth was put on the faces in the wrong order (e.g. the decomposePar z0 bug).
The surface is built like solar_illumination.py's illumination.vtp, so the two compare face by face.
Needs pyvista.
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pyvista as pv


def load_wall_surface(foam_file: Path, patch_names: list) -> pv.PolyData:
    """Wall patches merged patch by patch, faces in polyMesh order (as in solar_illumination.py)."""
    reader = pv.OpenFOAMReader(str(foam_file))
    reader.disable_all_cell_arrays()                # geometry only
    reader.disable_all_patch_arrays()
    for p in patch_names:
        reader.enable_patch_array(f"patch/{p}")
    reader.cell_to_point_creation = False
    bnd = reader.read()["boundary"]
    surfs = []
    for i, p in enumerate(patch_names):
        try:
            s = bnd[p].extract_surface(algorithm=None)   # newer pyvista
        except TypeError:
            s = bnd[p].extract_surface()                 # older pyvista: no 'algorithm'
        s.cell_data.clear()
        s.point_data.clear()
        s.cell_data["patchId"] = np.full(s.n_cells, i, np.int32)
        surfs.append(s)
    return surfs[0] if len(surfs) == 1 else pv.merge(surfs, merge_points=False)


def patch_block(data: bytes, patch: str) -> bytes:
    """Bytes of boundaryField { <patch> { ... } } (brace-matched, binary-safe)."""
    bf = data.index(b"boundaryField")
    m = re.compile(rb"\n\s*" + patch.encode() + rb"\s*\n\s*\{").search(data, bf)
    if not m:
        raise KeyError(patch)
    i, depth = m.end(), 1
    while depth:
        c = data[i:i + 1]
        if c == b"(":                     # skip list contents (binary may hold braces)
            i = skip_list(data, i)
            continue
        depth += {b"{": 1, b"}": -1}.get(c, 0)
        i += 1
    return data[m.end():i - 1]


def skip_list(data: bytes, i: int) -> int:
    """i at '(' of a list; return index after its ')'. Handles binary N( ... )."""
    m = re.search(rb"(\d+)\s*$", data[max(0, i - 20):i])
    j = data.index(b")", i)
    if m:                               # possible binary list N(<raw>)
        n = int(m.group(1))
        if data[i + 1 + 8 * n:i + 2 + 8 * n] == b")":
            return i + 2 + 8 * n
    return j + 1


def read_entry(block: bytes, key: str):
    """Scalar list (or uniform value) of 'key' in a patch block; None if absent."""
    m = re.search(rb"(?m)^\s*" + key.encode() + rb"\s+(uniform|nonuniform)", block)
    if not m:
        return None
    rest = block[m.end():]
    if m.group(1) == b"uniform":
        return float(re.match(rb"\s*([-+0-9.eE]+)", rest).group(1))
    m2 = re.match(rb"\s*List<scalar>\s*(\d+)\s*\(", rest)
    n, s = int(m2.group(1)), m2.end()
    if rest[s + 8 * n:s + 8 * n + 1] == b")" and n > 0:     # binary
        try:
            return np.frombuffer(rest[s:s + 8 * n], "<f8").copy()
        except ValueError:
            pass
    return np.array(rest[s:rest.index(b")", s)].split(), float)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", type=Path)
    ap.add_argument("time")
    ap.add_argument("--patches", nargs="+", default=["Terrain", "Buildings"])
    a = ap.parse_args()
    case, time, PATCHES = a.case, a.time, a.patches
    (case / "case.foam").touch()
    surf = load_wall_surface(case / "case.foam", PATCHES)
    pid = np.asarray(surf.cell_data["patchId"])
    arrays = {}
    for fld, keys in (("T", ["refValue", "valueFraction", "value"]), ("nut", ["z0", "value"]),
                      ("epsilon", ["z0"]), ("alphat", ["value"])):
        if not (case / time / fld).exists():
            continue
        data = (case / time / fld).read_bytes()
        for i, p in enumerate(PATCHES):
            blk = patch_block(data, p)
            for k in keys:
                v = read_entry(blk, k)
                if v is None:
                    continue
                name = f"{fld}_{k}"
                arr = arrays.setdefault(name, np.full(surf.n_cells, np.nan))
                sel = pid == i
                if np.ndim(v) and len(v) != sel.sum():
                    sys.exit(f"{fld}/{p}/{k}: {len(v)} values, {sel.sum()} faces")
                arr[sel] = v
    for name, pattern in (("file_z0", "constant/wallLists/{p}.z0"),
                          ("file_solar_z0", "solarIllumination/T_solar/{p}.z0"),
                          ("file_Tsolair", "solarIllumination/T_solar/{p}.Tsolair")):
        for i, p in enumerate(PATCHES):
            f = case / pattern.format(p=p)
            if f.exists():
                v = read_entry(b"x " + f.read_bytes(), "x")
                arr = arrays.setdefault(name, np.full(surf.n_cells, np.nan))
                arr[pid == i] = v
    for name, arr in arrays.items():
        surf.cell_data[name] = arr
        print(f"{name:20s} min {np.nanmin(arr):10.4g}  max {np.nanmax(arr):10.4g}")
    out = case / f"check_{time}.vtp"
    surf.save(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
