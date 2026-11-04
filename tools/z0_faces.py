#!/usr/bin/env python3
"""Per-face z0 list for a wall patch, in polyMesh face order (the recommended route).

    python3 tools/z0_faces.py <case> [--patch Terrain] [--time <t>]

Reads the z0 map written by make_z0.py (constant/boundaryData/<patch>/{points,0/z0}) and the mesh
(<time>/polyMesh if given or present, else constant/polyMesh; ascii or binary), gives every face of the
patch the value of the nearest map point (vertex-mean face centre) and writes

    constant/wallLists/<patch>.z0       nonuniform List<scalar> N ( ... )

for  z0  #include "$FOAM_CASE/constant/wallLists/<patch>.z0";  in 0/nut and 0/epsilon.
Rerun after every remesh (the list belongs to one mesh). In parallel, run tools/decompose_z0.py
after decomposePar: OpenFOAM v2412 decomposePar does not split this list (README).
"""
import argparse
import re
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def foam_body(path: Path):
    data = path.read_bytes()
    binary = re.search(rb"format\s+binary", data[:3000]) is not None
    return data[data.index(b"}", data.index(b"FoamFile")) + 1:], binary


def read_lists(path: Path, dtypes):
    """Consecutive top-level lists N(...) of a polyMesh file; dtypes per list (binary)."""
    body, binary = foam_body(path)
    out, pos = [], 0
    for dt in dtypes:
        m = re.compile(rb"(\d+)\s*\(").search(body, pos)
        n, s = int(m.group(1)), m.end()
        if binary:
            size = np.dtype(dt[0]).itemsize * dt[1] * n
            out.append(np.frombuffer(body[s:s + size], dt[0]).reshape(n, dt[1]) if dt[1] > 1
                       else np.frombuffer(body[s:s + size], dt[0]).copy())
            pos = s + size + 1
        else:
            e = s
            depth = 1
            while depth:                                     # matching ')', items may be (x y z)
                c = body[e:e + 1]
                depth += {b"(": 1, b")": -1}.get(c, 0)
                e += 1
            vals = np.array(body[s:e - 1].replace(b"(", b" ").replace(b")", b" ").split(), dt[0])
            out.append(vals.reshape(n, dt[1]) if dt[1] > 1 else vals)
            pos = e
    return out


def patch_range(boundary: Path, patch: str):
    txt = boundary.read_text(errors="replace")
    m = re.search(r"\n\s*" + patch + r"\s*\{([^}]*)\}", txt)
    if not m:
        raise SystemExit(f"patch {patch} not in {boundary}")
    return (int(re.search(r"startFace\s+(\d+)", m.group(1)).group(1)),
            int(re.search(r"nFaces\s+(\d+)", m.group(1)).group(1)))


def read_foam_list(path: Path, ncomp=1):
    text = path.read_text()
    body = text[text.index("(", text.index("}") if "FoamFile" in text else 0) + 1:text.rindex(")")]
    vals = np.array(body.replace("(", " ").replace(")", " ").split(), float)
    return vals.reshape(-1, ncomp) if ncomp > 1 else vals


def mesh_dir(case: Path, time):
    if time is not None:
        return case / time / "polyMesh"
    times = [(float(d.name), d / "polyMesh") for d in case.iterdir()
             if re.fullmatch(r"[0-9.eE+-]+", d.name) and (d / "polyMesh" / "boundary").exists()]
    return max(times)[1] if times else case / "constant" / "polyMesh"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", type=Path)
    ap.add_argument("--patch", default="Terrain")
    ap.add_argument("--time", default=None, help="time folder holding the mesh (default: latest, else constant)")
    ap.add_argument("--field", default="z0", help="map field in boundaryData/<patch>/0 (z0, albedo, ...)")
    a = ap.parse_args()

    poly = mesh_dir(a.case, a.time)
    start, n = patch_range(poly / "boundary", a.patch)
    (pts,) = read_lists(poly / "points", [("<f8", 3)])
    offsets, labels = read_lists(poly / "faces", [("<i4", 1), ("<i4", 1)])
    # face centres of the patch: mean of the face vertices
    o0, o1 = offsets[start:start + n], offsets[start + 1:start + n + 1]
    counts = o1 - o0
    idx = np.repeat(o0, counts) + (np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts))
    owner = np.repeat(np.arange(n), counts)
    centres = np.zeros((n, 3))
    np.add.at(centres, owner, pts[labels[idx]])
    centres /= counts[:, None]

    bd = a.case / "constant" / "boundaryData" / a.patch
    mp = read_foam_list(bd / "points", 3)
    vals = read_foam_list(bd / "0" / a.field)
    d, near = cKDTree(mp[:, :2]).query(centres[:, :2])
    out = a.case / "constant" / "wallLists" / f"{a.patch}.{a.field}"
    out.parent.mkdir(parents=True, exist_ok=True)
    v = vals[near]
    with open(out, "w") as f:
        f.write(f"nonuniform List<scalar> \n{n}\n(\n")
        f.write("\n".join(f"{x:.6g}" for x in v))
        f.write("\n)\n")
    print(f"mesh {poly.relative_to(a.case)}: {a.patch} {n} faces, map {len(mp)} points, "
          f"max distance to map point {d.max():.1f} m")
    print(f"{a.field}: min {v.min():.4g}, max {v.max():.4g}, face-mean {v.mean():.4g} -> {out}")


if __name__ == "__main__":
    main()
