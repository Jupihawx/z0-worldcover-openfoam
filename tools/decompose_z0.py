#!/usr/bin/env python3
"""Write the correct per-face z0 into processor*/0/{nut,epsilon} after decomposePar.

    python3 tools/decompose_z0.py <case> <patch> [<patch> ...] [--list <path pattern>]

decomposePar (OpenFOAM v2412) copies the z0 list of atmNutkWallFunction /
atmEpsilonWallFunction to every processor without splitting it, so each
processor gets the wrong values. This takes the full list (polyMesh face order)
and, through processor*/constant/polyMesh/faceProcAddressing, writes each
processor its own slice, in ascii or binary like the field file.

--list: path of the full list relative to the case, {patch} is replaced.
  default constant/wallLists/{patch}.z0      (tools/z0_faces.py)
  solar repo: solarIllumination/T_solar/{patch}.z0   (solar_illumination.py --material-map)
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np

FIELDS = ["nut", "epsilon", "omega"]   # those present in processor*/0 are patched


def patch_range(boundary: Path, patch: str):
    txt = boundary.read_text(errors="replace")
    m = re.search(r"\n\s*" + patch + r"\s*\{([^}]*)\}", txt)
    if not m:
        return None
    start = int(re.search(r"startFace\s+(\d+)", m.group(1)).group(1))
    n = int(re.search(r"nFaces\s+(\d+)", m.group(1)).group(1))
    return start, n


def is_binary(data: bytes) -> bool:
    return re.search(rb"format\s+binary", data[:2000]) is not None


def read_list(data: bytes, dtype: str):
    """First top-level list after the FoamFile header: ascii or binary N(raw)."""
    body = data[data.index(b"}", data.index(b"FoamFile")) + 1:]
    m = re.search(rb"\n(\d+)\s*\(", body)
    n, s = int(m.group(1)), m.end()
    if is_binary(data):
        return np.frombuffer(body[s:s + np.dtype(dtype).itemsize * n], dtype).copy()
    return np.array(body[s:body.index(b")", s)].split(), dtype)


def skip_list(data: bytes, i: int) -> int:
    """i at '('; index after the matching ')' (binary N(raw) aware)."""
    m = re.search(rb"(\d+)\s*$", data[max(0, i - 20):i])
    if m:
        n = int(m.group(1))
        if data[i + 1 + 8 * n:i + 2 + 8 * n] == b")":
            return i + 2 + 8 * n
    return data.index(b")", i) + 1


def replace_z0(data: bytes, patch: str, z0: np.ndarray) -> bytes:
    bf = data.index(b"boundaryField")
    m = re.compile(rb"\n\s*" + patch.encode() + rb"\s*\n\s*\{").search(data, bf)
    k = re.compile(rb"\n\s*z0\s").search(data, m.end())
    i = k.end()
    while data[i:i + 1] != b";":                  # end of the z0 entry
        i = skip_list(data, i) if data[i:i + 1] == b"(" else i + 1
    if is_binary(data):
        new = (f"nonuniform List<scalar> {len(z0)}(".encode()
               + z0.astype("<f8").tobytes() + b")")
    else:
        new = (f"nonuniform List<scalar> \n{len(z0)}\n(\n".encode()
               + "\n".join(f"{v:.6g}" for v in z0).encode() + b"\n)\n")
    return data[:k.end()] + new + data[i:]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", type=Path)
    ap.add_argument("patches", nargs="+")
    ap.add_argument("--list", default="constant/wallLists/{patch}.z0")
    a = ap.parse_args()
    case = a.case
    procs = sorted(case.glob("processor*"), key=lambda p: int(p.name[9:]))
    if not procs:
        sys.exit("no processor* directories: run decomposePar first")
    for patch in a.patches:
        g0, gn = patch_range(case / "constant/polyMesh/boundary", patch)
        z0_all = np.array((case / a.list.format(patch=patch))
                          .read_text().split("(", 1)[1].split(")")[0].split(), float)
        if len(z0_all) != gn:
            sys.exit(f"{patch}.z0 has {len(z0_all)} values, mesh patch has {gn} faces")
        done = 0
        for proc in procs:
            r = patch_range(proc / "constant/polyMesh/boundary", patch)
            if r is None or r[1] == 0:
                continue
            p0, pn = r
            fpa = read_list((proc / "constant/polyMesh/faceProcAddressing").read_bytes(), "<i4")
            glob = np.abs(fpa[p0:p0 + pn]) - 1 - g0      # +-(global face + 1)
            if glob.min() < 0 or glob.max() >= gn:
                sys.exit(f"{proc.name}: {patch} faces map outside the global patch")
            z0 = z0_all[glob]
            for fld in FIELDS:
                f = proc / "0" / fld
                if not f.exists():
                    continue
                f.write_bytes(replace_z0(f.read_bytes(), patch, z0))
            done += pn
        print(f"{patch}: wrote z0 for {done} of {gn} faces in {len(procs)} processors "
              f"(min {z0_all.min():.3g}, max {z0_all.max():.3g})")
        if done != gn:
            sys.exit("face count mismatch")


if __name__ == "__main__":
    main()
