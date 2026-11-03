#!/bin/bash
# TU Delft example, neutral ABL (simpleFoam, OpenFOAM v2412 openfoam.com):
#   geo/Buildings.stl (OSM via Blender + Blosm) -> blockMesh + snappyHexMesh -> z0 from ESA WorldCover
#   -> per-face z0 list -> solve
#
# Usage:  bash Allrun.sh            (from anywhere; source OpenFOAM's etc/bashrc first)
#   PYTHON=/path/to/venv/bin/python   Python with requirements.txt installed (default: python3)
#   NPROCS=4                          MPI ranks for the solver (1 = serial)
set -e
cd "$(dirname "$0")"
REPO=$(cd ../.. && pwd)
PY=${PYTHON:-python3}
nprocs=${NPROCS:-4}
mkdir -p logs

# 1. Geometry: geo/Buildings.stl and its projection origin geo/site.json (shipped with the example)
LAT0=$($PY -c "import json; print(json.load(open('geo/site.json'))['lat0'])")
LON0=$($PY -c "import json; print(json.load(open('geo/site.json'))['lon0'])")

# 2. Mesh: box in the wind frame (blockMesh), rotated about z so x points downwind, + buildings (snappy)
#    flowDir = (-sin dirFrom, -cos dirFrom): rotation angle = atan2(-cos dirFrom, -sin dirFrom)
DIRFROM=$(foamDictionary -entry dirFrom -value 0/include/ABLConditions)
ANGLE=$($PY -c "import math; d=math.radians($DIRFROM); print(round(math.degrees(math.atan2(-math.cos(d), -math.sin(d))), 4))")
mkdir -p constant/triSurface
cp geo/Buildings.stl constant/triSurface/
surfaceFeatureExtract > logs/surfaceFeatureExtract.log
blockMesh > logs/blockMesh.log
transformPoints -rotate-angle "((0 0 1) $ANGLE)" > logs/transformPoints.log
snappyHexMesh -overwrite > logs/snappyHexMesh.log
checkMesh > logs/checkMesh.log || echo "checkMesh reported problems, see logs/checkMesh.log"

# 3. z0 (and albedo, sensibleFraction) from WorldCover -> constant/boundaryData/Terrain, geo/*.png
#    --radius: the rotated box reaches ~2300 m from the centre (blockMesh.log has the unrotated box)
$PY "$REPO/make_z0.py" --case . --lat0 "$LAT0" --lon0 "$LON0" --dx 10 --radius 2400 | tee logs/make_z0.log
#    per-face list in polyMesh order -> constant/wallLists/Terrain.z0 (read by 0/nut and 0/epsilon)
$PY "$REPO/tools/z0_faces.py" . --patch Terrain | tee logs/z0_faces.log

# 4. Solve
if [ "$nprocs" -gt 1 ]; then
    sed -i "s/^numberOfSubdomains.*/numberOfSubdomains $nprocs;/" system/decomposeParDict
    rm -rf processor*
    decomposePar -force > logs/decomposePar.log
    # decomposePar (v2412) copies the z0 list unsplit: write each processor its own slice
    $PY "$REPO/tools/decompose_z0.py" . Terrain | tee logs/decompose_z0.log
    mpirun -np $nprocs simpleFoam -parallel > logs/simpleFoam.log
    reconstructPar -latestTime > logs/reconstructPar.log
else
    simpleFoam > logs/simpleFoam.log
fi
echo "Done. Open the case in ParaView (touch case.foam); probes in postProcessing/probes/"
