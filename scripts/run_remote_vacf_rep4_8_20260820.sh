#!/bin/bash
#PBS -N vacf_L5_allorg
#PBS -q i8cpu
#PBS -l select=1:ncpus=32:mpiprocs=32:ompthreads=1
#PBS -l walltime=12:00:00
#PBS -j oe
set -euo pipefail

ROOT="/lustre/home/users/ewu/vb_gcmc/MD/stage_vacf_tail_8_8_L2L5_extra5rep_weaknh_nomom_1ns_10fs_20260820_v2"
ANALYZER="$ROOT/scripts/analyze_vacf_tail.py"
OUT="$ROOT/vacf_analysis_allorigins_20260821/outputs"
mkdir -p "$OUT"
command -v python3
python3 -c 'import numpy; print("numpy", numpy.__version__)'

analyze_length() {
  local L="$1"
  for REP in 6 7 8; do
    local CASE="VACF_8_8_L${L}_extra_rep${REP}_1ns_10fs"
    local DUMP="$ROOT/${L}L/rep${REP}/${CASE}.oxygen_id_z_vz_10fs_1ns.dump"
    python3 "$ANALYZER" --dump "$DUMP" --timestep-ps 0.0005 --max-lag-ps 100 \
      --nblocks 0 --case-id "$CASE" --out "$OUT/${CASE}.csv"
  done
}
# The first all-origin attempt completed 2L--4L reps4--8 and 5L reps4--5.
# L5 requires a larger per-job memory allocation, so the remaining three
# cases run serially on this dedicated 32-CPU allocation.
analyze_length 5
python3 - <<'PY'
from pathlib import Path
import json
out=Path('/lustre/home/users/ewu/vb_gcmc/MD/stage_vacf_tail_8_8_L2L5_extra5rep_weaknh_nomom_1ns_10fs_20260820_v2/vacf_analysis_allorigins_20260821/outputs')
csv=sorted(out.glob('*.csv')); meta=sorted(out.glob('*.json'))
if len(csv)!=20 or len(meta)!=20:
    raise SystemExit(f'incomplete outputs: csv={len(csv)} json={len(meta)}')
(out.parent/'COMPLETE.json').write_text(json.dumps({'csv':len(csv),'json':len(meta),'definition':'all-origin full-trajectory lab and instantaneous oxygen-COM-subtracted VACF; max lag 100 ps; no block averaging'},indent=2))
PY
