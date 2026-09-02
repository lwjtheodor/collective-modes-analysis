#!/bin/bash
#PBS -N cjjT200
#PBS -q i8cpu
#PBS -l select=1:ncpus=2:mpiprocs=2:ompthreads=1
#PBS -l walltime=12:00:00
#PBS -j oe
set -euo pipefail
cd "$PBS_O_WORKDIR"
P="$(command -v python3 || command -v python)"
ROOT=/lustre/home/users/ewu/vb_gcmc/MD/transverse_velocity_5L_10fs_weakNH_nomom_4chirality_20260808/stage_transverse_cjj_kmin_5L_200ps_20260818
"$P" analyze_transverse_cjj_kmin_200ps.py --root "$ROOT" --out "$ROOT/analysis_kmin_cjj_200ps"
echo "CJJ analysis finished successfully."
