#!/bin/bash
#PBS -N Lz88L5w
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1
#PBS -l walltime=04:00:00
#PBS -j oe

set -euo pipefail
ROOT=/lustre/home/users/ewu/vb_gcmc/MD/stage_vacf_tail_8_8_L2L5_4rep_weaknh_nomom_10ns_100fs_20260821
WORK=${ROOT}/analysis_fullwater_Lz_20260824
PY=/lustre/home/users/ewu/.conda/envs/HB_analysis/bin/python
mkdir -p "${WORK}/fullwater_Lz_88L5_4rep"
cd "${WORK}"
"${PY}" audit_fullwater_orbital_Lz.py \
  --dumps \
  "${ROOT}/5L/rep1/LONG_WEAKNH_NOMOM_8_8_L5_rep1_100fs_10ns.water_100fs_10ns.dump" \
  "${ROOT}/5L/rep2/LONG_WEAKNH_NOMOM_8_8_L5_rep2_100fs_10ns.water_100fs_10ns.dump" \
  "${ROOT}/5L/rep3/LONG_WEAKNH_NOMOM_8_8_L5_rep3_100fs_10ns.water_100fs_10ns.dump" \
  "${ROOT}/5L/rep4/LONG_WEAKNH_NOMOM_8_8_L5_rep4_100fs_10ns.water_100fs_10ns.dump" \
  --output "${WORK}/fullwater_Lz_88L5_4rep" --stride 10
printf 'CCFEP full-water Lz audit finished successfully.\n' > "${WORK}/fullwater_Lz_88L5_4rep/PBS_FINISHED.txt"
