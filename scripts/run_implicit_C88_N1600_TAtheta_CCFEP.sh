#!/bin/bash
#PBS -N imp88th
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1
#PBS -l walltime=16:00:00
#PBS -j oe

set -euo pipefail
ROOT=/lustre/home/users/ewu/vb_gcmc/MD/N1600_weakNH_6ns
WORK=${ROOT}/analysis_TAtheta_linewidth_20260824
PY=/lustre/home/users/ewu/.conda/envs/HB_analysis/bin/python
mkdir -p "${WORK}"
cd "${WORK}"
for rep in 1 2 3 4; do
  "${PY}" analyze_implicitCNT_TAtheta_linewidth.py \
    --dump "${ROOT}/rep${rep}/production_100fs.dump" \
    --output "${WORK}/per_replica/rep${rep}" \
    --nmax 20 --nperseg 16384 --omega-fit-max 0.20 --primary-nmax 10 \
    --temperature-K 330 --case-label "Implicit CNT (8,8), N1600, rep${rep}"
done
printf 'C88 N1600 implicit-CNT TA_theta per-replica analysis finished successfully.\n' > "${WORK}/PBS_FINISHED.txt"
