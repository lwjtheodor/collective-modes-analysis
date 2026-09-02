# Canonical collective-dynamics commands

The entry point is:

```powershell
py scripts/collective_modes_cli.py <command> ...
```

On CCFEP, upload only `scripts/collective_modes/` and
`scripts/collective_modes_cli.py` into the analysis result/stage directory,
then invoke it with:

```bash
/lustre/home/users/ewu/.conda/envs/HB_analysis/bin/python collective_modes_cli.py <command> ...
```

Every command writes readable CSV plus `metadata.json`.  Its metadata must be
copied into the authority package under
`results/collective_mode_response/<topic>/<date>/`; completion still requires
the result-package README, source provenance and QA gates in `AGENT.md`.

## Required declaration and automatic detection

`audit` samples each dump header and reports fields, atom types, molecule IDs,
xyz/velocity/image availability, fluid-representation inference and supported observables in
`dump_capabilities.csv`.  It also reports a wall-model inference with a
confidence/basis field.

It is deliberately **not** allowed to infer implicit versus explicit physics
from a water-only dump: an explicit CNT can be absent from the output.  A
site-only dump also cannot prove oxygen-only water versus argon. A water-only
file therefore becomes `ambiguous_water_only` unless `--wall-model` is declared
or declared `--cnt-types` are actually present; a site-only fluid becomes
`single_site_fluid_ambiguous` unless `--fluid-kind water|oxygen_only|argon` is
declared. Filename hints
such as `implicit` are recorded as non-authoritative hints only.

For an implicit cylindrical analysis, `--rcnt-A` must be the CNT field radius
from protocol metadata, never a density-peak estimate.  For flexible explicit
CNT, the present command intentionally refuses a wall-relative calculation
until a CNT-frame extractor with actual CNT atoms/velocities is supplied.

## Input-capability contract

| Input fields | Commands enabled |
|---|---|
| `id,type,z[,iz],vz` | axial ISF, axial current, z-VACF, VACF-MSD-alpha |
| `id,type,x,y,z,vx,vy,vz` | all above plus cylindrical `Jr/Jtheta/L/Tinplane/Tr`, ordered cross kernels, r/theta VACF |
| `id,mol,type,x,y,z[,ix,iy,iz],vx,vy,vz` | additionally supports molecular/rotational plugins after their separate gate |
| CNT atom `x,y,z,vx,vy,vz` with declared CNT types | required later for flexible-CNT wall-relative frame |

## Commands and principal CSV products

| Command | Main products |
|---|---|
| `audit` | `dump_capabilities.csv` |
| `isf` | `isf_per_replica.csv`, `isf_ensemble_mean_sem.csv`: `F_total`, `F_self`, `F_distinct` for every `(n,m,lag)` |
| `current` | `current_per_replica.csv`, ensemble CJJ, ordered cross table, periodogram table; separate `Jz`, `Jr`, `Jtheta`, `L`, `Tinplane`, `Tr` |
| `vacf` | per-rep/ensemble VACF and `msd_alpha_from_vacf_*`; lab/selected-COM/wall-relative is explicit |
| `fit-current` | `current_mode_DHO_parameters.csv`; a row for each `(n,m)` retaining `Gamma`, `omega`, `a`, `b`, standard errors and `R2` |
| `construct` | per-mode and summed `W*Fs*Phi_J`, optional direct-VACF residual table; external static `W(n,m)` only, no amplitude fit |
| `plot` | a minimal data-readable PNG with zero line and columns named from the CSV |

The retained current model is

\[
\Phi_J(k_z,m;t)=\exp[-\Gamma(k_z,m)t]
 [a(k_z,m)\cos\omega(k_z,m)t+b(k_z,m)\sin\omega(k_z,m)t],
\]

and constructibility is written without a fitted global amplitude as

\[
C_{vv}^{\rm construct}(t)=\sum_{n,m}W(n,m)F_s(n,m;t)\Phi_J(n,m;t).
\]

`fit-current` does **not** fit a universal `Gamma(q)` or `omega(q)` law; its
CSV preserves the pointwise ​`kz`, `m/Rcnt`, `q`, `Gamma`, `omega`, `a`, and `b`
relation for a later protocol-specific fit.  This avoids pooling C88/C99,
explicit/implicit, or cadence-distinct data into a synthetic dispersion law.

## Minimal sequence

```powershell
# 1. Establish data capabilities before choosing a command.
py scripts/collective_modes_cli.py audit `
  --case-id C88_L5_rep1 --dumps H:\path\water.dump `
  --fluid-types 3 --timestep-ps 0.0005 --wall-model explicit_fixed `
  --axis-source box_center --rcnt-A 4.07 --output H:\out\audit

# 2. Same declared profile; m=0 remains axial and m>0 adds circumference.
py scripts/collective_modes_cli.py current `
  --case-id C88_L5 --dumps H:\path\rep1.dump H:\path\rep2.dump H:\path\rep3.dump `
  --fluid-types 3 --timestep-ps 0.0005 --wall-model explicit_fixed `
  --axis-source box_center --rcnt-A 4.07 --n 1:20 --m 0:4 `
  --max-lag-ps 100 --output H:\out\current

# 3. Match the same n/m grid for Fs; then construct only from measured W.
py scripts/collective_modes_cli.py isf ... --n 1:20 --m 0:4 --max-lag-ps 100 --output H:\out\isf
py scripts/collective_modes_cli.py construct `
  --current-csv H:\out\current\current_per_replica.csv `
  --current-channel L --isf-csv H:\out\isf\isf_per_replica.csv `
  --weights-csv H:\out\static\weights_nm.csv --vacf-csv H:\out\vacf\vacf_per_replica.csv `
  --output H:\out\construct
```

The static weight file has exactly `n,m,weight` columns.  A missing mode or
lag is not imputed.  The `constructibility_vs_direct_vacf.csv` residual is the
primary diagnostic; good CJJ alone is not a demonstration of VACF closure.
