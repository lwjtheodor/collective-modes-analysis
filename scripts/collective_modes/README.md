# Canonical collective-dynamics commands

The entry point is:

```powershell
py scripts/collective_modes_cli.py <command> ...
```

Every command writes readable CSV plus `metadata.json`.  Its metadata must be
copied into the authority result package. Completion still requires the
result-package README, source provenance, and QA gates in `AGENT.md`.

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

For a cylindrical analysis, the thin-shell mode-projection radius is computed
from the selected fluid itself:

\[
R_{\rm mode}=\langle r_{\rm selected}\rangle_{\rm all\ selected\ frames}.
\]

Thus an oxygen/water selection uses the ensemble O radial mean, and an argon
selection uses the ensemble Ar radial mean. CNT geometric radius and field
radius are never substituted for this quantity. For flexible explicit
CNT, the present command intentionally refuses a wall-relative calculation
until a CNT-frame extractor with actual CNT atoms/velocities is supplied.

## Input-capability contract

| Input fields | Commands enabled |
|---|---|
| `id,type,z[,iz],vz` | axial ISF, axial current, z-VACF, VACF-MSD-alpha |
| `id,type,x,y,z,vx,vy,vz` | all above plus cylindrical `Jr/Jtheta/L/Tinplane/Tr`, ordered cross kernels, r/theta VACF |
| `id,mol,type,x,y,z[,ix,iy,iz],vx,vy,vz` | additionally supports molecular/rotational plugins after their separate gate |
| CNT atom `x,y,z,vx,vy,vz` with declared CNT types | required later for flexible-CNT wall-relative frame |

## Replica and segment declaration

`--dumps a.dump b.dump` is retained only as shorthand for two independent
replicas. It must not be used for a continuation trajectory. Declare ordered
segments explicitly instead:

```powershell
--replica rep1=water_0_200ps.dump,water_200_1000ps.dump `
--replica rep2=water_rep2_0_1000ps.dump
```

or provide `--trajectory-manifest case.json`:

```json
{
  "replicas": [
    {"replica_id": "rep1", "segments": ["water_0_200ps.dump", "water_200_1000ps.dump"]},
    {"replica_id": "rep2", "segments": ["water_rep2_0_1000ps.dump"]}
  ]
}
```

The order in each `segments` array is physical time order. The reader rejects
field/selected-ID/box changes and nonmonotonic segment order; an identical
restart-boundary timestep is de-duplicated. The joined frames must then have
one uniform cadence, verified against `--timestep-ps`/`--dt-ps`.

## Commands and principal CSV products

| Command | Main products |
|---|---|
| `audit` | `dump_capabilities.csv` |
| `isf` | `isf_per_replica.csv`, `isf_ensemble_mean_sem.csv`: `F_total`, `F_self`, `F_distinct` for every `(n,m,lag)` |
| `current` | `current_per_replica.csv`, ensemble CJJ, ordered cross table, periodogram table; explicit `CJJ_extensive`, `CJJ_per_particle`, `CJJ_normalized`, `n_particles`; separate `Jz`, `Jr`, `Jtheta`, `L`, `Tinplane`, `Tr` |
| `vacf` | one native-cadence layer: per-rep/ensemble VACF and `msd_alpha_from_vacf_*`; lab/selected-COM/wall-relative is explicit |
| `vacf-stitch` | joins separately estimated cadence layers by declared physical-lag windows, then writes nonuniform-grid VACF-MSD-alpha products |
| `fit-current` | per-replica and ensemble-SEM mode-fit tables; `damped_carrier` is distinct from constrained `dho_physical` |
| `construct` | per-replica then ensemble `W*Fs*Phi_J`, optional matched direct-VACF residual tables; external static `W(n,m)` only, no amplitude fit |
| `plot` | a minimal data-readable PNG with zero line and columns named from the CSV |

The cylindrical Fourier phase is always

\[
\exp[-i(k_z z+m\theta)].
\]

The integer `m`, rather than `m/R`, appears in the angular phase.  The
thin-shell projection convention uses
\(k_\theta=m/R_{\rm mode}\) only for \(q\) and L/T projection. The current
metadata records that convention. Particles at the cylinder axis are rejected
for cylindrical observables because \(\theta,\mathbf e_r,\mathbf e_\theta\)
are undefined there.

`CJJ_extensive` is the all-origin fluctuating current ACF. Its explicitly
stored per-particle counterpart is

\[
C_{JJ}^{\rm per-particle}=C_{JJ}^{\rm extensive}/N,
\]

and `CJJ_normalized=CJJ_extensive/CJJ0_extensive`. Absolute CJJ, static
weights, and spectral intensity comparisons must name which normalization is
used.

The unconstrained retained carrier model is

\[
\Phi_J(k_z,m;t)=\exp[-\Gamma(k_z,m)t]
 [a(k_z,m)\cos\omega(k_z,m)t+b(k_z,m)\sin\omega(k_z,m)t],
\]

and constructibility is written without a fitted global amplitude as

\[
C_{vv}^{\rm construct}(t)=\sum_{n,m}W(n,m)F_s(n,m;t)\Phi_J(n,m;t).
\]

`fit-current --model damped_carrier` does **not** fit a universal `Gamma(q)` or `omega(q)` law; its
CSV preserves the pointwise `kz`, `m/R_mode`, `q`, `Gamma`, `omega`, `a`, and `b`
relation for a later protocol-specific fit.  This avoids pooling C88/C99,
explicit/implicit, or cadence-distinct data into a synthetic dispersion law.

`damped_carrier` is an intermediate-time descriptive fit and may carry a free
phase. It must not be called a strict DHO. For a normalized correlation that
includes zero time, use `--model dho_physical --fit-min-ps 0`; it constrains
\(C(0)=1\) and \(C'(0)=0\) through \(a=1,b=\Gamma/\omega\). Each replica is
fit independently; the authoritative uncertainty in the ensemble table is
replica SEM, not the nonlinear-fit covariance.

## Minimal sequence

```powershell
# 1. Establish data capabilities before choosing a command.
py scripts/collective_modes_cli.py audit `
  --case-id C88_L5_rep1 --dumps H:\path\water.dump `
  --fluid-types 3 --timestep-ps 0.0005 --wall-model explicit_fixed `
  --axis-source box_center --output H:\out\audit

# 2. Same declared profile; m=0 remains axial and m>0 adds circumference.
py scripts/collective_modes_cli.py current `
  --case-id C88_L5 --dumps H:\path\rep1.dump H:\path\rep2.dump H:\path\rep3.dump `
  --fluid-types 3 --timestep-ps 0.0005 --wall-model explicit_fixed `
  --axis-source box_center --n 1:20 --m 0:4 `
  --max-lag-ps 100 --output H:\out\current

# 3. Match the same n/m grid for Fs; then construct only from measured W.
py scripts/collective_modes_cli.py isf ... --n 1:20 --m 0:4 --max-lag-ps 100 --output H:\out\isf
py scripts/collective_modes_cli.py construct `
  --current-csv H:\out\current\current_per_replica.csv `
  --current-channel L --isf-csv H:\out\isf\isf_per_replica.csv `
  --weights-csv H:\out\static\weights_nm.csv --vacf-csv H:\out\vacf\vacf_per_replica.csv `
  --output H:\out\construct
```

The static weight file has exactly `n,m,weight` columns. A missing mode or lag
is not imputed. `construct` joins only identical
`(case_id,replica,n,m,lag_ps)` rows and rejects duplicate keys. It writes
`constructibility_sum_per_replica.csv` before calculating the replica mean/SEM
in `constructibility_sum_ensemble_mean_sem.csv`. The matched direct-VACF
residual is the primary diagnostic; good CJJ alone is not a demonstration of
VACF closure.

## Memory behaviour

The dump reader itself is iterator-based. `current` performs both its radius
pass and its current pass by streaming selected frames: it retains only compact
complex mode time series, not copies of atom tables. `isf` and `vacf` remain
correctness-first materialized paths at this revision because their all-origin
particle correlations require a full time history. A dedicated memmap
correlator is therefore required before routine multi-GB ISF/VACF production;
do not silently truncate or merge protocol-distinct files to work around RAM.

## Multirate VACF--MSD--alpha workflow

Do not place 1-fs, 10-fs, and 100-fs dumps into one `--replica`: a single
`vacf` call intentionally requires uniform cadence. Run `vacf` once per native
cadence layer, with identical fluid selection, component, and velocity-frame
definition. Then pass their `vacf_per_replica.csv` files to `vacf-stitch`:

```json
{
  "layers": [
    {"layer_id": "1fs", "csv": "H:/out/vacf_1fs/vacf_per_replica.csv", "lag_min_ps": 0.0, "lag_max_ps": 50.0},
    {"layer_id": "10fs", "csv": "H:/out/vacf_10fs/vacf_per_replica.csv", "lag_min_ps": 50.0, "lag_max_ps": 500.0, "include_lag_min": false},
    {"layer_id": "100fs", "csv": "H:/out/vacf_100fs/vacf_per_replica.csv", "lag_min_ps": 500.0, "lag_max_ps": 5000.0, "include_lag_min": false}
  ]
}
```

```powershell
py scripts/collective_modes_cli.py vacf-stitch `
  --layer-manifest H:\out\vacf_layers.json `
  --output H:\out\vacf_stitched
```

Windows are selected as `[min,max]` for a first layer and `(min,max]` when
`include_lag_min:false` is declared. They must cover lag zero, be strictly
nonoverlapping after selection, and contain the same `case_id`, replica,
component, velocity-frame convention, and stationary physical protocol. The
stitcher does not interpolate or average overlaps. It records `layer_id` and
source CSV for every output lag, and integrates the joined VACF with the
nonuniform-grid trapezoid rule before calculating MSD and alpha.
