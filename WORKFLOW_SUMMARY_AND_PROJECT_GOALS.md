# Collective-dynamics project: workflow summary and goals

**Project root:** `H:/gcmc_explore/translational_anomaly/02_isf_collective_modes`  
**Canonical private repository:** <https://github.com/lwjtheodor/collective-modes-analysis>  
**Status date:** 2026-09-07

## 1. Project objective

This project establishes an auditable, reusable analysis system for confined-fluid
collective dynamics in CNT-like cylinders.  The system must support water,
oxygen-only water representations, and argon; implicit CNT, explicit fixed CNT,
and eventually explicit flexible CNT; distinct box lengths, chiralities, and
sampling cadences.

The target is not a collection of case-specific plotting scripts.  It is a
traceable chain from a declared simulation protocol and raw dump to compact,
reusable scientific products with explicit definitions, provenance, QA, and
physical limits.

The principal observables are:

\[
F(k,m,t),\quad F_s(k,m,t),\quad F_d(k,m,t),
\]

\[
C_{JJ}(k,m,t),\quad C_{AB}(k,m,t),\quad S_{JJ}(k,m,\omega),
\]

\[
C_{vv}(t),\quad \mathrm{MSD}_{\rm VACF}(t),\quad \alpha(t),
\]

and the explicitly testable constructibility expression

\[
C_{vv}^{\rm construct}(t)=\sum_{n,m}W(n,m)F_s(n,m,t)\Phi_J(n,m,t).
\]

## 2. End-to-end dump lifecycle

The complete research workflow has two deliberately separated layers.

```text
simulation planning and execution layer
    -> dump production and manifest writing
    -> completion / raw-integrity / protocol gates
    -> compact-result recovery and asset registration
    -> canonical scientific-analysis layer
```

### 2.1 Simulation planning and task preparation

For every new production or continuation task, record before submission:

- case identity: chirality, geometry, `Lz`, temperature, pressure/density,
  molecule/atom count, and fluid representation;
- wall protocol: implicit, explicit fixed, or explicit flexible CNT, including
  force-field / field parameters;
- dynamical protocol: integrator timestep, thermostat/barostat/NVE state,
  relaxation time, momentum/COM treatment, and intended analysis velocity frame;
- replica provenance: independent configurational starts versus velocity seeds
  from a common parent configuration;
- dump plan: selected atom types, output fields, cadence layers, intended
  duration, dump boundaries, and restart/continuation order;
- expected scientific product and the required source cadence/window.

The simulation-side task should write a machine-readable
`trajectory_manifest.json` next to its dumps.  A canonical manifest must carry:

```text
case_id, physical protocol, fluid types, CNT/wall model,
integration timestep, dump interval, atom fields, box geometry,
replica ID/provenance, ordered segments, cadence layers, lag windows,
source paths and optional hashes.
```

Simulation submission, SSH, PBS, `jsub`, remote copying, and scheduler
monitoring intentionally live outside this analysis repository.  This project
consumes their declared outputs; it does not reintroduce cluster-control scripts
under `scripts/`.

### 2.2 Submission, completion, and raw-dump qualification

Submission acceptance, a visible directory, or a nonzero file size is never a
scientific completion condition.  A dump becomes **analysis-ready** only after:

1. normal simulation termination is demonstrated by the final log;
2. expected dump files exist at their declared paths;
3. frame count, first/last timestep, cadence, atom fields, particle IDs, and
   box protocol are audited;
4. any continuation boundary is checked for duplicate-timestep content identity;
5. protocol and replica provenance are registered in the manifest;
6. required compact logs/metadata are recovered locally and listed in the asset
   registry.

For a duplicate timestep at a segment hand-off, canonical analysis allows
de-duplication only when dump fields, selected IDs, box bounds, and all selected
particle records agree.  Same timestep with different coordinates, velocities,
or other fields is a fatal input error.

### 2.3 Dump recovery and local asset organization

Do not routinely copy multi-GB or multi-10-GB raw dumps to the local archive.
Keep raw sources at their declared CCFEP/Wisteria absolute paths and recover:

- compact CSV/NPZ/HDF5 scientific products;
- `metadata.json`, manifest, source hash/size/frame audit, and normal job log;
- figures only after their underlying compact tables are preserved;
- a short README describing protocol, definitions, and known limits.

Organize durable results by physical case/protocol and observable, not by an
execution timestamp.  Timestamps belong in run receipts and archival subfolders,
not in the primary scientific identity.

## 3. Canonical post-processing workflow

The current reusable entry point is:

```powershell
py scripts/collective_modes_cli.py <command> ...
```

### 3.1 Commands currently available

| Command | Function | Primary output |
|---|---|---|
| `audit` | dump fields, types, capability and protocol checks | `dump_capabilities.csv` |
| `radial-qa` | selected O/Ar radial distribution and thin-shell QA | `radial_mode_thin_shell_qa.csv` |
| `isf` | total/self/distinct axial or cylindrical ISF | `isf_*csv` |
| `current` | CJJ, L/T projection, ordered cross kernels, periodogram | `current_*csv` |
| `vacf` | one native-cadence peculiar VACF and VACF-integrated MSD/alpha | `vacf_*csv`, `msd_alpha_*csv` |
| `vacf-stitch` | cadence-layer VACF stitching on a nonuniform lag grid | `vacf_stitched_*csv` |
| `construct` | static-weighted `W Fs Phi_J` closure | `constructibility_*csv` |
| `fit-current` | per-replica current-kernel carrier/DHO fits | `current_mode_fit_*csv` |
| `plot` | minimal readable CSV plot | PNG |

### 3.2 Reader and input contract

The unified reader supports axial minimal dumps (`id,type,z,vz`) and full
three-dimensional dumps (`id,type,x,y,z,vx,vy,vz`), with image flags where
available.  Axial and cylindrical capability are never inferred merely from a
filename.

Ordered continuation files are declared either as repeated `--replica`
arguments or a JSON trajectory manifest.  Each replica must have a uniform
cadence; segment fields, selected IDs, box bounds, and timestep ordering must
remain compatible.  Different physical protocols and different native cadences
are not silently merged into one raw time series.

### 3.3 Cylindrical geometry and radial QA

The cylindrical Fourier phase is

\[
\exp[-i(k_z z+m\theta)].
\]

The angular phase always uses integer `m`, not `m/R`.  The fluid-selected mode
radius is

\[
R_{\rm mode}=\langle r_{\rm selected}\rangle,
\]

where selected O atoms define a water/oxygen mode and selected Ar atoms define
an argon mode.  CNT geometric radius or field radius must not replace this
quantity.

`m/R_mode` is used only for

\[
k_\theta=m/R_{\rm mode},\qquad q=\sqrt{k_z^2+k_\theta^2},
\]

and L/T projection.  `radial-qa` records radial mean, standard deviation,
quantiles, central-90-percent width, and `sigma_r/R_mode`.  Its thin-shell
pass/fail is an operational projection QA, not proof of a complete radial
eigenmode treatment.

### 3.4 CJJ, ISF, and cylindrical modes

For current modes, retain separately:

\[
J_z,\quad J_r,\quad J_\theta,\quad J_L,\quad J_{T,\theta},\quad J_{T,r}=J_r.
\]

The in-plane projections are

\[
J_L=\frac{k_zJ_z+k_\theta J_\theta}{q},\qquad
J_{T,\theta}=\frac{-k_\theta J_z+k_zJ_\theta}{q}.
\]

`Jr` and `Jtheta` are not averaged into a generic transverse channel.  Ordered
cross kernels `C_AB` and `C_BA` are retained separately.

CJJ output explicitly distinguishes extensive, per-particle, and normalized
forms:

\[
C_{JJ}^{\rm per-particle}=C_{JJ}^{\rm extensive}/N,
\qquad
C_{JJ}^{\rm normalized}=C_{JJ}^{\rm extensive}/C_{JJ}^{\rm extensive}(0).
\]

ISF output retains `F_total`, `F_self`, and `F_distinct=F_total-F_self` for
every declared `(n,m,lag)`.

### 3.5 VACF, MSD, alpha, and multirate stitching

Each `vacf` invocation handles exactly one native uniform cadence.  It must not
be given 1-fs, 10-fs, and 100-fs dumps as one raw replica.

For a multirate trajectory:

```text
native 1-fs VACF  -> short-time layer
native 10-fs VACF -> intermediate layer
native 100-fs VACF -> long-time layer
                         ↓
                declared physical-lag windows
                         ↓
              `vacf-stitch` nonuniform time grid
                         ↓
               VACF-integrated MSD and alpha(t)
```

The stitcher does not interpolate or silently average overlap regions.  It
requires a lag-zero point, non-overlapping selected lag points, and a common
component/physical protocol.  Integration uses the nonuniform-grid trapezoid
rule.

### 3.6 Constructibility and fits

The principal closure test is

\[
C_{vv}^{\rm construct}(t)=\sum_{n,m}W(n,m)F_s(n,m,t)\Phi_J(n,m,t).
\]

Static `W(n,m)` is external measured input, not a fitted global amplitude.
Current, ISF, and direct VACF join strictly within matching case/replica/mode/
lag keys.  Ensemble mean and SEM are calculated only after the per-replica
calculation.

The retained current carrier is

\[
\Phi_J=e^{-\Gamma t}[a\cos(\omega t)+b\sin(\omega t)].
\]

Fits are per replica; replica SEM, not a single nonlinear-fit covariance, is
the default uncertainty summary.

## 4. Current asset map

### 4.1 Canonical code and governance

```text
scripts/collective_modes/
tests/test_collective_modes_smoke.py
governance/MAINLINE.md
governance/CORRECTNESS_HOTFIX_2026-09-02.md
governance/OPERATIONS_SEPARATION.md
governance/inventory/2026-09-02/
```

The baseline inventory recorded 11,521 files, 1,579 scripts, 210 exact-script
duplicate clusters, 172 same-basename clusters, and 110 distinct CCFEP paths.
Those counts describe the historical asset estate, not the endorsed canonical
script set.

### 4.2 Local raw and compact assets

```text
assets/library/
assets/raw_oxygen_zvz_100fs_L2L10_weakNH_2026-08-14/
results/collective_mode_response/
remote_fetch/
```

The audited local estate contained approximately 7.87 GB under `assets/`,
58.02 GB under `results/`, and 2.05 GB under `remote_fetch/`.  Important result
families include C88/C99 matched-k CJJ, C88 transverse/helical/self-distinct,
multirate VACF/alpha, static vertex, and CJJ--VACF closure packages.

### 4.3 C88 N1600 authority and raw input

Local compact authority package:

```text
remote_fetch/implicit_C88_static_vertex_VACF_20260831/N1600_weakNH_6ns/
```

Important compact products:

```text
kernel_curves_per_replica_and_ensemble.csv
reconstruction_metrics_per_replica.csv
static_weights_per_replica.csv
metadata.json
analysis_1385918.ccpbs1.log
SUCCESS.txt
```

Verified remote raw inputs:

```text
/lustre/home/users/ewu/vb_gcmc/MD/N1600_weakNH_6ns/rep{1,2,3,4}/production_100fs.dump
```

For `rep1`, the dump was confirmed present on 2026-09-07 with size
20,430,588,691 bytes.  The authority protocol is implicit `(8,8)`, 350 K,
weak-NH, 100-fs dump over 6 ns, 60,001 frames, 1,600 selected oxygen atoms,
and four velocity-seed replicas.

Verified multirate C88 N1600 sources:

```text
/lustre/home/users/ewu/vb_gcmc/MD/stage_implicit_C88_N1600_LT_dualcadence_20260831/rep{1,2,3,4}/fine_1fs_first100ps.dump
/lustre/home/users/ewu/vb_gcmc/MD/stage_implicit_C88_N1600_LT_dualcadence_20260831/rep{1,2,3,4}/medium_10fs_full1ns.dump
/lustre/home/users/ewu/vb_gcmc/MD/N1600_weakNH_6ns/rep{1,2,3,4}/production_100fs.dump
```

## 5. C88 N1600 canonical-regression status

Canonical regression staging is present remotely:

```text
/lustre/home/users/ewu/vb_gcmc/MD/c88_n1600_canonical_regression_20260902/
```

It contains the canonical package, CLI wrapper, and `run_canonical_rep.pbs`.
The planned `rep1` task runs dump audit, radial QA, axial `m=0` CJJ for
`n=1:10`, and peculiar z-VACF through 250 ps.

As of 2026-09-07 it contains no `output/` directory and no `SUCCESS.txt`.
The earlier `jsub` attempt was rejected with `Reason: Disk`; no C88 canonical
regression job ID has been accepted.  Thus the code and input are ready, but
no raw-dump numerical regression result exists yet.

## 6. Open priorities

1. Resolve the CCFEP Disk submission gate, run C88 N1600 `rep1`, and validate
   normal PBS termination plus compact canonical outputs before expanding to
   `rep2`--`rep4`.
2. Compare canonical and historical C88 N1600 normalized CJJ/VACF numerically,
   with explicit velocity-frame and normalization columns.
3. Produce real-data multirate overlap and hand-off continuity QA for the
   1-fs/10-fs/100-fs C88 layers.
4. Implement memmap/blockwise particle-history correlation for very large ISF
   and VACF inputs; current-mode analysis is already streaming, ISF/VACF are
   not yet fully OOM-safe.
5. Introduce a versioned simulation-side `trajectory_manifest.json` generator
   and make the analysis CLI consume its full protocol/cadence/replica contract.

## 7. Evidence and interpretation rules

- A job submission, a queue row, directory existence, or a finite dump size is
  not a completed scientific result.
- A formal result requires source protocol, completed log, input QA, compact
  numerical products, metadata, and stated limits.
- Compare boxes at equal physical `k`, not merely equal integer mode index.
- Keep cadence layers, wall/thermostat protocols, velocity frames, and replica
  provenance visible in all comparisons.
- Treat velocity-seed spread from a shared parent as conditional variation, not
  automatically as independent-configurational uncertainty.
- Do not infer radial eigenmode completeness from a thin-shell pass, or infer
  high-frequency linewidths from low-cadence spectra.
