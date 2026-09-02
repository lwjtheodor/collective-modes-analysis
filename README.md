# ISF and collective modes

## Canonical analysis entry and Git sharing boundary

New work must enter through the unified collective-modes package rather than a
copied historical analysis script.  Start from
[`EXPOSED_COLLECTIVE_MODES_ASSETS.md`](EXPOSED_COLLECTIVE_MODES_ASSETS.md), then
use [`scripts/collective_modes/README.md`](scripts/collective_modes/README.md)
for the command contract.  The maintained source-and-provenance Git workflow,
including the private-first cross-host sharing policy, is in
[`governance/GIT_PUBLISHING_AND_CROSS_HOST.md`](governance/GIT_PUBLISHING_AND_CROSS_HOST.md).

**Question:** Is the finite-size response controlled by physical wave number rather than a box-label artefact?

**Asset entry point:** [`assets.md`](assets.md) is the project-level catalog for
all citable analysis assets, their authoritative reproducible archives, and
their protocol/qualification boundaries.

Required decomposition: total `F`, self `Fs`, and distinct `Fd`; compare the same physical `k` across lengths before interpreting a new `k_min` mode.  Keep KWW first-branch fits separate from oscillatory/current-mode metrics and from damping estimates.

## Existing findings

- Equal physical wave numbers represented by different mode indices in 2L/3L were found to have closely matching density/current behaviour, supporting a `k`-controlled rather than box-label-controlled description.
- Longer boxes admit lower `k_min` modes with slower distinct/current dynamics.  This is the core evidence compatible with an infrared cutoff.
- `F_s` and `F_d` must remain separate: a stretched or oscillatory total ISF does not by itself identify tagged-particle non-Gaussianity or collective compensation.

## Open questions

- Extend matched-physical-k `F/F_s/F_d` consistency from the early 2L/3L demonstration through all reliable 1L–5L representations.
- Determine KWW exponent and relaxation-time dependence from molecular spacing to long wavelengths with block uncertainty.
- Resolve whether apparent ISF/current oscillations are propagative peaks, finite-window spectral leakage, or damped non-propagative relaxation.

## Cross-chirality peculiar VACF and strict Cvv-ODE alpha (staged 2026-08-14)

The `(7,7)`, `(9,9)`, and `(17,0)` weak-NH/no-global-momentum-removal matrix
is being analysed with the same oxygen-level framework used for the completed
`(8,8)` calculation: instantaneous oxygen-COM subtraction, all-origin axial
peculiar VACF, and the strict reconstruction
\(\alpha_z(t)=tI(t)/J(t)\), where \(I=\int_0^t C_{vv}\) and
\(J=\int_0^t I\).  Inputs are 1L--5L (two replicas, 10 fs / 1 ns) and 10L
(three replicas, 100 fs / 10 ns), with analysis deliberately bounded to
0--200 ps.  The staged source and protocol README are in
`stage_crosschirality_vacf_alpha_2L10_weaknh_20260814/`; completed compact
CSV/JSON/figure assets will be archived under
`results/collective_mode_response/crosschirality_vacf_cvv_alpha_ode/2026-08-14/`.
If the current assets fail the SEM/shape gate, the prepared economical
resubmission is 1L/2L/4L/8L x four independent configuration starts x three
chiralities, with 2.5-ns production and oxygen `id,type,z,vz` every 100 fs.
The explicit conditional protocol and acceptance gate are in
`stage_crosschirality_vacf_alpha_2L10_weaknh_20260814/resubmission_plan_1L2L4L8L.md`.
For `(8,8)`, retain the completed 2L/4L middle points and add 1L and 8L with
eight starts each; 16L is an adaptive sentinel rather than a first-round case.
The dedicated resubmission uses a dual-rate oxygen-only output: 100 fs over
the full 2.5-ns production for 400-ps tail statistics, plus a concurrent
1-ns/10-fs stationary burst for high-resolution VACF and axial velocity DOS.

## Canonical scripts

The only supported scientific-analysis entry point is
`scripts/collective_modes_cli.py`, backed by the reusable package in
`scripts/collective_modes/`.  It covers dump auditing, axial/cylindrical
ISF, L/T/Tr current modes, VACF--MSD--alpha, mode fitting, constructibility,
and CSV-first plotting across case, chirality, box length, and declared
protocol parameters.

The former case-specific `analyze_*`, `aggregate_*`, `compare_*`,
`rebuild_*`, and `plot_*` scripts are intentionally not part of this
repository's current tree. They are historical provenance only; see
`governance/LEGACY_SCRIPT_ARCHIVE.md` and Git commit `9e06eae` if a past result
requires forensic reconstruction.

## Definitions and why this branch is the infrared test

With \(\rho_k(t)=\sum_j\exp[ikz_j(t)]\), the standard axial ISF decomposition is

\[
F(k,t)=N^{-1}\langle\rho_k(t)\rho_{-k}(0)\rangle,
\quad F_s(k,t)=N^{-1}\sum_i\langle e^{ik\Delta z_i(t)}\rangle,
\quad F_d=F-F_s.
\]

`F_s` asks whether a tagged molecule has decorrelated; `F_d` records the
other-molecule density rearrangement needed to sustain the collective pattern.
The longitudinal current is
\[
J_k(t)=\sum_j v_{z,j}(t)e^{ikz_j(t)},\qquad
C_J(k,t)=N^{-1}\langle J_k(t)J_{-k}(0)\rangle.
\]
Continuity, \(\dot\rho_k=-ikJ_k\), explains why a smooth density correlator
can coexist with a sign-changing current correlator.  Compare equal physical
\(k\), e.g. \(n=1\) in 1L with \(n=2\) in 2L, before discussing a new lowest
mode.  `C_J(k,0)` is retained in absolute units when comparing modal weight;
normalization is used only to compare shapes.

## Materialized evidence assets

- `assets/library/isf_current/matched_k_isf_total_self_distinct_1L5L.png`: total/self/distinct
  fixed-\(k\) comparison.
- `assets/library/isf_current/matched_k_self_isf_1L5L.png`: tagged relaxation alone.
- `assets/library/isf_current/current_mode_k_dependence_5L.png` and
  `assets/library/isf_current/absolute_current_rebound_time_vs_wavelength_5L.png`: absolute
  current-mode weights and operational first-rebound timings.
- The `5L_n1n10_*.png/.pdf` files are the wide wavevector scan.  Their rebound
  time is an operational shape measure, not a validated damping rate \(\Gamma\).

## Updated `(8,8)` 10L all-mode current-lobe result (2026-08-13)

For the 10L weak-NH/no-global-momentum-removal trajectory (100 fs / 10 ns;
three replicas; instantaneous water-COM axial-velocity subtraction), the first
complete negative lobe of the **C(0)-normalized** axial-current ACF was
extracted independently for modes \(n=1\ldots10\).  The mode wavelength is
\(\lambda=2\pi/k=L_z/n\), and the lobe is bounded by its first down-crossing
and first subsequent up-crossing.

The raw normalized lobe area is exceptionally well represented, over the
observed \(\lambda=10.084\)--100.840 nm range, by the empirical linear fit

\[
A_-(\lambda)=-0.851+0.143\lambda\ \mathrm{ps},\qquad R^2=0.99969.
\]

The nonzero negative intercept is **not** a physical negative-area prediction.
It marks finite-wavelength crossover: the fitted line must not be extrapolated
below the observed range.  In particular, direct division by wavelength gives
\(A_-/\lambda=0.143-0.851/\lambda\), so this ratio is not expected to be
constant even though the raw area is nearly linear.  At the two ends of the
sampled set, \(A_-/\lambda\) is 0.06748 (\(n=10\)) and 0.13441 ps nm\(^{-1}\)
(\(n=1\)).

**Stable panel-a source-data archive:**
`assets/library/cjj/cjj_88_10L_raw_source_area_lambda_linearity_audit.csv`.
It is derived directly from all 30 `CJJ_alln.json` mode-summary values (10
modes x 3 replicas), and records raw area, replica SEM, wavelength, linear
prediction, residual and \(A_-/\lambda\).  The associated static figure,
vector files and reproducible script are:

- `assets/library/cjj/cjj_88_10L_raw_source_area_lambda_linearity_audit_nature.png`
- `assets/library/cjj/cjj_88_10L_raw_source_area_lambda_linearity_audit_nature.pdf`
- `assets/library/cjj/cjj_88_10L_raw_source_area_lambda_linearity_audit_nature.svg`
- the historical generator recorded in `governance/LEGACY_SCRIPT_ARCHIVE.md`

The 5L-versus-10L matched-physical-\(k\) all-mode check remains available in
`assets/library/cjj/cjj_88_5L_10L_matched_k_first_lobe_ratios.csv`.  It should be used
before attributing the finite-wavelength correction to the box label itself.

## Updated `(8,8)` 2L--10L collapse (2026-08-10)

The latest completed dump assets were reanalysed using the validated de-COM,
all-origin MSD and the historical trailing-one-decade estimator.  In the
5--200 ps analysis interval,
\[
\alpha_z(t)=\operatorname{OLS}_{[t/10,t]}[\log M_z\ \mathrm{vs.}\ \log t].
\]
For 2L, 3L, 4L, 5L, and 10L, the mean valley locations are 24.56, 38.50,
48.20, 55.99, and 109.90 ps, or 12.28, 12.83, 12.05, 11.20, and 10.99 ps per
relative length.  Thus the timing remains close to linear through 10L.

The raw \(\alpha_z\) curves already align well against \(t/L_{\rm rel}\).
After only amplitude normalization of \(1-\alpha_z\), a pairwise shape scan
gives the descriptive optimum \(\beta_{\rm shape}=0.93\), with a
replica-bootstrap 68% interval 0.91--0.96.  This supports a near-ballistic
finite-size clock, but is not a thermodynamic-limit transport exponent: the
MSD amplitude and the final \(L\to\infty\) fate remain separate questions.

Use `assets/library/msd_alpha_88/8_8_MSD_alpha_collapse_2L10L_20260810.png` and
`assets/library/msd_alpha_88/8_8_alpha_shape_beta_scan_2L10L_20260810.png`; numerical tables and
the exact plotting script are in
`results/collective_mode_response/finite_size_collapse_2L10L/2026-08-10/`.
The key protocol limitation is explicit: 2L--5L use 10-fs/1-ns windows,
whereas 10L uses a 100-fs/10-ns window; all comparisons are restricted to
5--200 ps.
# Asset reference (2026-08-06)

| Asset | Coverage / cadence | Use boundary |
|---|---|---|
| `TRAJ-88-HF-2L5L` | `(8,8)`, 2L--5L, O 10/100 fs, case-specific duration | total/self/distinct ISF and current modes; audit duration before fitting attenuation |
| 5L molecular-current pilot | `(8,8)` 5L, 3 rep, first 320 ps at 10 fs | PSD resolution `0.003125 ps^-1`; use only with the documented window and replicate SEM |
| Cross-chirality HF comparison | `(7,7)/(15,0)`, 2L, O 10 fs for 300 ps | comparison of resolved short-time modes, not a length-series estimator |
| Weak-NH control | `(7,7)/(9,9)/(17,0)`, 1L--5L x2, O 10 fs/water 1 ps | protocol sensitivity only; see master index for provenance and last-audit status |

Paths and status definitions: [`dump_asset_inventory_20260806.md`](../shared/metadata/dump_asset_inventory_20260806.md).
