# Script canonicalization status — baseline 2026-09-02

## Non-negotiable interpretation of the first Git baseline

The existing 107 files below `scripts/` are now versioned as **source
candidates**, not as endorsed canonical implementations.  This distinction is
deliberate: the project contains many physically non-interchangeable
protocols, and the present inventory establishes that neither a filename nor a
Git path establishes scientific equivalence.

## Promotion gate

A candidate may be promoted to `scripts/<domain>/` only after a short review
record specifies all of the following.

1. Input schema and LAMMPS fields; O/water/atom selection.
2. Position/velocity reference frame and COM subtraction.
3. Observable formula and normalization (`C(0)` handling included).
4. `Lz`, discrete indices and physical `k` / cylindrical `m/Rcnt` convention.
5. Sampling cadence, total usable window, FFT/correlation/windowing details.
6. Replica provenance (velocity seed versus independent configuration) and
   uncertainty estimator.
7. A compact archived fixture plus expected numerical/shape QA.
8. Exact result-package and CCFEP input paths; a statement of what claims the
   implementation is *not* valid for.

## First review order

| Order | Domain target | Baseline candidates to inspect first | Reason |
|---|---|---|---|
| 1 | `scripts/longitudinal/` | `rebuild_axial_isf.py`, `rebuild_collective_isf_kseries.py`, `analyze_88_5L_full_dispersion.py`, `analyze_matched_k_5L_10L_protocolmatched.py` | defines the central matched-physical-`k` current/ISF chain |
| 2 | `scripts/cylindrical/` | `analyze_full_kz_static_longitudinal.py`, `analyze_kz_transverse_exploratory.py`, `analyze_implicitC88_transverse_SJJ.py` | must preserve separate axial/radial/circumferential and helical definitions |
| 3 | `scripts/self_transport/` | `analyze_vacf_tail.py`, `analyze_vacf_integral_lockin.py`, `analyze_windowed_alpha_from_allorigin_vacf.py` | connects velocity definition to MSD/alpha and lobe evidence |
| 4 | `scripts/qa/` | `audit_dump_assets.py`, `audit_fullwater_orbital_Lz.py`, `compact_tabular_assets.py` | provides the field/cadence/provenance gates required before any reconstruction |

The remaining candidates stay versioned and searchable, but should be marked
`historical`, `exploratory`, or `superseded` only after their input/output
contracts have been compared.  Do not infer status solely from modification
date or from duplicate bytes in an archive copy.
