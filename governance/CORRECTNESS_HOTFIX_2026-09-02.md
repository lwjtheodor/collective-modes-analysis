# Correctness hotfix: scientific-analysis mainline

This record addresses the review of the initial canonical implementation. No
formal dataset was analysed with the pre-hotfix code after this correction.

## Corrected before formal use

1. Cylindrical current modes now use the single-valued phase
   \(\exp[-i(k_z z+m\theta)]\). `m/R_mode` is limited to q and L/T projection.
2. Current CSVs distinguish extensive, per-particle, and normalized CJJ and
   record `n_particles`; no prior `CJJ_raw` meaning is silently retained.
3. Constructibility joins and sums within each `(case_id, replica)` before
   calculating an ensemble mean/SEM. Duplicate keys are fatal errors.
4. Mode fitting is per independent replica. A free-phase intermediate-time
   `damped_carrier` is explicitly distinct from constrained `dho_physical`.
5. `m=0` ISF/current/z-VACF have a real `id,type,z,vz` path and no longer
   require transverse fields.
6. Cylindrical observables reject axis-singular selected particles instead of
   silently assigning an arbitrary basis.

## Tests added

- projection reductions at `m=0` and `n=0`;
- Fourier phase independence from `R_mode`;
- axial-only dump execution;
- two-replica construct isolation and duplicate-key rejection.

## Still deliberately not claimed as complete

- a case manifest that distinguishes independent replicas from restart
  segments;
- memory-bounded streaming/memmap implementations for very large trajectories;
- an exact radial-cylinder continuity current rather than the declared
  thin-shell `m/R_mode` projection;
- a systematic lab-versus-selected-COM low-k sensitivity analysis;
- signed-frequency / complex quadrature spectral estimators.

These remain separate implementation work, not properties to infer from the
present CSV interface.
