# Collective-dynamics governance

This directory separates scientific provenance from day-to-day staging.  It is
the control layer for `02_isf_collective_modes`, not a replacement for the
authoritative result packages already required by `AGENT.md`.

## Repository boundary

Git tracks: project contracts, reusable-source candidates and subsequently
canonical source under `scripts/`, small templates, registry records, and reviewed inventory snapshots.  It does
not track trajectories, restart files, scheduler output, fetched raw data, or
large result trees.  Those stay at their existing local/CCFEP locations and
are referenced by exact path plus size/hash/frame audit in result manifests.

An analysis result can therefore be scientifically citable only when its
archive package and `assets.md` gates are met; a Git commit proves source
history, not run completion.

## Four data/code states

| State | Location/pattern | Git role | Scientific role |
|---|---|---|---|
| Candidate / canonical reusable code | `scripts/` then `scripts/<domain>/` | tracked | all baseline files are candidates until reviewed; only a promoted domain file may be extended for new work |
| Execution snapshot | `stage_*`, `remote_fetch/`, `heartbeat_fetch/` | ignored | immutable record of a submitted or recovered run; never silently treated as canonical |
| Authoritative analysis package | `results/collective_mode_response/<topic>/<date>/` | ignored | result-level README, compact tables, QA, figures and dump provenance |
| Display entry | `assets/library/<topic>/` | ignored | curated copies/links pointing back to the authoritative package |

## Required procedure for consolidating a script

1. Use the inventory to find all same-name and exact-content variants.
2. Compare input fields, COM/peculiar-velocity definition, cadence/window,
   wave-vector convention, replica semantics, estimator and output schema.
3. Select or write a reviewed canonical implementation in `scripts/`; give it
   a small protocol header and a deterministic CLI.
4. Add a `SCRIPT_PROVENANCE.yaml` in each new result package with the canonical
   relative path, Git commit, input asset IDs/absolute paths and output schema.
5. Mark the older snapshot as `superseded`, `historical`, or `exploratory` in
   its own README/manifest.  Never delete it during consolidation.

## Inventory cadence

Run the inventory before any consolidation wave and after each reviewed wave:

```powershell
py scripts/maintenance/inventory_collective_modes.py --root . --output governance/inventory/YYYY-MM-DD
```

It does not read trajectory contents or hash large files.  It records all file
paths/sizes, hashes source-sized files, detects exact script copies and
same-basename variants, and extracts literal CCFEP source paths.  Large dump
integrity remains the responsibility of the field-level dump registry defined
in `AGENT.md`.
