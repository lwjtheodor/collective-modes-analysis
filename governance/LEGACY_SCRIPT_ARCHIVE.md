# Legacy script archive: excluded from the scientific mainline

## Decision

On 2026-09-02, the current Git tree was deliberately reduced from 110 tracked
files under `scripts/` to 10 reusable mainline/maintenance files. The 100
removed files were named for particular chiralities, box lengths, cadence,
replica sets, fitting episodes, result archives, or presentation figures.
They are not generic interfaces and must not remain adjacent to the canonical
analysis entry point.

The removed sources were physically preserved, with their original relative
paths, in the local history archive:

```text
H:/gcmc_explore/translational_anomaly/collective_modes_legacy_analysis_archive/2026-09-02/
```

They also remain inspectable in the immutable initial repository snapshot
`9e06eae` (`chore: establish collective-modes analysis baseline`). Neither
location makes them supported current analysis code.

## Current allowed `scripts/` tree

```text
scripts/collective_modes_cli.py
scripts/collective_modes/
scripts/maintenance/inventory_collective_modes.py
```

The mainline package must express case/chirality/box-length/protocol choices as
declared command arguments and metadata, never by generating a new source file
whose name encodes a particular case.

## Promotion rule

Before a historical algorithm can return to the current tree, it must satisfy
all of the following:

1. Its physical definition is not already supplied by `audit`, `isf`,
   `current`, `vacf`, `fit-current`, `construct`, or `plot`.
2. The implementation has a case-independent API and does not embed a
   chirality, box length, file path, job identifier, or figure destination.
3. It accepts explicit input metadata and emits documented CSV/metadata.
4. It has a synthetic test and a real-dump regression against a named legacy
   reference product.
5. The addition is reviewed as a focused Git change, rather than copied from a
   historical script.
