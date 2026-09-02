# Exposed collective-modes analysis assets

## Canonical entry

The sole entry point for new collective-dynamics analysis is:

```powershell
py scripts/collective_modes_cli.py <command> ...
```

Do not clone a historical `analyze_*`, `rebuild_*`, or `plot_*` script for a
new case before checking whether this entry point already covers it.

## Asset map

| Asset | Canonical path | Purpose |
|---|---|---|
| CLI entry point | `scripts/collective_modes_cli.py` | Stable public command for all mainline analyses. |
| Shared reader | `scripts/collective_modes/dump.py` | Streaming LAMMPS dump reading and bounded protocol inference. |
| Observable schema | `scripts/collective_modes/schema.py` | Input-capability contract and safe refusal of unsupported physics. |
| ISF implementation | `scripts/collective_modes/commands.py` (`isf`) | Total, self, and distinct cylindrical/axial ISF. |
| Current implementation | `scripts/collective_modes/commands.py` (`current`) | L/T/Tr current correlations, cross terms, and spectra. |
| Self transport | `scripts/collective_modes/commands.py` (`vacf`) | VACF, VACF-integrated MSD, and alpha(t). |
| Mode fit and closure | `scripts/collective_modes/commands.py` (`fit-current`, `construct`) | DHO parameters and `W*Fs*Phi_J` closure tables. |
| Output contract | `scripts/collective_modes/output.py` | CSV plus `metadata.json` emitted by every command. |
| Minimal plotting | `scripts/collective_modes/commands.py` (`plot`) | CSV-only figure template; it never recomputes physics. |
| API and examples | `scripts/collective_modes/README.md` | Equations, capabilities, and command examples. |
| Regression gate | `tests/test_collective_modes_smoke.py` | Synthetic end-to-end test. |
| Mainline policy | `governance/MAINLINE.md` | Definitions and permitted consolidation boundaries. |
| Historical audit | `governance/SCRIPT_ARCHITECTURE_AUDIT_2026-09-02.md` | Legacy script classification and duplicate evidence. |
| Source inventory | `governance/inventory/2026-09-02/` | Immutable initial source inventory. |

## Required execution order

```text
audit → isf/current/vacf → fit-current/construct → plot
```

Run `audit` before choosing an observable. Record the declared wall model,
fluid kind, atom types, radius, axis convention, time step, and input path in
each generated `metadata.json`. Use identical `(n,m)` and lag grids whenever
ISF, current, and VACF products are joined.

## Non-negotiable protocol boundaries

- An explicit/implicit CNT model is not inferred solely from a water-only dump.
- A single-site dump is not automatically oxygen-only water or argon.
- `m=0` is axial; `m>0` requires xyz. Cylindrical currents require xyz
  velocities.
- Self/distinct ISF, `CLT`/`CTL`, velocity reference frames, and sampling
  cadences remain separate columns/products.
- A flexible explicit CNT wall frame requires actual CNT positions and
  velocities; the code will not invent one.

## Local regression gate

```powershell
py -3 -m unittest tests/test_collective_modes_smoke.py -v
git diff --check
git status --short
```

## Execution boundary

This package performs post-processing and scientific analysis only. It does
not submit jobs, poll queues, open remote connections, copy files, discover
remote paths, or encode a cluster-specific Python environment. Run it in the
environment that is local to the data, and record the input manifest, command,
environment, output checksums, and completion log in the result package.

## Sharing boundary

This repository owns reusable source, tests, and small governance metadata.
Raw trajectories, LAMMPS restarts, HDF5/NumPy arrays, and bulk results remain
in host-local/archive storage and are excluded by `.gitignore`.  See
`governance/GIT_PUBLISHING_AND_CROSS_HOST.md` before creating a Git remote.
