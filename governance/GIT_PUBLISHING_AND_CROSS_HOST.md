# Git publishing and cross-host sharing policy

## Scope

This is a **source-and-provenance repository**. It may be cloned on the
workstation, trusted data-host environments, and collaborator machines to run
the same analysis code. It does not version raw dumps or large generated
results, and it contains no job-submission or remote-control implementation.

`.gitignore` excludes trajectories, restarts, HDF5/NumPy arrays,
assets/results/stage directories, runtime logs, and Codex state. This is a
boundary, not a substitute for reviewing the staged diff before every push.

## Visibility decision

Start with a **private** remote repository. Do not make it public until every
tracked file has been reviewed for unpublished results, manuscript-sensitive
claims, cluster/account information, access tokens, credentials, third-party
redistribution restrictions, and protected data.

Never commit credentials. Use an SSH agent/key or the Git host credential
manager; never put a token in a remote URL.

## First publication workflow

First set a deliberate author identity; this repository currently has none.

```powershell
git config --local user.name "<your chosen author name>"
git config --local user.email "<your chosen author email>"
git diff --cached --check
git status --short
git commit -m "chore: establish collective-modes analysis baseline"
```

After creating an empty private repository on the chosen Git host:

```powershell
git branch -M main
git remote add origin <private-repository-url>
git push -u origin main
```

## Cross-host operating model

```text
private Git remote
      ├── workstation: edit, test, review, commit, push
      ├── data-host environment A: pull a reviewed commit and execute near data
      └── data-host environment B: pull the same commit and execute near data
```

Each host retains its large data in host-local campaign/archive storage. Every
run writes a compact manifest or result README containing the input/output
absolute paths, commit SHA, command line, Python environment, hashes, and QA
status. Commit the compact provenance after review, never the bulk payload.

Recommended sequence:

1. Workstation: implement, test, and commit a bounded change.
2. Push the reviewed commit to `main` or `feature/<topic>`.
3. On a data-host environment: `git pull --ff-only`, record
   `git rev-parse HEAD`, and run against its local data.
4. Archive outputs where they belong; transfer only compact reviewed products
   when appropriate.
5. Add/update the compact provenance manifest and push it after review.

## Version rules

- `main`: tested, reviewable analysis mainline.
- `feature/<topic>`: bounded implementation or method change.
- `archive/<date>-<topic>`: immutable historical source reference only, never
  a substitute for data storage.
- Tag `v0.1.0` only after real-dump regression against selected legacy cases,
  not merely after a synthetic smoke test.
- Every result package records the commit SHA; every definition-changing commit
  names the affected observable/schema/protocol.

## Required pre-push gate

```powershell
git diff --cached --check
py -3 -m unittest tests/test_collective_modes_smoke.py -v
git status --short
git ls-files | Select-String -Pattern '\.(dump|lammpstrj|dcd|xtc|trr|restart|h5|hdf5|npy|npz)$'
```

The final command must return no tracked large trajectories or states. Review
`git diff --cached` manually for sensitive absolute paths, job records, and
assets whose publication status is unclear.

## What Git does not solve

Git gives source history, reversible small-file changes, cross-host code
delivery, and an unambiguous code version for each analysis. It does not
replace trajectory archiving, backups, file checksums, data repositories, or
scientific qualification. Do not use Git LFS by default for multi-GB dumps;
that makes cross-cluster data ownership and recovery fragile.
