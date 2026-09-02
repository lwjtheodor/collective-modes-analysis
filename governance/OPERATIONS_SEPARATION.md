# Operational separation

This repository is deliberately limited to scientific post-processing,
observable definitions, tests, small governance metadata, and compact
provenance. It must not contain executable job submission, scheduler polling,
remote shell/SSH/SCP, remote storage probes, or automated data-discovery code.

On 2026-09-02, the seven tracked operational utilities present in the initial
baseline were moved to a local operations archive outside this repository.
They are retained there only as historical execution evidence and are not part
of the post-processing mainline or its Git sharing contract.

Future changes must preserve this boundary: analysis code accepts explicit
local input paths and writes reproducible local CSV/metadata products; an
independent operations layer is responsible for obtaining data and executing
that code in an appropriate environment.
