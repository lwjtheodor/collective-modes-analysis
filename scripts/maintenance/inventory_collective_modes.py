#!/usr/bin/env python3
"""Read-only inventory for the collective-dynamics project.

The inventory is intentionally metadata-first: all files are listed, while
only code/small text assets are hashed.  It therefore scales to multi-GiB dump
archives without turning an organizational scan into a data migration.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_SUFFIXES = {".py", ".sh", ".ps1", ".pbs", ".bash", ".ipynb"}
TEXT_SUFFIXES = SCRIPT_SUFFIXES | {".md", ".yaml", ".yml", ".toml", ".json", ".txt"}
REMOTE_RE = re.compile(r"/lustre/home/users/ewu(?:/[A-Za-z0-9_.+@%=-]+)+")
RAW_SUFFIXES = {".dump", ".lammpstrj", ".dcd", ".xtc", ".trr", ".restart"}
SKIP_PARTS = {".git", ".codex", "__pycache__"}


def rel_or_dot(path: Path, root: Path) -> str:
    value = path.relative_to(root).as_posix()
    return value or "."


def role_for(relative: str, suffix: str) -> str:
    first = relative.split("/", 1)[0]
    if suffix in RAW_SUFFIXES:
        return "raw_trajectory_or_restart"
    if first == "scripts":
        return "canonical_candidate_script"
    if first == "results":
        return "analysis_archive_or_derived_asset"
    if first == "assets":
        return "display_asset"
    if first in {"remote_fetch", "heartbeat_fetch"}:
        return "fetched_execution_snapshot"
    if first.startswith("stage_") or first in {"staging", "lowfreq", "highfreq", "high_protocol_1ps"}:
        return "stage_or_exploratory_snapshot"
    if suffix in SCRIPT_SUFFIXES:
        return "unclassified_script"
    return "project_support_or_other"


def sha256_if_small(path: Path, limit_bytes: int) -> str:
    if path.stat().st_size > limit_bytes:
        return ""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-hash-mib", type=float, default=16.0)
    args = parser.parse_args()

    root = args.root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    limit_bytes = int(args.max_hash_mib * 1024 * 1024)
    files: list[dict[str, object]] = []
    scripts: list[dict[str, object]] = []
    remote_rows: list[dict[str, object]] = []
    hash_groups: dict[str, list[str]] = defaultdict(list)
    basename_groups: dict[str, list[str]] = defaultdict(list)
    top_level: Counter[str] = Counter()
    role_counts: Counter[str] = Counter()
    byte_counts: Counter[str] = Counter()

    for path in root.rglob("*"):
        if not path.is_file() or any(part in SKIP_PARTS for part in path.parts):
            continue
        try:
            stat = path.stat()
        except OSError:
            continue
        relative = rel_or_dot(path, root)
        if relative.startswith(rel_or_dot(output, root) + "/"):
            continue
        suffix = path.suffix.lower()
        top = relative.split("/", 1)[0]
        role = role_for(relative, suffix)
        small_hash = sha256_if_small(path, limit_bytes) if suffix in TEXT_SUFFIXES else ""
        record = {
            "relative_path": relative,
            "bytes": stat.st_size,
            "modified_utc": datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
            "suffix": suffix or "[none]",
            "top_level": top,
            "role": role,
            "sha256_if_small": small_hash,
        }
        files.append(record)
        top_level[top] += 1
        role_counts[role] += 1
        byte_counts[top] += stat.st_size

        if suffix in SCRIPT_SUFFIXES:
            text = read_text(path)
            remote_paths = sorted(set(REMOTE_RE.findall(text)))
            script_record = record | {
                "basename": path.name,
                "remote_path_count": len(remote_paths),
                "remote_paths": " | ".join(remote_paths),
            }
            scripts.append(script_record)
            basename_groups[path.name].append(relative)
            if small_hash:
                hash_groups[small_hash].append(relative)
            for remote_path in remote_paths:
                remote_rows.append({"remote_path": remote_path, "referenced_by": relative})

    files.sort(key=lambda row: str(row["relative_path"]))
    scripts.sort(key=lambda row: str(row["relative_path"]))
    duplicate_rows = [
        {"sha256": digest, "copy_count": len(paths), "paths": " | ".join(sorted(paths))}
        for digest, paths in hash_groups.items() if len(paths) > 1
    ]
    duplicate_rows.sort(key=lambda row: (-int(row["copy_count"]), str(row["paths"])))
    collision_rows = [
        {"basename": name, "variant_count": len(paths), "paths": " | ".join(sorted(paths))}
        for name, paths in basename_groups.items() if len(paths) > 1
    ]
    collision_rows.sort(key=lambda row: (-int(row["variant_count"]), str(row["basename"])))
    unique_remote: dict[str, list[str]] = defaultdict(list)
    for row in remote_rows:
        unique_remote[str(row["remote_path"])].append(str(row["referenced_by"]))
    remote_summary = [
        {"remote_path": path, "reference_count": len(refs), "referenced_by": " | ".join(sorted(refs))}
        for path, refs in unique_remote.items()
    ]
    remote_summary.sort(key=lambda row: (-int(row["reference_count"]), str(row["remote_path"])))

    write_csv(output / "file_inventory.csv", files, list(files[0]) if files else [])
    write_csv(output / "script_inventory.csv", scripts, list(scripts[0]) if scripts else [])
    write_csv(output / "exact_script_duplicates.csv", duplicate_rows, ["sha256", "copy_count", "paths"])
    write_csv(output / "script_basename_collisions.csv", collision_rows, ["basename", "variant_count", "paths"])
    write_csv(output / "remote_path_references.csv", remote_summary, ["remote_path", "reference_count", "referenced_by"])

    summary = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(root),
        "hash_limit_mib": args.max_hash_mib,
        "file_count": len(files),
        "script_count": len(scripts),
        "exact_script_duplicate_clusters": len(duplicate_rows),
        "same_basename_script_clusters": len(collision_rows),
        "unique_ccfep_paths_referenced": len(remote_summary),
        "top_level": {key: {"files": top_level[key], "bytes": byte_counts[key]} for key in sorted(top_level)},
        "roles": dict(sorted(role_counts.items())),
    }
    (output / "inventory_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (output / "README.md").write_text(
        "# Read-only inventory snapshot\n\n"
        f"Generated: `{summary['generated_utc']}`  \n"
        f"Root: `{root}`  \n\n"
        f"- Files: **{summary['file_count']}**\n"
        f"- Scripts: **{summary['script_count']}**\n"
        f"- Exact script-copy clusters: **{summary['exact_script_duplicate_clusters']}**\n"
        f"- Same-basename variant clusters: **{summary['same_basename_script_clusters']}**\n"
        f"- Distinct literal CCFEP paths in scripts: **{summary['unique_ccfep_paths_referenced']}**\n\n"
        "`exact_script_duplicates.csv` is safe deduplication evidence only for identical small script bytes.\n"
        "`script_basename_collisions.csv` is a review queue, not proof of equivalence.\n"
        "`remote_path_references.csv` is a compact remote-inspection scope; it does not establish that a path exists today.\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
