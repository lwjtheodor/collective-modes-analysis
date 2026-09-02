#!/usr/bin/env python3
"""Create a field-aware inventory for LAMMPS dump assets.

The scanner reads at most 128 KiB from each dump, enough to identify the
``ITEM: ATOMS`` schema without loading trajectories.  It does not alter inputs.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


def first_match(pattern: str, text: str) -> str:
    match = re.search(pattern, text, flags=re.I)
    return match.group(1) if match else ""


def atom_fields(path: Path) -> str:
    try:
        with path.open("rb") as handle:
            head = handle.read(131_072).decode("utf-8", errors="replace")
        for line in head.splitlines():
            if line.startswith("ITEM: ATOMS"):
                return " ".join(line.split()[2:])
    except OSError:
        return "<unreadable>"
    return "<header-not-found>"


def classify(path: Path, root: Path) -> dict[str, str]:
    name = path.name
    full = str(path)
    upper = full.upper()
    # Prefer the filename and immediate case directory: campaign names such as
    # ``77_88`` are collection labels, not a physical CNT chirality.
    specific = "_".join((path.name, path.parent.name, path.parent.parent.name)).upper()
    chirality = first_match(r"(?:^|_)(\d+_\d+)(?:_|$)", specific)
    if not chirality:
        compact = first_match(r"(?:^|_)(7{2}|8{2}|9{2}|17_0|15_0|13_0)(?:_|$)", specific)
        chirality = {"77": "7_7", "88": "8_8", "99": "9_9"}.get(compact, compact)
    length = first_match(r"(?:^|[_\\/])L(\d+)(?:[_\\/]|$)", specific)
    if not length:
        length = first_match(r"(?:^|[_\\/])(\d+)L(?:[_\\/]|$)", specific)
    if not length:
        length = first_match(r"(?:^|[_\\/])(\d+)XL(?:[_\\/]|$)", specific)
    cadence = first_match(r"(?:_|\\b)(\d+(?:P\d+)?(?:FS|PS))(?:_|\\b|\.)", upper)
    duration = first_match(r"(?:_|\\b)(\d+(?:P\d+)?NS)(?:_|\\b|\.)", upper)
    if not duration:
        duration = first_match(r"FIRST(\d+(?:P\d+)?(?:NS|PS))", upper)
    rh = first_match(r"RH(\d+)", upper)
    temperature = first_match(r"T(\d+)K", upper)
    protocol = ",".join(tag for tag in ("WEAKNH", "WEAKNVT", "NOMOM", "MOM_OFF", "NVE", "NVT", "RETHERM", "CNTCSVR", "GCMC") if tag in upper)
    species = "oxygen" if "OXYGEN" in upper else ("water" if "WATER" in upper else "mixed/unspecified")
    return {
        "relative_path": str(path.relative_to(root)), "campaign": path.relative_to(root).parts[0],
        "filename": name, "chirality": chirality or "unparsed", "length_L": length or "unparsed",
        "cadence": cadence.lower() or "unparsed", "duration": duration.lower() or "unparsed",
        "temperature_K": temperature or "unparsed", "RH_percent": rh or "unparsed",
        "protocol_tags": protocol or "unparsed", "species_hint": species,
        "atom_fields": atom_fields(path), "bytes": str(path.stat().st_size),
        "local_status": "archive candidate; remote equivalence requires path+bytes/hash receipt",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    rows = [classify(path, root) for path in sorted(root.rglob("*.dump"))]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "local_dump_file_index.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [])
        writer.writeheader()
        writer.writerows(rows)

    groups: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row["campaign"], row["chirality"], row["length_L"], row["cadence"], row["duration"], row["atom_fields"])].append(row)
    md = ["# Local F-drive dump inventory", "", f"Scanned `{root}`: **{len(rows)}** `.dump` files.", "", "## Protocol/field matrix", "", "| campaign | chirality | L | cadence | duration | files | fields |", "|---|---|---:|---|---|---:|---|"]
    for key, group in sorted(groups.items()):
        campaign, chirality, length, cadence, duration, fields = key
        compact_fields = fields.replace("|", "\\|")
        md.append(f"| `{campaign}` | {chirality} | {length} | {cadence} | {duration} | {len(group)} | `{compact_fields}` |")
    md += ["", "## Field availability", "", "| atom fields | dump files |", "|---|---:|"]
    for fields, count in Counter(row["atom_fields"] for row in rows).most_common():
        md.append(f"| `{fields.replace('|', '\\|')}` | {count} |")
    md += ["", "`local_dump_file_index.csv` is the authoritative per-file list. `unparsed` means the filename did not encode that field; it is not evidence that the underlying simulation lacks the parameter."]
    (args.output_dir / "README.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"files={len(rows)} csv={csv_path} matrix={args.output_dir / 'README.md'}")


if __name__ == "__main__":
    main()
