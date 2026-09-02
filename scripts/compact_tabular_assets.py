#!/usr/bin/env python3
"""Convert numeric, row-oriented CSV analysis assets into compact HDF5 tables.

The source CSV is never changed.  Each HDF5 file contains a ``table`` group
with one compressed dataset per column, plus enough attributes to reconstruct
a pandas DataFrame.  This avoids a Parquet dependency while keeping ordinary
Python use simple (``h5py`` + ``pandas``/``numpy``).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def classify(value: str) -> str:
    try:
        int(value)
        return "int"
    except ValueError:
        try:
            float(value)
            return "float"
        except ValueError:
            return "str"


def infer_schema(source: Path) -> tuple[list[str], dict[str, str]]:
    with source.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames or []
        schema = {name: "int" for name in names}
        for i, row in enumerate(reader):
            for name in names:
                kind = classify(row[name])
                if schema[name] == "int" and kind == "float":
                    schema[name] = "float"
                elif kind == "str":
                    schema[name] = "str"
            if i >= 4095:
                break
    return names, schema


def csv_row_count(source: Path) -> int:
    with source.open("rb") as handle:
        return max(0, sum(block.count(b"\n") for block in iter(lambda: handle.read(8 * 1024 * 1024), b"")) - 1)


def convert(source: Path, target: Path, *, float64: bool, compression: int) -> dict:
    columns, schema = infer_schema(source)
    rows = csv_row_count(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    float_dtype = np.float64 if float64 else np.float32
    string_dtype = h5py.string_dtype(encoding="utf-8")
    datasets = {}

    with h5py.File(target, "w") as h5:
        table = h5.create_group("table")
        for name in columns:
            kind = schema[name]
            dtype = np.int32 if kind == "int" else (float_dtype if kind == "float" else string_dtype)
            datasets[name] = table.create_dataset(
                name, shape=(rows,), dtype=dtype, compression="gzip",
                compression_opts=compression, shuffle=(kind != "str"),
            )
        # Buffered column writes are orders of magnitude faster than assigning
        # one HDF5 scalar at a time for multi-million-row correlation tables.
        with source.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            start = 0
            buffer = {name: [] for name in columns}
            for row in reader:
                for name in columns:
                    buffer[name].append(row[name])
                if len(buffer[columns[0]]) == 100_000:
                    stop = start + len(buffer[columns[0]])
                    for name in columns:
                        kind = schema[name]
                        if kind == "int":
                            values = np.asarray(buffer[name], dtype=np.int32)
                        elif kind == "float":
                            values = np.asarray(buffer[name], dtype=float_dtype)
                        else:
                            values = np.asarray(buffer[name], dtype=object)
                        datasets[name][start:stop] = values
                        buffer[name].clear()
                    start = stop
            if buffer[columns[0]]:
                stop = start + len(buffer[columns[0]])
                for name in columns:
                    kind = schema[name]
                    if kind == "int":
                        values = np.asarray(buffer[name], dtype=np.int32)
                    elif kind == "float":
                        values = np.asarray(buffer[name], dtype=float_dtype)
                    else:
                        values = np.asarray(buffer[name], dtype=object)
                    datasets[name][start:stop] = values
        table.attrs["columns_json"] = json.dumps(columns)
        table.attrs["schema_json"] = json.dumps(schema, sort_keys=True)
        table.attrs["rows"] = rows
        h5.attrs["format"] = "analysis-hdf5-columnar-v1"
        h5.attrs["source_filename"] = source.name
        h5.attrs["source_sha256"] = file_sha256(source)
        h5.attrs["float_storage"] = "float64" if float64 else "float32"
        h5.attrs["reader"] = "h5py: pandas.DataFrame({key: ds[:] for key, ds in h5['table'].items()})"

    return {
        "source": source.name,
        "target": target.name,
        "rows": rows,
        "columns": columns,
        "schema": schema,
        "source_bytes": source.stat().st_size,
        "target_bytes": target.stat().st_size,
        "source_sha256": file_sha256(source),
        "target_sha256": file_sha256(target),
        "float_storage": "float64" if float64 else "float32",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, nargs="+", help="CSV source file(s)")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--float64", action="store_true", help="preserve float64 instead of compact float32")
    parser.add_argument("--compression", type=int, default=4, choices=range(1, 10))
    args = parser.parse_args()
    records = []
    for source in args.source:
        target = args.output_dir / f"{source.stem}.h5"
        records.append(convert(source, target, float64=args.float64, compression=args.compression))
        print(json.dumps(records[-1], ensure_ascii=False))
    manifest = args.output_dir / "manifest.json"
    manifest.write_text(json.dumps(records, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"manifest={manifest}")


if __name__ == "__main__":
    main()
