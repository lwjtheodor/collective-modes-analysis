"""One streaming LAMMPS text-dump reader with explicit capability detection."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np


@dataclass
class DumpSchema:
    fields: tuple[str, ...]
    atom_types: tuple[int, ...]
    has_molecule_ids: bool
    has_positions: bool
    has_velocities: bool
    has_image_flags: bool
    inferred_content: str
    confidence: str


@dataclass
class Frame:
    timestep: int
    bounds: np.ndarray
    fields: tuple[str, ...]
    values: np.ndarray

    @property
    def columns(self) -> dict[str, int]:
        return {name: index for index, name in enumerate(self.fields)}

    @property
    def box_lengths(self) -> np.ndarray:
        return self.bounds[:, 1] - self.bounds[:, 0]

    @property
    def box_center(self) -> np.ndarray:
        return self.bounds.mean(axis=1)

    def select_types(self, types: tuple[int, ...] | None) -> "Frame":
        if types is None:
            return self
        columns = self.columns
        if "type" not in columns:
            raise ValueError("type selection requested but dump has no type field")
        mask = np.isin(self.values[:, columns["type"]].astype(int), types)
        return Frame(self.timestep, self.bounds, self.fields, self.values[mask])

    def column(self, name: str) -> np.ndarray:
        try:
            return self.values[:, self.columns[name]]
        except KeyError as exc:
            raise ValueError(f"Dump frame at step {self.timestep} lacks field {name!r}") from exc


def _read_frame(handle, path: Path) -> Frame | None:
    marker = handle.readline()
    if not marker:
        return None
    if marker.strip() != "ITEM: TIMESTEP":
        raise ValueError(f"{path}: expected ITEM: TIMESTEP, found {marker!r}")
    timestep = int(handle.readline())
    if handle.readline().strip() != "ITEM: NUMBER OF ATOMS":
        raise ValueError(f"{path}: missing NUMBER OF ATOMS after step {timestep}")
    count = int(handle.readline())
    if not handle.readline().startswith("ITEM: BOX BOUNDS"):
        raise ValueError(f"{path}: missing BOX BOUNDS after step {timestep}")
    bounds = np.array([list(map(float, handle.readline().split()[:2])) for _ in range(3)], dtype=float)
    atom_marker = handle.readline().split()
    if atom_marker[:2] != ["ITEM:", "ATOMS"]:
        raise ValueError(f"{path}: missing ATOMS header after step {timestep}")
    fields = tuple(atom_marker[2:])
    rows = [handle.readline() for _ in range(count)]
    values = np.fromstring(" ".join(rows), sep=" ", dtype=float)
    expected = count * len(fields)
    if values.size != expected:
        raise ValueError(f"{path}: malformed atom block at step {timestep}; expected {expected} values, got {values.size}")
    return Frame(timestep, bounds, fields, values.reshape(count, len(fields)))


def iter_frames(path: Path) -> Iterator[Frame]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        while (frame := _read_frame(handle, path)) is not None:
            yield frame


def inspect_dump(path: Path, max_frames: int = 2) -> DumpSchema:
    frames: list[Frame] = []
    for frame in iter_frames(path):
        frames.append(frame)
        if len(frames) >= max_frames:
            break
    if not frames:
        raise ValueError(f"{path}: contains no readable LAMMPS frames")
    fields = frames[0].fields
    if any(frame.fields != fields for frame in frames[1:]):
        raise ValueError(f"{path}: atom fields change between sampled frames")
    types: Counter[int] = Counter()
    if "type" in fields:
        index = fields.index("type")
        for frame in frames:
            types.update(frame.values[:, index].astype(int).tolist())
    positions = {"x", "y", "z"}.issubset(fields)
    velocities = {"vx", "vy", "vz"}.issubset(fields)
    images = {"ix", "iy", "iz"}.issubset(fields)
    molecule_ids = "mol" in fields
    if positions and velocities and molecule_ids:
        content, confidence = "full_water_or_molecular_fluid", "capability_only"
    elif positions and velocities:
        content, confidence = "site_resolved_3d", "capability_only"
    elif {"z", "vz"}.issubset(fields):
        content, confidence = "axial_site_resolved", "capability_only"
    else:
        content, confidence = "partial_or_unknown", "capability_only"
    return DumpSchema(fields, tuple(sorted(types)), molecule_ids, positions, velocities, images, content, confidence)


def infer_protocol_hint(path: Path, declared_wall_model: str, declared_cnt_types: tuple[int, ...], schema: DumpSchema) -> tuple[str, str]:
    """Return a bounded wall-model hint; never manufacture physics from water fields.

    A water-only dump normally cannot distinguish an implicit analytic cylinder
    from an explicit CNT whose atoms were omitted at output time.  Such a case
    is deliberately reported as ambiguous rather than guessed.
    """
    if declared_wall_model != "unknown":
        return declared_wall_model, "manifest_declared"
    if declared_cnt_types and set(declared_cnt_types) & set(schema.atom_types):
        return "explicit_fixed_or_flexible", "declared_cnt_types_present_in_dump"
    lowered = str(path).lower()
    if "implicit" in lowered:
        return "implicit_hint", "filename_hint_non_authoritative"
    if "explicit" in lowered or "fullwater" in lowered:
        return "explicit_hint", "filename_hint_non_authoritative"
    return "ambiguous_water_only", "dump_fields_do_not_identify_wall_model"


def infer_fluid_hint(path: Path, declared_fluid_kind: str, schema: DumpSchema) -> tuple[str, str]:
    """Classify representation, distinguishing verified facts from filename hints."""
    if declared_fluid_kind != "auto":
        return declared_fluid_kind, "manifest_declared"
    lowered = str(path).lower()
    if "argon" in lowered or "_ar" in lowered:
        return "argon_hint", "filename_hint_non_authoritative"
    if schema.has_molecule_ids:
        return "molecular_fluid", "mol_field_present; species_not_proven"
    if schema.has_positions and schema.has_velocities:
        return "single_site_fluid_ambiguous", "site_fields_only; could be oxygen-only water or monatomic fluid"
    return "unknown_fluid", "insufficient_fields"


def validate_uniform_timestep(frames: list[Frame], timestep_ps: float | None) -> tuple[float | None, int | None]:
    if len(frames) < 2:
        return None, None
    steps = np.array([frame.timestep for frame in frames], dtype=np.int64)
    deltas = np.diff(steps)
    if not np.all(deltas == deltas[0]):
        raise ValueError(f"nonuniform dump-step cadence: {deltas[:8].tolist()}")
    interval = int(deltas[0])
    return (interval * timestep_ps if timestep_ps is not None else None), interval
