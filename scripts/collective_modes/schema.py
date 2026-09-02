"""Data contracts and capability gates for collective-dynamics analysis."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Literal


WallModel = Literal["implicit", "explicit_fixed", "explicit_flexible", "unknown"]
AxisSource = Literal["box_center", "fixed", "cnt_atoms", "unknown"]


@dataclass(frozen=True)
class ObservableRequirement:
    name: str
    fields: frozenset[str]
    needs_molecule_ids: bool = False
    needs_cnt_atoms: bool = False


REQUIREMENTS = {
    "audit": ObservableRequirement("audit", frozenset()),
    "axial_isf": ObservableRequirement("axial_isf", frozenset({"id", "z"})),
    "cylindrical_isf": ObservableRequirement("cylindrical_isf", frozenset({"id", "x", "y", "z"})),
    "axial_current": ObservableRequirement("axial_current", frozenset({"id", "z", "vz"})),
    "cylindrical_current": ObservableRequirement("cylindrical_current", frozenset({"id", "x", "y", "z", "vx", "vy", "vz"})),
    "vacf_z": ObservableRequirement("vacf_z", frozenset({"id", "vz"})),
    "vacf_cylindrical": ObservableRequirement("vacf_cylindrical", frozenset({"id", "x", "y", "vx", "vy"})),
    "molecular_rotation": ObservableRequirement("molecular_rotation", frozenset({"id", "mol", "x", "y", "z", "vx", "vy", "vz"}), True),
    "wall_relative": ObservableRequirement("wall_relative", frozenset({"id", "x", "y", "z", "vx", "vy", "vz"}), False, True),
}


@dataclass
class CaseProfile:
    """Explicit protocol metadata.  Auto-detection fills capabilities, not physics."""

    case_id: str
    dump_paths: list[Path]
    wall_model: WallModel = "unknown"
    axis_source: AxisSource = "unknown"
    rcnt_A: float | None = None
    oxygen_type: int | None = None
    fluid_types: tuple[int, ...] | None = None
    cnt_types: tuple[int, ...] = ()
    integration_timestep_ps: float | None = None
    velocity_frame_default: str = "selected_com"
    protocol_label: str = "unspecified"
    fluid_kind: str = "auto"
    source_locators: list[str] = field(default_factory=list)

    @property
    def selected_types(self) -> tuple[int, ...] | None:
        if self.fluid_types is not None:
            return self.fluid_types
        return (self.oxygen_type,) if self.oxygen_type is not None else None

    def require(self, observable: str, available_fields: Iterable[str], has_cnt_atoms: bool = False) -> None:
        if observable not in REQUIREMENTS:
            raise ValueError(f"Unknown observable requirement: {observable}")
        requirement = REQUIREMENTS[observable]
        missing = requirement.fields - set(available_fields)
        if missing:
            raise ValueError(
                f"{self.case_id}: cannot construct {observable}; missing dump fields {sorted(missing)}"
            )
        if requirement.needs_cnt_atoms and not has_cnt_atoms:
            raise ValueError(
                f"{self.case_id}: {observable} needs declared CNT atom types and CNT coordinates/velocities; "
                "a water-only dump cannot provide a flexible-wall frame"
            )
        if observable.startswith("cylindrical") and self.axis_source == "unknown":
            raise ValueError(
                f"{self.case_id}: cylindrical observable needs --axis-source and, for analytic CNT, --rcnt-A"
            )
        if self.wall_model == "implicit" and observable == "wall_relative":
            raise ValueError(f"{self.case_id}: implicit CNT has no CNT atom velocities for wall-relative VACF")
