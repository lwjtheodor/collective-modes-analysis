"""Protocol-gated collective-dynamics analysis for CNT-confined fluids.

This package is the canonical successor to the project's scattered dump readers.
It is intentionally metadata-first: source data are never identified solely by
file name, and unsupported observables fail before numerical processing.
"""

from .schema import CaseProfile, ObservableRequirement

__all__ = ["CaseProfile", "ObservableRequirement"]
