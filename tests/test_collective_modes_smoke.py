"""Synthetic end-to-end smoke test for canonical collective-mode commands."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "scripts" / "collective_modes_cli.py"


def synthetic_dump(path: Path) -> None:
    lines: list[str] = []
    for step in range(0, 80, 10):
        lines += ["ITEM: TIMESTEP\n", f"{step}\n", "ITEM: NUMBER OF ATOMS\n", "3\n", "ITEM: BOX BOUNDS pp pp pp\n", "0 10\n", "0 10\n", "0 20\n", "ITEM: ATOMS id mol type x y z vx vy vz ix iy iz\n"]
        shift = step / 80.0
        lines += [f"1 1 1 {5+shift:.6f} 5.0 {1+shift:.6f} 0.10 0.02 0.20 0 0 0\n", f"2 2 1 5.0 {6+shift:.6f} {7+shift:.6f} -0.10 0.03 -0.20 0 0 0\n", f"3 3 2 1.0 1.0 1.0 0.0 0.0 0.0 0 0 0\n"]
    path.write_text("".join(lines), encoding="utf-8")


class CanonicalCommandsSmokeTest(unittest.TestCase):
    def run_cli(self, *args: str) -> None:
        result = subprocess.run([sys.executable, str(CLI), *args], cwd=ROOT, capture_output=True, text=True)
        if result.returncode:
            self.fail(f"command failed: {' '.join(args)}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}")

    def test_audit_isf_current_vacf_construct_and_plot(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            temp_path = Path(temp); dump = temp_path / "water.dump"; synthetic_dump(dump)
            common = ["--case-id", "synthetic", "--dumps", str(dump), "--fluid-types", "1", "--timestep-ps", "0.001", "--wall-model", "implicit", "--axis-source", "box_center", "--rcnt-A", "4.0", "--max-frames", "8"]
            audit_dir = temp_path / "audit"; self.run_cli("audit", *common, "--output", str(audit_dir))
            with (audit_dir / "dump_capabilities.csv").open() as handle:
                row = next(csv.DictReader(handle)); self.assertEqual(row["supports_cylindrical_current"], "True")
            ambiguous_dir = temp_path / "ambiguous"; self.run_cli("audit", "--case-id", "unknown-wall", "--dumps", str(dump), "--fluid-types", "1", "--output", str(ambiguous_dir))
            with (ambiguous_dir / "dump_capabilities.csv").open() as handle:
                row = next(csv.DictReader(handle)); self.assertEqual(row["wall_model_inferred"], "ambiguous_water_only")
            isf_dir = temp_path / "isf"; self.run_cli("isf", *common, "--output", str(isf_dir), "--n", "1", "--m", "0:1", "--max-lag-ps", "0.06")
            current_dir = temp_path / "current"; self.run_cli("current", *common, "--output", str(current_dir), "--n", "1", "--m", "0:1", "--max-lag-ps", "0.06")
            vacf_dir = temp_path / "vacf"; self.run_cli("vacf", *common, "--output", str(vacf_dir), "--component", "z", "--max-lag-ps", "0.06")
            weights = temp_path / "weights.csv"; weights.write_text("n,m,weight\n1,0,0.2\n1,1,0.1\n", encoding="utf-8")
            construct_dir = temp_path / "construct"; self.run_cli("construct", "--current-csv", str(current_dir / "current_per_replica.csv"), "--isf-csv", str(isf_dir / "isf_per_replica.csv"), "--weights-csv", str(weights), "--vacf-csv", str(vacf_dir / "vacf_per_replica.csv"), "--output", str(construct_dir))
            fit_dir = temp_path / "fit"; self.run_cli("fit-current", "--current-csv", str(current_dir / "current_per_replica.csv"), "--output", str(fit_dir), "--fit-max-ps", "0.06")
            plot = temp_path / "plot.png"; self.run_cli("plot", "--csv", str(construct_dir / "constructibility_sum.csv"), "--x", "lag_ps", "--y", "construct_sum_WFsPhi", "--output", str(plot))
            for expected in [isf_dir / "isf_ensemble_mean_sem.csv", current_dir / "current_cross_ordered_per_replica.csv", vacf_dir / "msd_alpha_from_vacf_ensemble_mean_sem.csv", construct_dir / "constructibility_sum.csv", construct_dir / "constructibility_vs_direct_vacf.csv", fit_dir / "current_mode_DHO_parameters.csv", plot]:
                self.assertTrue(expected.is_file(), expected)


if __name__ == "__main__":
    unittest.main()
