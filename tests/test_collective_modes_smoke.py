"""Synthetic end-to-end smoke test for canonical collective-mode commands."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from scripts.collective_modes.core import cylindrical_currents


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "scripts" / "collective_modes_cli.py"


def synthetic_dump(path: Path) -> None:
    lines: list[str] = []
    for step in range(0, 80, 10):
        lines += ["ITEM: TIMESTEP\n", f"{step}\n", "ITEM: NUMBER OF ATOMS\n", "3\n", "ITEM: BOX BOUNDS pp pp pp\n", "0 10\n", "0 10\n", "0 20\n", "ITEM: ATOMS id mol type x y z vx vy vz ix iy iz\n"]
        shift = step / 80.0
        lines += [f"1 1 1 {6+shift:.6f} 5.0 {1+shift:.6f} 0.10 0.02 0.20 0 0 0\n", f"2 2 1 5.0 {6+shift:.6f} {7+shift:.6f} -0.10 0.03 -0.20 0 0 0\n", f"3 3 2 1.0 1.0 1.0 0.0 0.0 0.0 0 0 0\n"]
    path.write_text("".join(lines), encoding="utf-8")


def axial_dump(path: Path, steps=range(0, 80, 10)) -> None:
    lines: list[str] = []
    for step in steps:
        lines += ["ITEM: TIMESTEP\n", f"{step}\n", "ITEM: NUMBER OF ATOMS\n", "2\n", "ITEM: BOX BOUNDS pp pp pp\n", "0 10\n", "0 10\n", "0 20\n", "ITEM: ATOMS id type z vz\n", f"1 1 {1+step/100:.6f} 0.20\n", f"2 1 {7-step/100:.6f} -0.20\n"]
    path.write_text("".join(lines), encoding="utf-8")


class CanonicalCommandsSmokeTest(unittest.TestCase):
    def run_cli(self, *args: str) -> None:
        result = subprocess.run([sys.executable, str(CLI), *args], cwd=ROOT, capture_output=True, text=True)
        if result.returncode:
            self.fail(f"command failed: {' '.join(args)}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}")

    def test_audit_isf_current_vacf_construct_and_plot(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            temp_path = Path(temp); dump = temp_path / "water.dump"; synthetic_dump(dump)
            common = ["--case-id", "synthetic", "--dumps", str(dump), "--fluid-types", "1", "--timestep-ps", "0.001", "--wall-model", "implicit", "--axis-source", "box_center", "--max-frames", "8"]
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
            plot = temp_path / "plot.png"; self.run_cli("plot", "--csv", str(construct_dir / "constructibility_sum_ensemble_mean_sem.csv"), "--x", "lag_ps", "--y", "construct_sum_WFsPhi_mean", "--output", str(plot))
            for expected in [isf_dir / "isf_ensemble_mean_sem.csv", current_dir / "current_cross_ordered_per_replica.csv", vacf_dir / "msd_alpha_from_vacf_ensemble_mean_sem.csv", construct_dir / "constructibility_sum_per_replica.csv", construct_dir / "constructibility_sum_ensemble_mean_sem.csv", construct_dir / "constructibility_vs_direct_vacf_per_replica.csv", fit_dir / "current_mode_fit_per_replica.csv", plot]:
                self.assertTrue(expected.is_file(), expected)

    def test_axial_minimal_field_dump_runs_all_axial_observables(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            temp_path = Path(temp); dump = temp_path / "oxygen_z_vz.dump"; axial_dump(dump)
            common = ["--case-id", "axial", "--dumps", str(dump), "--fluid-types", "1", "--timestep-ps", "0.001", "--max-frames", "8"]
            self.run_cli("isf", *common, "--output", str(temp_path / "isf"), "--n", "1", "--m", "0", "--max-lag-ps", "0.06")
            self.run_cli("current", *common, "--output", str(temp_path / "current"), "--n", "1", "--m", "0", "--max-lag-ps", "0.06")
            self.run_cli("vacf", *common, "--output", str(temp_path / "vacf"), "--component", "z", "--max-lag-ps", "0.06")
            with (temp_path / "current" / "current_per_replica.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual({row["channel"] for row in rows}, {"Jz", "L"})
            self.assertIn("CJJ_per_particle", rows[0])

    def test_ordered_segments_deduplicate_restart_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            temp_path = Path(temp)
            first = temp_path / "segment_001.dump"; axial_dump(first, range(0, 50, 10))
            second = temp_path / "segment_002.dump"; axial_dump(second, range(40, 80, 10))
            output = temp_path / "current"
            self.run_cli("current", "--case-id", "joined", "--replica", f"rep1={first},{second}", "--fluid-types", "1", "--timestep-ps", "0.001", "--output", str(output), "--n", "1", "--m", "0", "--max-lag-ps", "0.06")
            with (output / "current_per_replica.csv").open() as handle:
                row0 = next(csv.DictReader(handle))
            self.assertEqual(row0["replica"], "rep1")
            self.assertEqual(row0["n_time_origins"], "8")
            bad = temp_path / "bad_boundary.dump"; axial_dump(bad, range(40, 80, 10))
            bad.write_text(bad.read_text(encoding="utf-8").replace("1 1 1.400000 0.20", "1 1 1.401000 0.20", 1), encoding="utf-8")
            result = subprocess.run([sys.executable, str(CLI), "vacf", "--case-id", "joined", "--replica", f"rep1={first},{bad}", "--fluid-types", "1", "--timestep-ps", "0.001", "--output", str(temp_path / "bad"), "--component", "z", "--max-lag-ps", "0.06"], cwd=ROOT, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("differs from the preceding segment", result.stderr)

    def test_vacf_stitch_keeps_native_cadence_layers_and_nonuniform_integral(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            temp_path = Path(temp)
            fine = temp_path / "vacf_1fs.csv"; coarse = temp_path / "vacf_20fs.csv"
            header = "case_id,replica,component,lag_ps,VACF,n_time_origins\n"
            fine.write_text(header + "case,rep1,z,0.0,1.0,100\ncase,rep1,z,0.01,0.9,90\ncase,rep1,z,0.02,0.8,80\n", encoding="utf-8")
            coarse.write_text(header + "case,rep1,z,0.0,1.0,100\ncase,rep1,z,0.02,0.8,80\ncase,rep1,z,0.04,0.6,60\ncase,rep1,z,0.06,0.4,40\n", encoding="utf-8")
            manifest = temp_path / "layers.json"
            manifest.write_text('{"layers":[{"layer_id":"1fs","csv":"' + str(fine).replace("\\", "\\\\") + '","lag_min_ps":0.0,"lag_max_ps":0.02},{"layer_id":"20fs","csv":"' + str(coarse).replace("\\", "\\\\") + '","lag_min_ps":0.02,"lag_max_ps":0.06,"include_lag_min":false}]}', encoding="utf-8")
            output = temp_path / "stitched"
            self.run_cli("vacf-stitch", "--layer-manifest", str(manifest), "--output", str(output))
            with (output / "vacf_stitched_per_replica.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual([float(row["lag_ps"]) for row in rows], [0.0, 0.01, 0.02, 0.04, 0.06])
            self.assertEqual([row["layer_id"] for row in rows], ["1fs", "1fs", "1fs", "20fs", "20fs"])
            with (output / "msd_alpha_from_stitched_vacf_per_replica.csv").open() as handle:
                derived = list(csv.DictReader(handle))
            self.assertTrue(np.isfinite(float(derived[-1]["alpha_from_VACF"])))

    def test_cylindrical_phase_and_projection_identities(self) -> None:
        xyz = np.asarray([[6.0, 5.0, 2.0], [5.0, 6.0, 7.0]])
        velocity = np.asarray([[0.3, 1.0, 2.0], [-0.7, 2.0, -1.0]])
        n = np.asarray([1]); m0 = np.asarray([0]); m1 = np.asarray([1])
        mode_a = cylindrical_currents(xyz, velocity, 20.0, np.asarray([5.0, 5.0]), 4.0, n, m1)
        mode_b = cylindrical_currents(xyz, velocity, 20.0, np.asarray([5.0, 5.0]), 9.0, n, m1)
        theta = np.arctan2(xyz[:, 1] - 5.0, xyz[:, 0] - 5.0)
        expected = np.sum(velocity[:, 2] * np.exp(-1j * ((2*np.pi/20.0) * xyz[:, 2] + theta)))
        self.assertTrue(np.allclose(mode_a["Jz"][0, 0], expected), "cylindrical current must use integer m*theta phase")
        self.assertTrue(np.allclose(mode_a["Jz"], mode_b["Jz"]), "R_mode must not alter m*theta Fourier phase")
        mzero = cylindrical_currents(xyz, velocity, 20.0, np.asarray([5.0, 5.0]), 4.0, n, m0)
        self.assertTrue(np.allclose(mzero["L"], mzero["Jz"]))
        self.assertTrue(np.allclose(mzero["Tinplane"], mzero["Jtheta"]))
        nzero = cylindrical_currents(xyz, velocity, 20.0, np.asarray([5.0, 5.0]), 4.0, np.asarray([0]), m1)
        self.assertTrue(np.allclose(nzero["L"], nzero["Jtheta"]))

    def test_construct_matches_replicas_without_cross_mixing(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            temp_path = Path(temp)
            current = temp_path / "current.csv"; isf = temp_path / "isf.csv"; weights = temp_path / "weights.csv"
            current.write_text("case_id,replica,channel,n,m,lag_ps,CJJ_normalized\ncase,1,L,1,0,0.0,2.0\ncase,2,L,1,0,0.0,5.0\n", encoding="utf-8")
            isf.write_text("case_id,replica,n,m,lag_ps,F_self\ncase,1,1,0,0.0,3.0\ncase,2,1,0,0.0,7.0\n", encoding="utf-8")
            weights.write_text("n,m,weight\n1,0,0.1\n", encoding="utf-8")
            output = temp_path / "construct"
            self.run_cli("construct", "--current-csv", str(current), "--isf-csv", str(isf), "--weights-csv", str(weights), "--output", str(output))
            with (output / "constructibility_sum_per_replica.csv").open() as handle:
                values = {row["replica"]: float(row["construct_sum_WFsPhi"]) for row in csv.DictReader(handle)}
            self.assertAlmostEqual(values["1"], 0.6)
            self.assertAlmostEqual(values["2"], 3.5)
            duplicate = temp_path / "duplicate_isf.csv"
            duplicate.write_text(isf.read_text(encoding="utf-8") + "case,2,1,0,0.0,7.0\n", encoding="utf-8")
            result = subprocess.run([sys.executable, str(CLI), "construct", "--current-csv", str(current), "--isf-csv", str(duplicate), "--weights-csv", str(weights), "--output", str(temp_path / "duplicate_out")], cwd=ROOT, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("duplicate key", result.stderr)


if __name__ == "__main__":
    unittest.main()
