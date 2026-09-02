#!/usr/bin/env python3
"""Copy an initial, whole-frame window from a text LAMMPS dump.

This intentionally keeps the original atom columns and only trims in time, so
the resulting file remains a provenance-preserving input for local diagnostics.
"""
import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--frames", required=True, type=int)
    args = parser.parse_args()
    if args.frames < 2:
        raise ValueError("--frames must be at least 2")

    copied = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.input.open("r", errors="replace") as src, args.output.open("w", newline="") as dst:
        while copied < args.frames:
            marker = src.readline()
            if not marker:
                break
            if marker != "ITEM: TIMESTEP\n":
                raise ValueError(f"unexpected dump marker before frame {copied}: {marker!r}")
            frame = [marker]
            step = src.readline()
            natom_marker = src.readline()
            natom_line = src.readline()
            bounds_marker = src.readline()
            if not (step and natom_marker == "ITEM: NUMBER OF ATOMS\n" and natom_line and bounds_marker.startswith("ITEM: BOX BOUNDS")):
                raise ValueError(f"truncated frame {copied}")
            frame.extend([step, natom_marker, natom_line, bounds_marker])
            for _ in range(3):
                line = src.readline()
                if not line: raise ValueError(f"truncated bounds in frame {copied}")
                frame.append(line)
            atoms_marker = src.readline()
            if not atoms_marker.startswith("ITEM: ATOMS"):
                raise ValueError(f"missing atom header in frame {copied}")
            frame.append(atoms_marker)
            natom = int(natom_line)
            for _ in range(natom):
                line = src.readline()
                if not line: raise ValueError(f"truncated atoms in frame {copied}")
                frame.append(line)
            dst.writelines(frame)
            copied += 1
    if copied != args.frames:
        raise ValueError(f"only found {copied} of requested {args.frames} frames")
    print(f"copied_frames={copied} output={args.output}")


if __name__ == "__main__":
    main()
