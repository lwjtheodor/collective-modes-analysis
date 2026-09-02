"""Command-line interface for the canonical collective-dynamics mainline."""

from __future__ import annotations

import argparse
from pathlib import Path

from . import commands


def _case_arguments(parser: argparse.ArgumentParser, cylindrical: bool = False) -> None:
    parser.add_argument("--case-id", required=True); parser.add_argument("--dumps", nargs="+", help="compatibility shorthand: each file is one independent replica")
    parser.add_argument("--replica", action="append", help="repeat ID=segment1,segment2 in physical time order; distinct IDs are independent replicas")
    parser.add_argument("--trajectory-manifest", type=Path, help="JSON manifest with ordered replicas/segments; mutually exclusive with --dumps/--replica")
    parser.add_argument("--output", required=True, type=Path); parser.add_argument("--protocol", default="unspecified")
    parser.add_argument("--fluid-kind", choices=["auto", "water", "oxygen_only", "argon"], default="auto")
    parser.add_argument("--wall-model", choices=["implicit", "explicit_fixed", "explicit_flexible", "unknown"], default="unknown")
    parser.add_argument("--axis-source", choices=["box_center", "fixed", "cnt_atoms", "unknown"], default="box_center" if cylindrical else "unknown")
    parser.add_argument("--axis-xy", nargs=2, type=float)
    parser.add_argument("--oxygen-type", type=int); parser.add_argument("--fluid-types", nargs="*", type=int); parser.add_argument("--cnt-types", nargs="*", type=int)
    parser.add_argument("--timestep-ps", type=float); parser.add_argument("--dt-ps", type=float); parser.add_argument("--max-frames", type=int, default=0)
    parser.add_argument("--velocity-frame", choices=["lab", "selected_com", "wall_relative"], default="selected_com")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="collective-modes", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("audit", help="detect dump fields/capabilities without inferring wall physics"); _case_arguments(p); p.set_defaults(func=commands.audit)
    p = sub.add_parser("isf", help="total/self/distinct cylindrical density ISF"); _case_arguments(p, cylindrical=True); p.add_argument("--n", default="1:5"); p.add_argument("--m", default="0"); p.add_argument("--max-lag-ps", type=float, required=True); p.set_defaults(func=commands.isf)
    p = sub.add_parser("current", help="cylindrical current ACF, ordered cross kernels and spectra"); _case_arguments(p, cylindrical=True); p.add_argument("--n", default="1:5"); p.add_argument("--m", default="0"); p.add_argument("--max-lag-ps", type=float, required=True); p.set_defaults(func=commands.current)
    p = sub.add_parser("vacf", help="VACF, consistent VACF-MSD and alpha"); _case_arguments(p, cylindrical=True); p.add_argument("--component", choices=["z", "r", "theta"], default="z"); p.add_argument("--max-lag-ps", type=float, required=True); p.set_defaults(func=commands.vacf)
    p = sub.add_parser("vacf-stitch", help="stitch separately estimated native-cadence VACFs, then integrate MSD/alpha on the nonuniform lag grid"); p.add_argument("--layer-manifest", type=Path, required=True); p.add_argument("--vacf-column", default="VACF"); p.add_argument("--output", required=True, type=Path); p.set_defaults(func=commands.vacf_stitch)
    p = sub.add_parser("construct", help="build no-free-amplitude sum W Fs Phi_J from canonical CSV products"); p.add_argument("--current-csv", required=True, type=Path); p.add_argument("--isf-csv", required=True, type=Path); p.add_argument("--weights-csv", required=True, type=Path); p.add_argument("--vacf-csv", type=Path); p.add_argument("--output", required=True, type=Path); p.add_argument("--current-channel", default="L"); p.add_argument("--current-column", default="CJJ_normalized"); p.add_argument("--self-column", default="F_self"); p.add_argument("--vacf-column", default="VACF_normalized"); p.set_defaults(func=commands.construct)
    p = sub.add_parser("fit-current", help="per-replica current-kernel fit followed by replica SEM"); p.add_argument("--current-csv", required=True, type=Path); p.add_argument("--output", required=True, type=Path); p.add_argument("--channel", default="L"); p.add_argument("--column", default="CJJ_normalized"); p.add_argument("--model", choices=["damped_carrier", "dho_physical"], default="damped_carrier"); p.add_argument("--fit-min-ps", type=float, default=0.0); p.add_argument("--fit-max-ps", type=float, required=True); p.set_defaults(func=commands.fit_current)
    p = sub.add_parser("plot", help="minimal readable CSV plot template"); p.add_argument("--csv", required=True, type=Path); p.add_argument("--x", required=True); p.add_argument("--y", required=True); p.add_argument("--group", default="case_id"); p.add_argument("--output", required=True, type=Path); p.set_defaults(func=commands.plot)
    return parser


def main() -> None:
    args = build_parser().parse_args(); args.func(args)


if __name__ == "__main__":
    main()
