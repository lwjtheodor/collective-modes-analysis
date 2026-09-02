"""Read-only CCFEP storage preflight for a prospective 20L VACF campaign."""
from __future__ import annotations

import sys

sys.path.insert(0, r"E:\ssh_scp")
from cluster_control.core.ssh_common import execute_command


def main() -> None:
    command = " ; ".join(
        [
            "df -h /lustre/home/users/ewu",
            "lfs quota -u ewu /lustre",
            "find /lustre/home/users/ewu/vb_gcmc/MD/stage_vacf_tail_8_8_L10_8rep_continue10ns_100fs_tdamp1ns_20260823 -name '*_tdamp1ns_10ns_100fs.dump' -printf '%s %p\\n' | sort -n",
        ]
    )
    out, err, code = execute_command("ccfep", command, timeout=120, get_pty=True)
    print(out)
    print(err, file=sys.stderr)
    raise SystemExit(code)


if __name__ == "__main__":
    main()
