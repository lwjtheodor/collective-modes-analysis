from __future__ import annotations

import sys

sys.path.insert(0, r"E:\ssh_scp")
from cluster_control.core.ssh_common import execute_command


def main() -> None:
    cmd = "find /lustre/home/users/ewu/vb_gcmc/MD -maxdepth 4 -type f \\( -name '*8_8*L20*.restart' -o -name '*L20*8_8*.restart' -o -name '*8_8*20L*.restart' \\) -printf '%s %p\\n' | sort -n | tail -100"
    out, err, code = execute_command("ccfep", cmd, timeout=120, get_pty=True)
    print(out)
    print(err, file=sys.stderr)
    raise SystemExit(code)


if __name__ == "__main__":
    main()
