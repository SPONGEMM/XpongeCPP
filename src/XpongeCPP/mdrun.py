"""Legacy-compatible SPONGE runner helpers."""

from __future__ import annotations

import os
import sys
from pathlib import Path

from ._compat.imports import Xopen, Xprint

_THIS_PATH = Path(__file__).resolve().parent


def _runtime_state_dir():
    root = os.environ.get("XPONGECPP_STATE_DIR")
    if root:
        return Path(root)
    xdg = os.environ.get("XDG_CACHE_HOME")
    if xdg:
        return Path(xdg) / "XpongeCPP"
    return Path.home() / ".cache" / "XpongeCPP"


_BIN_PATH_FILE = _runtime_state_dir() / "BIN_PATH.dat"


def _read_bin_path():
    if not _BIN_PATH_FILE.exists():
        _BIN_PATH_FILE.parent.mkdir(parents=True, exist_ok=True)
        with Xopen(str(_BIN_PATH_FILE), "w") as handle:
            handle.write("../../bin")
    with Xopen(str(_BIN_PATH_FILE), "r") as handle:
        that_path = handle.read().strip()
    if not os.path.isabs(that_path):
        that_path = os.path.join(str(_THIS_PATH), that_path)
    return that_path


def run(args):
    """Run SPONGE in the legacy ``Xponge.mdrun`` calling style."""

    if isinstance(args, str):
        args = args.split()
        args.insert(0, sys.argv[0])

    if len(args) < 2 or args[1] in ("-h", "--help"):
        Xprint(
            """ mdrun: run the SPONGE md simulation
        Usage:
            mdrun, mdrun -h, mdrun --help: see this help
            mdrun -set BIN_PATH: set the SPONGE path to BIN_PATH
                                 BIN_PATH can be an absolute path or a relative path to this module file
            mdrun -reset: reset the SPONGE path to ../../bin
            mdrun SPONGE*:  run SPONGE"""
        )
        sys.exit()

    if args[1] == "-set":
        _BIN_PATH_FILE.parent.mkdir(parents=True, exist_ok=True)
        with Xopen(str(_BIN_PATH_FILE), "w") as handle:
            handle.write(args[2])
        sys.exit()
    elif args[1] == "-reset":
        _BIN_PATH_FILE.parent.mkdir(parents=True, exist_ok=True)
        with Xopen(str(_BIN_PATH_FILE), "w") as handle:
            handle.write("../../bin")
        sys.exit()

    that_path = _read_bin_path()
    cmd = os.path.join(that_path, args[1])
    if not (os.path.exists(cmd) or os.path.exists(cmd + ".exe")):
        probe = os.system(args[1] + f" -v > {os.devnull}")
        if probe != 0:
            Xprint(
                "No MD Engine found.\n"
                f"  There is no executable program named '{args[1]}' in '{that_path}' or PATH\n"
                "  Maybe you need to use Xponge.mdrun -set SPONGE_PATH to set the path to MD Engine, "
                "or add the path to your environment variables",
            )
            sys.exit(1)
        cmd = args[1]
    if len(args) > 2:
        cmd += " " + " ".join(args[2:]) + f" < {os.devnull}"
    return os.system(cmd)


__all__ = ["run"]
