"""Command-line entry point for XpongeCPP."""

from __future__ import annotations

import argparse

from . import __version__
from .io_bundle.cli import add_legacy_to_bundle_parser


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="XpongeCPP")
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=__version__,
    )
    subparsers = parser.add_subparsers(dest="command")
    add_legacy_to_bundle_parser(subparsers)
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
