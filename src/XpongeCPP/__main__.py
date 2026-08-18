"""Command-line entry point for XpongeCPP."""

from __future__ import annotations

import argparse

from . import __version__
from .io_bundle.cli import add_legacy_to_bundle_parser


def _add_legacy_test_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "test",
        help="test the basic function of Xponge",
    )
    parser.add_argument(
        "-p",
        "--purpose",
        metavar="programmatic",
        default="programmatic",
        choices=("programmatic", "academic"),
        help="select the programmatic or academic test profile",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        metavar="INFO",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        help="set the test output verbosity",
    )
    parser.add_argument("-d", "--do", metavar="todo", default="base")
    parser.add_argument("-f", "--file", metavar="file")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="XpongeCPP")
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=__version__,
        help="show the version of Xponge",
    )
    subparsers = parser.add_subparsers(dest="command")
    _add_legacy_test_parser(subparsers)
    add_legacy_to_bundle_parser(subparsers)
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
