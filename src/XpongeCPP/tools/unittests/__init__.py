"""Minimal legacy-compatible test helper namespace.

This preserves the import shape of ``Xponge.tools.unittests`` for migration and
basic scripting, without attempting to replicate the full historical test runner
environment from upstream Xponge.
"""

from __future__ import annotations

import logging
import unittest

from ...helper import GlobalSetting, Xdict, Xopen, Xprint, source

warnings = __import__("warnings")
warnings.filterwarnings("ignore")

CATEGORY = Xdict(
    {
        "0": "base",
        "1": "building",
        "2": "forcefield_loading",
        "3": "forcefield_using",
        "4": "MD_efficiency",
        "5": "MD_function",
        "6": "MD_thermodynamics",
        "7": "MD_kinetics",
        "8": "enhancing_sampling",
        "9": "workflow",
        "100": "application",
    },
    not_found_message="{} is not a valid unittest category",
)


class XpongeTestRunner(unittest.TextTestRunner):
    """Compatibility wrapper matching the historical Xponge test runner shape."""

    def run(self, test):
        result = self._makeResult()
        test(result)
        if result.errors:
            for error in result.errors:
                Xprint(error[1], "ERROR")
        if result.failures:
            for error in result.failures:
                Xprint(error[1], "ERROR")
        return result


def mytest(args):
    """Minimal compatibility entrypoint for historical ``Xponge test`` hooks."""

    if hasattr(args, "verbose"):
        GlobalSetting.logger.setLevel(args.verbose)
    if hasattr(args, "purpose"):
        GlobalSetting.purpose = args.purpose
    Xprint("XpongeCPP minimal unittest compatibility shim is active.")
    return None


__all__ = [
    "CATEGORY",
    "GlobalSetting",
    "Xdict",
    "Xopen",
    "XpongeTestRunner",
    "Xprint",
    "logging",
    "mytest",
    "source",
]
