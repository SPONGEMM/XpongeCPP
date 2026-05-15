"""First-wave legacy-compatible analysis helpers."""

from __future__ import annotations

import numpy as np

from . import wham
from .md_analysis import XpongeMoleculeReader, mda
from .sasa import SASA


class MdoutReader:
    """Legacy-compatible reader for SPONGE ``mdout`` text files."""

    def __init__(self, filename):
        with open(filename, encoding="utf-8") as handle:
            self.content = handle.readline().split()
            words = handle.read().replace("****", "nan").split()
        width = len(self.content)
        self.data = np.array(
            [[float(words[i + j]) for j in range(width)] for i in range(0, len(words), width)]
        )
        self.content_index = {self.content[i]: i for i in range(width)}

    def __getattribute__(self, attr):
        if attr not in ("content", "content_index", "data") and attr in object.__getattribute__(self, "content_index"):
            return object.__getattribute__(self, "data")[:, object.__getattribute__(self, "content_index")[attr]]
        return object.__getattribute__(self, attr)


__all__ = ["MdoutReader", "SASA", "XpongeMoleculeReader", "mda", "wham"]
