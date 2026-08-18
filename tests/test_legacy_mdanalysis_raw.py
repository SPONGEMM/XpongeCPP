from __future__ import annotations

import pytest


pytest.importorskip("MDAnalysis")


def test_sponge_input_reader_exposes_origin_format_contract():
    from Xponge.analysis import md_analysis as xmda

    assert xmda.SpongeInputReader.format == "SPONGE_MASS"
    assert xmda.SpongeInputReader._format_hint("system_mass.txt") is True
    assert xmda.SpongeInputReader._format_hint("system_charge.txt") is False

