import importlib


CORE_SCOPE_LEGACY_IMPORTS = [
    "Xponge",
    "Xponge.analysis",
    "Xponge.analysis.md_analysis",
    "Xponge.analysis.sasa",
    "Xponge.analysis.wham",
    "Xponge.assign",
    "Xponge.build",
    "Xponge.forcefield.amber.bsc1",
    "Xponge.forcefield.amber.ff14sb",
    "Xponge.forcefield.amber.ff19sb",
    "Xponge.forcefield.amber.gaff",
    "Xponge.forcefield.amber.ol3",
    "Xponge.forcefield.amber.rsff2c",
    "Xponge.forcefield.amber.tip3p",
    "Xponge.forcefield.amber.tip4pew",
    "Xponge.forcefield.base.angle_base",
    "Xponge.forcefield.base.cmap_base",
    "Xponge.forcefield.base.charge_base",
    "Xponge.forcefield.base.dihedral_base",
    "Xponge.forcefield.base.lj_base",
    "Xponge.forcefield.base.mass_base",
    "Xponge.forcefield.base.nb14_base",
    "Xponge.forcefield.base.virtual_atom_base",
    "Xponge.forcefield.charmm.charmm27",
    "Xponge.forcefield.charmm.charmm36",
    "Xponge.forcefield.charmm.tip3p_charmm",
    "Xponge.forcefield.opls.oplsaam",
    "Xponge.forcefield.special.gb",
    "Xponge.forcefield.special.fep",
    "Xponge.forcefield.special.min",
    "Xponge.forcefield.sw.mw",
    "Xponge.helper",
    "Xponge.helper.cv",
    "Xponge.helper.file",
    "Xponge.helper.gromacs",
    "Xponge.helper.math",
    "Xponge.helper.namespace",
    "Xponge.load",
    "Xponge.mdrun",
    "Xponge.process",
    "Xponge.tools",
    "Xponge.tools.unittests",
]


def test_core_scope_legacy_import_matrix():
    failures = []
    for module_name in CORE_SCOPE_LEGACY_IMPORTS:
        try:
            importlib.import_module(module_name)
        except Exception as exc:  # pragma: no cover - this is the failure path we want reported
            failures.append((module_name, repr(exc)))
    assert not failures, f"legacy imports failed: {failures}"
