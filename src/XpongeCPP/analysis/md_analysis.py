"""Legacy-compatible MDAnalysis integration surface."""

from __future__ import annotations

import os.path
import time
from collections.abc import Iterable

import numpy as np

from .. import ResidueType
from .._compat.imports import Xopen
from ..helper.math import guess_element_from_mass

try:  # pragma: no cover - optional dependency branch
    import MDAnalysis as mda
    from MDAnalysis.coordinates import H5MD, base
    from MDAnalysis.core import topologyattrs
    from MDAnalysis.core.topology import Topology
    from MDAnalysis.lib import util
    from MDAnalysis.topology.base import TopologyReaderBase
    from MDAnalysis.topology.guessers import guess_masses
except ModuleNotFoundError:
    mda = None


def _missing_mdanalysis():
    raise ModuleNotFoundError("'MDAnalysis' package needed. Maybe you need 'pip install MDAnalysis'")


if mda is None:
    class XpongeMoleculeReader:  # pragma: no cover - exercised only on missing dependency path
        def __init__(self, *_args, **_kwargs):
            _missing_mdanalysis()


    __all__ = ["XpongeMoleculeReader", "mda"]
else:
    # pylint: disable=abstract-method, arguments-differ, protected-access, unused-argument
    class SpongeNoneReader(base.ReaderBase):
        def __init__(self, _, n_atoms, **kwargs):
            super().__init__(_, **kwargs)
            self._n_atoms = n_atoms

        @property
        def n_atoms(self):
            return self._n_atoms

        @property
        def n_frames(self):
            return 0

        def close(self):
            return


    class SpongeInputReader(TopologyReaderBase):
        def parse(self, **kwargs):
            attrs = [topologyattrs.Segids(np.array(["SYSTEM"], dtype=object))]
            has_names = False
            has_type_names = False
            self.filename = self.filename.replace("_mass.txt", "")
            if os.path.exists(self.filename + "_atom_name.txt"):
                with util.openany(self.filename + "_atom_name.txt") as fm:
                    natoms = int(fm.readline())
                    names = [line.strip() for line in fm]
                    has_names = True
                    attrs.append(topologyattrs.Atomnames(names))
            if os.path.exists(self.filename + "_atom_type_name.txt"):
                with util.openany(self.filename + "_atom_type_name.txt") as fm:
                    natoms = int(fm.readline())
                    names = [line.strip() for line in fm]
                    has_type_names = True
                    attrs.append(topologyattrs.Atomtypes(names))
            if os.path.exists(self.filename + "_mass.txt"):
                with util.openany(self.filename + "_mass.txt") as fm:
                    natoms = int(fm.readline())
                    masses = [float(line.strip()) for line in fm]
                    atom_names = [guess_element_from_mass(mass) for mass in masses]
                    attrs.append(topologyattrs.Masses(masses))
                    if not has_names:
                        attrs.append(topologyattrs.Atomnames(atom_names, guessed=True))
                    if not has_type_names:
                        attrs.append(topologyattrs.Atomtypes(atom_names, guessed=True))
                    attrs.append(topologyattrs.Elements(atom_names, guessed=True))
            if os.path.exists(self.filename + "_charge.txt"):
                with util.openany(self.filename + "_charge.txt") as fm:
                    natoms = int(fm.readline())
                    charges = [float(line.strip()) / 18.2223 for line in fm]
                    attrs.append(topologyattrs.Charges(charges))
            resid = np.zeros(natoms, dtype=np.int32)
            nres = 1
            if os.path.exists(self.filename + "_residue.txt"):
                with util.openany(self.filename + "_residue.txt") as fm:
                    natoms, nres = fm.readline().split()
                    natoms, nres = int(natoms), int(nres)
                    resid = np.zeros(natoms, dtype=np.int32)
                    count = 0
                    for i, line in enumerate(fm):
                        res_length = int(line.strip())
                        resid[count : count + res_length] = i
                        count += res_length
            if os.path.exists(self.filename + "_resname.txt"):
                with util.openany(self.filename + "_resname.txt") as fm:
                    nres = int(fm.readline())
                    resname = [line.strip() for line in fm]
                    attrs.append(topologyattrs.Resnames(resname))
            attrs.append(topologyattrs.Resids(np.arange(nres) + 1))
            attrs.append(topologyattrs.Atomids(np.arange(natoms) + 1))
            attrs.append(topologyattrs.Resnums(np.arange(nres) + 1))
            if os.path.exists(self.filename + "_bond.txt"):
                with util.openany(self.filename + "_bond.txt") as fm:
                    fm.readline()
                    bonds = [[int(words) for words in line.split()[:2]] for line in fm]
                    attrs.append(topologyattrs.Bonds(bonds))
            self._n_atoms = natoms
            return Topology(natoms, nres, 1, attrs, resid, None)


    class XpongeMoleculeReader(base.ReaderBase):
        def __init__(self, filename, **kwargs):
            self.molecule = filename
            self.molecule.get_atoms()
            super().__init__(filename, **kwargs)
            self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
            self.ts.positions = np.array([[getattr(atom, i) for i in "xyz"] for atom in self.molecule.atoms])

        @property
        def n_atoms(self):
            return len(self.molecule.atoms)

        @property
        def n_frames(self):
            return 1

        def parse(self, **kwargs):
            attrs = [topologyattrs.Segids(np.array(["SYSTEM"], dtype=object))]
            molecule = self.molecule
            natoms = len(molecule.atoms)
            nres = len(molecule.residues)
            masses = [atom.mass for atom in molecule.atoms]
            elements = [guess_element_from_mass(mass) for mass in masses]
            attrs.append(topologyattrs.Masses(masses))
            attrs.append(topologyattrs.Elements(elements, guessed=True))
            attrs.append(topologyattrs.Atomnames([atom.name for atom in molecule.atoms]))
            attrs.append(topologyattrs.Atomtypes([atom.type.name for atom in molecule.atoms]))
            attrs.append(topologyattrs.Charges([atom.charge for atom in molecule.atoms]))
            attrs.append(topologyattrs.Resids(np.arange(nres) + 1))
            attrs.append(topologyattrs.Atomids(np.arange(natoms) + 1))
            attrs.append(topologyattrs.Resnums(np.arange(nres) + 1))
            attrs.append(topologyattrs.Resnames([residue.name for residue in molecule.residues]))
            resid = np.zeros(natoms, dtype=np.int32)
            count = 0
            for i, res in enumerate(molecule.residues):
                resid[count : count + len(res.atoms)] = i
                count += len(res.atoms)
            bonds = [
                [molecule.atom_index[self._t2a(residue, ai)], molecule.atom_index[self._t2a(residue, aj)]]
                for residue in molecule.residues
                for ai, bondi in residue.type.connectivity.items()
                for aj in bondi
            ]
            bonds.extend([[int(atom1), int(atom2)] for atom1, atom2 in molecule.residue_links])
            attrs.append(topologyattrs.Bonds(bonds))
            return Topology(natoms, nres, 1, attrs, resid, None)

        @staticmethod
        def _t2a(residue, atom):
            return residue.name2atom(residue.type.atom2name(atom))

        def _reopen(self):
            self.ts.frame = -1

        def _read_next_timestep(self):
            if self.ts.frame == 0:
                raise EOFError
            self.ts.frame += 1
            return self.ts


    __all__ = ["XpongeMoleculeReader", "mda", "SpongeInputReader", "SpongeNoneReader"]
