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

if mda is not None:
    from .bundle_mdanalysis import (  # noqa: F401
        BundleTopologyParser,
        SPONGEH5MDReader,
        SpongeH5MDReader,
        load_bundle_universe,
        register_mdanalysis_formats,
        validate_bundle_pair,
    )


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
        format = "SPONGE_MASS"

        @staticmethod
        def _format_hint(thing):
            return isinstance(thing, str) and thing.endswith("_mass.txt")

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


    class SpongeTrajectoryReader(base.ReaderBase):
        """Read a legacy SPONGE float32 ``.dat`` trajectory and its box data."""

        format = "SPONGE_TRAJ"

        @staticmethod
        def _format_hint(thing):
            return isinstance(thing, str) and thing.endswith(".dat")

        def __init__(self, dat_file_name, n_atoms, **kwargs):
            super().__init__(dat_file_name, **kwargs)
            box = kwargs.get("box", None)
            if isinstance(box, str):
                self.boxname = box
                self.box = None
                self._get_box_offset()
            elif box is None:
                raise TypeError(
                    "box should be provided for the sponge trajectory file "
                    f"{dat_file_name}"
                )
            else:
                self.boxname = None
                self.box = box
            self._n_atoms = n_atoms
            self._n_frames = os.path.getsize(dat_file_name) // 12 // self.n_atoms
            self.trajfile = None
            self.boxfile = None
            self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
            self._read_next_timestep()

        @property
        def n_frames(self):
            return self._n_frames

        @property
        def n_atoms(self):
            return self._n_atoms

        @classmethod
        def with_arguments(cls, **kwargs):
            class SpongeTrajectoryReaderWithArguments(cls):
                def __init__(self, dat_file_name, n_atoms, **kwargs_):
                    kwargs_.update(kwargs)
                    super().__init__(dat_file_name, n_atoms, **kwargs_)

            return SpongeTrajectoryReaderWithArguments

        def close(self):
            if self.trajfile is not None:
                self.trajfile.close()
                self.trajfile = None
            if self.boxfile is not None:
                self.boxfile.close()
                self.boxfile = None

        def open_trajectory(self):
            self.trajfile = util.anyopen(self.filename, "rb")
            if self.box is None:
                self.boxfile = util.anyopen(self.boxname)
            self.ts.frame = -1
            return self.trajfile, self.boxfile

        def _reopen(self):
            self.close()
            self.open_trajectory()

        def _read_frame(self, frame):
            if self.trajfile is None:
                self.open_trajectory()
            if self.boxfile is not None:
                self.boxfile.seek(self._offsets[frame])
            self.trajfile.seek(self.n_atoms * 12 * frame)
            self.ts.frame = frame - 1
            return self._read_next_timestep()

        def _read_next_timestep(self):
            ts = self.ts
            if self.trajfile is None:
                self.open_trajectory()
            payload = self.trajfile.read(12 * self.n_atoms)
            if not payload:
                raise EOFError
            ts.positions = np.frombuffer(payload, dtype=np.float32).reshape(
                self.n_atoms, 3
            )
            if self.box is not None:
                ts.dimensions = self.box
            else:
                ts.dimensions = list(map(float, self.boxfile.readline().split()))
            ts.frame += 1
            return ts

        def _get_box_offset(self):
            self._offsets = [0]
            with util.openany(self.boxname) as handle:
                line = handle.readline()
                while line:
                    self._offsets.append(handle.tell())
                    line = handle.readline()
            self._offsets.pop()


    class SpongeTrajectoryWriter:
        """Write a legacy SPONGE float32 ``.dat`` trajectory and ``.box`` sidecar."""

        def __init__(self, filename, write_box=True, **_kwargs):
            self.write_box = write_box
            if not filename.endswith(".dat"):
                raise ValueError(
                    "the name of the SPONGE trajectory file should end with '.dat'"
                )
            self.datname = filename
            self.boxname = filename[::-1].replace(".dat"[::-1], ".box"[::-1], 1)[::-1]
            self.datfile = None
            self.boxfile = None

        def __enter__(self):
            self.open()
            return self

        def __exit__(self, exc_type, exc_val, exc_tb):
            self.close()

        def open(self):
            self.datfile = Xopen(self.datname, "wb")
            if self.write_box:
                self.boxfile = Xopen(self.boxname, "w")

        def close(self):
            if self.datfile is not None:
                self.datfile.close()
                self.datfile = None
            if self.boxfile is not None:
                self.boxfile.close()
                self.boxfile = None

        def write(self, universe):
            if isinstance(universe, mda.Universe):
                ts = universe.coord
            elif isinstance(universe, mda.AtomGroup):
                ts = universe.ts
            else:
                raise TypeError(
                    f"u should be Universe or AtomGroup, but {type(universe)} got"
                )
            self.datfile.write(ts.positions.astype(np.float32).tobytes())
            if (
                self.write_box
                and hasattr(ts, "dimensions")
                and isinstance(ts.dimensions, Iterable)
                and len(ts.dimensions) == 6
            ):
                self.boxfile.write(" ".join(f"{item:.6f}" for item in ts.dimensions) + "\n")


    class SpongeCoordinateReader(base.ReaderBase):
        """Read a legacy single-frame ``*_coordinate.txt`` SPONGE coordinate file."""

        format = "SPONGE_CRD"

        @staticmethod
        def _format_hint(thing):
            return isinstance(thing, str) and thing.endswith("_coordinate.txt")

        def __init__(self, file_name, n_atoms, **kwargs):
            super().__init__(file_name, **kwargs)
            self._n_atoms = n_atoms
            self._n_frames = 1
            self.file = None
            self.start = 0
            self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
            self._read_next_timestep()

        @property
        def n_frames(self):
            return self._n_frames

        @property
        def n_atoms(self):
            return self._n_atoms

        def close(self):
            if self.file is not None:
                self.file.close()
                self.file = None

        def open_file(self):
            self.file = util.anyopen(self.filename, "r")
            self.file.readline()
            self.start = self.file.tell()

        def _reopen(self):
            self.close()
            self.open_file()

        def _read_frame(self, frame):
            if self.file is None:
                self.open_file()
            self.file.seek(self.start)
            self.ts.frame = frame - 1
            return self._read_next_timestep()

        def _read_next_timestep(self):
            ts = self.ts
            if self.file is None:
                self.open_file()
            if self.file.tell() != self.start:
                raise EOFError
            ts.positions = np.loadtxt(self.file, max_rows=self.n_atoms)
            ts.dimensions = list(map(float, self.file.readline().split()))
            ts.frame += 1
            return ts


    class SpongeCoordinateWriter:
        """Write a legacy single-frame ``*_coordinate.txt`` SPONGE coordinate file."""

        def __init__(self, file_name, n_atoms=None):
            self.filename = file_name
            self.file = None
            self.n_atoms = n_atoms

        def __enter__(self):
            self.open()
            return self

        def __exit__(self, exc_type, exc_val, exc_tb):
            self.close()

        def open(self):
            self.file = Xopen(self.filename, "w")

        def close(self):
            if self.file is not None:
                self.file.close()
                self.file = None

        def write(self, universe):
            if isinstance(universe, mda.Universe):
                ts = universe.coord
            elif isinstance(universe, mda.AtomGroup):
                ts = universe.ts
            else:
                raise TypeError(
                    f"u should be Universe or AtomGroup, but {type(universe)} got"
                )
            if self.n_atoms is None:
                self.n_atoms = len(ts.positions)
            lines = [str(self.n_atoms)]
            lines.extend(
                f"{coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}"
                for coord in ts.positions[: self.n_atoms]
            )
            if hasattr(ts, "dimensions"):
                lines.append(" ".join(f"{item:.6f}" for item in ts.dimensions))
            else:
                lines.append("999 999 999 90 90 90")
            self.file.write("\n".join(lines) + "\n")


    mda._SINGLEFRAME_WRITERS["SPONGE_CRD"] = SpongeCoordinateWriter
    mda._SINGLEFRAME_WRITERS["TXT"] = SpongeCoordinateWriter
    mda._MULTIFRAME_WRITERS["SPONGE_TRAJ"] = SpongeTrajectoryWriter
    mda._MULTIFRAME_WRITERS["DAT"] = SpongeTrajectoryWriter


    __all__ = [
        "BundleTopologyParser",
        "SPONGEH5MDReader",
        "SpongeH5MDReader",
        "SpongeInputReader",
        "SpongeNoneReader",
        "SpongeTrajectoryReader",
        "SpongeTrajectoryWriter",
        "SpongeCoordinateReader",
        "SpongeCoordinateWriter",
        "XpongeMoleculeReader",
        "load_bundle_universe",
        "mda",
        "register_mdanalysis_formats",
        "validate_bundle_pair",
    ]
