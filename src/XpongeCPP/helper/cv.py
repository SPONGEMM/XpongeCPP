"""Minimal legacy-compatible CV helpers for XpongeCPP.

This first-wave compatibility module focuses on the parts of ``Xponge.helper.cv``
that are used by core legacy workflows to generate ``cv.txt`` style SPONGE input.
It intentionally avoids the heavier MDAnalysis dependency from legacy Xponge.
"""

from pathlib import Path
from collections.abc import Iterable

import numpy as np

from .. import Atom
from .._compat.imports import Xdict


class _CVVirtualAtom:
    def __str__(self):
        return self.name

    def to_string(self, folder):
        raise NotImplementedError


class _CV:
    def __str__(self):
        return self.name

    def to_string(self, folder):
        raise NotImplementedError


class _CVBias:
    def to_string(self, folder):
        raise NotImplementedError


class _Center(_CVVirtualAtom):
    def __init__(self, name, atom, weight):
        self.name = name
        self.atom = atom
        self.weight = weight

    def to_string(self, folder):
        prefix = Path(folder)
        if len(self.atom) > 10:
            with open(prefix / f"{self.name}_atom.txt", "w", encoding="utf-8") as f:
                f.write("\n".join(map(str, self.atom)))
            with open(prefix / f"{self.name}_weight.txt", "w", encoding="utf-8") as f:
                f.write("\n".join(map(str, self.weight)))
            atom_line = f"atom_in_file = {prefix / (self.name + '_atom.txt')}"
            weight_line = f"weight_in_file = {prefix / (self.name + '_weight.txt')}"
        else:
            atom_line = "atom = " + " ".join(map(str, self.atom))
            weight_line = "weight = " + " ".join(map(str, self.weight))
        return f"""{self.name}
{{
    vatom_type = center
    {atom_line}
    {weight_line}
}}
"""


class _COM(_CVVirtualAtom):
    def __init__(self, name, atom):
        self.name = name
        self.atom = atom

    def to_string(self, folder):
        prefix = Path(folder)
        if len(self.atom) > 10:
            with open(prefix / f"{self.name}_atom.txt", "w", encoding="utf-8") as f:
                f.write("\n".join(map(str, self.atom)))
            atom_line = f"atom_in_file = {prefix / (self.name + '_atom.txt')}"
        else:
            atom_line = "atom = " + " ".join(map(str, self.atom))
        return f"""{self.name}
{{
    vatom_type = center_of_mass
    {atom_line}
}}
"""


class _Position(_CV):
    def __init__(self, name, atom, xyz, scaled):
        self.name = name
        self.atom = atom
        self.xyz = xyz
        if not scaled:
            self.scaled = ""
        else:
            self.scaled = "scaled_"
            self.period = 1

    def to_string(self, folder):
        del folder
        return f"""{self.name}
{{
    CV_type = {self.scaled}position_{self.xyz}
    atom = {self.atom}
}}
"""


class _Dihedral(_CV):
    period = 2 * np.pi

    def __init__(self, name, atom1, atom2, atom3, atom4):
        self.name = name
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

    def to_string(self, folder):
        del folder
        return f"""{self.name}
{{
    CV_type = dihedral
    atom = {self.atom1} {self.atom2} {self.atom3} {self.atom4}
}}
"""


class _RMSD(_CV):
    def __init__(self, name, atom, coordinate):
        self.name = name
        self.atom = atom
        self.coordinate = coordinate

    def to_string(self, folder):
        prefix = Path(folder)
        with open(prefix / f"{self.name}_atom.txt", "w", encoding="utf-8") as f:
            f.write("\n".join(map(str, self.atom)))
        np.savetxt(prefix / f"{self.name}_coordinate.txt", self.coordinate)
        return f"""{self.name}
{{
    CV_type = rmsd
    atom_in_file = {prefix / (self.name + '_atom.txt')}
    coordinate_in_file = {prefix / (self.name + '_coordinate.txt')}
}}
"""


class _PrintCV(_CVBias):
    def __init__(self, cv):
        self.cv = list(cv) if isinstance(cv, Iterable) and not isinstance(cv, (str, _CV)) else [cv]

    def to_string(self, folder):
        del folder
        return f"""print
{{
    CV = {" ".join(str(cv) for cv in self.cv)}
}}
"""


class _SteerCV(_CVBias):
    def __init__(self, cv, weight):
        if isinstance(cv, Iterable) and not isinstance(cv, (str, _CV)):
            self.cv = list(cv)
            self.weight = list(weight)
        else:
            self.cv = [cv]
            self.weight = [weight]

    def to_string(self, folder):
        del folder
        return f"""steer
{{
    CV = {" ".join(str(cv) for cv in self.cv)}
    weight = {" ".join(str(weight) for weight in self.weight)}
}}
"""


class _RestrainCV(_CVBias):
    def __init__(self, cv, weight, reference):
        if isinstance(cv, Iterable) and not isinstance(cv, (str, _CV)):
            self.cv = list(cv)
            self.weight = list(weight)
            self.reference = list(reference)
        else:
            self.cv = [cv]
            self.weight = [weight]
            self.reference = [reference]
        self.start_step = []
        self.max_step = []
        self.reduce_step = []
        self.stop_step = []

    def to_string(self, folder):
        del folder
        need_period = any(hasattr(cv, "period") for cv in self.cv)
        period_line = ""
        if need_period:
            period_line = "    period = " + " ".join(str(getattr(cv, "period", 0)) for cv in self.cv) + "\n"
        return f"""restrain
{{
    CV = {" ".join(str(cv) for cv in self.cv)}
    weight = {" ".join(str(weight) for weight in self.weight)}
    reference = {" ".join(str(reference) for reference in self.reference)}
    start_step = {" ".join(str(v) for v in self.start_step)}
    max_step = {" ".join(str(v) for v in self.max_step)}
    reduce_step = {" ".join(str(v) for v in self.reduce_step)}
    stop_step = {" ".join(str(v) for v in self.stop_step)}
{period_line}}}
"""


class _Meta1D(_CVBias):
    def __init__(self, cv, **kwargs):
        self.cv = cv
        self.kwargs = kwargs

    def to_string(self, folder):
        del folder
        lines = [f"    {key} = {value}\n" for key, value in self.kwargs.items()]
        return f"""metad
{{
    CV = {self.cv}
{"".join(lines)}}}
"""


class CVSystem:
    """Minimal compatibility version of ``Xponge.helper.cv.CVSystem``."""

    __slots__ = ["molecule", "virtual_atom", "cv", "bias", "names", "_association"]
    BIASES = {"print", "restrain", "steer", "meta1d"}
    SPONGE_NAMES = {"bond", "angle", "dihedral"}
    _SOLVENT_RESIDUES = {"WAT", "H2O", "HOH", "SOL"}
    _BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}

    def __init__(self, molecule):
        self._association = Xdict(not_found_message="No name {} found in the system")
        self.molecule = molecule
        self.virtual_atom = Xdict(not_found_message="No virtual atom named {} found in the system")
        self.cv = Xdict(not_found_message="No collected variable named {} found in the system")
        self.bias = Xdict(not_found_message="No bias named {} found in the system")
        self.names = Xdict(not_found_message="No name {} found in the system")
        self.names.fromkeys(self.BIASES, [])
        self.names.fromkeys(self.SPONGE_NAMES, [])

    def _all_atoms(self):
        return list(self.molecule.atoms)

    def _atom_index(self, atom):
        return int(self.molecule.atom_index[atom])

    def _protein_atoms(self):
        atoms = []
        for residue in self.molecule.residues:
            if residue.name not in self._SOLVENT_RESIDUES:
                atoms.extend(residue.atoms)
        return atoms

    def _select_atoms(self, select):
        select = " ".join(str(select).lower().split())
        if select == "protein":
            return self._protein_atoms()
        if select == "protein and backbone":
            return [atom for atom in self._protein_atoms() if atom.name.upper() in self._BACKBONE_ATOMS]
        raise NotImplementedError(f"Unsupported legacy CV selection: {select!r}")

    def get_atom_index(self, atom):
        if atom in self.virtual_atom:
            return atom
        return self._atom_index(atom)

    def add_center(self, name, select, weight=None):
        if name in self.names:
            raise ValueError(f"{name} has been defined in the name system")
        atoms = [self._atom_index(atom) for atom in self._select_atoms(select)]
        if weight is None:
            weight = [1.0 / len(atoms)] * len(atoms)
            self.virtual_atom[name] = _Center(name, atoms, weight)
        elif weight == "mass":
            self.virtual_atom[name] = _COM(name, atoms)
        else:
            self.virtual_atom[name] = _Center(name, atoms, weight)
        self.names[name] = self.virtual_atom[name]
        self._association[name] = []

    def add_cv_position(self, name, atom, xyz, scaled):
        if name in self.names:
            raise ValueError(f"{name} has been defined in the name system")
        if xyz not in ("x", "y", "z"):
            raise ValueError(f"xyz should be one of 'x', 'y' or 'z', but {xyz} is given")
        self.cv[name] = _Position(name, self.get_atom_index(atom), xyz, scaled)
        self.names[name] = self.cv[name]
        self._association[name] = []
        if atom in self.virtual_atom:
            self._association[atom].append(self.names[name])
        elif not isinstance(atom, Atom):
            raise TypeError(f"atom should be an Xponge.Atom or a name of virtual atom, but {atom} is given")

    def add_cv_dihedral(self, name, atom1, atom2, atom3, atom4):
        if name in self.names:
            raise ValueError(f"{name} has been defined in the name system")
        self.cv[name] = _Dihedral(
            name,
            self.get_atom_index(atom1),
            self.get_atom_index(atom2),
            self.get_atom_index(atom3),
            self.get_atom_index(atom4),
        )
        self.names[name] = self.cv[name]
        self._association[name] = []
        for atom in [atom1, atom2, atom3, atom4]:
            if atom in self.virtual_atom:
                self._association[atom].append(self.names[name])
            elif not isinstance(atom, Atom):
                raise TypeError(f"atom should be an Xponge.Atom or a name of virtual atom, but {atom} is given")

    def add_cv_rmsd(self, name, select):
        if name in self.names:
            raise ValueError(f"{name} has been defined in the name system")
        atoms = self._select_atoms(select)
        atom_indices = [self._atom_index(atom) for atom in atoms]
        coordinates = np.array([[atom.x, atom.y, atom.z] for atom in atoms], dtype=float)
        self.cv[name] = _RMSD(name, atom_indices, coordinates)
        self.names[name] = self.cv[name]
        self._association[name] = []

    def print(self, name):
        if name not in self.cv:
            raise ValueError(f"{name} is not a valid CV")
        if "print" not in self.bias:
            self.bias["print"] = _PrintCV(self.cv[name])
            self.names["print"] = self.bias["print"]
        else:
            self.bias["print"].cv.append(self.cv[name])
        self._association[name].append(self.bias["print"])

    def steer(self, name, weight):
        if name not in self.cv:
            raise ValueError(f"{name} is not a valid CV")
        if "steer" not in self.bias:
            self.bias["steer"] = _SteerCV(self.cv[name], weight)
            self.names["steer"] = self.bias["steer"]
        else:
            self.bias["steer"].cv.append(self.cv[name])
            self.bias["steer"].weight.append(weight)

    def restrain(self, name, weight, reference, start_step=0, max_step=0, reduce_step=0, stop_step=0):
        if name not in self.cv:
            raise ValueError(f"{name} is not a valid CV")
        if "restrain" not in self.bias:
            self.bias["restrain"] = _RestrainCV(self.cv[name], weight, reference)
            self.names["restrain"] = self.bias["restrain"]
        else:
            self.bias["restrain"].cv.append(self.cv[name])
            self.bias["restrain"].weight.append(weight)
            self.bias["restrain"].reference.append(reference)
        self.bias["restrain"].start_step.append(start_step)
        self.bias["restrain"].max_step.append(max_step)
        self.bias["restrain"].reduce_step.append(reduce_step)
        self.bias["restrain"].stop_step.append(stop_step)

    def meta1d(self, name, CV_grid, CV_minimal, CV_maximum, CV_sigma, welltemp_factor, height):
        if name not in self.cv:
            raise ValueError(f"{name} is not a valid CV")
        if "metad" not in self.bias:
            kwargs = {
                "CV_grid": CV_grid,
                "CV_minimal": CV_minimal,
                "CV_maximum": CV_maximum,
                "welltemp_factor": welltemp_factor,
                "height": height,
                "CV_sigma": CV_sigma,
            }
            if hasattr(self.cv[name], "period"):
                kwargs["CV_period"] = self.cv[name].period
            self.bias["metad"] = _Meta1D(name, **kwargs)
            self.names["metad"] = self.bias["metad"]

    def output(self, filename, folder="."):
        with open(filename, "w", encoding="utf-8") as f:
            if self.virtual_atom:
                f.write("#" * 30 + "\n#definition of virtual atoms\n" + "#" * 30 + "\n")
                for virtual_atom in self.virtual_atom.values():
                    f.write(virtual_atom.to_string(folder))
            if self.cv:
                f.write("#" * 30 + "\n#definition of collected variables\n" + "#" * 30 + "\n")
                for cv in self.cv.values():
                    f.write(cv.to_string(folder))
            if self.bias:
                f.write("#" * 30 + "\n#definition of bias\n" + "#" * 30 + "\n")
                for bias in self.bias.values():
                    f.write(bias.to_string(folder))


CVSystem.Add_Center = CVSystem.add_center
CVSystem.AddCenter = CVSystem.add_center
CVSystem.addCenter = CVSystem.add_center
CVSystem.Add_CV_Position = CVSystem.add_cv_position
CVSystem.AddCVPosition = CVSystem.add_cv_position
CVSystem.addCvPosition = CVSystem.add_cv_position
CVSystem.Add_CV_Dihedral = CVSystem.add_cv_dihedral
CVSystem.AddCVDihedral = CVSystem.add_cv_dihedral
CVSystem.addCvDihedral = CVSystem.add_cv_dihedral
CVSystem.Add_CV_RMSD = CVSystem.add_cv_rmsd
CVSystem.AddCVRMSD = CVSystem.add_cv_rmsd
CVSystem.addCvRMSD = CVSystem.add_cv_rmsd
CVSystem.Print = CVSystem.print
CVSystem.Steer = CVSystem.steer
CVSystem.Restrain = CVSystem.restrain
CVSystem.Meta1D = CVSystem.meta1d
CVSystem.Meta_1D = CVSystem.meta1d
CVSystem.Output = CVSystem.output


__all__ = ["CVSystem"]
