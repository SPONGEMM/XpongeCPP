"""Centralized legacy Assign compatibility patches."""

from __future__ import annotations

from collections import OrderedDict
from importlib.resources import files

from .._core import Assign, ResidueType
from ..assign import AssignRule

_core_determine_connectivity = Assign.determine_connectivity
_core_determine_bond_order = Assign.determine_bond_order
_core_determine_atom_type = Assign.determine_atom_type
_core_save_as_mol2 = Assign.save_as_mol2
_core_set_atom_type = Assign.set_atom_type
_core_to_residuetype = Assign.to_residuetype

_CONNECTIVITY_RADII = {
    "H": 0.35,
    "C": 0.73,
    "N": 0.66,
    "O": 0.69,
    "F": 0.68,
    "P": 1.06,
    "S": 1.02,
    "Cl": 0.99,
    "Br": 1.14,
    "I": 1.33,
}


def _tpacm4_tables():
    base = files("XpongeCPP.data.assign.tpacm4")
    return (base.joinpath("ATOMTYPE.dat").read_text(), base.joinpath("CHARGE.dat").read_text())


def _assignment_to_rdkit_mol(assignment):
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise ImportError("RDKit is required for Gasteiger charge calculation") from exc

    mol = Chem.RWMol()
    for atom_index, element in enumerate(assignment.atoms):
        atom = Chem.Atom(element)
        formal_charges = getattr(assignment, "formal_charges", [])
        if atom_index < len(formal_charges):
            atom.SetFormalCharge(int(formal_charges[atom_index]))
        mol.AddAtom(atom)
    bond_type = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
    for atom1, bonded_atoms in enumerate(assignment.bonds):
        for atom2, order in bonded_atoms.items():
            if atom1 < atom2:
                mol.AddBond(atom1, int(atom2), bond_type.get(int(order), Chem.BondType.UNSPECIFIED))
    rdmol = mol.GetMol()
    conf = Chem.Conformer(assignment.atom_count)
    for atom in range(assignment.atom_count):
        x, y, z = assignment.coordinates[atom]
        conf.SetAtomPosition(atom, (x, y, z))
    rdmol.AddConformer(conf)
    flags = Chem.rdmolops.SanitizeFlags
    Chem.SanitizeMol(rdmol, flags.SANITIZE_ALL & ~flags.SANITIZE_PROPERTIES)
    return rdmol


def _assign_calculate_charge(self, method, **parameters):
    method = method.upper()
    if method == "TPACM4":
        atom_type_table, charge_table = _tpacm4_tables()
        charge = parameters.get("charge", int(round(sum(self.formal_charges))))
        self._calculate_tpacm4(atom_type_table, charge_table, int(charge))
        return None
    if method == "GASTEIGER":
        try:
            from rdkit.Chem import rdPartialCharges
        except ImportError as exc:
            raise ImportError("RDKit is required for Gasteiger charge calculation") from exc
        rdmol = _assignment_to_rdkit_mol(self)
        rdPartialCharges.ComputeGasteigerCharges(rdmol)
        self.set_charges([float(atom.GetProp("_GasteigerCharge")) for atom in rdmol.GetAtoms()])
        return None
    if method == "RESP":
        from ..assign import resp

        self.set_charges(
            resp.resp_fit(
                self,
                basis=parameters.get("basis", "6-31g*"),
                opt=parameters.get("opt", False),
                charge=parameters.get("charge", int(round(sum(self.formal_charges)))),
                spin=parameters.get("spin", 0),
                extra_equivalence=parameters.get("extra_equivalence", []),
                grid_density=parameters.get("grid_density", 6),
                grid_cell_layer=parameters.get("grid_cell_layer", 4),
                a1=parameters.get("a1", 0.0005),
                a2=parameters.get("a2", 0.001),
                two_stage=parameters.get("two_stage", True),
                only_esp=parameters.get("only_esp", False),
                radius=parameters.get("radius", None),
            )
        )
        return None
    raise ValueError("methods should be one of the following: 'RESP', 'GASTEIGER', 'TPACM4' (case-insensitive)")


def _first_neighbor(assignment, atom):
    return next(iter(assignment.bonds[atom]))


def _is_carboxylic_acid_hydrogen(assignment, atom):
    if assignment.atoms[atom] != "H":
        return False
    oxygen = _first_neighbor(assignment, atom)
    if assignment.atoms[oxygen] != "O":
        return False
    carbon = None
    for neighbor in assignment.bonds[oxygen]:
        if neighbor != atom:
            carbon = neighbor
            break
    if carbon is None or not assignment.atoms[carbon] == "C" or len(assignment.bonds[carbon]) != 3:
        return False
    for neighbor, order in assignment.bonds[carbon].items():
        if neighbor != oxygen and assignment.atoms[neighbor] == "O" and int(order) == 2:
            return True
    return False


def _is_carboxylate_oxygen(assignment, atom):
    if assignment.atoms[atom] != "O" or assignment.formal_charges[atom] != -1:
        return False
    carbon = _first_neighbor(assignment, atom)
    for neighbor, order in assignment.bonds[carbon].items():
        if neighbor != atom and assignment.atoms[neighbor] == "O" and int(order) == 2:
            return True
    return False


def _is_phenol_hydrogen(assignment, atom):
    if assignment.atoms[atom] != "H":
        return False
    oxygen = _first_neighbor(assignment, atom)
    if assignment.atoms[oxygen] != "O":
        return False
    for neighbor in assignment.bonds[oxygen]:
        if neighbor != atom and assignment.has_atom_marker(neighbor, "AR0"):
            return True
    return False


def _is_phenolate_oxygen(assignment, atom):
    if assignment.atoms[atom] != "O" or assignment.formal_charges[atom] != -1:
        return False
    return assignment.has_atom_marker(_first_neighbor(assignment, atom), "AR0")


def _hydrogen_position_for(assignment, atom):
    x = y = z = 0.0
    coord = assignment.coordinates[atom]
    count = 0
    for neighbor in assignment.bonds[atom]:
        ncoord = assignment.coordinates[neighbor]
        dx = coord[0] - ncoord[0]
        dy = coord[1] - ncoord[1]
        dz = coord[2] - ncoord[2]
        norm = (dx * dx + dy * dy + dz * dz) ** 0.5 or 1.0
        x += coord[0] + dx / norm
        y += coord[1] + dy / norm
        z += coord[2] + dz / norm
        count += 1
    if count == 0:
        return coord
    return x / count, y / count, z / count


def _assign_set_ph(self, ph):
    self.kekulize()
    to_delete = []
    to_add = []
    for atom in range(self.atom_count):
        if _is_carboxylic_acid_hydrogen(self, atom) and 4.0 < ph:
            to_delete.append(atom)
        elif _is_phenol_hydrogen(self, atom) and 10.0 < ph:
            to_delete.append(atom)
        elif self.atoms[atom] == "H" and self.atoms[_first_neighbor(self, atom)] == "O" and 15.9 < ph:
            to_delete.append(atom)
        elif _is_carboxylate_oxygen(self, atom) and 4.0 > ph:
            to_add.append(atom)
        elif _is_phenolate_oxygen(self, atom) and 10.0 > ph:
            to_add.append(atom)
        elif self.atoms[atom] == "O" and self.formal_charges[atom] == -1 and 15.9 > ph:
            to_add.append(atom)

    for atom in to_add:
        x, y, z = _hydrogen_position_for(self, atom)
        self.add_atom("H", x, y, z)
        self.add_bond(self.atom_count - 1, atom, 1)
    for atom in sorted(set(to_delete), reverse=True):
        self.delete_atom(atom)
    self.determine_bond_order(True, None)
    return int(round(sum(self.formal_charges)))


def _assign_set_atom_type_compat(self, atom, atom_type):
    if hasattr(atom_type, "name"):
        atom_type = atom_type.name
    return _core_set_atom_type(self, atom, str(atom_type))


def _assign_set_atom_types_compat(self, atom_types):
    atom_types = list(atom_types)
    if len(atom_types) != self.atom_count:
        raise ValueError("atom_types length should match assignment atom_count")
    for atom_index, atom_type in enumerate(atom_types):
        _assign_set_atom_type_compat(self, atom_index, atom_type)
    return None


def _assign_to_residuetype_compat(self, name, charge=None):
    try:
        return _core_to_residuetype(self, name)
    except ValueError as exc:
        if "duplicate atom name" not in str(exc):
            raise

    residue_type = ResidueType(name)
    if charge is None:
        charges = list(getattr(self, "charges", []) or [])
        if not charges:
            charges = [0.0] * self.atom_count
    else:
        charges = list(charge)

    counts = {}
    names = []
    for index, element in enumerate(self.atoms):
        atom_name = self.names[index] or element
        count = counts.get(atom_name, 0)
        counts[atom_name] = count + 1
        if count:
            atom_name = f"{atom_name}{count}"
        names.append(atom_name)

    atom_types = list(self.atom_types)
    for index, atom_name in enumerate(names):
        residue_type.add_atom(
            atom_name,
            atom_types[index],
            self.coordinates[index][0],
            self.coordinates[index][1],
            self.coordinates[index][2],
        )
        residue_type.atoms[-1].charge = charges[index]

    for atom1, bonded in enumerate(self.bonds):
        for atom2 in bonded:
            if atom1 < int(atom2):
                residue_type.add_connectivity(names[atom1], names[int(atom2)])
    return residue_type


def _distance(coord1, coord2):
    dx = coord1[0] - coord2[0]
    dy = coord1[1] - coord2[1]
    dz = coord1[2] - coord2[2]
    return (dx * dx + dy * dy + dz * dz) ** 0.5


def _assign_determine_connectivity(self, simple_cutoff=None, tolerance=1.0):
    if simple_cutoff is not None:
        cutoff = float(simple_cutoff)
        for i, coord_i in enumerate(self.coordinates):
            for j in range(i + 1, self.atom_count):
                if _distance(coord_i, self.coordinates[j]) < cutoff:
                    self.add_bond(i, j, -1)
        return None
    for i, coord_i in enumerate(self.coordinates):
        for j in range(i + 1, self.atom_count):
            distance = _distance(coord_i, self.coordinates[j])
            if distance == 0:
                continue
            rij = _CONNECTIVITY_RADII.get(self.atoms[i], 1.25) + _CONNECTIVITY_RADII.get(self.atoms[j], 1.25)
            if distance <= 1.5:
                factor = 1 - 0.15
            elif distance <= 1.9:
                factor = 1 - 0.11
            elif distance <= 2.05:
                factor = 1 - 0.09
            else:
                factor = 1 - 0.08
            factor /= tolerance
            ratio = rij / distance
            if factor < ratio < 2:
                self.add_bond(i, j, -1)
    return None


def _normalise_penalty_scores(assign, penalty_scores):
    if penalty_scores is None:
        return []
    if len(penalty_scores) != assign.atom_count:
        raise ValueError("penalty_scores should have one entry per atom")
    result = []
    for entry in penalty_scores:
        if isinstance(entry, dict):
            items = list(entry.items())
        else:
            items = list(entry)
        if not items:
            raise ValueError("penalty_scores entries should not be empty")
        result.append([(int(valence), int(penalty)) for valence, penalty in items])
    return result


def _assign_determine_bond_order(self, max_step=2000, max_stat=20000, penalty_scores=None,
                                 check_formal_charge=True, total_charge=None, extra_criteria=None):
    if isinstance(max_step, bool):
        check_formal_charge = max_step
        total_charge = max_stat
        max_step = 2000
        max_stat = 20000
    if penalty_scores is None and extra_criteria is None and max_step == 2000 and max_stat == 20000:
        return _core_determine_bond_order(self, check_formal_charge, total_charge)
    penalties = _normalise_penalty_scores(self, penalty_scores)
    return self._determine_bond_order_custom(
        check_formal_charge,
        total_charge,
        max_step,
        max_stat,
        penalties,
        extra_criteria,
    )


def _phmodel_type(self, atom):
    if _is_phenolate_oxygen(self, atom):
        return "B-phenol"
    if _is_phenol_hydrogen(self, atom):
        return "A-phenol"
    if _is_carboxylate_oxygen(self, atom):
        return "B-carboxylic"
    if _is_carboxylic_acid_hydrogen(self, atom):
        return "A-carboxylic"
    if self.atoms[atom] == "O" and self.formal_charges[atom] == -1:
        return "B-alcohol"
    if self.atoms[atom] == "H" and self.bonds[atom] and self.atoms[_first_neighbor(self, atom)] == "O":
        return "A-alcohol"
    return "N"


def _assign_determine_atom_type_from_rule(assign, rule):
    if not rule.built:
        rule.rules = OrderedDict(sorted(rule.rules.items(), key=lambda item: rule.priority[item[0]]))
        rule.built = True
    backup = list(assign.atom_types)
    if rule.pre_action:
        rule.pre_action(assign)
    assigned_types = []
    for atom in range(assign.atom_count):
        for atom_type, judge in rule.rules.items():
            if judge(atom, assign):
                assigned_types.append(atom_type)
                assign.set_atom_type(atom, atom_type)
                break
        else:
            raise KeyError(f"No atom type found for assignment {assign.name} of atom #{atom}")
    if rule.post_action:
        rule.post_action(assign)
    if rule.pure_string:
        for atom, atom_type in enumerate(backup):
            assign.set_atom_type(atom, atom_type)
        return assigned_types
    return None


def _assign_determine_atom_type(self, rule):
    if str(rule).lower() == "phmodel":
        self.kekulize()
        return [_phmodel_type(self, atom) for atom in range(self.atom_count)]
    if isinstance(rule, AssignRule):
        return _assign_determine_atom_type_from_rule(self, rule)
    if str(rule).lower() not in {"gaff", "gaff2", "sybyl"} and str(rule) in AssignRule.all:
        return _assign_determine_atom_type_from_rule(self, AssignRule.all[str(rule)])
    return _core_determine_atom_type(self, rule)


def _assign_determine_equal_atoms(self):
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise ImportError("RDKit is required for determine_equal_atoms") from exc
    smiles = []
    for atom in range(self.atom_count):
        mol = _assignment_to_rdkit_mol(self)
        mol.GetAtoms()[atom].SetIsotope(1)
        smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True))
    group = {i: i for i in range(self.atom_count)}
    for i in range(self.atom_count):
        if group[i] == i:
            for j in range(i + 1, self.atom_count):
                if smiles[j] == smiles[i]:
                    group[j] = i
    out = []
    realmap = {}
    for atom in range(self.atom_count):
        if group[atom] == atom:
            realmap[atom] = len(out)
            out.append([atom])
        else:
            out[realmap[group[atom]]].append(atom)
    return [atoms for atoms in out if len(atoms) > 1]


def _assign_uff_optimize(self):
    try:
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise ImportError("RDKit is required for uff_optimize") from exc
    mol = _assignment_to_rdkit_mol(self)
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    for atom in range(self.atom_count):
        pos = conf.GetAtomPosition(atom)
        self.set_coordinate(atom, pos.x, pos.y, pos.z)


def _assign_uff_energy(self):
    try:
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise ImportError("RDKit is required for uff_energy") from exc
    mol = _assignment_to_rdkit_mol(self)
    ff = AllChem.UFFGetMoleculeForceField(mol)
    return ff.CalcEnergy()


def _assign_save_as_mol2(self, filename, residue_name="MOL", atomtype="sybyl"):
    if atomtype and str(atomtype).lower() != "sybyl":
        raise ValueError(f"save_as_mol2 only supports atomtype='sybyl'; got {atomtype!r}")
    if atomtype:
        self.determine_atom_type("sybyl")
    return _core_save_as_mol2(self, filename, residue_name)


def install_legacy_assign_patches():
    Assign.calculate_charge = _assign_calculate_charge
    Assign.Calculate_Charge = _assign_calculate_charge
    Assign.set_ph = _assign_set_ph
    Assign.Set_PH = _assign_set_ph
    Assign.set_atom_type = _assign_set_atom_type_compat
    Assign.to_residuetype = _assign_to_residuetype_compat
    Assign.determine_connectivity = _assign_determine_connectivity
    Assign.Determine_Connectivity = _assign_determine_connectivity
    Assign.determine_bond_order = _assign_determine_bond_order
    Assign.Determine_Bond_Order = _assign_determine_bond_order
    Assign.determine_atom_type = _assign_determine_atom_type
    Assign.Determine_Atom_Type = _assign_determine_atom_type
    Assign.determine_equal_atoms = _assign_determine_equal_atoms
    Assign.Determine_Equal_Atoms = _assign_determine_equal_atoms
    Assign.Determine_Ring_And_Bond_Type = Assign.determine_ring_and_bond_type
    Assign.uff_optimize = _assign_uff_optimize
    Assign.UFF_Optimize = _assign_uff_optimize
    Assign.uff_energy = property(_assign_uff_energy)
    Assign.save_as_mol2 = _assign_save_as_mol2
    Assign.Save_As_Mol2 = _assign_save_as_mol2
    Assign.Save_As_PDB = Assign.save_as_pdb
    Assign.Set_Charge = Assign.set_charge
    Assign.Set_Charges = Assign.set_charges
    Assign.Set_Coordinate = Assign.set_coordinate
    Assign.Set_Formal_Charge = Assign.set_formal_charge
    Assign.Set_Atom_Type = Assign.set_atom_type
    Assign.To_ResidueType = Assign.to_residuetype
    Assign.Add_Bond_Marker = Assign.add_bond_marker
    Assign.Has_Bond_Marker = Assign.has_bond_marker
    Assign.atom_types = property(Assign.atom_types.fget, _assign_set_atom_types_compat)


__all__ = ["install_legacy_assign_patches"]
