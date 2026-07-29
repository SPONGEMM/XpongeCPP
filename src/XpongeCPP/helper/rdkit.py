"""Legacy-compatible RDKit helper interfaces."""

from __future__ import annotations

try:
    from rdkit import Chem
except ModuleNotFoundError as exc:
    raise ModuleNotFoundError(
        "'rdkit' package needed. Maybe you need 'conda install -c rdkit rdkit'"
    ) from exc

from .._core import Assign


def assign_to_rdmol(assign, ignore_bond_type=False):
    """Convert an ``Assign`` object to an RDKit molecule."""
    mol_a = Chem.RWMol()
    for atom in assign.atoms:
        mol_a.AddAtom(Chem.Atom(atom))
    for atom, bonds in enumerate(assign.bonds):
        for aton, order in bonds.items():
            if aton < atom:
                continue
            if ignore_bond_type or order == -1:
                temp_bond = Chem.BondType.UNSPECIFIED
            elif order == 1:
                temp_bond = Chem.BondType.SINGLE
            elif order == 2:
                temp_bond = Chem.BondType.DOUBLE
            elif order == 3:
                temp_bond = Chem.BondType.TRIPLE
            else:
                raise NotImplementedError(f"Unsupported bond order {order}")
            mol_a.AddBond(atom, aton, temp_bond)
    conf = Chem.Conformer(assign.atom_count)
    for i in range(assign.atom_count):
        conf.SetAtomPosition(i, assign.coordinates[i])
    mol = mol_a.GetMol()
    mol.AddConformer(conf)
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetFormalCharge(assign.formal_charges[i])
    flags = Chem.rdmolops.SanitizeFlags
    Chem.SanitizeMol(mol, flags.SANITIZE_ALL & ~flags.SANITIZE_PROPERTIES)
    return mol


def rdmol_to_assign(rdmol):
    """Convert an RDKit molecule to ``Assign`` with upstream-compatible bond handling."""
    assign = Assign()
    kekulized_bonds = {}

    for atom in rdmol.GetAtoms():
        pos = rdmol.GetConformer().GetAtomPosition(atom.GetIdx())
        assign.add_atom(atom.GetSymbol(), pos.x, pos.y, pos.z)
        assign.formal_charges[-1] = atom.GetFormalCharge()

    if any(bond.GetBondType() == Chem.BondType.AROMATIC for bond in rdmol.GetBonds()):
        try:
            kekulized = Chem.Mol(rdmol)
            Chem.Kekulize(kekulized, clearAromaticFlags=True)
        except Exception:
            kekulized_bonds = {}
        else:
            for bond in kekulized.GetBonds():
                key = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                kekulized_bonds[key] = bond.GetBondType()
                kekulized_bonds[(key[1], key[0])] = bond.GetBondType()

    has_unknown_bond = False
    for bond in rdmol.GetBonds():
        temp_bond = bond.GetBondType()
        if temp_bond == Chem.BondType.AROMATIC and kekulized_bonds:
            temp_bond = kekulized_bonds[(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())]
        if temp_bond == Chem.BondType.UNSPECIFIED:
            order = -1
        elif temp_bond == Chem.BondType.SINGLE:
            order = 1
        elif temp_bond == Chem.BondType.DOUBLE:
            order = 2
        elif temp_bond == Chem.BondType.TRIPLE:
            order = 3
        elif temp_bond == Chem.BondType.AROMATIC:
            order = -1
            has_unknown_bond = True
        else:
            raise NotImplementedError(f"Unknown bond type {temp_bond}")
        assign.add_bond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), order)
    assign.determine_ring_and_bond_type()
    if has_unknown_bond:
        assign.determine_bond_order()
    return assign


def insert_atom_type_to_rdmol(mol, res, assign, atom_type_dict=None):
    """Insert residue atom-type information onto an RDKit molecule via isotopes."""
    if atom_type_dict is None:
        atom_type_dict = {}
    for i, atom in enumerate(mol.GetAtoms()):
        atom_type = res.name2atom(assign.names[i]).type.name
        if atom_type not in atom_type_dict:
            atom_type_dict[atom_type] = len(atom_type_dict)
        atom.SetIsotope(atom_type_dict[atom_type])


def find_equal_atoms(assign):
    """Return chemically equivalent atom groups using RDKit canonical SMILES."""
    mols = []
    canon_smiles = []
    for i in range(len(assign.atoms)):
        mols.append(assign_to_rdmol(assign))
        mols[-1].GetAtoms()[i].SetIsotope(1)
        canon_smiles.append(Chem.MolToSmiles(mols[-1], isomericSmiles=True))
    group = {i: i for i in range(len(assign.atoms))}
    for i in range(len(assign.atoms)):
        if group[i] == i:
            for j in range(i + 1, len(assign.atoms)):
                if canon_smiles[j] == canon_smiles[i]:
                    group[j] = i
    ret = []
    realmap = {}
    for i in group:
        if group[i] == i:
            ret.append([i])
            realmap[i] = len(realmap)
        else:
            ret[realmap[group[i]]].append(i)
    return [group_atoms for group_atoms in ret if len(group_atoms) > 1]
