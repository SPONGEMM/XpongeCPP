"""FEP softcore helpers compatible with the common Xponge workflow."""


def prepare_lj_soft_core(molecule, lj_type_b_by_atom, subsys=1):
    """Set B-state LJ types and subsystem ids, then enable softcore export."""
    for atom_index, lj_type_b in lj_type_b_by_atom.items():
        atom = _atom_by_index(molecule, atom_index)
        atom.lj_type_b = lj_type_b
        atom.subsys = subsys
    molecule.enable_lj_soft_core()
    return molecule


def set_lj_type_b(molecule, atom_index, lj_type_b):
    atom = _atom_by_index(molecule, atom_index)
    atom.lj_type_b = lj_type_b
    return molecule


def set_subsys(molecule, atom_index, subsys):
    atom = _atom_by_index(molecule, atom_index)
    atom.subsys = int(subsys)
    molecule.enable_subsys_division()
    return molecule


def enable_lj_soft_core(molecule):
    molecule.enable_lj_soft_core()
    return molecule


def add_soft_bond_from_a(molecule, atom1, atom2, k, b):
    molecule.add_bond_soft(atom1, atom2, k, b, 0)
    return molecule


def add_soft_bond_from_b(molecule, atom1, atom2, k, b):
    molecule.add_bond_soft(atom1, atom2, k, b, 1)
    return molecule


def _atom_by_index(molecule, atom_index):
    target = int(atom_index)
    if target < 0 or target >= molecule.atom_count:
        raise IndexError("atom index out of range")
    seen = 0
    for residue in molecule.residues:
        for atom in residue.atoms:
            if seen == target:
                return atom
            seen += 1
    raise IndexError("atom index out of range")


Prepare_LJ_Soft_Core = prepare_lj_soft_core
Set_LJ_Type_B = set_lj_type_b
Set_Subsys = set_subsys
Enable_LJ_Soft_Core = enable_lj_soft_core
Add_Soft_Bond_From_A = add_soft_bond_from_a
Add_Soft_Bond_From_B = add_soft_bond_from_b
