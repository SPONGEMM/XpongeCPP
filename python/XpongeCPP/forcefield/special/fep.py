"""FEP softcore helpers compatible with the common Xponge workflow."""

from ... import (
    merge_dual_topology as _merge_dual_topology,
    merge_force_field as _merge_force_field,
    register_amber_lj_parameter,
)

register_amber_lj_parameter("ZERO_LJ_ATOM", "ZERO_LJ_ATOM", 0.0, 0.0)


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


def merge_dual_topology(molecule, residue, residue_b, assign_a=None, assign_b=None, tmcs=60,
                        image_path=None, similarity_limit=0, imcs=None):
    """Build Xponge-style A/B dual topology using an explicit or RDKit-derived atom map."""
    residue_index = int(getattr(residue, "index", residue))
    match_b_to_a = _match_map_from_inputs(assign_a, assign_b, tmcs, image_path, similarity_limit, imcs)
    return _merge_dual_topology(molecule, residue_index, residue_b, match_b_to_a)


def merge_force_field(molecule_a, molecule_b, default_lambda, specific_lambda=None, intra_fep=False):
    """Merge A/B force-field states at a lambda value."""
    del intra_fep
    return _merge_force_field(molecule_a, molecule_b, float(default_lambda), specific_lambda or {})


def get_free_molecule(molecule, perturbating_residues, intra_fep=False):
    """Return a copy whose selected residues are decoupled with zero LJ and zero charge."""
    del intra_fep
    if not isinstance(perturbating_residues, (list, tuple, set)):
        perturbating_residues = [perturbating_residues]
    residue_indices = {int(getattr(residue, "index", residue)) for residue in perturbating_residues}
    free = molecule.copy()
    for residue_index, residue in enumerate(free.residues):
        if residue_index not in residue_indices:
            continue
        for atom in residue.atoms:
            atom.charge = 0.0
            atom.type = "ZERO_LJ_ATOM"
            atom.subsys = 1
    return free


def intramolecule_nb_to_nb14(molecule, perturbating_residues):
    """Compatibility hook; XpongeCPP builds exclusions from topology during export."""
    del perturbating_residues
    return molecule


def save_soft_core_lj(molecule=None):
    if molecule is not None:
        molecule.enable_lj_soft_core()
    return molecule


def save_hard_core_lj(molecule=None):
    if molecule is not None:
        molecule.enable_lj_soft_core(False)
        molecule.enable_subsys_division(False)
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


def _match_map_from_inputs(assign_a, assign_b, tmcs, image_path, similarity_limit, imcs):
    if isinstance(assign_a, dict):
        return {int(key): int(value) for key, value in assign_a.items()}
    if imcs is not None:
        return _match_map_from_imcs(imcs)
    if assign_a is None or assign_b is None:
        raise ValueError("merge_dual_topology requires an explicit match map or Assign objects")
    try:
        from rdkit.Chem import AllChem, Draw, RemoveHs, SanitizeMol, rdmolops
        from rdkit.Chem import rdFMCS as MCS
    except ImportError as exc:
        raise ImportError("RDKit is required to derive an FEP dual-topology atom map") from exc

    rdmol_a = _assign_to_rdkit(assign_a)
    rdmol_b = _assign_to_rdkit(assign_b)
    flags = rdmolops.SanitizeFlags
    rdmol_a_no_h = RemoveHs(rdmol_a, sanitize=False)
    SanitizeMol(rdmol_a_no_h, flags.SANITIZE_ALL & ~flags.SANITIZE_PROPERTIES)
    rdmol_b_no_h = RemoveHs(rdmol_b, sanitize=False)
    SanitizeMol(rdmol_b_no_h, flags.SANITIZE_ALL & ~flags.SANITIZE_PROPERTIES)
    result = MCS.FindMCS(
        [rdmol_a_no_h, rdmol_b_no_h],
        completeRingsOnly=True,
        bondCompare=MCS.BondCompare.CompareOrderExact,
        timeout=int(tmcs),
    )
    if result.queryMol is None:
        raise OSError("No common substructure found")
    match_a = rdmol_a.GetSubstructMatch(result.queryMol)
    a_hydrogens = {i: _hydrogen_count(rdmol_a, atom_id) for i, atom_id in enumerate(match_a)}
    match_b = sorted(
        rdmol_b.GetSubstructMatches(result.queryMol, uniquify=False),
        key=lambda match: sum(min(_hydrogen_count(rdmol_b, atom_id), a_hydrogens[i]) for i, atom_id in enumerate(match)),
    )[-1]
    match_b_to_a = {int(match_b[i]): int(match_a[i]) for i in range(len(match_a))}
    extra_map = {}
    for b_idx, a_idx in match_b_to_a.items():
        a_atom = rdmol_a.GetAtomWithIdx(a_idx)
        a_h = [neighbor.GetIdx() for neighbor in a_atom.GetNeighbors() if neighbor.GetAtomicNum() == 1]
        b_atom = rdmol_b.GetAtomWithIdx(b_idx)
        for neighbor in b_atom.GetNeighbors():
            if not a_h:
                break
            if neighbor.GetAtomicNum() == 1:
                extra_map[int(neighbor.GetIdx())] = int(a_h.pop())
    match_b_to_a.update(extra_map)
    tanimoto = len(match_b_to_a) / (len(assign_a.atoms) + len(assign_b.atoms) - len(match_b_to_a))
    if tanimoto <= similarity_limit:
        raise ValueError(f"similarity (={tanimoto}) should be greater than {similarity_limit}")
    if image_path:
        AllChem.Compute2DCoords(rdmol_a)
        AllChem.Compute2DCoords(rdmol_b)
        img = Draw.MolsToGridImage(
            [rdmol_a, rdmol_b],
            molsPerRow=1,
            subImgSize=(1200, 600),
            highlightAtomLists=[tuple(match_b_to_a.values()), tuple(match_b_to_a.keys())],
        )
        img.save(image_path)
    return match_b_to_a


def _match_map_from_imcs(path):
    match_a = None
    match_b = None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("mcs atoms in"):
                continue
            left, right = line.split(":", 1)
            atoms = tuple(int(item.strip()) for item in right.replace("(", "").replace(")", "").split(",") if item.strip())
            if "r1" in left:
                match_a = atoms
            elif "r2" in left:
                match_b = atoms
    if match_a is None or match_b is None or len(match_a) != len(match_b):
        raise ValueError(f"the format of the mcs file {path} is not right")
    return {int(match_b[i]): int(match_a[i]) for i in range(len(match_a))}


def _assign_to_rdkit(assign):
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise ImportError("RDKit is required to convert Assign to an RDKit molecule") from exc
    rw_mol = Chem.RWMol()
    for element in assign.atoms:
        rw_mol.AddAtom(Chem.Atom(element))
    bond_types = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
        -1: Chem.BondType.AROMATIC,
    }
    for atom_i, neighbors in enumerate(assign.bonds):
        for atom_j, order in neighbors.items():
            if int(atom_j) <= atom_i:
                continue
            rw_mol.AddBond(atom_i, int(atom_j), bond_types.get(int(order), Chem.BondType.SINGLE))
    return rw_mol.GetMol()


def _hydrogen_count(rdmol, atom_id):
    return rdmol.GetAtomWithIdx(int(atom_id)).GetTotalNumHs(includeNeighbors=True)


Prepare_LJ_Soft_Core = prepare_lj_soft_core
Merge_Dual_Topology = merge_dual_topology
Merge_Force_Field = merge_force_field
Get_Free_Molecule = get_free_molecule
Intramolecule_NB_To_NB14 = intramolecule_nb_to_nb14
Save_Soft_Core_LJ = save_soft_core_lj
Save_Hard_Core_LJ = save_hard_core_lj
Set_LJ_Type_B = set_lj_type_b
Set_Subsys = set_subsys
Enable_LJ_Soft_Core = enable_lj_soft_core
Add_Soft_Bond_From_A = add_soft_bond_from_a
Add_Soft_Bond_From_B = add_soft_bond_from_b
