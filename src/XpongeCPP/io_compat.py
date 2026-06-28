"""Compatibility readers layered on top of the core assignment API."""

import math
import re

from ._core import Assign


def _hy36_decode(width, field):
    text = field.strip()
    if not text:
        return None
    if text[0] in "+-" or text[0].isdigit():
        return int(text)
    digits_upper = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    digits_lower = "0123456789abcdefghijklmnopqrstuvwxyz"
    lower = field[0].islower()
    digits = digits_lower if lower else digits_upper
    value = 0
    for char in field:
        value = value * 36 + digits.index(char)
    power = 36 ** (width - 1)
    if lower:
        return value + 16 * power + 10 ** width
    return value - 10 * power + 10 ** width


def _guess_element(atom_name, explicit=""):
    explicit = explicit.strip()
    if explicit:
        return explicit[0].upper() + explicit[1:].lower()
    text = atom_name.strip()
    text = "".join(char for char in text if not char.isdigit())
    if not text:
        return ""
    if len(text) >= 2 and text[:2].capitalize() in {"Cl", "Br", "Na", "Mg", "Ca", "Zn", "Fe"}:
        return text[:2].capitalize()
    return text[0].upper()


def _read_text_source(source):
    if hasattr(source, "read"):
        return source.read()
    if isinstance(source, str) and ("\n" in source or source.lstrip().startswith("data_")):
        return source
    with open(source, encoding="utf-8") as handle:
        return handle.read()


def _rdmol_to_assignment(rdmol):
    from .helper.rdkit import rdmol_to_assign

    return rdmol_to_assign(rdmol)


def get_assignment_from_pdb(file, only_residue="", bond_tolerance=1.0, total_charge=None):
    text = _read_text_source(file)
    assignment = Assign()
    serial_to_atom = {}
    conect = []
    has_conect = False
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20].strip()
            if only_residue and resname != only_residue:
                continue
            serial = _hy36_decode(5, line[6:11])
            if serial is None:
                continue
            atom_name = line[12:16].strip()
            element = _guess_element(atom_name, line[76:78] if len(line) >= 78 else "")
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            assignment.name = resname
            serial_to_atom[serial] = assignment.atom_count
            assignment.add_atom(element, x, y, z, atom_name)
        elif line.startswith("CONECT"):
            has_conect = True
            atom = _hy36_decode(5, line[6:11])
            if atom is None:
                continue
            for start in range(11, min(len(line), 31), 5):
                field = line[start:start + 5]
                if not field.strip():
                    continue
                bonded_atom = _hy36_decode(5, field)
                if bonded_atom is not None:
                    conect.append((atom, bonded_atom))
    if assignment.atom_count == 0:
        raise OSError("The input is not a pdb file")
    for atom, bonded_atom in conect:
        if atom in serial_to_atom and bonded_atom in serial_to_atom:
            atom1 = serial_to_atom[atom]
            atom2 = serial_to_atom[bonded_atom]
            if atom1 < atom2:
                assignment.add_bond(atom1, atom2, 1)
    if not has_conect:
        assignment.determine_connectivity(tolerance=bond_tolerance)
    assignment.determine_bond_order(total_charge=total_charge)
    return assignment


def get_assignment_from_xyz(file, bond_tolerance=1.0, total_charge=None):
    text = _read_text_source(file)
    lines = text.splitlines()
    if not lines:
        raise OSError("The input is not a xyz file")
    atom_numbers = int(lines[0].strip())
    assignment = Assign()
    assignment.name = lines[1].strip()
    for index, line in enumerate(lines[2:2 + atom_numbers]):
        atom_name, x, y, z = line.split()[:4]
        assignment.add_atom(atom_name, float(x), float(y), float(z), f"{atom_name}{index + 1}")
    assignment.determine_connectivity(tolerance=bond_tolerance)
    assignment.determine_bond_order(total_charge=total_charge)
    return assignment


def get_assignment_from_smiles(smiles):
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise ImportError("RDKit is required for get_assignment_from_smiles") from exc

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    return _rdmol_to_assignment(mol)


def get_assignment_from_pubchem(parameter, keyword="name", total_charge=None, **kwargs):
    try:
        import pubchempy as pcp
    except ImportError as exc:
        raise ImportError("PubChemPy is required for get_assignment_from_pubchem") from exc

    compounds = pcp.get_compounds(parameter, keyword, record_type="3d")
    if not compounds:
        try:
            raise pcp.NotFoundError
        except AttributeError as exc:
            raise ValueError(f"PubChem query returned no compounds: {parameter}") from exc
    if len(compounds) != 1:
        raise NotImplementedError
    compound = compounds[0]
    assignment = Assign(str(getattr(compound, "cid", "PUBCHEM")))
    for atom in compound.atoms:
        assignment.add_atom(atom.element, atom.x, atom.y, atom.z)
    for bond in compound.bonds:
        assignment.add_bond(bond.aid1 - 1, bond.aid2 - 1, bond.order)
    assignment.determine_ring_and_bond_type()
    return assignment


def _cif_float(value):
    value = value.strip("'\"")
    if "(" in value:
        value = value.split("(", 1)[0]
    return float(value)


def _cif_loops(text):
    loops = []
    lines = text.splitlines()
    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue
        index += 1
        headers = []
        while index < len(lines) and lines[index].strip().startswith("_"):
            headers.append(lines[index].strip())
            index += 1
        rows = []
        while index < len(lines):
            line = lines[index].strip()
            if not line or line.startswith("#"):
                index += 1
                continue
            if line == "loop_" or line.startswith("_") or line.startswith("data_"):
                break
            parts = line.split()
            if len(parts) >= len(headers):
                rows.append(dict(zip(headers, parts)))
            index += 1
        loops.append((headers, rows))
    return loops


def _basis_vectors_from_length_and_angle(la, lb, lc, alpha, beta, gamma):
    alpha = math.radians(alpha)
    beta = math.radians(beta)
    gamma = math.radians(gamma)
    ax, ay, az = la, 0.0, 0.0
    bx, by, bz = lb * math.cos(gamma), lb * math.sin(gamma), 0.0
    cx = lc * math.cos(beta)
    cy = lc * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
    cz2 = lc * lc - cx * cx - cy * cy
    cz = math.sqrt(max(cz2, 0.0))
    return [(ax, ay, az), (bx, by, bz), (cx, cy, cz)]


def _parse_cif_symops(text, lattice_info):
    match = re.search(r"(_symmetry_equiv_pos_as_xyz|_space_group_symop_operation_xyz)\s+(.+?)(?!_)\n(_\S+|loop_\S*)", text, flags=re.DOTALL)
    if not match:
        return
    symops = match.group(2).replace("'", "")
    if set(symops) - set("+-,xyz0123456789\n /"):
        raise ValueError("the symmetry operator strings can only be simple math expression of x, y, z")
    symops = symops.replace("x", "1").replace("y", "1").replace("z", "1").strip()
    lattice_info["basis_position"] = [
        [eval(op) for op in symop.split(",") if op]  # pylint: disable=eval-used
        for symop in symops.split("\n")
    ]


def get_assignment_from_cif(source, total_charge=0, orthogonal_threshold=None, keep_cell_angle=True):
    text = _read_text_source(source)
    data_name = "CIF"
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("data_"):
            data_name = stripped[5:]
            break
    values = {}
    for raw in text.splitlines():
        parts = raw.strip().split(None, 1)
        if len(parts) == 2 and parts[0].startswith("_cell_"):
            values[parts[0]] = parts[1]

    lattice_info = {"scale": 1, "style": "custom"}
    basis = None
    if "_cell_length_a" in values:
        lengths = [
            _cif_float(values["_cell_length_a"]),
            _cif_float(values["_cell_length_b"]),
            _cif_float(values["_cell_length_c"]),
        ]
        lattice_info["cell_length"] = lengths
        angles = [
            _cif_float(values["_cell_angle_alpha"]),
            _cif_float(values["_cell_angle_beta"]),
            _cif_float(values["_cell_angle_gamma"]),
        ]
        if not keep_cell_angle:
            angles = [90, 90, 90]
        elif orthogonal_threshold is not None:
            angles = [90 if abs(angle - 90) < orthogonal_threshold else angle for angle in angles]
        lattice_info["cell_angle"] = angles
        basis = _basis_vectors_from_length_and_angle(*lengths, *angles)
    _parse_cif_symops(text, lattice_info)

    rows = []
    bond_rows = []
    for headers, loop_rows in _cif_loops(text):
        if "_atom_site_type_symbol" in headers:
            rows = loop_rows
        if "_geom_bond_atom_site_label_1" in headers:
            bond_rows = loop_rows
    if not rows:
        raise ValueError("CIF input does not contain a simple _atom_site loop")

    assignment = Assign(data_name)
    name_to_atom = {}
    for row in rows:
        element = row.get("_atom_site_type_symbol", row.get("_atom_site.type_symbol", "")).strip("'\"")
        name = row.get("_atom_site_label", row.get("_atom_site_label_atom_id", element)).strip("'\"")
        if "_atom_site_Cartn_x" in row or "_atom_site.Cartn_x" in row:
            x = _cif_float(row.get("_atom_site_Cartn_x", row.get("_atom_site.Cartn_x", "0")))
            y = _cif_float(row.get("_atom_site_Cartn_y", row.get("_atom_site.Cartn_y", "0")))
            z = _cif_float(row.get("_atom_site_Cartn_z", row.get("_atom_site.Cartn_z", "0")))
        elif "_atom_site_fract_x" in row:
            if basis is None:
                raise ValueError("fractional CIF coordinates require cell information")
            fx = _cif_float(row["_atom_site_fract_x"])
            fy = _cif_float(row["_atom_site_fract_y"])
            fz = _cif_float(row["_atom_site_fract_z"])
            x = fx * basis[0][0] + fy * basis[1][0] + fz * basis[2][0]
            y = fx * basis[0][1] + fy * basis[1][1] + fz * basis[2][1]
            z = fx * basis[0][2] + fy * basis[1][2] + fz * basis[2][2]
        else:
            raise ValueError("There is no atom position found in CIF input")
        name_to_atom[name] = assignment.atom_count
        assignment.add_atom(element, x, y, z, name)
    if bond_rows:
        for row in bond_rows:
            atom1 = row["_geom_bond_atom_site_label_1"].strip("'\"")
            atom2 = row["_geom_bond_atom_site_label_2"].strip("'\"")
            if atom1 in name_to_atom and atom2 in name_to_atom:
                assignment.add_bond(name_to_atom[atom1], name_to_atom[atom2], -1)
    else:
        assignment.determine_connectivity(1.2)
    assignment.determine_bond_order(True, total_charge)
    return assignment, lattice_info
