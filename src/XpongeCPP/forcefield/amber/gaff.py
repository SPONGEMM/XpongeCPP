"""Amber GAFF compatibility helpers."""

from ._forcefield_family import activate_forcefield_family

activate_forcefield_family("small_molecule", "gaff")

import os
from importlib import import_module
from tempfile import TemporaryDirectory

from ... import (
    AtomType,
    Molecule,
    Save_Mol2,
    get_template_molecule,
    has_template,
    implemented_gaff_assign_types,
    molecule_from_residuetype,
    register_amber_parmdat_file,
)
from . import data_path
from . import load_parameters_from_frcmod

register_amber_parmdat_file(str(data_path("gaff.dat")))
AtomType.New_From_String("\n".join(implemented_gaff_assign_types()))


def _import_xpongelib():
    try:
        xlib = import_module("XpongeLib")
    except ImportError as exc:
        raise ImportError(
            "gaff.parmchk2_gaff requires the external XpongeLib runtime "
            "(install mokda-xpongelib / XpongeLib first)"
        ) from exc
    if not hasattr(xlib, "_parmchk2"):
        raise ImportError("Installed XpongeLib does not expose _parmchk2")
    return xlib


def _coerce_parmchk2_input(ifname):
    if isinstance(ifname, (str, os.PathLike)):
        return str(ifname), None
    tempdir = TemporaryDirectory()
    tempfile = os.path.join(tempdir.name, "temp.mol2")
    if isinstance(ifname, Molecule):
        Save_Mol2(ifname, tempfile)
        return tempfile, tempdir
    if hasattr(ifname, "atom_count") and hasattr(ifname, "bond_count") and hasattr(ifname, "name"):
        Save_Mol2(molecule_from_residuetype(ifname), tempfile)
        return tempfile, tempdir
    if hasattr(ifname, "name") and has_template(ifname.name):
        Save_Mol2(get_template_molecule(ifname.name), tempfile)
        return tempfile, tempdir
    tempdir.cleanup()
    raise TypeError("parmchk2_gaff expects a mol2 path or a template-like molecule object")


def _is_gaff_like_type(atom_type):
    return any(ch.islower() for ch in atom_type)


def _filter_mixed_gaff_mol2(mol2_path):
    sections = {}
    current = None
    with open(mol2_path, encoding="utf-8", errors="ignore") as handle:
        for line in handle.read().splitlines():
            if line.startswith("@<TRIPOS>"):
                current = line.strip()
                sections[current] = []
                continue
            if current is not None:
                sections[current].append(line)
    atom_records = []
    for line in sections.get("@<TRIPOS>ATOM", []):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) < 6:
            return mol2_path, None
        atom_records.append(fields)
    if not atom_records:
        return mol2_path, None
    keep_old_ids = {int(fields[0]) for fields in atom_records if _is_gaff_like_type(fields[5])}
    if not keep_old_ids or len(keep_old_ids) == len(atom_records):
        return mol2_path, None
    tempdir = TemporaryDirectory()
    filtered_path = os.path.join(tempdir.name, "gaff_subset.mol2")
    old_to_new = {old_id: index + 1 for index, old_id in enumerate(sorted(keep_old_ids))}
    sub_old_to_new = {}
    filtered_atoms = []
    for fields in atom_records:
        old_id = int(fields[0])
        if old_id not in old_to_new:
            continue
        fields = list(fields)
        fields[0] = str(old_to_new[old_id])
        sub_id = int(fields[6]) if len(fields) > 6 else 1
        sub_old_to_new.setdefault(sub_id, len(sub_old_to_new) + 1)
        if len(fields) > 6:
            fields[6] = str(sub_old_to_new[sub_id])
        filtered_atoms.append(fields)
    filtered_bonds = []
    for line in sections.get("@<TRIPOS>BOND", []):
        fields = line.split()
        if len(fields) < 4:
            continue
        old_a = int(fields[1])
        old_b = int(fields[2])
        if old_a not in old_to_new or old_b not in old_to_new:
            continue
        fields = list(fields)
        fields[0] = str(len(filtered_bonds) + 1)
        fields[1] = str(old_to_new[old_a])
        fields[2] = str(old_to_new[old_b])
        filtered_bonds.append(fields)
    filtered_substructures = []
    for old_sub_id, new_sub_id in sorted(sub_old_to_new.items(), key=lambda item: item[1]):
        source = None
        for line in sections.get("@<TRIPOS>SUBSTRUCTURE", []):
            fields = line.split()
            if fields and int(fields[0]) == old_sub_id:
                source = list(fields)
                break
        if source is None:
            source = [str(old_sub_id), "MOL", "1", "TEMP", "0", "****", "****", "0", "ROOT"]
        source[0] = str(new_sub_id)
        if len(source) > 2:
            source[2] = "1"
        filtered_substructures.append(source)
    with open(filtered_path, "w", encoding="utf-8") as handle:
        handle.write("@<TRIPOS>MOLECULE\n")
        handle.write("GAFF_SUBSET\n")
        handle.write(f"{len(filtered_atoms):6d}{len(filtered_bonds):6d}{len(filtered_substructures):6d}     0     1\n")
        handle.write("SMALL\nUSER_CHARGES\n")
        handle.write("@<TRIPOS>ATOM\n")
        for fields in filtered_atoms:
            handle.write(" ".join(fields) + "\n")
        handle.write("@<TRIPOS>BOND\n")
        for fields in filtered_bonds:
            handle.write(" ".join(fields) + "\n")
        handle.write("@<TRIPOS>SUBSTRUCTURE\n")
        for fields in filtered_substructures:
            handle.write(" ".join(fields) + "\n")
    return filtered_path, tempdir


def parmchk2_gaff(ifname, ofname, direct_load=True, keep=True):
    """Generate frcmod parameters with legacy Xponge-compatible semantics."""
    xlib = _import_xpongelib()
    mol2_path, tempdir = _coerce_parmchk2_input(ifname)
    filtered_tempdir = None
    try:
        mol2_path, filtered_tempdir = _filter_mixed_gaff_mol2(mol2_path)
        datapath = os.path.dirname(xlib.__file__)
        xlib._parmchk2(mol2_path, "mol2", str(ofname), datapath, 0, 1, 1)
        if direct_load:
            load_parameters_from_frcmod(ofname, prefix=False)
        if not keep:
            os.remove(ofname)
    finally:
        if filtered_tempdir is not None:
            filtered_tempdir.cleanup()
        if tempdir is not None:
            tempdir.cleanup()
