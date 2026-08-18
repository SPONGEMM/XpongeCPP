"""Amber frcmod text parser compatible with Xponge 1.7b8."""

from __future__ import annotations

import re


def _nb14_row(line, atoms):
    number = r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
    scee = re.search(r"SCEE\s*=\s*" + number, line)
    scnb = re.search(r"SCNB\s*=\s*" + number, line)
    if not scee and not scnb:
        return ""
    electrostatic = 1.0 / float(scee.group(1)) if scee else 1.0 / 1.2
    lennard_jones = 1.0 / float(scnb.group(1)) if scnb else 1.0 / 2.0
    return f"{'-'.join(atoms)} {lennard_jones} {electrostatic}\n"


def _cmap_row(line, cmap, temporary, cmap_flag):
    if line.startswith("%FLAG"):
        if "CMAP_COUNT" in line:
            if temporary:
                for residue in temporary["residues"]:
                    cmap[f"C-N-{residue}@XC-C-N"] = {
                        "resolution": temporary["info"]["resolution"],
                        "parameters": temporary["info"]["parameters"],
                    }
            temporary = {
                "residues": [],
                "info": {
                    "resolution": 24,
                    "count": int(line.split()[-1]),
                    "parameters": [],
                },
            }
            cmap_flag = "CMAP_COUNT"
        elif "CMAP_RESOLUTION" in line:
            temporary["info"]["resolution"] = int(line.split()[-1])
            cmap_flag = "CMAP_RESOLUTION"
        elif "CMAP_RESLIST" in line:
            cmap_flag = "CMAP_RESLIST"
        elif "CMAP_TITLE" in line:
            cmap_flag = "CMAP_TITLE"
        elif "CMAP_PARAMETER" in line:
            cmap_flag = "CMAP_PARAMETER"
    elif cmap_flag == "CMAP_RESLIST":
        temporary["residues"].extend(line.split())
    elif cmap_flag == "CMAP_PARAMETER":
        temporary["info"]["parameters"].extend(float(value) for value in line.split())
    return temporary, cmap_flag


def _atoms_and_words(line, width, previous=None):
    atom_field = line[:width]
    if atom_field.strip():
        return [word.strip() for word in atom_field.split("-")], line[width:].split()
    if previous is None:
        raise ValueError("Amber parameter continuation line has no preceding atom types")
    return list(previous), line[width:].split()


def parse_frcmod(filename, nbtype="RE", include_nb14=False):
    """Return Xponge-compatible parameter text blocks from an Amber frcmod."""
    atom_types = {}
    bonds = ["name  k[kcal/mol·A^-2]    b[A]\n"]
    angles = ["name  k[kcal/mol·rad^-2]    b[degree]\n"]
    propers = ["name  k[kcal/mol]    phi0[degree]    periodicity    reset\n"]
    nb14s = ["name    kLJ    kee\n"]
    impropers = ["name  k[kcal/mol]    phi0[degree]    periodicity\n"]
    cmap = {}
    cmap_flag = None
    temporary_cmap = {"residues": []}
    last_dihedral_atoms = None
    reset = 1
    flag = None

    if nbtype == "SK":
        raise NotImplementedError
    if nbtype == "AC":
        ljs = ["name A[kcal/mol·A^-12]   B[kcal/mol·A^-6]\n"]
    elif nbtype == "RE":
        ljs = ["name rmin[A]   epsilon[kcal/mol]\n"]
    else:
        raise ValueError(f"Unsupported Amber frcmod nonbonded type: {nbtype}")

    with open(filename, encoding="utf-8") as frcmod:
        next(frcmod, None)
        for line in frcmod:
            if not line.strip():
                continue
            words = line.split()
            if flag != "CMAP" and len(words) == 1:
                flag = line.strip()
                if flag[:4] == "DIHE":
                    last_dihedral_atoms = None
                    reset = 1
            elif flag and flag[:4] == "MASS":
                atom_types[words[0]] = words[1]
            elif flag and flag[:4] == "BOND":
                atoms, words = _atoms_and_words(line, 5)
                bonds.append(
                    "-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\n"
                )
            elif flag and flag[:4] == "ANGL":
                atoms, words = _atoms_and_words(line, 8)
                angles.append(
                    "-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\n"
                )
            elif flag and flag[:4] == "DIHE":
                atoms, words = _atoms_and_words(line, 11, last_dihedral_atoms)
                last_dihedral_atoms = atoms
                propers.append(
                    "-".join(atoms)
                    + "\t"
                    + str(float(words[1]) / int(words[0]))
                    + "\t"
                    + words[2]
                    + "\t"
                    + str(abs(int(float(words[3]))))
                    + "\t"
                    + str(reset)
                    + "\n"
                )
                nb14s.append(_nb14_row(line, atoms))
                reset = 0 if int(float(words[3])) < 0 else 1
            elif flag and flag[:4] == "IMPR":
                atoms, words = _atoms_and_words(line, 11)
                impropers.append(
                    "-".join(atoms)
                    + "\t"
                    + words[0]
                    + "\t"
                    + words[1]
                    + "\t"
                    + str(int(float(words[2])))
                    + "\n"
                )
            elif flag and flag[:4] == "NONB":
                ljs.append(
                    words[0] + "-" + words[0] + "\t" + words[1] + "\t" + words[2] + "\n"
                )
            elif flag and flag[:4] == "CMAP":
                temporary_cmap, cmap_flag = _cmap_row(
                    line, cmap, temporary_cmap, cmap_flag
                )

    for residue in temporary_cmap["residues"]:
        cmap[f"C-N-{residue}@XC-C-N"] = {
            "resolution": temporary_cmap["info"]["resolution"],
            "parameters": temporary_cmap["info"]["parameters"],
        }
    atoms = ["name  mass  LJtype\n"]
    atoms.extend(
        atom + "\t" + mass + "\t" + atom + "\n"
        for atom, mass in atom_types.items()
    )
    values = [
        "".join(atoms),
        "".join(bonds),
        "".join(angles),
        "".join(propers),
        "".join(impropers),
        "".join(ljs),
    ]
    if include_nb14:
        values.append("".join(nb14s))
    values.append(cmap)
    return values


def _read_harmonic_rows(stream, rows, width):
    for line in stream:
        if not line.strip():
            break
        atoms, words = _atoms_and_words(line, width)
        rows.append(
            "-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\n"
        )
    return rows


def parse_parmdat(filename):
    """Return Xponge-compatible parameter blocks from an Amber parmdat."""
    with open(filename, encoding="utf-8") as parmdat:
        next(parmdat, None)
        atom_types = {}
        lj_types = {}
        for line in parmdat:
            if not line.strip():
                break
            words = line.split()
            atom_types[words[0]] = words[1]
            lj_types[words[0]] = words[0]

        next(parmdat, None)
        bonds = _read_harmonic_rows(
            parmdat,
            ["name  k[kcal/mol·A^-2]    b[A]\n"],
            5,
        )
        angles = _read_harmonic_rows(
            parmdat,
            ["name  k[kcal/mol·rad^-2]    b[degree]\n"],
            8,
        )

        reset = 1
        propers = ["name  k[kcal/mol]    phi0[degree]    periodicity    reset\n"]
        nb14s = ["name    kLJ     kee\n"]
        atoms = None
        for line in parmdat:
            if not line.strip():
                break
            atoms, words = _atoms_and_words(line, 11, atoms)
            nb14s.append(_nb14_row(line, atoms))
            propers.append(
                "-".join(atoms)
                + "\t"
                + str(float(words[1]) / int(words[0]))
                + "\t"
                + words[2]
                + "\t"
                + str(abs(int(float(words[3]))))
                + "\t"
                + str(reset)
                + "\n"
            )
            reset = 0 if int(float(words[3])) < 0 else 1

        impropers = ["name  k[kcal/mol]    phi0[degree]    periodicity\n"]
        for line in parmdat:
            if not line.strip():
                break
            improper_atoms, words = _atoms_and_words(line, 11)
            impropers.append(
                "-".join(improper_atoms)
                + "\t"
                + words[0]
                + "\t"
                + words[1]
                + "\t"
                + str(int(float(words[2])))
                + "\n"
            )

        next(parmdat, None)
        next(parmdat, None)
        for line in parmdat:
            if not line.strip():
                break
            aliases = line.split()
            primary = aliases.pop(0)
            for alias in aliases:
                lj_types[alias] = primary

        mode = next(parmdat).split()[1]
        if mode == "SK":
            raise NotImplementedError
        if mode == "AC":
            ljs = ["name A[kcal/mol·A^-12]   B[kcal/mol·A^-6]\n"]
        elif mode == "RE":
            ljs = ["name rmin[A]   epsilon[kcal/mol]\n"]
        else:
            raise ValueError(f"Unsupported Amber parmdat nonbonded type: {mode}")
        for line in parmdat:
            if not line.strip():
                break
            words = line.split()
            ljs.append(
                words[0] + "-" + words[0] + "\t" + words[1] + "\t" + words[2] + "\n"
            )

    atoms = ["name  mass  LJtype\n"]
    atoms.extend(
        atom + "\t" + mass + "\t" + lj_types[atom] + "\n"
        for atom, mass in atom_types.items()
    )
    return [
        "".join(atoms),
        "".join(bonds),
        "".join(angles),
        "".join(propers),
        "".join(impropers),
        "".join(ljs),
        "".join(nb14s),
    ]
