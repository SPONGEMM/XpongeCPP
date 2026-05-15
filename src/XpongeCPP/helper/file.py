"""Legacy helper.file compatibility helpers."""

from __future__ import annotations

import io
import re
import sys
from importlib import import_module
from pathlib import Path


def import_python_script(path):
    path = path if isinstance(path, Path) else Path(path)
    sys.path.append(str(path.parent))
    if path.suffix != ".py":
        raise TypeError(f"{path} should be a python script")
    return import_module(path.stem)


def file_filter(infile, outfile, reg_exp, replace_dict):
    if not isinstance(reg_exp, list):
        raise TypeError("reg_exp should be a list of regular expressions")
    if not isinstance(replace_dict, dict):
        raise TypeError("replace_dict should be a dict of regular expressions and the replacement")
    if not isinstance(infile, io.IOBase):
        infile = open(infile, "r", encoding="utf-8")
    lines = ""
    with infile as f:
        for line in infile:
            for keyword in reg_exp:
                if not isinstance(keyword, str):
                    raise TypeError("reg_exp should be a list of regular expressions")
                if re.match(keyword, line):
                    for reg, rep in replace_dict.items():
                        line = re.sub(reg, rep, line)
                    lines += line
                    break
    if not isinstance(outfile, io.IOBase):
        with open(outfile, "w", encoding="utf-8") as f:
            f.write(lines)
    else:
        outfile.write(lines)


def pdb_filter(infile, outfile, heads, hetero_residues, chains=None, rename_ions=None):
    if not isinstance(heads, list):
        raise TypeError("heads should be a list")
    if not isinstance(hetero_residues, list):
        raise TypeError("hetero_residues should be a list")
    if rename_ions is None:
        rename_ions = {}
    if not isinstance(rename_ions, dict):
        raise TypeError("replace_dict should be a dict of regular expressions and the replacement")
    replace_dict = {}
    for a, b in rename_ions.items():
        if len(a) == 1:
            aname = f"{a}   | {a}  |  {a} |   {a}"
            rname = f"{a}  | {a} |  {a}"
        elif len(a) == 2:
            aname = f"{a}  | {a} |  {a}"
            rname = f" {a}|{a} "
        elif len(a) == 3:
            aname = f"{a} | {a}"
            rname = a
        else:
            raise ValueError("The ion name in a pdb file should not be longer than 3 characters")
        replace_dict[f"(^HETATM [ 0-9]{{4}} )({aname})(.)({rname})"] = r"\g<1>" + f"{b:4s}" + r"\g<3>" + f"{b:3s}"
    reg_exp = []
    for head in heads:
        if head == "ATOM" and chains is not None:
            for chain in chains:
                reg_exp.append("^ATOM.{17}%s" % chain)
        elif head == "SEQRES" and chains is not None:
            for chain in chains:
                reg_exp.append("^SEQRES.{5}%s" % chain)
        elif head == "TER" and chains is not None:
            reg_exp.append(r"^TER\s*$")
            for chain in chains:
                reg_exp.append("^TER.{18}%s" % chain)
        else:
            reg_exp.append(f"^{head}")
    for hetres in hetero_residues:
        if len(hetres) == 1:
            hetres = f"{hetres}  | {hetres} |  {hetres}"
        elif len(hetres) == 2:
            hetres = f" {hetres}|{hetres} "
        reg_exp.append("^HETATM.{11}%s" % (hetres))
    file_filter(infile, outfile, reg_exp, replace_dict)


__all__ = ["file_filter", "import_python_script", "pdb_filter"]
