"""Xponge-compatible GROMACS topology loaders."""

from io import StringIO
import os


class GlobalSetting:
    GMXIncludePaths = []
    _gmx_bonded_type_parsers = {}
    PDBProteinResidueNames = set()
    HISMap = {}

    @classmethod
    def Add_GMX_Include_Path(cls, path):
        cls.GMXIncludePaths.append(path)

    @classmethod
    def Set_GMX_Bonded_Type_Parser(cls, force_name, func, parser):
        cls._gmx_bonded_type_parsers[(force_name, int(func))] = parser

    @classmethod
    def get_gmx_bonded_type_parser(cls, force_name, func):
        return cls._gmx_bonded_type_parsers.get((force_name, int(func)), lambda words, mol, stat: None)

    @classmethod
    def Add_PDB_Residue_Name_Mapping(cls, place, pdb_name, real_name):
        from . import register_pdb_residue_name_mapping

        register_pdb_residue_name_mapping(place, pdb_name, real_name)

    @classmethod
    def Add_PDB_Residue_Alias_Mapping(cls, pdb_name, real_name):
        from . import register_pdb_residue_alias_mapping

        register_pdb_residue_alias_mapping(pdb_name, real_name)

    @classmethod
    def Add_HIS_Mapping(cls, residue_name, hid, hie, hip):
        from . import register_his_mapping

        register_his_mapping(residue_name, hid, hie, hip)
        cls.HISMap[residue_name] = {"HID": hid, "HIE": hie, "HIP": hip}


class GromacsTopologyIterator:
    def __init__(self, filename=None, macros=None):
        self.files = []
        self.filenames = []
        self.stack = []
        self.flag = ""
        self.macro_define_stat = []
        self.defined_macros = dict(macros or {})
        if filename:
            self._add_iterator_file(filename)

    def __iter__(self):
        self.flag = ""
        self.macro_define_stat = []
        self.stack = []
        return self

    def _read_raw_line(self):
        while self.files:
            line = self.files[-1].readline()
            if line:
                return line
            self.files[-1].close()
            self.files.pop()
            self.filenames.pop()
        return None

    def __next__(self):
        while True:
            line = self._read_raw_line()
            if line is None:
                raise StopIteration
            line = self._line_preprocess(line)
            if line is None:
                continue
            if line[0] == "#":
                self._line_define(line)
                continue
            if self.macro_define_stat and not self.macro_define_stat[-1]:
                continue
            if "[" in line and "]" in line:
                self.flag = line.strip()[1:-1].strip()
                continue
            for macro, tobecome in self.defined_macros.items():
                line = line.replace(macro, tobecome)
            return line

    def _add_iterator_file(self, filename):
        if self.files:
            filename = os.path.abspath(os.path.join(os.path.dirname(self.filenames[-1]), filename.replace('"', "")))
        else:
            filename = os.path.abspath(filename.replace('"', ""))
        if not os.path.exists(filename):
            base_name = os.path.basename(filename)
            for include_dir in GlobalSetting.GMXIncludePaths:
                candidate = os.path.abspath(os.path.join(include_dir, base_name))
                if os.path.exists(candidate):
                    filename = candidate
                    break
        self.files.append(open(filename, encoding="utf-8"))
        self.filenames.append(filename)

    def _line_preprocess(self, line):
        line = line.strip()
        comment = line.find(";")
        if comment >= 0:
            line = line[:comment]
        while line and line[-1] == "\\":
            extra = self._read_raw_line()
            if extra is None:
                line = line[:-1]
                break
            extra = extra.strip()
            comment = extra.find(";")
            if comment >= 0:
                extra = extra[:comment]
            line = line[:-1] + " " + extra
        if not line:
            return None
        return line

    def _line_define(self, line):
        words = line.split()
        if self.macro_define_stat and not self.macro_define_stat[-1] and words[0] not in {
            "#ifdef", "#ifndef", "#else", "#endif"
        }:
            return None
        if words[0] == "#ifdef":
            macro = words[1]
            self.macro_define_stat.append(bool((not self.macro_define_stat or self.macro_define_stat[-1]) and macro in self.defined_macros))
        elif words[0] == "#ifndef":
            macro = words[1]
            self.macro_define_stat.append(bool((not self.macro_define_stat or self.macro_define_stat[-1]) and macro not in self.defined_macros))
        elif words[0] == "#else":
            if len(self.macro_define_stat) <= 1 or self.macro_define_stat[-2]:
                self.macro_define_stat[-1] = not self.macro_define_stat[-1]
        elif words[0] == "#endif":
            self.macro_define_stat.pop()
        elif words[0] == "#define":
            self.defined_macros[words[1]] = line[line.find(words[1]) + len(words[1]):].strip() if len(words) > 2 else ""
        elif words[0] == "#include":
            self._add_iterator_file(words[1])
        elif words[0] == "#undef":
            self.defined_macros.pop(words[1])
        elif words[0] == "#error":
            raise AssertionError(line)
        return None


def _ffitp_dihedrals(line, buffers):
    words = line.split()
    func = words[4]
    if func == "1":
        buffers["dihedrals"].append("-".join(words[:4]) + " " + " ".join(words[5:]) + " 0\n")
    elif func == "2":
        temp = [words[1], words[2], words[0], words[3]]
        temp2 = [words[1], words[2], words[3], words[0]]
        if words[0][0] == "O" or words[3] == "C" or words[3] == "CN3T":
            temp = temp2
        buffers["impropers"].append("-".join(temp) + f" {float(words[5])} {float(words[6]) / 2}\n")
    elif func == "3":
        buffers["RB_dihedrals"].append("-".join(words[:4]) + " " + " ".join(words[5:]) + "\n")
    elif func == "4":
        buffers["periodic_impropers"].append("-".join(words[:4]) + " " + " ".join(words[5:]) + "\n")
    elif func == "9":
        for i in range(5, len(words), 20):
            buffers["dihedrals"].append("-".join(words[:4]) + " " + " ".join(words[i:i + 3]) + " 0\n")
    else:
        raise NotImplementedError(f"Unsupported dihedral function type {func} for line:\n{line}")


def load_ffitp(filename, macros=None):
    iterator = GromacsTopologyIterator(filename, macros)
    output = {"cmaps": {}, "bond_type_names": {}}
    buffers = {
        "nb14": ["name  kLJ  kee\n"],
        "atomtypes": ["name mass charge[e] LJtype\n"],
        "bonds": ["name b[nm] k[kJ/mol·nm^-2]\n"],
        "angles": ["name b[degree] k[kJ/mol·rad^-2]\n"],
        "Urey-Bradley": ["name b[degree] k[kJ/mol·rad^-2] r13[nm] kUB[kJ/mol·nm^-2]\n"],
        "dihedrals": ["name phi0[degree] k[kJ/mol] periodicity  reset\n"],
        "periodic_impropers": ["name phi0[degree] k[kJ/mol] periodicity\n"],
        "impropers": ["name phi0[degree] k[kJ/mol·rad^-2]\n"],
        "RB_dihedrals": ["name c0[kJ/mol] c1[kJ/mol] c2[kJ/mol] c3[kJ/mol] c4[kJ/mol] c5[kJ/mol]\n"],
    }
    fudge_lj = 1.0
    fudge_qq = 1.0
    for line in iterator:
        if not iterator.flag:
            continue
        words = line.split()
        if iterator.flag == "defaults":
            assert int(words[0]) == 1, "SPONGE Only supports Lennard-Jones now"
            if int(words[1]) == 1:
                buffers["LJ"] = ["name A[kJ/mol·nm^6] B[kJ/mol·nm^12]\n"]
                buffers["nb14_extra"] = ["name A[kJ/mol·nm^6] B[kJ/mol·nm^12] kee\n"]
            else:
                buffers["LJ"] = ["name sigma[nm] epsilon[kJ/mol] \n"]
                buffers["nb14_extra"] = ["name sigma[nm] epsilon[kJ/mol] kee\n"]
            fudge_lj = float(words[3]) if len(words) > 3 else 1.0
            fudge_qq = float(words[4]) if len(words) > 4 else 1.0
            if len(words) > 2 and words[2] == "yes":
                buffers["nb14"].append(f"X-X {fudge_lj} {fudge_qq}\n")
        elif iterator.flag == "atomtypes":
            offset = 0
            if len(words) == 8:
                output["bond_type_names"][words[0]] = words.pop(1)
            elif len(words) == 6:
                offset = 1
            buffers["atomtypes"].append(f"{words[0]} {float(words[2 - offset])} {float(words[3 - offset])} {words[0]}\n")
            buffers.setdefault("LJ", ["name sigma[nm] epsilon[kJ/mol] \n"]).append(
                f"{words[0]}-{words[0]} {float(words[5 - offset])} {float(words[6 - offset])}\n"
            )
        elif iterator.flag == "pairtypes":
            if len(words) <= 3:
                buffers["nb14"].append(f"{words[0]}-{words[1]} {fudge_lj} {fudge_qq}\n")
            elif words[2] == "1":
                buffers.setdefault("nb14_extra", ["name sigma[nm] epsilon[kJ/mol] kee\n"]).append(
                    f"{words[0]}-{words[1]} {float(words[3])} {float(words[4])} {fudge_qq}\n"
                )
                buffers["nb14"].append(f"{words[0]}-{words[1]} 0 0\n")
            elif words[2] == "2":
                raise NotImplementedError
        elif iterator.flag == "bondtypes":
            if words[2] != "1":
                raise NotImplementedError
            buffers["bonds"].append(f"{words[0]}-{words[1]} {float(words[3])} {float(words[4]) / 2}\n")
        elif iterator.flag == "angletypes":
            if words[3] == "1":
                buffers["angles"].append("-".join(words[:3]) + f" {float(words[4])} {float(words[5]) / 2}\n")
            elif words[3] == "5":
                buffers["Urey-Bradley"].append(
                    "-".join(words[:3]) + f" {float(words[4])} {float(words[5]) / 2} {float(words[6])} {float(words[7]) / 2}\n"
                )
            else:
                raise NotImplementedError
        elif iterator.flag == "dihedraltypes":
            _ffitp_dihedrals(line, buffers)
        elif iterator.flag == "cmaptypes":
            output["cmaps"]["-".join(words[:5])] = {
                "resolution": int(words[7]),
                "parameters": [float(word) / 4.184 for word in words[8:]],
            }
        elif iterator.flag == "nonbond_params":
            buffers.setdefault("LJ", ["name sigma[nm] epsilon[kJ/mol] \n"]).append(
                f"{words[0]}-{words[1]} {float(words[3])} {float(words[4])}\n"
            )
    for key, value in buffers.items():
        output[key] = "".join(value)
    return output


def _molitp_find_tail_residue(filename, macros, water_replace):
    heads = {}
    tails = {}
    system_molecules = set()
    current = None
    iterator = GromacsTopologyIterator(filename, macros)
    for line in iterator:
        if iterator.flag == "moleculetype":
            current = line.split()[0]
            tails[current] = -999999
        elif iterator.flag == "atoms":
            resnr = int(line.split()[2])
            tails[current] = max(tails[current], resnr)
            heads.setdefault(current, resnr)
        elif iterator.flag == "molecules":
            molname = line.split()[0]
            if not (molname == "SOL" and water_replace):
                system_molecules.add(molname)
    return heads, tails, system_molecules


def _molitp_to_mol2(name, atoms, bonds):
    residue_ids = {}
    for atom in atoms:
        residue_ids.setdefault(atom["resnr"], len(residue_ids) + 1)
    out = [
        "@<TRIPOS>MOLECULE",
        name,
        f"{len(atoms)} {len(bonds)} {len(residue_ids)}",
        "SMALL",
        "USER_CHARGES",
        "@<TRIPOS>ATOM",
    ]
    for atom in atoms:
        out.append(
            f"{atom['nr']} {atom['atom_name']} 0.0 0.0 0.0 {atom['atom_type']} "
            f"{residue_ids[atom['resnr']]} {atom['resname']} {atom['charge']}"
        )
    out.append("@<TRIPOS>BOND")
    for index, (atom1, atom2) in enumerate(bonds, 1):
        out.append(f"{index} {atom1} {atom2} 1")
    return "\n".join(out) + "\n"


def _molitp_next_generated_residue_name(base_name, generated_variants, has_template):
    candidate = base_name
    suffix = 1
    while candidate in generated_variants or has_template(candidate):
        candidate = f"{base_name}_{suffix}"
        suffix += 1
    return candidate


def _molitp_residue_signature_key(signature):
    return tuple(sorted((atom_name, atom_type, round(charge, 8)) for atom_name, (atom_type, charge) in signature.items()))


def _molitp_apply_template_variants(mol_specs, molecule_order):
    from . import get_template_molecule, has_template

    template_signatures = {}
    generated_variants = {}
    signature_to_name = {}

    def get_signature(name):
        if name in generated_variants:
            return generated_variants[name]
        if name not in template_signatures:
            template = get_template_molecule(name)
            template_signatures[name] = {
                atom.name: (atom.type, atom.charge)
                for atom in template.residues[0].atoms
            }
        return template_signatures[name]

    for mol_name in molecule_order:
        residues = {}
        residue_order = []
        for atom in mol_specs[mol_name]["atoms"]:
            residues.setdefault(atom["resnr"], []).append(atom)
            if atom["resnr"] not in residue_order:
                residue_order.append(atom["resnr"])

        for resnr in residue_order:
            residue_atoms = residues[resnr]
            base_name = residue_atoms[0]["resname"]
            if not has_template(base_name):
                continue

            current_name = base_name
            current_signature = dict(get_signature(base_name))
            for atom in residue_atoms:
                atom_name = atom["atom_name"]
                desired = (atom["atom_type"], atom["charge"])
                if current_signature.get(atom_name) == desired:
                    continue

                compatible_name = None
                compatible_signature = None
                for candidate_name, candidate_signature in generated_variants.items():
                    if candidate_name != base_name and not candidate_name.startswith(base_name + "_"):
                        continue
                    if candidate_signature.get(atom_name) != desired:
                        continue
                    compatible_name = candidate_name
                    compatible_signature = candidate_signature
                    break

                if compatible_name is not None:
                    current_name = compatible_name
                    current_signature = dict(compatible_signature)
                    continue

                new_name = _molitp_next_generated_residue_name(base_name, generated_variants, has_template)
                new_signature = dict(current_signature)
                new_signature[atom_name] = desired
                generated_variants[new_name] = new_signature
                signature_to_name.setdefault(_molitp_residue_signature_key(new_signature), new_name)
                current_name = new_name
                current_signature = new_signature

            final_signature_key = _molitp_residue_signature_key(current_signature)
            current_name = signature_to_name.get(final_signature_key, current_name)
            signature_to_name.setdefault(final_signature_key, current_name)
            for atom in residue_atoms:
                atom["resname"] = current_name


def load_molitp(filename, water_replace=True, head_prefix="N", tail_prefix="C", macros=None):
    from . import get_template_molecule, load_mol2

    heads, tails, system_molecules = _molitp_find_tail_residue(filename, macros, water_replace)
    mol_specs = {}
    molecule_order = []
    current_name = None
    skip = False
    system_name = None
    system_counts = []
    iterator = GromacsTopologyIterator(filename, macros)
    for line in iterator:
        words = line.split()
        if iterator.flag == "moleculetype":
            current_name = words[0]
            skip = bool(water_replace and current_name == "SOL") or current_name not in system_molecules
            if not skip:
                mol_specs[current_name] = {
                    "atoms": [],
                    "bonds": [],
                    "pairs": [],
                    "angles": [],
                    "dihedrals": [],
                    "virtual_sites2": [],
                    "virtual_sites3": [],
                    "virtual_sites4": [],
                }
                molecule_order.append(current_name)
        elif iterator.flag == "atoms" and not skip:
            nr = int(words[0])
            resnr = int(words[2])
            resname = words[3]
            if tails[current_name] != heads[current_name] and resnr == heads[current_name]:
                resname = head_prefix + resname
            if tails[current_name] != heads[current_name] and resnr == tails[current_name]:
                resname = tail_prefix + resname
            mol_specs[current_name]["atoms"].append(
                {
                    "nr": nr,
                    "atom_type": words[1],
                    "resnr": resnr,
                    "resname": resname,
                    "atom_name": words[4],
                    "charge": float(words[6]),
                }
            )
        elif iterator.flag == "bonds" and not skip:
            mol_specs[current_name]["bonds"].append((int(words[0]), int(words[1])))
            mol_specs[current_name]["bonds"].append(tuple(words))
        elif iterator.flag in {
            "pairs", "angles", "dihedrals", "virtual_sites2", "virtual_sites3", "virtual_sites4"
        } and not skip:
            mol_specs[current_name][iterator.flag].append(tuple(words))
        elif iterator.flag == "system":
            system_name = line.strip()
        elif iterator.flag == "molecules":
            system_counts.append((words[0], int(words[1])))

    _molitp_apply_template_variants(mol_specs, molecule_order)

    mols = {}
    for name in molecule_order:
        spec = mol_specs[name]
        bond_pairs = [bond for bond in spec["bonds"] if len(bond) == 2]
        mol = load_mol2(StringIO(_molitp_to_mol2(name, spec["atoms"], bond_pairs)))
        atom_views = []
        for residue in mol.residues:
            atom_views.extend(residue.atoms)
        stat = {index: atom for index, atom in enumerate(atom_views, 1)}

        for words in spec["bonds"]:
            if len(words) != 2:
                parser = GlobalSetting.get_gmx_bonded_type_parser("bond", int(words[2]))
                parser(words, mol, stat)
        for section, force_name, func_index in (
            ("pairs", "pair", 2),
            ("angles", "angle", 3),
            ("dihedrals", "dihedral", 4),
            ("virtual_sites2", "virtual_site2", 3),
            ("virtual_sites3", "virtual_site3", 4),
            ("virtual_sites4", "virtual_site4", 5),
        ):
            for words in spec[section]:
                parser = GlobalSetting.get_gmx_bonded_type_parser(force_name, int(words[func_index]))
                parser(words, mol, stat)
        mols[name] = mol
    if water_replace and "SOL" in [name for name, _ in system_counts]:
        try:
            mols["SOL"] = get_template_molecule("WAT")
        except RuntimeError:
            pass

    system = None
    if system_name is not None:
        from . import Molecule

        system = Molecule(system_name)
        for name, count in system_counts:
            if name not in mols:
                continue
            for _ in range(count):
                system.add_molecule(mols[name])
    return system, mols
