"""Register packaged Martini 3.0.0 compatibility data."""

from ... import GlobalSetting, load_ffitp, register_amber_lj_parameter
from . import data_path


def _register_same_type_lj_from_ffitp(filename):
    output = load_ffitp(str(filename))
    registered = {}
    for line in output["LJ"].splitlines()[1:]:
        words = line.split()
        if len(words) != 3:
            continue
        atom_pair, sigma, epsilon = words
        atom_type1, atom_type2 = atom_pair.split("-", 1)
        if atom_type1 == atom_type2:
            registered[atom_type1] = (float(epsilon), float(sigma) / 2.0)
    for atom_type, (epsilon, rmin) in registered.items():
        register_amber_lj_parameter(atom_type, atom_type, epsilon, rmin)


MARTINI300_DIR = data_path("martini300")

if str(MARTINI300_DIR) not in GlobalSetting.GMXIncludePaths:
    GlobalSetting.Add_GMX_Include_Path(str(MARTINI300_DIR))

_register_same_type_lj_from_ffitp(MARTINI300_DIR / "martini_v3.0.0.itp")
