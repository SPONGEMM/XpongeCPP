"""Register RSFF2C residue-specific CMAP parameters."""

from ... import register_amber_cmap_parameter
from . import data_path


def _register_rsff2c_cmaps():
    text = data_path("RSFF2C.dat").read_text()
    chunks = text.strip().replace("%FLAG", "%FORMAT(8(F9.5))").split("%FORMAT(8(F9.5))")[1:]
    names = [line.split("COMMENT")[1].split() for line in chunks[::2]]
    templates = {
        "phi/psi": "C-N-{res}@CA-C-N",
        "chi/phi": "G-CB-{res}@CA-N-C",
        "chi/psi": "G-CB-{res}@CA-C-N",
    }
    for name, parameter_block in zip(names, chunks[1::2]):
        key = templates[name[1]].format(res=name[0])
        parameters = [float(word) for word in parameter_block.split()]
        register_amber_cmap_parameter(key, 24, parameters)


_register_rsff2c_cmaps()
