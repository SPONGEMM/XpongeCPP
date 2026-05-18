"""Public entrypoints for the unified legacy compatibility layer."""

from .assign import install_legacy_assign_patches
from .bootstrap import install_legacy_bootstrap
from .molecule import install_molecule_io_methods
from .runtime import get_legacy_residue_links_override, install_legacy_runtime_patches
from .symbols import install_template_globals


def enable_legacy_namespace(namespace=None, template_names=None, overwrite=False):
    """Enable common legacy Xponge conveniences in one call."""
    install_legacy_bootstrap(namespace=namespace)
    return install_template_globals(namespace=namespace, template_names=template_names, overwrite=overwrite)


Enable_Legacy_Namespace = enable_legacy_namespace
Get_Legacy_Residue_Links_Override = get_legacy_residue_links_override
Install_Legacy_Assign_Patches = install_legacy_assign_patches
Install_Legacy_Runtime_Patches = install_legacy_runtime_patches
Install_Template_Globals = install_template_globals
Install_Molecule_IO_Methods = install_molecule_io_methods
Install_Legacy_Bootstrap = install_legacy_bootstrap
