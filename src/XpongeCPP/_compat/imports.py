"""Legacy import-path and package-alias helpers.

This module is the centralized place for thin compatibility shims.  Legacy-facing
modules should prefer importing helpers from here instead of hand-writing
re-export logic in many files.
"""

from __future__ import annotations

import logging
import os
import stat
import sys
from functools import update_wrapper
from inspect import currentframe
from importlib import import_module
from typing import Iterable


class Xdict(dict):
    """Minimal legacy-compatible dict with lazy/key-error hooks."""

    def __init__(self, *args, **kwargs):
        self.not_found_method = kwargs.pop("not_found_method", None)
        self.not_found_message = kwargs.pop("not_found_message", None)
        super().__init__(*args, **kwargs)
        self.id = hash(id(self))

    def __getitem__(self, key):
        toget = self.get(key, self.id)
        if toget != self.id:
            return toget
        if self.not_found_method:
            value = self.not_found_method(key)
            self[key] = value
            return self[key]
        if self.not_found_message:
            raise KeyError(self.not_found_message.format(key))
        raise KeyError(f"{key}")


def xopen(filename, flag, mode=None, encoding=None):
    """Open files with legacy Xponge defaults."""

    if encoding is None and "b" not in flag:
        encoding = "UTF-8"
    if mode is None:
        mode = stat.S_IRWXO | stat.S_IRWXG | stat.S_IRWXU
    if flag in ("w", "wb"):
        real_flags = os.O_RDWR | os.O_CREAT | os.O_TRUNC
    elif flag in ("r", "rb"):
        real_flags = os.O_RDONLY
    else:
        raise NotImplementedError(flag)
    fd = os.open(filename, real_flags, mode)
    return os.fdopen(fd, flag, encoding=encoding)


def xprint(to_print, verbose=logging.INFO):
    """Legacy print helper.

    XpongeCPP does not currently carry the full logger surface from Xponge, so
    the compatibility path falls back to plain printing.
    """

    del verbose
    print(to_print)


Xprint = xprint
Xpri = xprint
Xopen = xopen


def generate_new_bonded_force_type(type_name, atoms, properties, is_compulsory, is_multiple=False):
    """Create a minimal legacy-compatible bonded force type class.

    This first-wave compatibility version focuses on the legacy decorators and
    registry shape needed by packaged force-field Python modules.
    """

    topology_like = [int(i) for i in str(atoms).split("-")]

    class BondedForceType:
        _name = str(type_name)
        topology_like = topology_like
        compulsory = bool(is_compulsory)
        multiple = bool(is_multiple)
        atom_numbers = len(topology_like)
        far = max(topology_like) if topology_like else 0
        _parameters = {"name": str, **dict(properties)}
        _types = Xdict(not_found_message="Bonded Force Type {} not found. Did you import the proper force field?")
        _types_different_name = Xdict(
            not_found_message="Bonded Force Type {} not found. Did you import the proper force field?"
        )

        @classmethod
        def Same_Force(cls, atom_list):
            if isinstance(atom_list, str):
                atom_list_temp = [atom.strip() for atom in atom_list.split("-")]
                return [atom_list, "-".join(atom_list_temp[::-1])]
            return [atom_list, atom_list[::-1]]

        @classmethod
        def Set_Same_Force_Function(cls, func):
            update_wrapper(func, cls.Same_Force)
            cls.Same_Force = classmethod(func)
            return func

        @classmethod
        def Get_Type_Name(cls, atoms_):
            return "-".join([atom.type.name for atom in atoms_])

        @classmethod
        def Type_Name_Getter(cls, func):
            update_wrapper(func, cls.Get_Type_Name)
            cls.Get_Type_Name = func
            return func

        @classmethod
        def Add_Property(cls, extra_properties):
            cls._parameters.update(dict(extra_properties))
            return cls

        @classmethod
        def New_From_String(cls, text):
            del text
            return cls

        @classmethod
        def New_From_Dict(cls, mapping):
            cls._types.update(dict(mapping))
            return cls

    from ..gromacs import GlobalSetting

    GlobalSetting.BondedForces.append(BondedForceType)
    GlobalSetting.BondedForcesMap[getattr(BondedForceType, "_name")] = BondedForceType
    if BondedForceType.far > GlobalSetting.farthest_bonded_force:
        GlobalSetting.farthest_bonded_force = BondedForceType.far
    return BondedForceType


def generate_new_pairwise_force_type(type_name, properties):
    """Create a minimal legacy-compatible pairwise force type class."""

    class PairwiseForceType:
        _name = str(type_name)
        _parameters = {"name": str, **dict(properties)}
        _types = Xdict(not_found_message="{} not found. Did you import the proper force field?")

        @classmethod
        def Add_Property(cls, extra_properties):
            cls._parameters.update(dict(extra_properties))
            return cls

    return PairwiseForceType


Generate_New_Bonded_Force_Type = generate_new_bonded_force_type
Generate_New_Pairwise_Force_Type = generate_new_pairwise_force_type


def debug(mode=True):
    """Compatibility debug toggle placeholder."""

    return bool(mode)


def set_global_alternative_names(*args, **kwargs):
    """Legacy namespace helper placeholder."""

    del args, kwargs
    return None


def get_main_namespace_injection_policy():
    """Return the current legacy ``__main__`` injection policy.

    The default remains strict legacy compatibility so old scripts keep working
    unchanged.  A dedicated helper makes the policy explicit and auditable, and
    also leaves room for future narrowing once script dependence is better
    understood.
    """

    raw = os.environ.get("XPONGECPP_LEGACY_MAIN_NAMESPACE", "strict")
    value = str(raw).strip().lower()
    if value in {"0", "false", "off", "disabled", "disable", "none"}:
        return "disabled"
    return "strict"


def should_install_main_namespace_exports():
    """Whether legacy names should be mirrored into ``__main__``."""

    return get_main_namespace_injection_policy() != "disabled"


def set_real_global_variable(name, value, namespace=None):
    """Install a legacy global symbol into a namespace."""

    if namespace is None:
        namespace = sys.modules.get("__main__").__dict__
    namespace[name] = value
    return value


def remove_real_global_variable(name, namespace=None):
    """Remove a legacy global symbol from a namespace if it exists."""

    if namespace is None:
        namespace = sys.modules.get("__main__").__dict__
    namespace.pop(name, None)
    return None


def source(module, into_global=True, reload_module=False):
    """Import a module relative to the caller and optionally merge globals."""

    caller_globals = currentframe().f_back.f_globals
    module_obj = import_module(module, package=caller_globals["__name__"])
    if reload_module:
        from importlib import reload

        module_obj = reload(module_obj)
    if into_global:
        for key, value in module_obj.__dict__.items():
            if not key.startswith("_"):
                caller_globals[key] = value
    return module_obj


def reexport_module(target: str, namespace: dict, public: Iterable[str] | None = None):
    """Populate *namespace* with public names from *target*."""

    module = import_module(target)
    if public is None:
        public = getattr(module, "__all__", None)
        if public is None:
            public = [name for name in dir(module) if not name.startswith("_")]
    public = list(public)
    for name in public:
        namespace[name] = getattr(module, name)
    namespace.setdefault("__all__", public)
    namespace.setdefault("__legacy_target__", target)
    return module


def copy_public_attributes(source, namespace: dict, *, skip: Iterable[str] | None = None):
    """Copy public attributes from *source* into *namespace*.

    This is used by legacy package-name shims such as ``src/Xponge/__init__.py``
    so that the shims do not each carry hand-written loops for public attribute
    forwarding.
    """

    skip = set(skip or ())
    for name in dir(source):
        if name.startswith("_") or name in skip:
            continue
        namespace.setdefault(name, getattr(source, name))
    return namespace


def install_main_namespace_exports(source_namespace: dict, target_namespace: dict | None = None):
    """Mirror public names from a shim namespace into ``__main__``.

    This preserves the existing bare-name legacy script behavior while
    centralizing the high-intrusion namespace injection logic in one place for
    later auditing or gating.
    """

    if not should_install_main_namespace_exports():
        return None
    if target_namespace is None:
        main_module = sys.modules.get("__main__")
        if main_module is None:  # pragma: no cover - defensive only
            return None
        target_namespace = main_module.__dict__
    for name, value in list(source_namespace.items()):
        if name.startswith("_"):
            continue
        target_namespace.setdefault(name, value)
    return target_namespace


def extend_package_path(namespace: dict, target_package: str):
    """Make a legacy package also search the target package's submodule path."""

    module = import_module(target_package)
    legacy_path = list(namespace.get("__path__", []))
    target_path = list(getattr(module, "__path__", []))
    for path in target_path:
        if path not in legacy_path:
            legacy_path.append(path)
    namespace["__path__"] = legacy_path
    namespace.setdefault("__legacy_target__", target_package)
    return module
