"""Canonical hash and HDF5 helpers for SPONGE v2 input bundles."""

from __future__ import annotations

import hashlib
from typing import Any

import numpy as np


_CANONICAL_NUMERIC_DTYPES = {
    "float32",
    "float64",
    "int8",
    "uint8",
    "int32",
    "uint32",
    "int64",
    "uint64",
}


def canonical_dataset_hash(
    bundle_file: str,
    datasets: dict[str, Any],
    *,
    path_prefixes: tuple[str, ...] = (),
) -> str:
    """Hash logical datasets using SPONGE v2 canonical encoding."""

    digest = hashlib.sha256(bundle_file.encode("utf-8"))
    for dataset_path, value in sorted(datasets.items()):
        if path_prefixes and not dataset_path.startswith(path_prefixes):
            continue
        array = np.asarray(value)
        is_string = array.dtype.kind in {"O", "U", "S"}
        dtype_name = "object" if is_string else array.dtype.name
        if dtype_name == "bool":
            array = array.astype(np.uint8, copy=False)
            dtype_name = "uint8"
        if not is_string and dtype_name not in _CANONICAL_NUMERIC_DTYPES:
            raise TypeError(
                f"unsupported canonical dtype {array.dtype} at {dataset_path}"
            )
        digest.update(b"\0")
        digest.update(dataset_path.encode("utf-8"))
        digest.update(b"\0")
        digest.update(dtype_name.encode("ascii"))
        digest.update(b"\0")
        digest.update(repr(array.shape).encode("ascii"))
        digest.update(b"\0")
        if is_string:
            values = (
                item.decode("utf-8")
                if isinstance(item, (bytes, np.bytes_))
                else str(item)
                for item in array.reshape(-1)
            )
            digest.update("\0".join(values).encode("utf-8"))
        else:
            native_dtype = np.dtype(dtype_name)
            digest.update(
                np.ascontiguousarray(array, dtype=native_dtype).tobytes()
            )
    return "sha256:" + digest.hexdigest()


def write_dataset(handle, path: str, value, *, string: bool = False) -> None:
    """Create one dataset and its parent groups."""

    import h5py

    parent, _, name = path.rpartition("/")
    group = handle.require_group(parent or "/")
    data = np.asarray(value, dtype=object) if string else value
    dtype = h5py.string_dtype("utf-8") if string else None
    group.create_dataset(name, data=data, dtype=dtype)


def write_string(handle, path: str, value: str) -> None:
    write_dataset(handle, path, value, string=True)
