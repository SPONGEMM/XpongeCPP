"""SPONGE bundled-input discovery, validation, and read APIs."""

from .case import (
    BundleCase,
    bundle_case_from_prefix,
    parse_mdin_text,
    scan_bundle_case,
)
from .errors import (
    BundleCapabilityError,
    BundleConflictError,
    BundleError,
    BundleExportError,
    BundlePathError,
    BundleSchemaError,
    BundleValidationError,
)
from .manifest import ManifestEntry, ReverseConversionManifest
from .reader import BundleReader
from .reverse_converter import (
    BundleToLegacyConverter,
    convert_bundle_to_legacy,
)

__all__ = [
    "BundleCase",
    "BundleCapabilityError",
    "BundleConflictError",
    "BundleError",
    "BundleExportError",
    "BundlePathError",
    "BundleReader",
    "BundleSchemaError",
    "BundleValidationError",
    "BundleToLegacyConverter",
    "ManifestEntry",
    "ReverseConversionManifest",
    "bundle_case_from_prefix",
    "convert_bundle_to_legacy",
    "parse_mdin_text",
    "scan_bundle_case",
]
