"""SPONGE bundled-input discovery, validation, and read APIs."""

from .case import (
    BundleCase,
    bundle_case_from_prefix,
    parse_mdin_text,
    scan_bundle_case,
)
from .bundle_builder import canonical_dataset_hash
from .converter import LegacyToBundleConverter, convert_legacy_to_bundle
from .errors import (
    BundleCapabilityError,
    BundleConflictError,
    BundleError,
    BundleExportError,
    BundlePathError,
    BundleSchemaError,
    BundleValidationError,
)
from .legacy_case import LegacyCase, scan_legacy_case
from .manifest import (
    ConversionManifest,
    ManifestEntry,
    ReverseConversionManifest,
)
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
    "ConversionManifest",
    "LegacyCase",
    "LegacyToBundleConverter",
    "ManifestEntry",
    "ReverseConversionManifest",
    "bundle_case_from_prefix",
    "canonical_dataset_hash",
    "convert_bundle_to_legacy",
    "convert_legacy_to_bundle",
    "parse_mdin_text",
    "scan_bundle_case",
    "scan_legacy_case",
]
