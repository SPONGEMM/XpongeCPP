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
    AmbiguousH5MDLayoutError,
    BundleCapabilityError,
    BundleConflictError,
    BundleError,
    BundleExportError,
    BundleMDAnalysisError,
    BundlePathError,
    BundleSchemaError,
    BundleTopologyError,
    BundleTrajectoryError,
    BundleUnitError,
    BundleValidationError,
    IncompleteBundleError,
    UnverifiedBundlePairError,
)
from .legacy_case import LegacyCase, scan_legacy_case
from .manifest import (
    ConversionManifest,
    ManifestEntry,
    ReverseConversionManifest,
)
from .output_writer import (
    LegacyOutputBundleWriter,
    LegacyOutputConversionError,
    convert_legacy_outputs_to_bundle,
)
from .reader import BundleReader
from .protocol import (
    ProtocolCVRestraint,
    ProtocolCollectiveVariable,
    ProtocolDistanceConstraints,
    ProtocolHardWall,
    ProtocolMetadynamics,
    ProtocolPositionalRestraint,
    ProtocolSITS,
    ProtocolSoftWall,
    ProtocolSteering,
    ProtocolVirtualAtom,
    SpongeProtocol,
)
from .reverse_converter import (
    BundleToLegacyConverter,
    convert_bundle_to_legacy,
)
from .saver import save_sponge_input_bundle

__all__ = [
    "BundleCase",
    "AmbiguousH5MDLayoutError",
    "BundleCapabilityError",
    "BundleConflictError",
    "BundleError",
    "BundleExportError",
    "BundleMDAnalysisError",
    "BundlePathError",
    "BundleReader",
    "BundleSchemaError",
    "BundleTopologyError",
    "BundleTrajectoryError",
    "BundleUnitError",
    "BundleValidationError",
    "IncompleteBundleError",
    "BundleToLegacyConverter",
    "ConversionManifest",
    "LegacyCase",
    "LegacyOutputBundleWriter",
    "LegacyOutputConversionError",
    "LegacyToBundleConverter",
    "ManifestEntry",
    "ProtocolCVRestraint",
    "ProtocolCollectiveVariable",
    "ProtocolDistanceConstraints",
    "ProtocolHardWall",
    "ProtocolMetadynamics",
    "ProtocolPositionalRestraint",
    "ProtocolSITS",
    "ProtocolSoftWall",
    "ProtocolSteering",
    "ProtocolVirtualAtom",
    "ReverseConversionManifest",
    "SpongeProtocol",
    "UnverifiedBundlePairError",
    "bundle_case_from_prefix",
    "canonical_dataset_hash",
    "convert_bundle_to_legacy",
    "convert_legacy_to_bundle",
    "convert_legacy_outputs_to_bundle",
    "parse_mdin_text",
    "scan_bundle_case",
    "scan_legacy_case",
    "save_sponge_input_bundle",
]
