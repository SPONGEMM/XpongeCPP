"""Errors raised while reading and validating bundled SPONGE artifacts."""


class BundleError(RuntimeError):
    """Base class for bundled input failures."""


class BundleSchemaError(BundleError):
    """Raised when an HDF5 artifact has an unsupported schema."""


class BundleValidationError(BundleError):
    """Raised when bundled artifacts are internally inconsistent."""


class BundlePathError(BundleError):
    """Raised when a bundle path escapes its allowed root."""


class BundleCapabilityError(BundleError):
    """Raised when a bundle cannot be represented by legacy inputs."""


class BundleConflictError(BundleError):
    """Raised when conversion would overwrite an existing output."""


class BundleExportError(BundleError):
    """Raised when typed bundle data cannot be exported safely."""


class BundleMDAnalysisError(BundleValidationError):
    """Base class for bundle-to-MDAnalysis integration failures."""


class BundleTopologyError(BundleMDAnalysisError):
    """Raised when a bundled topology cannot be exposed to MDAnalysis."""


class BundleTrajectoryError(BundleMDAnalysisError):
    """Raised when a SPONGE H5MD trajectory is invalid or unsupported."""


class BundleUnitError(BundleMDAnalysisError):
    """Raised when bundle units cannot be converted to MDAnalysis units."""


class UnverifiedBundlePairError(BundleMDAnalysisError):
    """Raised when topology/trajectory compatibility cannot be proven."""


class IncompleteBundleError(BundleMDAnalysisError):
    """Raised when a trajectory bundle is not finalized or is truncated."""


class AmbiguousH5MDLayoutError(BundleTrajectoryError):
    """Raised when an H5MD file matches conflicting SPONGE layouts."""
