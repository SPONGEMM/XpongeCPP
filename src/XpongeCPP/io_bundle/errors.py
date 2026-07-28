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
