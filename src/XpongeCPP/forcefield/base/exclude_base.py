"""Legacy exclude_base import shim."""


class Exclude:
    """Minimal compatibility shape for legacy exclusion configuration."""

    current = None

    def __init__(self, n=4):
        self.n = int(n)
        Exclude.current = self


__all__ = ["Exclude"]
