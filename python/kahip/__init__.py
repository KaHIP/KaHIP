"""
KaHIP - Karlsruhe High Quality Graph Partitioning Framework

Python bindings for the KaHIP graph partitioning library.
"""

try:
    from ._version import __version__
except ImportError:
    # Fallback version for development installs without setuptools-scm
    __version__ = "0.0.0+unknown"

try:
    from .kahip import kaffpa
except ImportError as e:
    raise ImportError(
        "Failed to import the KaHIP C++ extension. "
        "Make sure the package is properly installed."
    ) from e

__all__ = ["kaffpa", "__version__"]


# Partitioning mode constants
FAST = 0
ECO = 1
STRONG = 2
FASTSOCIAL = 3
ECOSOCIAL = 4
STRONGSOCIAL = 5
