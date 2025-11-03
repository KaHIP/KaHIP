"""
KaHIP - Karlsruhe High Quality Graph Partitioning Framework

Python bindings for the KaHIP graph partitioning library.
"""

__version__ = "3.19.0"

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
