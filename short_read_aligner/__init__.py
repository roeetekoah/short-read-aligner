"""A small educational ungapped short-read aligner."""

from .core import AlignContext, align_all, align_read, build_index, build_minimizer_index

__version__ = "0.1.0"
__all__ = ["AlignContext", "align_all", "align_read", "build_index", "build_minimizer_index"]
