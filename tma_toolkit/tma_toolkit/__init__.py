"""tma_toolkit -- TMA core metadata, selection and comparison for Xenium data (AnnData)."""

from .cohort import CoreSelection, TMACohort
from .metadata import attach_core_metadata, make_template, read_core_metadata, validate_core_metadata

__all__ = ["TMACohort", "CoreSelection", "make_template", "read_core_metadata",
           "validate_core_metadata", "attach_core_metadata"]
