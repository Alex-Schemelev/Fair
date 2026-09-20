"""Shared helpers for in silico PCR web app."""

from .paths import (
    REF_FILENAME,
    PRIMERS_FILENAME,
    CONSERVATION_FILENAME,
    LEGACY_REF_NAMES,
    LEGACY_PRIMERS_NAMES,
    LEGACY_CONS_NAMES,
    find_existing,
)
from .iupac import (
    IUPAC_FROM_BASES,
    consensus_iupac_from_row,
    IUPAC_COMP,
    IUPAC_SET,
    iupac_complement,
    mismatch_stats_vs_consensus,
)
from .align import simple_align
from .conservation import compute_conservation_for_primers, compute_consensus_qc_for_products
from .filters import apply_product_filters
from .annotation import (
    rebuild_annotation_products,
    stable_product_id,
    stable_primer_feature_id,
    extract_annotation_map,
)
from .session_io import (
    load_dimers_matrix_from_run,
    is_path_within_directory,
    safe_extract_zip,
    slim_session_payload,
    cleanup_old_runs,
)

__all__ = [
    "REF_FILENAME",
    "PRIMERS_FILENAME",
    "CONSERVATION_FILENAME",
    "LEGACY_REF_NAMES",
    "LEGACY_PRIMERS_NAMES",
    "LEGACY_CONS_NAMES",
    "find_existing",
    "IUPAC_FROM_BASES",
    "consensus_iupac_from_row",
    "IUPAC_COMP",
    "IUPAC_SET",
    "iupac_complement",
    "mismatch_stats_vs_consensus",
    "simple_align",
    "compute_conservation_for_primers",
    "compute_consensus_qc_for_products",
    "apply_product_filters",
    "rebuild_annotation_products",
    "stable_product_id",
    "stable_primer_feature_id",
    "extract_annotation_map",
    "load_dimers_matrix_from_run",
    "is_path_within_directory",
    "safe_extract_zip",
    "slim_session_payload",
    "cleanup_old_runs",
]
