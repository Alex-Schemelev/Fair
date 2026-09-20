import os

REF_FILENAME = "reference.gb"
PRIMERS_FILENAME = "primers.csv"
CONSERVATION_FILENAME = "conservation.tsv"
LEGACY_REF_NAMES = ("reference.gb", "K03455.1.gb")
LEGACY_PRIMERS_NAMES = ("primers.csv", "HIV_primers_table.csv")
LEGACY_CONS_NAMES = ("conservation.tsv", "S_results_table_nt.txt")


def find_existing(run_dir, candidates):
    """Return first existing path under run_dir among candidates, or None."""
    for name in candidates:
        path = os.path.join(run_dir, name)
        if os.path.isfile(path):
            return path
    return None
