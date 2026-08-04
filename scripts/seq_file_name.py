"""Parse and display metadata encoded in sequence filenames."""

import re
from pathlib import Path


ERROR_RATES_RE = re.compile(
    r"\.txt\.gz_(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)\.txt\.gz_"
)


def get_error_rates(seq_file: Path) -> tuple[float, float] | None:
    match = ERROR_RATES_RE.search(seq_file.name)
    if match is None:
        return None
    return float(match.group(1)), float(match.group(2))


def get_seq_label(seq_file: Path) -> str:
    rates = get_error_rates(seq_file)
    if rates is None:
        return "Baseline (no injected errors)"
    bit_flip, switch = rates
    return f"Bit flip {bit_flip:g}, switch {switch:g}"


def get_seq_sort_key(seq_file: Path) -> tuple[bool, float, float]:
    rates = get_error_rates(seq_file)
    if rates is None:
        return False, 0.0, 0.0
    return True, *rates
