"""
Utilities for handling motifs, PWMs, and annotations.
"""

import json
import logging
import math
from pathlib import Path

import numpy as np

try:
    import numba

    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

    # Fallback if numba is not available
    def numba_njit(*args, **kwargs):
        def decorator(func):
            return func

        return decorator


logger = logging.getLogger(__name__)


# FIMO p-value calculation functions
if HAS_NUMBA:

    @numba.njit("float64(float64, float64)", fastmath=True, cache=True)
    def logaddexp2(x, y):
        """Calculate the logaddexp in a numerically stable manner in base 2."""
        if x == float("-inf") and y == float("-inf"):
            return float("-inf")
        if x == float("inf") or y == float("inf"):
            return float("inf")
        vmax, vmin = max(x, y), min(x, y)
        return vmax + math.log2(math.pow(2, vmin - vmax) + 1)

    @numba.njit(cache=True)
    def _pwm_to_mapping(log_pwm, bin_size):
        """Calculate score <-> log p-value mappings using NRDB background frequencies."""
        n, motif_len = log_pwm.shape
        # NRDB background frequencies
        log_bg_A = math.log2(0.281774)
        log_bg_C = math.log2(0.222020)
        log_bg_G = math.log2(0.228876)
        log_bg_T = math.log2(0.267330)
        int_log_pwm = np.round(log_pwm / bin_size).astype(np.int32)
        smallest, largest = 9999999, -9999999
        log_pwm_min_csum, log_pwm_max_csum = 0, 0
        for i in range(motif_len):
            log_pwm_min = 9999999
            log_pwm_max = -9999999
            for j in range(n):
                log_pwm_min = min(log_pwm_min, int_log_pwm[j, i])
                log_pwm_max = max(log_pwm_max, int_log_pwm[j, i])
            log_pwm_min_csum += log_pwm_min
            log_pwm_max_csum += log_pwm_max
            smallest = min(smallest, log_pwm_min_csum)
            largest = max(largest, log_pwm_max_csum)
        largest += motif_len
        logpdf_size = largest - smallest + 1
        logpdf = np.empty(logpdf_size, dtype=np.float64)
        old_logpdf = np.full(logpdf_size, -np.inf, dtype=np.float64)
        # Map nucleotide index to background log frequency
        log_bg_values = [log_bg_A, log_bg_C, log_bg_G, log_bg_T]
        for i in range(n):
            idx = int_log_pwm[i, 0] - smallest
            old_logpdf[idx] = logaddexp2(old_logpdf[idx], log_bg_values[i])
        for i in range(1, motif_len):
            for j in range(logpdf_size):
                logpdf[j] = -np.inf
            for j in range(logpdf_size):
                x = old_logpdf[j]
                if x != -np.inf:
                    for k in range(n):
                        idx = j + int_log_pwm[k, i]
                        if 0 <= idx < logpdf_size:
                            logpdf[idx] = logaddexp2(logpdf[idx], log_bg_values[k] + x)
            for j in range(logpdf_size):
                old_logpdf[j] = logpdf[j]
        for i in range(logpdf_size - 2, -1, -1):
            logpdf[i] = logaddexp2(logpdf[i], logpdf[i + 1])
        return smallest, logpdf
else:

    def logaddexp2(x, y):
        """Calculate the logaddexp in a numerically stable manner in base 2."""
        if x == float("-inf") and y == float("-inf"):
            return float("-inf")
        if x == float("inf") or y == float("inf"):
            return float("inf")
        vmax, vmin = max(x, y), min(x, y)
        return vmax + math.log2(math.pow(2, vmin - vmax) + 1)

    def _pwm_to_mapping(log_pwm, bin_size):
        """Calculate score <-> log p-value mappings using NRDB background frequencies."""
        n, motif_len = log_pwm.shape
        # NRDB background frequencies
        log_bg_A = math.log2(0.281774)
        log_bg_C = math.log2(0.222020)
        log_bg_G = math.log2(0.228876)
        log_bg_T = math.log2(0.267330)
        int_log_pwm = np.round(log_pwm / bin_size).astype(np.int32)
        smallest, largest = 9999999, -9999999
        log_pwm_min_csum, log_pwm_max_csum = 0, 0
        for i in range(motif_len):
            log_pwm_min = 9999999
            log_pwm_max = -9999999
            for j in range(n):
                log_pwm_min = min(log_pwm_min, int_log_pwm[j, i])
                log_pwm_max = max(log_pwm_max, int_log_pwm[j, i])
            log_pwm_min_csum += log_pwm_min
            log_pwm_max_csum += log_pwm_max
            smallest = min(smallest, log_pwm_min_csum)
            largest = max(largest, log_pwm_max_csum)
        largest += motif_len
        logpdf_size = largest - smallest + 1
        logpdf = np.empty(logpdf_size, dtype=np.float64)
        old_logpdf = np.full(logpdf_size, -np.inf, dtype=np.float64)
        # Map nucleotide index to background log frequency
        log_bg_values = [log_bg_A, log_bg_C, log_bg_G, log_bg_T]
        for i in range(n):
            idx = int_log_pwm[i, 0] - smallest
            old_logpdf[idx] = logaddexp2(old_logpdf[idx], log_bg_values[i])
        for i in range(1, motif_len):
            for j in range(logpdf_size):
                logpdf[j] = -np.inf
            for j in range(logpdf_size):
                x = old_logpdf[j]
                if x != -np.inf:
                    for k in range(n):
                        idx = j + int_log_pwm[k, i]
                        if 0 <= idx < logpdf_size:
                            logpdf[idx] = logaddexp2(logpdf[idx], log_bg_values[k] + x)
            for j in range(logpdf_size):
                old_logpdf[j] = logpdf[j]
        for i in range(logpdf_size - 2, -1, -1):
            logpdf[i] = logaddexp2(logpdf[i], logpdf[i + 1])
        return smallest, logpdf


def trim_motif_padding(motif_log_odds: np.ndarray) -> np.ndarray:
    """Remove zero-padding positions from motif log-odds to get actual motif length."""
    valid_positions = []
    for pos in range(motif_log_odds.shape[1]):
        column = motif_log_odds[:, pos]
        is_all_zero = np.allclose(column, 0.0, atol=1e-6)
        is_all_same_negative = (
            np.allclose(column, column[0], atol=1e-6) and column[0] < -3.0
        )
        if not (is_all_zero or is_all_same_negative):
            valid_positions.append(pos)
    if not valid_positions:
        return motif_log_odds[:, :0]
    start_idx = valid_positions[0]
    end_idx = valid_positions[-1] + 1
    return motif_log_odds[:, start_idx:end_idx]


def pvalue_from_logpdf(
    scores: np.ndarray, smallest: int, logpdf: np.ndarray, bin_size: float
) -> np.ndarray:
    """
    Vectorized conversion of raw scores to p-values using a precomputed lookup table.
    """
    score_indices = (scores / bin_size).astype(int) - smallest
    # Clip indices to be within the bounds of the logpdf table
    score_indices = np.clip(score_indices, 0, len(logpdf) - 1)
    log_pvalues = logpdf[score_indices]
    return 2.0**log_pvalues


def compute_pvalue_mapping(
    pwm: np.ndarray, bin_size: float = 0.01
) -> tuple[int, np.ndarray, float]:
    """
    Compute score-to-pvalue mapping for a PWM.

    Args:
        pwm: Position weight matrix with shape (length, 4) in probability space
        bin_size: Bin size for discretizing scores (default: 0.01)

    Returns:
        Tuple of (smallest_score_index, log_pdf_array, bin_size)
    """
    # PWMs are stored as (length, 4), need to transpose to (4, length) for computation
    if pwm.shape[1] == 4:
        pwm = pwm.T  # Convert (length, 4) to (4, length)
    elif pwm.shape[0] == 4:
        # Already in (4, length) format
        pass
    else:
        raise ValueError(
            f"PWM must have shape (length, 4) or (4, length), got {pwm.shape}"
        )

    # Check if PWM is already in log-odds space (negative values common)
    # or probability space (all positive, typically sums to 1 per column)
    col_sums = pwm.sum(axis=0)
    is_probability = np.allclose(col_sums, 1.0, atol=0.1) and np.all(pwm >= -1e-6)

    if is_probability:
        # Convert probability PWM to log-odds space
        # Add small epsilon to avoid log(0)
        epsilon = 1e-10
        pwm_normalized = pwm + epsilon
        pwm_normalized = pwm_normalized / pwm_normalized.sum(axis=0, keepdims=True)

        # Convert to log2 space (log-odds) with background frequencies
        # Background: A=0.281774, C=0.222020, G=0.228876, T=0.267330
        bg_freqs = np.array([0.281774, 0.222020, 0.228876, 0.267330])[:, np.newaxis]
        log_pwm = np.log2(pwm_normalized / bg_freqs)
    else:
        # Assume already in log-odds space, but check if it's log10 or log2
        # If values are very negative (like -10 to -20), likely log10
        # If values are moderate (like -5 to -10), likely log2
        max_abs_val = np.abs(pwm).max()
        log_pwm = pwm * np.log10(2.0) if max_abs_val > 15 else pwm

    # Compute mapping
    smallest, logpdf = _pwm_to_mapping(log_pwm, bin_size)

    return smallest, logpdf, bin_size


def read_pwms(pwm_file):
    """Read HOCOMOCO PWM file and return motifs with their PWMs."""
    motifs = {}
    current_motif = None
    pwm_lines = []

    with Path(pwm_file).open() as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_motif and pwm_lines:
                    pwm = np.array(
                        [[float(x) for x in row.split("\t")] for row in pwm_lines]
                    )
                    motifs[current_motif] = pwm
                current_motif = line[1:]
                pwm_lines = []
            elif line and current_motif:
                pwm_lines.append(line)

    if current_motif and pwm_lines:
        pwm = np.array([[float(x) for x in row.split("\t")] for row in pwm_lines])
        motifs[current_motif] = pwm

    return motifs


def read_annotations(annotation_file):
    """Read HOCOMOCO annotation file and extract thresholds."""
    annotations = {}

    with Path(annotation_file).open() as f:
        for line in f:
            data = json.loads(line.strip())
            name = data["name"]
            thresholds = data["standard_thresholds"]
            annotations[name] = thresholds

    return annotations


def pad_motifs_to_same_length(motifs_dict, annotations_dict):
    """Pad all motifs to the same length for batch processing."""
    max_length = max(pwm.shape[0] for pwm in motifs_dict.values())

    motif_names = []
    padded_kernels = []
    original_lengths = []
    thresholds = {"p0.001": [], "p0.0005": [], "p0.0001": []}

    pad_value = 0.0

    for motif_name, pwm in motifs_dict.items():
        if motif_name not in annotations_dict:
            continue

        orig_len = pwm.shape[0]
        original_lengths.append(orig_len)

        total_pad = max_length - orig_len
        pad_left = total_pad // 2
        pad_right = total_pad - pad_left

        padded_pwm = np.pad(
            pwm, ((pad_left, pad_right), (0, 0)), "constant", constant_values=pad_value
        )

        padded_kernels.append(padded_pwm)
        motif_names.append(motif_name)

        for p_val in ["0.001", "0.0005", "0.0001"]:
            thresholds[f"p{p_val}"].append(annotations_dict[motif_name][p_val])

    pad_offsets = np.array(
        [(max_length - orig_len) // 2 for orig_len in original_lengths]
    )
    original_lengths = np.array(original_lengths)

    return (
        motif_names,
        np.array(padded_kernels),
        original_lengths,
        thresholds,
        max_length,
        pad_offsets,
    )


def filter_motifs(motifs_dict, annotations_dict, motif_selection=None):
    """
    Filter motifs based on selection criteria.

    Supports multiple selection formats:
    - Index ranges: "0-10" selects motifs 0 through 10
    - Comma-separated indices: "0,1,5,10" selects specific indices
    - Exact motif names: "CTCF.H13CORE.0.P.A"
    - Prefix matching: "CTCF" matches all motifs starting with "CTCF."
    - Mixed formats: "CTCF,TP53,0-5" combines prefix and index selection

    Args:
        motifs_dict (dict): Dictionary of motif name -> PWM matrix
        annotations_dict (dict): Dictionary of motif name -> annotation data
        motif_selection (str): Selection criteria string

    Returns:
        tuple: (filtered_motifs_dict, filtered_annotations_dict)
    """
    if motif_selection is None:
        return motifs_dict, annotations_dict

    all_motif_names = list(motifs_dict.keys())
    logger.info(f"Total available motifs: {len(all_motif_names)}")

    # Parse selection criteria
    if "-" in motif_selection and "," not in motif_selection:
        # Handle simple range like "0-10"
        start, end = motif_selection.split("-")
        selections = [f"{i}" for i in range(int(start), int(end) + 1)]
    else:
        # Handle comma-separated list that may include ranges
        selections = [s.strip() for s in motif_selection.split(",")]

    selected_motifs = []

    for selection in selections:
        # Handle range within comma-separated list
        if "-" in selection:
            try:
                start, end = selection.split("-")
                range_indices = [f"{i}" for i in range(int(start), int(end) + 1)]
                for idx_str in range_indices:
                    idx = int(idx_str)
                    if 0 <= idx < len(all_motif_names):
                        selected_motifs.append(all_motif_names[idx])
                        logger.info(f"Selected motif {idx}: {all_motif_names[idx]}")
                    else:
                        logger.warning(
                            f"Index {idx} out of range (0-{len(all_motif_names)-1})"
                        )
                continue
            except ValueError:
                logger.warning(f"Invalid range format: '{selection}'")
                continue

        # Handle numeric index
        try:
            idx = int(selection)
            if 0 <= idx < len(all_motif_names):
                selected_motifs.append(all_motif_names[idx])
                logger.info(f"Selected motif {idx}: {all_motif_names[idx]}")
            else:
                logger.warning(f"Index {idx} out of range (0-{len(all_motif_names)-1})")
            continue
        except ValueError:
            pass

        # Handle exact motif name match
        if selection in motifs_dict:
            selected_motifs.append(selection)
            logger.info(f"Selected motif by exact name: {selection}")
            continue

        # Handle prefix matching
        prefix_matches = []
        for motif_name in all_motif_names:
            if motif_name.startswith(selection + ".") or motif_name.startswith(
                selection + "_"
            ):
                prefix_matches.append(motif_name)

        if prefix_matches:
            selected_motifs.extend(prefix_matches)
            logger.info(
                f"Selected {len(prefix_matches)} motifs by prefix '{selection}': {prefix_matches}"
            )
            continue

        # Handle partial name matching (case-insensitive)
        partial_matches = []
        selection_upper = selection.upper()
        for motif_name in all_motif_names:
            if selection_upper in motif_name.upper():
                partial_matches.append(motif_name)

        if partial_matches:
            selected_motifs.extend(partial_matches)
            logger.info(
                f"Selected {len(partial_matches)} motifs by partial match '{selection}': {partial_matches}"
            )
            continue

        logger.warning(f"No motifs found matching '{selection}'")

    # Remove duplicates while preserving order
    selected_motifs = list(dict.fromkeys(selected_motifs))

    if not selected_motifs:
        logger.error("No valid motifs selected!")
        return motifs_dict, annotations_dict

    # Create filtered dictionaries
    filtered_motifs = {
        name: motifs_dict[name] for name in selected_motifs if name in motifs_dict
    }
    filtered_annotations = {
        name: annotations_dict[name]
        for name in selected_motifs
        if name in annotations_dict
    }

    logger.info(
        f"Filtered to {len(filtered_motifs)} motifs from {len(motifs_dict)} total"
    )
    if len(filtered_motifs) <= 10:
        logger.info(f"Selected motifs: {list(filtered_motifs.keys())}")
    else:
        logger.info(
            f"Selected motifs: {list(filtered_motifs.keys())[:5]} ... {list(filtered_motifs.keys())[-5:]}"
        )

    return filtered_motifs, filtered_annotations


def resolve_target_motif(target_motif_str, motif_names):
    """
    Resolve target motif string to motif index with advanced matching.

    Supports multiple matching strategies:
    - Exact name matching
    - Prefix matching (e.g., "FOXA1" matches "FOXA1.H13CORE.0.P.B")
    - Partial name matching (case-insensitive)

    Args:
        target_motif_str: Either motif name, partial name, or index as string
        motif_names: List of motif names in order

    Returns:
        int: Index of the target motif

    Raises:
        ValueError: If motif cannot be resolved
    """
    # Try to parse as integer index first
    try:
        idx = int(target_motif_str)
        if 0 <= idx < len(motif_names):
            logger.info(f"Target motif resolved by index {idx}: {motif_names[idx]}")
            return idx
        else:
            raise ValueError(
                f"Target motif index {idx} out of range (0-{len(motif_names)-1})"
            )
    except ValueError as e:
        if "out of range" in str(e):
            raise e

    # Try exact name match
    if target_motif_str in motif_names:
        idx = motif_names.index(target_motif_str)
        logger.info(
            f"Target motif resolved by exact name '{target_motif_str}': index {idx}"
        )
        return idx

    # Try prefix matching with preference for forward strand (non-RC)
    prefix_matches = []
    for i, name in enumerate(motif_names):
        if name.startswith(target_motif_str + ".") or name.startswith(
            target_motif_str + "_"
        ):
            prefix_matches.append((i, name))

    if len(prefix_matches) == 1:
        idx, name = prefix_matches[0]
        logger.info(
            f"Target motif resolved by prefix '{target_motif_str}': {name} (index {idx})"
        )
        return idx
    elif len(prefix_matches) > 1:
        # Prefer non-reverse complement matches
        non_rc_matches = [
            (i, name) for i, name in prefix_matches if not name.endswith("_RC")
        ]
        if len(non_rc_matches) == 1:
            idx, name = non_rc_matches[0]
            logger.info(
                f"Target motif resolved by prefix '{target_motif_str}': {name} (index {idx})"
            )
            return idx
        elif len(non_rc_matches) > 1:
            # If multiple non-RC matches, take the first one
            idx, name = non_rc_matches[0]
            logger.info(
                f"Target motif resolved by prefix '{target_motif_str}' (first of {len(non_rc_matches)} matches): {name} (index {idx})"
            )
            return idx
        else:
            # Only RC matches available, take the first one
            idx, name = prefix_matches[0]
            logger.info(
                f"Target motif resolved by prefix '{target_motif_str}' (RC match): {name} (index {idx})"
            )
            return idx

    # Try partial name match (case-insensitive)
    partial_matches = []
    target_motif_upper = target_motif_str.upper()
    for i, name in enumerate(motif_names):
        if target_motif_upper in name.upper():
            partial_matches.append((i, name))

    if len(partial_matches) == 1:
        idx, name = partial_matches[0]
        logger.info(
            f"Target motif resolved by partial match '{target_motif_str}': {name} (index {idx})"
        )
        return idx
    elif len(partial_matches) > 1:
        # Prefer non-reverse complement matches
        non_rc_matches = [
            (i, name) for i, name in partial_matches if not name.endswith("_RC")
        ]
        if len(non_rc_matches) == 1:
            idx, name = non_rc_matches[0]
            logger.info(
                f"Target motif resolved by partial match '{target_motif_str}': {name} (index {idx})"
            )
            return idx
        elif len(non_rc_matches) > 1:
            # Multiple matches - provide helpful error message
            match_names = [name for i, name in non_rc_matches]
            raise ValueError(
                f"Ambiguous target motif '{target_motif_str}'. Multiple non-RC matches found: {match_names[:5]}"
            )
        else:
            # Only RC matches
            match_names = [name for i, name in partial_matches]
            raise ValueError(
                f"Ambiguous target motif '{target_motif_str}'. Only RC matches found: {match_names[:5]}"
            )

    # No matches found
    raise ValueError(
        f"Target motif '{target_motif_str}' not found in {len(motif_names)} available motifs. "
        f"Try exact name, prefix (e.g., 'FOXA1' for 'FOXA1.H13CORE.0.P.B'), or motif index."
    )


def create_reverse_complement_motifs(motifs_dict, annotations_dict):
    """
    Create reverse complement motifs for strand-specific scanning.

    Args:
        motifs_dict: Dictionary of motif_name -> pwm
        annotations_dict: Dictionary of motif_name -> thresholds

    Returns:
        Extended dictionaries with reverse complement motifs
    """
    extended_motifs = motifs_dict.copy()
    extended_annotations = annotations_dict.copy()

    for motif_name, pwm in motifs_dict.items():
        if motif_name in annotations_dict:
            rc_pwm = np.flip(pwm, axis=0)
            rc_pwm = rc_pwm[:, [3, 2, 1, 0]]
            rc_name = f"{motif_name}_RC"
            extended_motifs[rc_name] = rc_pwm
            extended_annotations[rc_name] = annotations_dict[motif_name]

    logger.info(
        f"Created reverse complement motifs: {len(motifs_dict)} -> {len(extended_motifs)} total"
    )

    return extended_motifs, extended_annotations
