"""
Utilities for handling motifs, PWMs, and annotations.
"""

import json
import logging
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


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
