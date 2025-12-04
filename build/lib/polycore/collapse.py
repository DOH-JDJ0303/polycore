from typing import List, Tuple, Dict, Optional
import numpy as np, logging

def collapse_sequences(sequences: List[str], names: List[str]) -> Tuple[List[str], List[str], Dict[int, List[int]]]:
    """
    Collapse identical sequences.

    Returns:
        unique_sequences : list of unique sequence strings
        rep_names        : representative names
        idx_map          : dict {rep_idx -> [orig_indices]} mapping reps to their group
    """
    logger = logging.getLogger(__name__)

    seen = {}
    idx_map: Dict[int, List[int]] = {}
    unique_sequences = []
    rep_names = []

    # always keep the reference (index 0)
    unique_sequences.append(sequences[0])
    rep_names.append(names[0])
    seen[sequences[0]] = 0
    idx_map[0] = [0]

    for i, (seq, name) in enumerate(zip(sequences[1:], names[1:]), start=1):
        if seq in seen:
            rep_idx = seen[seq]
            idx_map[rep_idx].append(i)
        else:
            rep_idx = len(unique_sequences)
            seen[seq] = rep_idx
            unique_sequences.append(seq)
            rep_names.append(name)
            idx_map[rep_idx] = [i]

    # summary: report which originals collapsed into which rep
    for rep_idx, orig_idxs in idx_map.items():
        group_names = [names[j] for j in orig_idxs]
        if len(group_names) > 1:
            logger.info(f"Identical samples will be treated as one: {group_names} -> ")

    return unique_sequences, rep_names, idx_map

def expand_results(filtered_array: np.ndarray, filter_mask: np.ndarray, idx_map: Dict[int, List[int]], orig_names: List[str], keep_filtered: bool = True):
    """
    Expand filtered results back to original sample space.
    Works for both 1D and 2D arrays.
    """
    fi = 0  # Index into filtered_array
    expanded_rows = []
    expanded_names = []
    
    for mi, oi in idx_map.items():
        if filter_mask[mi]:
            val = filtered_array[fi]
            fi += 1
        else:
            if not keep_filtered:
                continue
            # Create appropriate null value based on array dimensions
            if filtered_array.ndim == 1:
                val = np.nan
            else:
                # For 2D, create a row of NaNs with same width as filtered_array
                val = np.full(filtered_array.shape[1], np.nan)
        
        # Expand to all original samples in this group
        for i in oi:
            expanded_rows.append(val)
            expanded_names.append(orig_names[i])
    
    return np.array(expanded_rows), expanded_names

def expand_distances(diffs, filter_mask, idx_map, orig_names):
    """
    Expand a rep-level distance matrix to original samples,
    but only for reps where filter_mask is True.
    """
    # 1) Ordered kept reps by rep index so fi/fj match rows/cols in `diffs`
    kept_reps = [r for r, keep in enumerate(filter_mask) if keep]

    # 2) Build expanded name list once (length == total originals among kept reps)
    expanded_names = [orig_names[i]
                      for rep in kept_reps
                      for i in idx_map[rep]]

    # 3) Allocate expanded distance matrix and fill by blocks
    sizes = [len(idx_map[rep]) for rep in kept_reps]
    n = sum(sizes)
    expanded = np.empty((n, n), dtype=diffs.dtype)

    r0 = 0
    for a, rep_a in enumerate(kept_reps):
        ga = idx_map[rep_a]        # original indices in this group
        ra = len(ga)
        c0 = 0
        for b, rep_b in enumerate(kept_reps):
            gb = idx_map[rep_b]
            rb = len(gb)
            expanded[r0:r0+ra, c0:c0+rb] = diffs[a, b]
            c0 += rb
        r0 += ra

    return expanded, expanded_names



def expand_vector(vector: np.ndarray, names: List[str], idx_map: Dict[int, List[int]]):
    expanded = []
    for rep_idx, val in enumerate(vector):
        for _ in idx_map[rep_idx]:
            expanded.append(val)
    return expanded