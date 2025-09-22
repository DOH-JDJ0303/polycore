from typing import List, Tuple, Dict, Optional
import numpy as np, logging

def collapse_sequences(sequences: List[str], names: List[str]):
    """
    Collapse identical sequences.
    Returns:
        unique_sequences: list of unique sequence strings
        rep_names: representative names
        idx_map: dict {rep_idx -> [orig_indices]} mapping reps to their group
    """
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

    return unique_sequences, rep_names, idx_map

def expand_results(array: np.ndarray, names: List[str], idx_map: Dict[int, List[int]], orig_names: List[str]):
    expanded_array = []
    expanded_names = []
    for rep_idx, rep in enumerate(names):
        for orig_idx in idx_map[rep_idx]:
            expanded_array.append(array[rep_idx])
            expanded_names.append(orig_names[orig_idx])
    return np.vstack(expanded_array), expanded_names


def expand_distances(diffs: np.ndarray, names: List[str], idx_map: Dict[int, List[int]], orig_names: List[str]):
    expanded_names = []
    for rep_idx, rep in enumerate(names):
        for orig_idx in idx_map[rep_idx]:
            expanded_names.append(orig_names[orig_idx])

    n = len(expanded_names)
    expanded_diffs = np.zeros((n, n), dtype=int)

    row_map = {}
    k = 0
    for rep_idx, rep in enumerate(names):
        for orig_idx in idx_map[rep_idx]:
            row_map[orig_idx] = k
            k += 1

    for i, rep_i in enumerate(names):
        for j, rep_j in enumerate(names):
            d = diffs[i, j]
            for orig_i in idx_map[i]:
                for orig_j in idx_map[j]:
                    expanded_diffs[row_map[orig_i], row_map[orig_j]] = d

    return expanded_diffs, expanded_names


def expand_vector(vector: np.ndarray, names: List[str], idx_map: Dict[int, List[int]]):
    expanded = []
    for rep_idx, val in enumerate(vector):
        for _ in idx_map[rep_idx]:
            expanded.append(val)
    return expanded