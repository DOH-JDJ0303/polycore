from typing import List, Tuple, Dict, Optional
import numpy as np, logging

def collapse_sequences(sequences: List[str], names: List[str]) -> Tuple[List[str], List[str], Dict[int, List[int]]]:
    """
    Collapse identical sequences.

    Returns:
        unique_sequences : list of unique sequence strings
        idx_map          : dict {rep_idx -> [orig_indices]} mapping reps to their group
    """
    logger = logging.getLogger(__name__)

    seen = {}
    idx_map: Dict[int, List[int]] = {}
    unique_sequences = []

    # always keep the reference (index 0)
    unique_sequences.append(sequences[0])
    seen[''.join(sequences[0])] = 0
    idx_map[0] = [0]

    for i, (seq_list, name) in enumerate(zip(sequences[1:], names[1:]), start=1):
        seq = ''.join(seq_list)
        if seq in seen:
            rep_idx = seen[seq]
            idx_map[rep_idx].append(i)
        else:
            rep_idx = len(unique_sequences)
            seen[seq] = rep_idx
            unique_sequences.append(seq)
            idx_map[rep_idx] = [i]

    # summary: report which originals collapsed into which rep
    for rep_idx, orig_idxs in idx_map.items():
        group_names = [names[j] for j in orig_idxs]
        if len(group_names) > 1:
            logger.info(f"Identical samples will be treated as one: {group_names} -> ")

    return unique_sequences, idx_map

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

def expand_stacks(names_orig, names_filt, idx_maps, stacks, ref_masks, core_masks, const_masks):
    full_stacks = {}
    core_stacks = {}
    for i, n in enumerate(names_orig):
        if n not in names_filt:
            continue
        idx_map = []
        for m in idx_maps:
            for k, v in m.items():
                if i in v:
                    idx_map.append(k)

        # full stack masked
        full_stacks_by_name = []
        core_stacks_by_name = []
        for j, stack in enumerate(stacks):
            idx = idx_map[j]

            ref_mask   = ref_masks[j]
            core_mask  = core_masks[j]
            const_mask = const_masks[j]
            full_stack = stack[idx].copy()

            full_stack[(ref_mask & core_mask)] = 'N'
            core_stack = full_stack[(~const_mask & ~ref_mask)]
            full_stacks_by_name.append(full_stack)
            core_stacks_by_name.append(core_stack)
        full_stacks[n] = full_stacks_by_name
        core_stacks[n] = core_stacks_by_name

    return full_stacks, core_stacks