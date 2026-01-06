import numpy as np, logging
from typing import List, Tuple, Dict, Optional
from collections import Counter
from .io_ops import create_plot, write_distances
from .distance import to_bits, calculate_distances
from .utils import auto_chunk_size
from .collapse import expand_distances

def read_bed(path):
    """
    Read the first 3 columns of a BED file.
    Returns a list of (chrom, start, end) tuples.
    """
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end = line.split()[:3]
            out.append((chrom, int(start), int(end)))
    return out

def apply_bed_mask(stack, bed):
    # Mask sites specified in BED (set to 'N' in all sequences)
    mask = None
    if bed is not None:
        mask = read_bed(bed)   # [(chrom, start, end), ...]
        logging.info(f"Loaded {len(bed_mask)} mask intervals from {bed}")
        if bed_mask:
            n_sites = stack.shape[1]
            mask_cols = np.zeros(n_sites, dtype=bool)

            for _, start, end in bed_mask:
                # Clip to alignment bounds, and be robust to non-int inputs
                start = max(0, int(start))
                end   = min(n_sites, int(end))
                if start < end:
                    mask_cols[start:end] = True

            n_masked = int(mask_cols.sum())
            if n_masked > 0:
                stack[:, mask_cols] = 'N'
                logging.info(f"Masked {n_masked} sites from BED regions")
                
                # Check if all sites are now masked
                if n_masked == n_sites:
                    raise ValueError("All sites were masked by BED regions")


import numpy as np
import logging
from collections import Counter

def find_core(stacks, idx_maps, ref_masks, names_orig, names_filt, gfs,
              threshold=1.0, progressive=True):
    """
    For each stack/contig:
      - Build a progressive order of originals (core_order) by gf, keeping ref first.
      - Convert that into a rep_order (rep indices, with duplicates) using idx_map.
      - Compute per-site weighted non-N fraction across the first t originals (t grows).
      - core_mask is True where fraction < threshold (i.e., NOT core).
      - Track progress stats in cfs_prog: {orig_index: [n_core_sites, n_valid_sites]} aggregated across contigs.
      - Return final core_mask per contig at the point where the last filtered sample is included
        (or the full progressive end if progressive=True and no filtered endpoint exists).
    """
    # originals sorted by gf (ref first)
    core_order = [0] + sorted(range(1, len(gfs)), key=lambda j: gfs[j], reverse=True)

    names_filt_set = set(names_filt)

    # store counts as mutable lists so we can sum
    cfs_prog = {}
    core_masks = []

    logging.info(f"Determining core genome (progressive={progressive})")

    for contig_i, stack in enumerate(stacks):
        idx_map = idx_maps[contig_i]     # rep -> list[orig_idx]
        ref_mask = ref_masks[contig_i]   # True where ref is invalid ('N') or masked

        # Build orig -> rep lookup once (fast)
        orig_to_rep = {}
        for rep, origs in idx_map.items():
            for o in origs:
                orig_to_rep[o] = rep

        # Build rep_order in the same sequence as core_order, including duplicates
        rep_order = []
        final_step = None  # step t (1-based) where we have included the last "filtered" original seen so far

        for t, orig_idx in enumerate(core_order, start=1):
            rep = orig_to_rep.get(orig_idx)
            if rep is None:
                continue
            rep_order.append(rep)

            # mark the latest step where a kept name is included
            if names_orig[orig_idx] in names_filt_set:
                final_step = t

        if not rep_order:
            raise ValueError(f"Contig {contig_i}: rep_order is empty (idx_map/core_order mismatch).")

        # If no filtered samples matched, fall back to full set
        if final_step is None:
            final_step = len(rep_order)

        # Decide which steps to evaluate
        if progressive:
            start, stop = 1, len(rep_order) + 1
        else:
            start, stop = final_step, final_step + 1

        rep_set = None
        valid_matrix = None  # (n_unique_reps, n_cols) int
        final_core_mask = None

        for t in range(start, stop):  # t = number of originals included so far
            reps = rep_order[:t]              # includes duplicates
            rep_counts = Counter(reps)

            # order-preserving unique reps
            rep_set_new = list(dict.fromkeys(reps))

            if rep_set is None or rep_set_new != rep_set:
                rep_set = rep_set_new
                # 1 where base is valid (not N), 0 otherwise
                valid_matrix = (stack[rep_set, :] != "N").astype(np.int32)

            weights = np.array([rep_counts[r] for r in rep_set], dtype=np.int32)  # (n_unique_reps,)

            # weighted fraction per site across t originals
            fractions = (valid_matrix * weights[:, None]).sum(axis=0) / t

            # core_mask True = NOT core (below threshold)
            core_mask = fractions < threshold

            # stats excluding ref-masked sites
            n_valid_sites = int((~ref_mask).sum())
            n_core_sites = int((~core_mask[~ref_mask]).sum())

            # attribute this step's stats to the original index included at step t
            orig_at_t = core_order[t - 1]
            if orig_at_t not in cfs_prog:
                cfs_prog[orig_at_t] = [0, 0]
            cfs_prog[orig_at_t][0] += n_core_sites
            cfs_prog[orig_at_t][1] += n_valid_sites

            if t == final_step:
                final_core_mask = core_mask

        if final_core_mask is None:
            # should not happen, but be defensive
            final_core_mask = core_mask

        core_masks.append(final_core_mask)

    if progressive:
        create_plot([c / t for k, (c, t) in cfs_prog.items()], [names_orig[i] for i in core_order])

    return core_masks, cfs_prog


def find_const(stacks, idx_filt, ref_masks, core_masks, ploidy):
    const_counts = {'A': 0, 'T':0, 'C':0, 'G':0}
    const_masks  = []
    diffs        = []
    for i, stack in enumerate(stacks):

        ref_mask  = ref_masks[i]
        core_mask = core_masks[i]

        stack_valid = stack[idx_filt, :]

        const_mask = np.sum(stack_valid[0] == stack_valid, axis=0) == sum(idx_filt)
        const_masks.append(const_mask)
            
        for k, v in Counter(stack_valid[0,:][const_mask]).items():
            if k in const_counts:
                const_counts[k] += v

    # Save base composition of constants
    filename = 'fconst.txt'
    with open(filename, 'w') as f:
        f.write(f"{const_counts.get('A',0)*ploidy},"
                f"{const_counts.get('C',0)*ploidy},"
                f"{const_counts.get('G',0)*ploidy},"
                f"{const_counts.get('T',0)*ploidy}\n")
    logging.info(f'Saved file -> {filename}')

    return const_masks

def find_diffs(orig_names, stacks, idx_maps, idx_filt, filt_map, ref_masks, core_masks, const_masks, ploidy, bit_lut):
    diffs = []
    names = None
    for i, stack in enumerate(stacks):

        idx_map     = idx_maps[i]
        ref_mask   = ref_masks[i]
        core_mask  = core_masks[i]
        const_mask = const_masks[i]

        stack_valid = stack[idx_filt, :]

        stack_vars = stack_valid[:,(~const_mask & ~core_mask & ~ref_mask)]
        bits = to_bits(stack_vars, bit_lut)
        chunk_size = auto_chunk_size(bits.shape[0])
        diff = calculate_distances(bits, ploidy, chunk_size)
        diff_expanded, names_expanded = expand_distances(diff, filt_map, idx_map, orig_names)
        if names is not None and names_expanded != names:
            raise ValueError("Names are wrong")
        else:
            names = names_expanded

        diffs.append(diff_expanded)

    for i, diff in enumerate(diffs):
        if i == 0:
            diff_sum = diff.copy()
        else:
            diff_sum += diff

    write_distances(names, diff_sum)
