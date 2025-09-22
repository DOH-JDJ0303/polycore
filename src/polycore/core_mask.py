import numpy as np, logging
from typing import List, Tuple, Dict, Optional
from collections import Counter
from .io_ops import create_plot

def filter_sequences(stack, names, valid_bases, min_gf=0.9):
    names = np.array(names)  # ensure numpy array
    # Filter invalid bases in reference
    ref = stack[0]
    stack_valid = stack[:, np.isin(ref, valid_bases)]
    logging.info(f"Removed {stack.shape[1] - stack_valid.shape[1]} invalid reference positions")
    # Mask invalid bases in all samples
    stack_valid[~np.isin(stack_valid, valid_bases)] = 'N'

    logging.info(f"Filtering {stack.shape[0]} sequences, min_gf={min_gf}")
    # Calculate genome fraction
    gf = np.sum(stack_valid != 'N', axis=1) / stack_valid.shape[1]
    keep = gf >= min_gf

    logging.info(f"Kept {np.sum(keep)}/{len(keep)} sequences")
    if np.any(~keep):
        logging.info(f"Sequences with genome fraction below {min_gf}: {names[~keep]}")

    return stack_valid[keep, :], names[keep], gf[keep]

def find_core(stack, names, gf, threshold=1.0, progressive=True):
    """Calculate progressive core genome fraction."""
    n_rows, n_cols = stack.shape
    if not progressive:
        logging.info('Determining core (non-progressive)')
        fractions = np.sum(stack != 'N', axis=0) / n_rows
        core_mask = fractions >= threshold
        core_fraction = np.sum(core_mask) / n_cols
        logging.info(f"Sites below min-cf ({threshold}): {np.sum(~core_mask)}")
        logging.info(f"Final core fraction: {core_fraction:.2f}")
        stack_core = stack[:, core_mask]
        return stack_core, names, [core_fraction] * len(names)

    logging.info('Determing soft-core (progressive):')
    sort_indices = np.argsort(gf[1:])[::-1] + 1
    sorted_stack = np.vstack([stack[0:1], stack[sort_indices]])
    sorted_names = np.concatenate([names[0:1], names[sort_indices]])

    cfs = []                          # progression trajectory (sorted order)
    per_sample_cf = {}                # sample -> core fraction at its addition
    for i in range(1, n_rows + 1):
        fractions = np.sum(sorted_stack[0:i] != 'N', axis=0) / i
        core_mask = fractions >= threshold
        core_fraction = np.sum(core_mask) / n_cols
        cfs.append(core_fraction)
        sample_name = sorted_names[i-1]
        per_sample_cf[sample_name] = core_fraction
        logging.info(f"  {i}/{n_rows}: {sample_name} ({core_fraction:.2f})")

    logging.info(f"Sites below min-cf ({threshold}): {np.sum(~core_mask)}")
    logging.info(f"Final core fraction: {core_fraction:.2f}")

    create_plot(cfs, sorted_names)

    stack_core_sorted = sorted_stack[:, core_mask]

    # re-sort rows back to the original input order
    combined_idx = np.concatenate(([0], sort_indices))
    pos_to_original = np.argsort(combined_idx)
    stack_core_original = stack_core_sorted[pos_to_original, :]

    # now produce cfs in the same order as names
    cfs_original = [per_sample_cf[n] for n in names]

    return stack_core_original, names, cfs_original

def find_const(stack, names, ploidy, min_pf, min_pn):
    logging.info("Finding constant / variable sites")
    ref = stack[0]
    samples = stack[1:]
    n_samples = samples.shape[0]

    # Count matches against ref, ignoring N's
    is_match = (samples == ref) & (samples != 'N')
    n_matches = np.sum(is_match, axis=0)

    # Count informative samples (exclude N's)
    informative = np.sum(samples != 'N', axis=0)

    # Find SNVs (at least one informative mismatch)
    var_mask = (n_matches < informative)
    logging.info(f"Found {np.sum(var_mask)} variants")

    # Apply frequency/number filters if requested
    if min_pf > 0 or min_pn > 0:
        frac_nonref = (informative - n_matches) / np.maximum(informative, 1)  # avoid div/0
        num_nonref = informative - n_matches
        var_mask = var_mask & (frac_nonref >= min_pf) & (num_nonref >= min_pn)
        logging.info(f"Filtered to {np.sum(var_mask)} variants (min-pf: {min_pf}, min-pn: {min_pn})")

    # Constant = everything else
    const_mask = ~var_mask
    logging.info(f"Remaining {np.sum(const_mask)} sites treated as constant")

    # Count constant bases (using ref row)
    const = Counter(ref[const_mask])
    vars = stack[:, var_mask]

    # Ensure vars has consistent shape even if empty
    if vars.shape[1] == 0:
        vars = np.empty((stack.shape[0], 0), dtype=stack.dtype)

    # Save base composition of constants
    filename = 'fconst.txt'
    with open(filename, 'w') as f:
        f.write(f"{const.get('A',0)*ploidy},"
                f"{const.get('C',0)*ploidy},"
                f"{const.get('G',0)*ploidy},"
                f"{const.get('T',0)*ploidy}\n")
    logging.info(f'Saved file -> {filename}')

    return vars