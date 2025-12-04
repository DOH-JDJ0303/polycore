import numpy as np, logging
from typing import List, Tuple, Dict, Optional
from collections import Counter
from .io_ops import create_plot

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

def filter_sequences(stack, names, valid_bases, min_gf=0.9, bed=None):
    """
    Filter and mask sequences.

    - stack: 2D numpy array of shape (n_seqs, n_sites), dtype='U1' (chars)
    - names: list of sequence names, length n_seqs
    - valid_bases: iterable of valid characters (e.g. ['A','C','G','T','-'])
    - min_gf: minimum genome fraction to keep a sequence
    - bed: optional path to BED file with regions to mask (0-based, half-open)

    Returns: (stack_valid, gf, keep_mask)
    
    Raises:
        ValueError: If all sites are masked or no sequences pass the filter
    """
    names = np.array(names)  # ensure numpy array

    # Filter invalid bases in reference (drop columns where ref is not in valid_bases)
    ref = stack[0]
    ref_valid_mask = np.isin(ref, valid_bases)
    stack_valid = stack[:, ref_valid_mask]
    logging.info(f"Removed {stack.shape[1] - stack_valid.shape[1]} invalid reference positions")
    
    # Check if all sites were removed
    if stack_valid.shape[1] == 0:
        raise ValueError("All sites were filtered out due to invalid reference bases")

    # Mask invalid bases in all samples
    stack_valid[~np.isin(stack_valid, valid_bases)] = 'N'

    logging.info(f"Filtering {stack_valid.shape[0]} sequences, min_gf={min_gf}")

    # Calculate genome fraction (fraction of non-N sites)
    gf = (np.sum(stack_valid != 'N', axis=1) / stack_valid.shape[1]).astype(float)
    keep = gf >= min_gf

    logging.info(f"Kept {np.sum(keep)}/{len(keep)} sequences")
    
    # Check if no sequences pass the filter
    if np.sum(keep) == 0:
        raise ValueError(f"No sequences passed the genome fraction filter (min_gf={min_gf}). "
                        f"Maximum genome fraction was {gf.max():.3f}")
    
    if np.any(~keep):
        logging.info(f"Sequences with genome fraction below {min_gf}: {names[~keep]}")

    # Mask sites specified in BED (set to 'N' in all sequences)
    mask = None
    if bed is not None:
        mask = read_bed(bed)   # [(chrom, start, end), ...]
        logging.info(f"Loaded {len(mask)} mask intervals from {bed}")
        if mask:
            n_sites = stack_valid.shape[1]
            mask_cols = np.zeros(n_sites, dtype=bool)

            for _, start, end in mask:
                # Clip to alignment bounds, and be robust to non-int inputs
                start = max(0, int(start))
                end   = min(n_sites, int(end))
                if start < end:
                    mask_cols[start:end] = True

            n_masked = int(mask_cols.sum())
            if n_masked > 0:
                stack_valid[:, mask_cols] = 'N'
                logging.info(f"Masked {n_masked} sites from BED regions")
                
                # Check if all sites are now masked
                if n_masked == n_sites:
                    raise ValueError("All sites were masked by BED regions")

    return stack_valid, gf, keep


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
        return stack_core, names, [np.nan] * len(names)

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