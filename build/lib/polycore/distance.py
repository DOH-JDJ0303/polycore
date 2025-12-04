from typing import List, Tuple, Dict, Optional
from .utils import IUPAC_BITS, ALLELES, POPCOUNT16, ambiguity_size
import numpy as np, logging

def to_bits(sequences: np.ndarray, lut) -> np.ndarray:
    logging.info('Converting variants to bits')
    up = np.char.upper(sequences)
    codes_obj = np.frompyfunc(ord, 1, 1)(up)
    codes_int = codes_obj.astype(np.int32)
    if np.any(codes_int > 127):
        bad_chars = np.unique(up[codes_int > 127])
        bad_desc = ", ".join(f"{repr(c)} (U+{ord(c):04X})" for c in bad_chars)
        raise ValueError(
            f"Non-ASCII character(s) detected in sequences: {bad_desc}. "
            "Please sanitize the FASTA to ASCII IUPAC letters (A,C,G,T and ambiguity codes) and '-' only."
        )
    codes = codes_int.astype(np.uint8)
    return lut[codes]

def build_match_table(ploidy: int) -> np.ndarray:
    """Table[m1, m2] = number of matching copies between two masks at this ploidy."""
    table = np.zeros((16, 16), dtype=np.uint8)
    # Precompute counts for all masks once
    cnt = [ _mask_to_counts(m, ploidy) for m in range(16) ]
    for m1 in range(16):
        c1 = cnt[m1]
        for m2 in range(16):
            c2 = cnt[m2]
            table[m1, m2] = np.minimum(c1, c2).sum()
    return table

def _mask_to_counts(mask: int, ploidy: int) -> np.ndarray:
    """Convert a 4-bit IUPAC mask (0..15) to a length-4 allele count vector
       under the fixed-composition interpretation:
       - k = popcount(mask)
       - if k == 0: [0,0,0,0] (unknown)
       - if k == 1: put all ploidy copies on that allele
       - if k == ploidy: 1 copy on each allele in the mask
       (Assumes k <= ploidy; otherwise treat as invalid upstream.)
    """
    k = POPCOUNT16[mask]
    counts = np.zeros(4, dtype=np.uint8)
    if k == 0:
        return counts
    idxs = [i for i, bit in enumerate(ALLELES) if mask & bit]
    if k == 1:
        counts[idxs[0]] = ploidy
    elif k == ploidy:
        for i in idxs:
            counts[i] = 1
    else:
        # If you *might* see k in (1, ploidy) but not equal to either, you need a rule.
        # For strict pipelines, treat as unknown (all zeros) or raise.
        # Here we'll mark as unknown to avoid bias:
        counts[:] = 0
    return counts

def calculate_distances(bits: np.ndarray, ploidy: int, chunk_size: int) -> np.ndarray:
    """
    bits: (n_samples, n_sites) uint8 bitmasks (0=unknown, A=1, C=2, G=4, T=8, combos via OR)
    """
    logging.info(f"Calculating pairwise distances in {chunk_size} bp chunks")
    n, L = bits.shape
    diffs = np.zeros((n, n), dtype=np.int64)
    match_table = build_match_table(ploidy)

    for start in range(0, L, chunk_size):
        end = min(start + chunk_size, L)
        chunk = bits[:, start:end]  # (n, w)
        logging.info(f"  Processing sites {start}-{end}")

        for i in range(n):
            bi = chunk[i]
            for j in range(i+1, n):
                bj = chunk[j]

                # only compare where both known
                both_known = (bi > 0) & (bj > 0)

                if not np.any(both_known):
                    continue

                # matches per site from lookup; vectorized gather
                matches = match_table[bi[both_known], bj[both_known]].astype(np.int64)

                # mismatches per site = ploidy - matches
                d = (ploidy - matches).sum()

                diffs[i, j] += d
                diffs[j, i] += d

    return diffs

def set_ploidy(sequences, bit_map, ploidy=None):
    # Auto-detect ploidy if not specified
    if ploidy is None:
        all_chars = set(''.join(sequences))
        ploidy = max((ambiguity_size(c) for c in all_chars if ambiguity_size(c) > 0), default=1)
        logging.info(f"Auto-detected ploidy: {ploidy}")
    else:
        logging.info(f"Using specified ploidy: {ploidy}")
        
    # Update accepted IUPAC characters based on ploidy
    # Characters not in IUPAC_BITS_FILT will be treated as zero
    bit_map_filt = { b: v for b, v in bit_map.items() if ambiguity_size(b) <= ploidy }

    # Precompute lookup table for fast conversion
    lut = np.zeros(128, dtype=np.uint8)
    for ch, bits in bit_map_filt.items():
        lut[ord(ch)] = lut[ord(ch.lower())] = bits
    return ploidy, bit_map_filt, lut

def create_stack(sequences, bit_map):
    logging.info(f"Creating array stack")
    valid_bases = list(bit_map.keys())
    logging.info(f"Valid bases: {valid_bases}")
    stack = np.vstack([np.array(list(s), dtype='U1') for s in sequences])
    return stack, valid_bases