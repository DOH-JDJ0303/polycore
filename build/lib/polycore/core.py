import logging
from collections import Counter
from typing import List, Optional

import numpy as np

from .io import create_plot, read_bed, write_distances, write_fconst
from .utils import (
    auto_chunk_size,
    collapse_sequences,
    to_bits,
    calculate_distances,
    expand_distances,
    ambiguity_size,
    MIXED_IUPAC,
    COUNT_DENOM
)

logger = logging.getLogger(__name__)


def set_ploidy(state, bit_map):
    """Set ploidy and filter bit map based on auto-detection or config."""
    cfg = state.cfg

    # Auto-detect ploidy if not specified
    if cfg.ploidy is None:
        all_chars = set("".join("".join(contigs) for contigs in state.sequences))
        ploidy = max(
            (ambiguity_size(ch) for ch in all_chars if ambiguity_size(ch) > 0),
            default=1,
        )
        logger.info(f"Auto-detected ploidy: {ploidy}")
    else:
        ploidy = cfg.ploidy
        logger.info(f"Using specified ploidy: {ploidy}")

    # Filter bit map and create ASCII lookup table
    bit_map_filt = {ch: bits for ch, bits in bit_map.items() if ambiguity_size(ch) <= ploidy}

    lut = np.zeros(128, dtype=np.uint8)
    for ch, bits in bit_map_filt.items():
        for char in (ch, ch.lower()):
            o = ord(char)
            if o < 128:
                lut[o] = bits

    dropped = sorted(set(bit_map) - set(bit_map_filt))
    if dropped:
        logger.info(f"Codes above ploidy {ploidy} will be treated as missing: {', '.join(dropped)}")

    state.ploidy = ploidy
    state.bit_map = bit_map_filt
    state.bit_lut = lut
    state.valid_bases = list(bit_map_filt.keys())


# -----------------------------------------------------------------------------
# BED masking
# -----------------------------------------------------------------------------


def gather_bed_mask(state, bed_file: Optional[str]) -> None:
    """Build and apply a BED mask per contig."""
    if state.contigs is None or state.stacks is None:
        raise ValueError("gather_bed_mask: state.contigs and state.stacks must be set")

    # Always create aligned masks (even if no BED)
    state.bed_mask = [np.zeros(st.shape[1], dtype=bool) for st in state.stacks]

    if not bed_file:
        logger.info("No BED mask provided; skipping BED masking.")
        return

    bed = read_bed(bed_file)
    logger.info(f"Loaded {len(bed)} mask interval(s) from {bed_file}")

    if not bed:
        logger.info("BED contained no intervals; nothing to mask.")
        return

    contig_to_idx = {c: i for i, c in enumerate(state.contigs)}
    unknown_chroms = set()

    # Fill masks
    for chrom, start, end in bed:
        idx = contig_to_idx.get(chrom)
        if idx is None:
            unknown_chroms.add(chrom)
            continue

        mask_cols = state.bed_mask[idx]
        n_sites = mask_cols.shape[0]

        s = max(0, int(start) - 1)
        e = min(n_sites, int(end) - 1)
        if s < e:
            mask_cols[s:e] = True

    if unknown_chroms:
        logger.warning(f"BED chrom(s) not found in contigs (skipped): {', '.join(sorted(unknown_chroms))}")

    # Apply masks and log per contig
    total = 0
    for contig, st, m in zip(state.contigs, state.stacks, state.bed_mask):
        n_sites = st.shape[1]
        n_masked = int(m.sum())
        total += n_masked

        if n_masked == 0:
            logger.info(f"BED mask: contig={contig} masked_sites=0")
            continue

        if n_masked == n_sites:
            raise ValueError(f"All sites were masked by BED regions for contig {contig!r}")

        pct = 100.0 * n_masked / n_sites
        logger.info(f"BED mask: contig={contig} masked_sites={n_masked} ({pct:.2f}%)")

    state.n_mask = total

    logger.info(f"BED masking complete. Total masked sites across contigs: {total}")


# -----------------------------------------------------------------------------
# Sample filtering / metrics
# -----------------------------------------------------------------------------


def calculate_genome_fraction(state) -> None:
    """Compute genome fraction per sample relative to the reference (index 0)."""
    if state.n_valid is None:
        raise ValueError("calculate_genome_fraction: state.n_valid is missing")

    if state.n_valid == 0:
        raise ValueError("calculate_genome_fraction: no valid sites")

    if not state.n_miss:
        raise ValueError("calculate_genome_fraction: state.n_miss is missing")

    # Reference is always index 0
    gfs = []
    gf_ref = [1.0]

    for v in state.n_miss[1:]:
        gfs.append((state.n_valid - float(v)) / state.n_valid)

    state.gfs = gf_ref + gfs

    logger.info(
        f"Genome fraction: min={np.min(gfs):.3f} "
        f"median={np.median(gfs):.3f} "
        f"max={np.max(gfs):.3f}"
    )



def filter_samples(state) -> None:
    """Filter samples by genome fraction (min_gf)."""
    cfg = state.cfg
    if state.gfs is None:
        raise ValueError("filter_samples: state.gfs not set (run calculate_genome_fraction first)")
    if state.names is None:
        raise ValueError("filter_samples: state.names not set")
    if state.idx_maps is None or state.stacks is None:
        raise ValueError("filter_samples: state.idx_maps/state.stacks not set")

    keep_orig = [gf > cfg.min_gf for gf in state.gfs]
    state.names_filt = [nm for nm, keep in zip(state.names, keep_orig) if keep]
    keep_set = set(state.names_filt)

    logger.info(f"Filtering samples by min_gf={cfg.min_gf:.3f}: kept={len(state.names_filt)}/{len(state.names)}")

    # Build per-contig representative keep masks
    filt_maps: List[List[bool]] = []
    for contig_i, idx_map in enumerate(state.idx_maps):
        keep_reps = [any(state.names[o] in keep_set for o in orig_idxs) for _, orig_idxs in idx_map.items()]
        filt_maps.append(keep_reps)

        n_kept = sum(keep_reps)
        logger.debug(f"Contig {contig_i}: kept {n_kept}/{len(keep_reps)} representative row(s)")

    state.filt_maps = filt_maps


# -----------------------------------------------------------------------------
# Core genome
# -----------------------------------------------------------------------------


def find_core(state) -> None:
    """Determine per-contig core masks and per-sample core fractions."""
    cfg = state.cfg
    if state.gfs is None or state.ref_masks is None:
        raise ValueError("find_core: requires state.gfs and state.ref_masks")
    if state.names is None or state.names_filt is None:
        raise ValueError("find_core: requires state.names and state.names_filt")
    if state.idx_maps is None or state.stacks is None:
        raise ValueError("find_core: requires state.idx_maps and state.stacks")

    # Sort originals by gf (ref first)
    state.core_order = [0] + sorted(range(1, len(state.gfs)), key=lambda j: state.gfs[j], reverse=True)
    names_filt_set = set(state.names_filt)

    logger.info(f"Determining core genome (min_cf={cfg.min_cf:.3f} progressive={cfg.progressive})")

    cfs_prog: dict[int, List[int]] = {}  # orig_idx -> [core_sites_sum, valid_sites_sum]
    core_masks: List[np.ndarray] = []

    for contig_i, stack in enumerate(state.stacks):
        idx_map = state.idx_maps[contig_i]
        ref_mask = state.ref_masks[contig_i]

        # Build orig -> rep lookup
        orig_to_rep = {}
        for rep, origs in idx_map.items():
            for o in origs:
                orig_to_rep[o] = rep

        rep_order = []
        final_step = None

        for t, orig_idx in enumerate(state.core_order):
            rep = orig_to_rep.get(orig_idx)
            rep_order.append(rep)
            if state.names[orig_idx] in names_filt_set:
                final_step = t + 1

        if not rep_order:
            raise ValueError(f"Contig {contig_i}: rep_order is empty (idx_map/core_order mismatch)")

        # Choose progressive steps
        t_values = range(1, len(rep_order) + 1) if cfg.progressive else range(final_step, final_step + 1)

        rep_set = None
        valid_matrix = None
        final_core_mask = None

        for t in t_values:
            reps = rep_order[:t]
            rep_counts = Counter(reps)

            rep_set_new = list(dict.fromkeys(reps))  # order-preserving unique
            if rep_set is None or rep_set_new != rep_set:
                rep_set = rep_set_new
                valid_matrix = (stack[rep_set, :] != "N").astype(np.int32)

            weights = np.array([rep_counts[r] for r in rep_set], dtype=np.int32)
            fractions = (valid_matrix * weights[:, None]).sum(axis=0) / t

            core_mask = fractions < cfg.min_cf  # True means NOT core
            n_valid_sites = int((~ref_mask).sum())
            n_core_sites = int((~core_mask[~ref_mask]).sum())

            orig_at_t = state.core_order[t - 1]
            cfs_prog.setdefault(orig_at_t, [0, 0])
            cfs_prog[orig_at_t][0] += n_core_sites
            cfs_prog[orig_at_t][1] += n_valid_sites

            if t == final_step:
                final_valid_sites = n_valid_sites
                final_core_sites  = n_core_sites
                final_core_mask = core_mask

        core_masks.append(final_core_mask)

        logger.info(
            f"Core summary: contig={state.contigs[contig_i]} "
            f"valid_sites={final_valid_sites} "
            f"core_sites={final_core_sites} "
            f"({(100 * final_core_sites / final_valid_sites):.1f}%)"
        )


    # Core fractions in the same order as core_order
    cfs = []
    labels = []
    for orig_idx in state.core_order:
        c, t = cfs_prog.get(orig_idx, (0, 1))
        cfs.append(c / t)
        labels.append(state.names[orig_idx])

    state.core_masks = core_masks
    state.cfs = cfs

    if cfg.progressive:
        create_plot(state.cfs, labels)


# -----------------------------------------------------------------------------
# Constant sites + diffs/distances
# -----------------------------------------------------------------------------


def find_const(state) -> None:
    """Find constant sites among filtered representatives per contig."""
    if state.filt_maps is None or state.ref_masks is None or state.core_masks is None:
        raise ValueError("find_const: requires filt_maps, ref_masks, core_masks")
    if state.stacks is None or state.names_filt is None:
        raise ValueError("find_const: requires stacks and names_filt")

    const_counts = {"A": 0, "T": 0, "C": 0, "G": 0}
    const_masks: List[np.ndarray] = []

    n_filt = len(state.names_filt)
    logger.info(f"Finding constant sites using {n_filt} filtered sample(s)")

    for contig_i, stack in enumerate(state.stacks):
        filt_map = np.array(state.filt_maps[contig_i], dtype=bool)
        stack_filt = stack[filt_map, :]

        if stack_filt.shape[0] != n_filt:
            logger.debug(
                f"Contig {contig_i}: filtered reps={stack_filt.shape[0]} (names_filt={n_filt}) "
                f"due to representative collapsing"
            )

        # Site is constant if all rows equal row0
        const_mask = (np.sum(stack_filt[0] == stack_filt, axis=0) == stack_filt.shape[0])
        const_masks.append(const_mask)

        # Count constant bases from reference row
        for base, count in Counter(stack_filt[0, :][const_mask]).items():
            if base in const_counts:
                const_counts[base] += int(count)

    write_fconst(const_counts, state.ploidy)
    state.const_masks = const_masks


def find_diffs(state) -> None:
    """Compute per-contig distance matrices from variable sites and sum across contigs."""
    if state.const_masks is None or state.core_masks is None or state.ref_masks is None:
        raise ValueError("find_diffs: requires const_masks, core_masks, ref_masks")
    if state.filt_maps is None or state.idx_maps is None:
        raise ValueError("find_diffs: requires filt_maps and idx_maps")
    if state.stacks is None or state.names is None:
        raise ValueError("find_diffs: requires stacks and names")
    if state.bit_lut is None or state.ploidy is None:
        raise ValueError("find_diffs: requires bit_lut and ploidy")
    if state.core_order is None:
        raise ValueError("find_diffs: requires state.core_order (run find_core first)")

    logger.info(f"Computing differences/distances across {len(state.stacks)} contig(s)")

    # Canonical order for summing
    name_to_idx = {n: i for i, n in enumerate(state.names_filt)}
    n_sites = len(state.names_filt)
    diff_sum = np.zeros((n_sites, n_sites), dtype=np.int64)

    for contig_i, stack in enumerate(state.stacks):
        idx_map = state.idx_maps[contig_i]
        filt_map = state.filt_maps[contig_i]
        ref_mask = state.ref_masks[contig_i]
        core_mask = state.core_masks[contig_i]
        const_mask = state.const_masks[contig_i]

        stack_filt = stack[np.array(filt_map, dtype=bool), :]

        var_cols = (~const_mask) & (~core_mask) & (~ref_mask)
        n_var = int(var_cols.sum())

        contig_name = state.contigs[contig_i] if state.contigs else contig_i
        logger.info(f"Contig {contig_name}: variable_sites={n_var}")

        if n_var == 0:
            bits = np.zeros((stack_filt.shape[0], 0), dtype=np.uint8)
        else:
            stack_vars = stack_filt[:, var_cols]
            bits = to_bits(stack_vars, state.bit_lut)

        chunk_size = state.cfg.chunk_size or auto_chunk_size(bits.shape[0])
        diff = calculate_distances(bits, state.ploidy, chunk_size)

        max_possible = COUNT_DENOM * state.ploidy * n_var
        if diff.max() > max_possible:
            raise ValueError(f"contig {contig_name}: distance {diff.max()} exceeds max {max_possible}")

        diff_expanded, names_expanded = expand_distances(diff, filt_map, idx_map, state.names, state.names_filt)

        # Map this contig's expanded matrix into canonical indices
        try:
            idx = np.array([name_to_idx[n] for n in names_expanded], dtype=np.int64)
        except KeyError as e:
            raise ValueError(f"find_diffs: expanded name not found in canonical names: {e}") from e

        # Add into the correct rows/cols
        diff_sum[np.ix_(idx, idx)] += diff_expanded.astype(np.int64, copy=False)

    state.dist = diff_sum / COUNT_DENOM
    
    write_distances(state.names_filt, state.dist)


# -----------------------------------------------------------------------------
# Stack building / valid sites / expansion
# -----------------------------------------------------------------------------


def build_stacks(state) -> None:
    """Build per-contig stacks (2D arrays) and representative idx_maps."""
    if state.sequences is None or state.names is None:
        raise ValueError("build_stacks: state.sequences and state.names must be set")
    if state.contigs is None:
        raise ValueError("build_stacks: state.contigs must be set")
    if state.valid_bases is None:
        raise ValueError("build_stacks: state.valid_bases must be set")

    stacks: List[np.ndarray] = []
    idx_maps: List[dict] = []

    n_contigs = len(state.sequences[0])
    if n_contigs != len(state.contigs):
        logger.warning(f"Contig count mismatch: sequences={n_contigs} contigs={len(state.contigs)}")

    valid_set = set(state.valid_bases)
    logger.info(f"Building stacks for {n_contigs} contig(s)")

    for contig_i in range(n_contigs):
        contig_name = state.contigs[contig_i] if contig_i < len(state.contigs) else str(contig_i)
        logger.info(f"Creating stack for contig={contig_name}")

        # Per-contig sequences
        contig_seqs = [sample[contig_i] for sample in state.sequences]

        # Representative collapsing map
        rep_seqs, idx_map = collapse_sequences(contig_seqs, state.names)
        idx_maps.append(idx_map)

        # Build stack from original contig_seqs
        stack = np.vstack([np.array(list("".join(rep_seqs[i])), dtype="U1") for i in idx_map.keys()])

        # Normalize invalid characters to 'N'
        mask_invalid = ~np.isin(stack, list(valid_set))
        if mask_invalid.any():
            n_bad = int(mask_invalid.sum())
            logger.debug(f"Contig {contig_name}: normalized {n_bad} invalid base(s) to 'N'")
            stack[mask_invalid] = "N"

        stacks.append(stack)

    state.stacks = stacks
    state.idx_maps = idx_maps


def find_valid_sites(state) -> None:
    """Determine ref_masks and n_valid per original sample across contigs."""
    if state.stacks is None or state.idx_maps is None or state.names is None:
        raise ValueError("find_valid_sites: requires stacks, idx_maps, names")

    # Apply BED mask if configured (in-place)
    if getattr(state.cfg, "bed_file", None):
        gather_bed_mask(state, state.cfg.bed_file)

    ref_masks: List[np.ndarray] = []
    n_miss  = [0 for _ in state.names]
    n_mix   = [0 for _ in state.names]
    n_total = 0
    n_valid = 0

    for contig_i, stack in enumerate(state.stacks):
        n_total += stack[0].size

        ref_mask = (stack[0] == "N") | (state.bed_mask[contig_i] if state.bed_mask is not None else False)
        ref_masks.append(ref_mask)

        n_valid += np.sum(~ref_mask)
        
        # missing and mixed counts per representative row, excluding ref-masked columns
        sub = stack[:, ~ref_mask]

        n_miss_rep  = (sub == "N").sum(axis=1)
        n_mixed_rep = np.isin(sub, MIXED_IUPAC).sum(axis=1)

        # Add rep counts back to original samples via idx_map
        idx_map = state.idx_maps[contig_i]
        for rep_idx, (miss, mix) in enumerate(zip(n_miss_rep, n_mixed_rep)):
            for orig_idx in idx_map[rep_idx]:
                n_miss[orig_idx] += int(miss)
                n_mix[orig_idx]  += int(mix)

        contig_name = state.contigs[contig_i] if state.contigs else contig_i
        logger.info(f"Valid sites: contig={contig_name} valid_cols={int((~ref_mask).sum())} (of {stack.shape[1]})")

    # Update missing sites in reference
    n_miss_ref = n_total - n_valid - state.n_mask
    n_miss[0] = n_miss_ref

    state.n_total    = n_total
    state.n_valid    = n_valid
    state.n_miss_ref = n_miss_ref
    state.n_miss     = n_miss
    state.n_mix      = n_mix
    state.ref_masks  = ref_masks
    
    logger.info(
        f"Missing sites per sample (summed contigs): "
        f"min={int(np.min(n_miss))} median={int(np.median(n_miss))} max={int(np.max(n_miss))}"
    )

def expand_stacks(state) -> None:
    """Expand stacks back to per-original-sample outputs."""
    if state.names is None or state.names_filt is None:
        raise ValueError("expand_stacks: requires names and names_filt")
    if state.idx_maps is None or state.stacks is None:
        raise ValueError("expand_stacks: requires idx_maps and stacks")
    if state.ref_masks is None or state.core_masks is None or state.const_masks is None:
        raise ValueError("expand_stacks: requires ref_masks, core_masks, const_masks")

    full_stacks: dict[str, List[np.ndarray]] = {}
    full_mask_stacks: dict[str, List[np.ndarray]] = {}
    core_mask_stacks: dict[str, List[np.ndarray]] = {}

    kept = 0
    names_filt_set = set(state.names_filt)

    for orig_i, name in enumerate(state.names):
        if name not in names_filt_set:
            continue
        kept += 1

        # Determine representative index per contig for this original sample
        rep_idxs: List[int] = []
        for idx_map in state.idx_maps:
            rep_for_orig = None
            for rep_idx, origs in idx_map.items():
                if orig_i in origs:
                    rep_for_orig = rep_idx
                    break
            if rep_for_orig is None:
                raise ValueError(f"expand_stacks: could not map original index {orig_i} to representative")
            rep_idxs.append(rep_for_orig)

        full_by_contig: List[np.ndarray] = []
        full_mask_by_contig: List[np.ndarray] = []
        core_mask_by_contig: List[np.ndarray] = []

        for contig_i, stack in enumerate(state.stacks):
            rep_idx = rep_idxs[contig_i]

            ref_mask = state.ref_masks[contig_i]
            core_mask = state.core_masks[contig_i]
            const_mask = state.const_masks[contig_i]

            full_row = stack[rep_idx].copy()
            
            full_mask_row = stack[rep_idx].copy()
            full_mask_row[ref_mask | core_mask] = "N"

            core_mask_row = full_mask_row[(~const_mask) & (~ref_mask)]

            full_by_contig.append(full_row)
            full_mask_by_contig.append(full_mask_row)
            core_mask_by_contig.append(core_mask_row)

        full_stacks[name] = full_by_contig
        full_mask_stacks[name] = full_mask_by_contig
        core_mask_stacks[name] = core_mask_by_contig

    state.full_stacks      = full_stacks
    state.full_mask_stacks = full_mask_stacks
    state.core_mask_stacks = core_mask_stacks

    logger.info(f"Expanded stacks for {kept} filtered sample(s)")
