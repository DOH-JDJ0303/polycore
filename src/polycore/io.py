from __future__ import annotations

import logging
import os
import re
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import screed

import plotly.graph_objects as go
import plotly.io as pio


logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Naming / sanitization
# -----------------------------------------------------------------------------

def get_fasta_name(filepath: str) -> str:
    """Return base filename without .gz and without extension."""
    basename = os.path.basename(filepath)
    if basename.endswith(".gz"):
        basename = os.path.splitext(basename)[0]
    return os.path.splitext(basename)[0]


def sanitize_contig_name(name: str) -> str:
    """Make names safe and filesystem-friendly."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", name.split()[0])


def sanitize_name(name: str) -> str:
    """Make names safe and filesystem-friendly."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


# -----------------------------------------------------------------------------
# FASTA / input
# -----------------------------------------------------------------------------

def load_sequences(state) -> None:
    """
    Load FASTA sequences into state.

    Populates:
      - state.sequences: list[list[str]] where outer is sample and inner is contig sequence string
      - state.names:     list[str] sample names (ref first)
      - state.contigs:   list[str] contig names in reference order

    Behavior:
      - cfg.split=True  -> each record becomes its own sample (sample name = record name, contig="contig")
      - cfg.split=False -> each file is one sample (sample name = file base name, contig = record name)
    """
    files = state.files or []
    cfg = state.cfg
    if not files:
        raise ValueError("load_sequences: state.files is empty")

    logger.info("Loading %d FASTA file(s)...", len(files))

    # seq_map[sample_name][rec_index] = {"contig": contig_name, "sequence": SEQ}
    seq_map: Dict[str, Dict[int, Dict[str, str]]] = defaultdict(dict)
    first_name: Optional[str] = None

    for file_i, filepath in enumerate(files, start=1):
        logger.debug("Reading FASTA %d/%d: %s", file_i, len(files), filepath)
        with screed.open(filepath) as handle:
            rec_count = 0
            for rec_i, rec in enumerate(handle):
                rec_count += 1
                if cfg.split:
                    sample_name = sanitize_contig_name(rec["name"])
                    contig_name = "Reference"
                    rec_index = 0
                else:
                    sample_name = get_fasta_name(filepath)
                    contig_name = sanitize_contig_name(rec["name"])
                    rec_index = rec_i

                if file_i == 1 and rec_i == 0:
                    first_name = sample_name

                if rec_index == 0 and sample_name in seq_map:
                    raise ValueError(f"{ 'Contig' if cfg.split else 'Files' } names must be unique: { sample_name }")

                seq_map[sample_name][rec_index] = {
                    "contig": contig_name,
                    "sequence": rec["sequence"].upper(),
                }

        if rec_count == 0:
            raise ValueError(f"FASTA file had no records: {filepath}")

    # Determine reference sample name
    if cfg.ref_by_name:
        ref_name = next((n for n in seq_map.keys() if n.lower() == "reference"), None)
    else:
        ref_name = first_name

    if not ref_name:
        raise ValueError("Reference could not be determined (no records read?)")

    logger.info("Using %s as reference.", ref_name)

    ref_records = seq_map[ref_name]
    if not ref_records:
        raise ValueError("Reference has no records")

    # Reference defines contig order/names and expected lengths
    ref_keys = sorted(ref_records.keys())
    contigs = [ref_records[i]["contig"] for i in ref_keys]
    ref_seqs = [ref_records[i]["sequence"] for i in ref_keys]

    names: List[str] = [ref_name]
    sequences: List[List[str]] = [ref_seqs]
    ref_n_contigs = len(ref_seqs)

    # Validate other samples against reference
    for sample_name, records in seq_map.items():
        if sample_name == ref_name:
            continue

        if len(records) != ref_n_contigs:
            raise ValueError(
                f"{sample_name} contains {len(records)} contigs, but reference has {ref_n_contigs}"
            )

        seq_list: List[str] = []
        for i in ref_keys:
            seq = records[i]["sequence"]
            ref_len = len(ref_records[i]["sequence"])
            if len(seq) != ref_len:
                raise ValueError(
                    f"Contig {i} in sample {sample_name} length {len(seq)} != reference length {ref_len}"
                )
            seq_list.append(seq)

        names.append(sample_name)
        sequences.append(seq_list)

    state.sequences = sequences
    state.names = names
    state.contigs = contigs

    logger.info("Loaded %d sample(s) with %d contig(s)", len(names), len(contigs))


# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------

def write_fconst(const_counts: Dict[str, int], ploidy: int, filename: str = "fconst.txt") -> None:
    """Write constant-base composition scaled by ploidy."""
    if ploidy is None:
        raise ValueError("write_fconst: ploidy is None")

    a = const_counts.get("A", 0) * ploidy
    c = const_counts.get("C", 0) * ploidy
    g = const_counts.get("G", 0) * ploidy
    t = const_counts.get("T", 0) * ploidy

    with open(filename, "w") as f:
        f.write(f"{a},{c},{g},{t}\n")

    logger.info("Saved file -> %s", filename)


def write_distances(names: List[str], diffs: np.ndarray,
                    wide_path: str = "dist_wide.csv",
                    long_path: str = "dist_long.csv") -> None:
    """Write distance matrix in wide and long formats."""
    if diffs.ndim != 2 or diffs.shape[0] != diffs.shape[1]:
        raise ValueError(f"write_distances: diffs must be square 2D; got shape {diffs.shape}")
    if len(names) != diffs.shape[0]:
        raise ValueError(f"write_distances: len(names)={len(names)} != diffs size={diffs.shape[0]}")

    # Wide format
    cols = np.array(["name"] + names, dtype=object)
    rows = np.array(names, dtype=object)
    out = np.column_stack((rows, diffs.astype(str)))
    out = np.vstack((cols, out))
    np.savetxt(wide_path, out, delimiter=",", fmt="%s")
    logger.info("Saved file -> %s", wide_path)

    # Long format
    n_pairs = 0
    with open(long_path, "w") as f:
        f.write("sample1,sample2,diff\n")
        n = len(names)
        for i in range(n):
            for j in range(i + 1, n):
                f.write(f"{names[i]},{names[j]},{int(diffs[i, j])}\n")
                n_pairs += 1

    logger.info("Saved file -> %s", long_path)
    logger.info("Distance matrices: %d pairwise comparison(s)", n_pairs)


def write_summary(state, path: str = "summary.csv") -> None:
    """
    Write per-sample summary (in state.core_order) with:
      - name
      - length: total usable alignment length across contigs (excluding ref-masked sites)
      - missing: length - n_valid, plus ref-masked sites added back to reference missing count
      - genome_fraction
      - core_fraction
      - variants: optional (proxy = sum of row distances) if state.dist exists and matches names
      - MEAN_SNP_DIST, MIN_SNP_DIST, MAX_SNP_DIST: pairwise SNP distance stats,
        excluding self-comparisons and the Reference sample

    Expects:
      - state.names, state.core_order, state.ref_masks, state.n_valid, state.gfs, state.cfs
      - optionally state.dist (square matrix aligned to expanded names order)
    """
    required = ["names", "core_order", "ref_masks", "n_valid", "gfs", "cfs"]
    for k in required:
        if getattr(state, k, None) is None:
            raise ValueError(f"write_summary: state.{k} is required but missing")
        
    n_miss_arr = np.asarray(state.n_miss, dtype=np.int64)
    n_mix_arr = np.asarray(state.n_mix, dtype=np.int64)
    gfs_arr = np.asarray(state.gfs, dtype=float)
    cfs_arr = np.asarray(state.cfs, dtype=float)

    if len(n_miss_arr) != len(state.names):
        raise ValueError("write_summary: state.n_miss length does not match state.names")
    if len(gfs_arr) != len(state.names):
        raise ValueError("write_summary: state.gfs length does not match state.names")

    # cfs is stored in core_order order (from find_core cleanup); enforce that
    if len(cfs_arr) != len(state.core_order):
        # If you instead store cfs per original sample, adjust logic here.
        raise ValueError(
            "write_summary: state.cfs should be aligned to state.core_order "
            f"(len(cfs)={len(cfs_arr)} len(core_order)={len(state.core_order)})"
        )

    variants = state.dist[0]

    # Build SNP distance stats per sample, excluding Reference and self-comparisons
    # state.dist is aligned to state.names_filt, so use names_filt here
    ref_label = "Reference"
    snp_stats = {}
    dist_matrix = state.dist
    names = state.names_filt

    for i, name in enumerate(names):
        if name == ref_label:
            continue
        values = [
            int(dist_matrix[i][j])
            for j, other in enumerate(names)
            if other != name and other != ref_label
        ]
        if values:
            snp_stats[name] = {
                "MEAN_SNP_DIST": round(sum(values) / len(values)),
                "MIN_SNP_DIST": min(values),
                "MAX_SNP_DIST": max(values),
            }

    # Emit rows in core_order
    lines = ["name,length,masked,missing,mixed,genome_fraction,core_fraction,included,variants,MEAN_SNP_DIST,MIN_SNP_DIST,MAX_SNP_DIST"]
    for out_i, orig_idx in enumerate(state.core_order):
        nm = state.names[orig_idx]
        gf = float(gfs_arr[orig_idx])
        cf = float(cfs_arr[out_i])  # cfs aligned to core_order
        miss = int(n_miss_arr[orig_idx])
        mix = int(n_mix_arr[orig_idx])
        incl = nm in state.names_filt
        var = str(variants[state.names_filt.index(nm)]) if incl else 'null'

        s = snp_stats.get(nm, {})
        mean_snp = s.get("MEAN_SNP_DIST", "")
        min_snp = s.get("MIN_SNP_DIST", "")
        max_snp = s.get("MAX_SNP_DIST", "")

        lines.append(
            ",".join([
                nm,
                str(state.n_total),
                str(state.n_mask),
                str(miss),
                str(mix),
                f"{gf:.6f}",
                f"{cf:.6f}",
                str(incl),
                var,
                str(mean_snp),
                str(min_snp),
                str(max_snp),
            ])
        )

    with open(path, "w", newline="") as f:
        f.write("\n".join(lines) + "\n")

    logger.info("Saved file -> %s", path)
    logger.info(
        f"Summary stats: n_total={state.n_total} n_valid={state.n_valid} n_mask={state.n_mask} n_miss_ref={state.n_miss_ref} samples={len(state.names)}"
    )


def create_plot(core_fractions: List[float], names: List[str], path: str = "core_fraction_plot.html") -> bool:
    """Create a progressive core fraction plot (HTML)."""
    try:
        if not core_fractions:
            logger.info("No progression data for plotting")
            return False

        logger.info("Creating progressive core plot -> %s", path)

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=list(range(1, len(core_fractions) + 1)),
            y=core_fractions,
            mode="lines+markers",
            text=names[:len(core_fractions)],
            hovertemplate="<b>%{text}</b><br>Order: %{x}<br>Core Fraction: %{y:.3f}<extra></extra>",
        ))
        fig.update_layout(
            title="Soft-core Genome Fraction vs Sample Addition",
            xaxis_title="Number of sequences included",
            yaxis_title="Soft-core Genome Fraction",
            width=900,
            height=520,
            showlegend=False,
        )
        pio.write_html(fig, path)
        logger.info("Core progression: %.3f -> %.3f", core_fractions[0], core_fractions[-1])
        return True

    except Exception:
        # Plotly is installed, so ImportError isn't the only failure mode; be robust.
        logger.exception("Progressive core plot failed; skipping.")
        return False


def stacks_to_fasta(stack_dict: Dict[str, List[np.ndarray]], filename: str) -> None:
    """Write concatenated stacks per sample to FASTA."""
    n_written = 0
    with open(filename, "w") as f:
        for name, stacks in stack_dict.items():
            seq = "".join(["".join(s.tolist()) for s in stacks])
            f.write(f">{name}\n{seq}\n")
            n_written += 1
    logger.info("Saved file -> %s (%d record(s))", filename, n_written)


def stacks_to_csv(state, path: str = "full.csv") -> None:
    """
    Write a per-site matrix with CHROM/POS/FILTER and each sample’s full-stack base.

    Filter exclusivity rules (highest priority first):
      1) bed_masked        -> ONLY filter reported
      2) ambiguous_reference -> ONLY filter reported
      3) below_min_cf
      4) constant

    Expects:
      - state.full_stacks: dict[name] -> list(per-contig arrays)
      - state.contigs, state.ref_masks, state.core_masks, state.const_masks
      - optionally state.bed_mask: list[bool array] aligned to contigs/stacks
    """
    if state.full_stacks is None:
        raise ValueError("stacks_to_csv: state.full_stacks is missing")
    if state.contigs is None or state.ref_masks is None or state.core_masks is None or state.const_masks is None:
        raise ValueError("stacks_to_csv: missing contigs and/or masks")

    has_bed = getattr(state, "bed_mask", None) is not None
    if has_bed and len(state.bed_mask) != len(state.contigs):
        raise ValueError(
            f"stacks_to_csv: state.bed_mask length ({len(state.bed_mask)}) "
            f"!= n_contigs ({len(state.contigs)})"
        )

    sample_names = list(state.full_stacks.keys())
    out = np.array(["CHROM", "POS", "FILTER"] + sample_names, dtype=object)

    n_bed_total = 0

    for contig_i, contig in enumerate(state.contigs):
        ref_mask = np.asarray(state.ref_masks[contig_i], dtype=bool)     # ambiguous reference
        core_mask = np.asarray(state.core_masks[contig_i], dtype=bool)   # below_min_cf mask (True => below)
        const_mask = np.asarray(state.const_masks[contig_i], dtype=bool) # constant

        bed_mask = None
        if has_bed:
            bed_mask = np.asarray(state.bed_mask[contig_i], dtype=bool)
            if bed_mask.shape[0] != ref_mask.shape[0]:
                raise ValueError(
                    f"stacks_to_csv: bed_mask length mismatch for contig {contig!r}: "
                    f"bed={bed_mask.shape[0]} ref_mask={ref_mask.shape[0]}"
                )
            n_bed = int(bed_mask.sum())
            n_bed_total += n_bed
            logger.info("CSV filters: contig=%s bed_masked_sites=%d", contig, n_bed)

        # Build the static per-site columns once, based on the first sample's stack
        first = True
        block = None

        for name in sample_names:
            stack = state.full_stacks[name][contig_i]

            if first:
                n_sites = int(stack.shape[0])
                chrom = np.array([contig] * n_sites, dtype=object)
                pos = np.arange(1, n_sites + 1, dtype=np.int64)

                # Start with empty filters
                filt = np.empty(n_sites, dtype=object)

                # Precompute boolean conditions
                is_bed = bed_mask if bed_mask is not None else np.zeros(n_sites, dtype=bool)
                is_amb = ref_mask
                is_below = core_mask
                is_const = const_mask

                # Apply exclusivity:
                # bed_masked wins everything
                # ambiguous_reference wins everything except bed
                only_bed = is_bed
                only_amb = (~only_bed) & is_amb
                remaining = (~only_bed) & (~only_amb)

                # For remaining sites, allow combinations of below_min_cf and constant
                below_ok = remaining & is_below
                const_ok = remaining & is_const

                # Fill FILTER strings
                # Defaults to "" (no filter)
                filt[:] = ""

                filt[only_bed] = "bed_masked"
                filt[only_amb] = "ambiguous_reference"

                # For the rest, join any non-empty flags
                # (vectorize by building arrays of strings and joining per index)
                below_f = np.where(below_ok, "below_min_cf", "")
                const_f = np.where(const_ok, "constant", "")

                for k in np.flatnonzero(remaining):
                    flags = []
                    if below_f[k]:
                        flags.append(below_f[k])
                    if const_f[k]:
                        flags.append(const_f[k])
                    filt[k] = ";".join(flags)

                block = np.vstack([chrom, pos.astype(object), filt])
                first = False

            block = np.vstack([block, stack.astype(object)])

        out = np.vstack([out, block.T])

    np.savetxt(path, out, delimiter=",", fmt="%s")
    logger.info("Saved file -> %s", path)
    if has_bed:
        logger.info("CSV filters: total bed_masked_sites=%d", n_bed_total)


# -----------------------------------------------------------------------------
# BED
# -----------------------------------------------------------------------------

def read_bed(path: str) -> List[Tuple[str, int, int]]:
    """
    Read the first 3 columns of a BED file.
    Returns a list of (chrom, start, end).
    """
    out: List[Tuple[str, int, int]] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end = line.split()[:3]
            out.append((chrom, int(start), int(end)))

    logger.info("Read BED -> %s (%d interval(s))", path, len(out))
    return out
