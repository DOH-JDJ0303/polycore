import sys, time, argparse, logging, numpy as np
from .utils import set_up_logging, IUPAC_BITS, auto_chunk_size
from .io_ops import load_sequences, write_distances, write_fasta_from_array, write_vcf_from_array, write_summary
from .collapse import collapse_sequences, expand_results, expand_distances, expand_vector
from .distance import set_ploidy, create_stack, to_bits, calculate_distances
from .core_mask import filter_sequences, find_core, find_const
from typing import List, Tuple, Dict, Optional


__version__ = "1.1"

def build_parser():
    p = argparse.ArgumentParser(description="PolyCore - Core genome analysis on polyploid organisms")
    p.add_argument("--ref", required=True, help="Reference FASTA file")
    p.add_argument("--sample", nargs="+", required=True, help="Sample sequences")
    p.add_argument("--min-gf", type=float, default=0.9, help="Minimum genome fraction per sample")
    p.add_argument("--min-cf", type=float, default=0.95, help="Minimum fraction with valid data per site")
    p.add_argument("--min-pf", type=float, default=0, help="Min fraction with alt per site (SNP vs SNV)")
    p.add_argument("--min-pn", type=float, default=0, help="Min # samples with alt per site (SNP vs SNV)")
    p.add_argument("--progressive", action='store_true')
    p.add_argument("--ploidy", type=int)
    p.add_argument("--chunk-size", type=int, help="Sites per chunk for pairwise diffs (controls memory)")
    p.add_argument("--version", action="version", version=__version__)
    return p

def main(argv=None):
    argv = argv or sys.argv[1:]
    set_up_logging()
    logger = logging.getLogger(__name__)
    parser = build_parser()
    args = parser.parse_args(argv)

    start = time.time()
    logger.info(f"PolyCore v{__version__} starting")

    try:
        files = [args.ref] + args.sample
        sequences, orig_names = load_sequences(files)

        # Collapse → bitmaps → stack
        sequences, names_rep, idx_map = collapse_sequences(sequences, orig_names)
        ploidy, bit_map, bit_lut = set_ploidy(sequences, IUPAC_BITS, args.ploidy)
        stack, valid_bases = create_stack(sequences, bit_map)

        # Filtering
        stack_valid, gf, filter_mask = filter_sequences(
            stack, np.array(names_rep), valid_bases, args.min_gf
        )
        stack_filt  = stack_valid[filter_mask, :]
        names_filt  = np.array(names_rep)[filter_mask]
        gf_filt     = gf[filter_mask]

        # Core only on kept reps
        core, names_core, cfs = find_core(
            stack_filt, names_filt, gf_filt,
            threshold=args.min_cf, progressive=args.progressive
        )

        # Vars on core set
        vars = find_const(core, names_core, ploidy, args.min_pf, args.min_pn)

        # Distances on core variants
        bits = to_bits(vars, bit_lut)
        chunk_size = args.chunk_size or auto_chunk_size(bits.shape[0])
        diffs = calculate_distances(bits, ploidy, chunk_size)

        # ---------------- Expansion strategy ----------------
        no_mask = np.full(stack_valid.shape[0], True, dtype=bool)
        stack_exp, names_exp = expand_results(stack_valid, no_mask, idx_map, orig_names)
        gf_exp, _ = expand_results(np.array(gf),  no_mask, idx_map,  orig_names)
        cfs_exp, _ = expand_results(np.array(cfs), filter_mask, idx_map, orig_names)
        core_exp, names_core_exp = expand_results(core, filter_mask, idx_map, orig_names, keep_filtered=False)
        vars_exp, _ = expand_results(vars, filter_mask, idx_map, orig_names, keep_filtered=False)
        diffs0_exp, _ = expand_results(diffs[0,:], filter_mask, idx_map, orig_names)
        diffs_exp, diffs_exp_names  = expand_distances(diffs, filter_mask, idx_map, orig_names)

        # ---------------- Write outputs ----------------
        # Distances/FASTA/VCF: only core set (no filtered samples)
        write_distances(diffs_exp_names, diffs_exp)
        write_fasta_from_array(core_exp, names_core_exp, "core.full.aln")
        write_fasta_from_array(vars_exp, names_core_exp, "core.aln")
        write_vcf_from_array(vars_exp, names_core_exp, "core.vcf")

        # Summary: all originals
        write_summary(names_exp, stack_exp, gf_exp, cfs_exp, diffs0_exp)

        logger.info(f"Done in {time.time()-start:.1f}s")

    except Exception:
        logger = logging.getLogger(__name__)
        logger.exception("Failed after %.1fs", time.time() - start)
        sys.exit(1)

