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

        sequences, rep_names, idx_map = collapse_sequences(sequences, orig_names)
        ploidy, bit_map, bit_lut = set_ploidy(sequences, IUPAC_BITS, args.ploidy)
        stack, valid_bases = create_stack(sequences, bit_map)

        stack, names_filt, gf = filter_sequences(stack, np.array(rep_names), valid_bases, args.min_gf)
        core, names_core, cfs = find_core(stack, names_filt, gf, threshold=args.min_cf, progressive=args.progressive)
        vars = find_const(core, names_core, ploidy, args.min_pf, args.min_pn)

        bits = to_bits(vars, bit_lut)
        chunk_size = args.chunk_size or auto_chunk_size(bits.shape[0])  # note: n_samples drives memory, not sites

        diffs = calculate_distances(bits, ploidy, chunk_size)

        # Expand back to original names
        stack_full, names_full = expand_results(stack, names_filt, idx_map, orig_names)
        core_full, names_full  = expand_results(core, names_core, idx_map, orig_names)
        vars_full, names_full  = expand_results(vars, names_core, idx_map, orig_names)
        diffs_full, names_full = expand_distances(diffs, names_core, idx_map, orig_names)
        gf_full  = expand_vector(gf,  names_filt, idx_map)
        cfs_full = expand_vector(cfs, names_core, idx_map)

        write_distances(names_full, diffs_full)
        write_summary(names_full, stack_full, gf_full, cfs_full, diffs_full[0, :].tolist())
        write_fasta_from_array(core_full, names_full, 'core.full.aln')
        write_fasta_from_array(vars_full, names_full, 'core.aln')
        write_vcf_from_array(vars_full, names_full, 'core.vcf')

        logger.info(f"Done in {time.time()-start:.1f}s")

    except Exception as e:
        logger = logging.getLogger(__name__)
        logger.error(f"Failed after {time.time()-start:.1f}s: {e}")
        sys.exit(1)
