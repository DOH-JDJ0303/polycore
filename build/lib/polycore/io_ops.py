from typing import List, Tuple, Dict, Optional
import numpy as np, logging, screed, os
import plotly.graph_objects as go
import plotly.io as pio
import re
from collections import defaultdict

def get_fasta_name(filepath):
    basename = os.path.basename(filepath)

    # Remove .gz if present
    if basename.endswith('.gz'):
        basename = os.path.splitext(basename)[0]

    # Remove fasta extension and return
    return os.path.splitext(basename)[0]


def sanitize_name(name: str) -> str:
    """Make contig names safe and filesystem-friendly."""
    return re.sub(r'[^A-Za-z0-9._-]', '_', name)


def load_sequences(files: List[str], split: bool = False, ref_by_name: bool = False) -> Dict[str, Dict[str, str]]:
    logging.info(f"Loading {len(files)} FASTA files...")

    seq_map = defaultdict(lambda: defaultdict(dict))
    first_name = None
    
    for i, filepath in enumerate(files, 1):
        with screed.open(filepath) as handle:
            for e, rec in enumerate(handle):

                if split:
                    name = sanitize_name(rec['name'])
                    contig = 'contig'
                else:
                    name   = get_fasta_name(filepath)
                    contig = sanitize_name(rec['name'])

                if i == 1 and e == 0:
                    first_name = name

                seq_map[name][e] = { 'contig': contig, 'sequence': rec['sequence'].upper() }

    # Check for duplicate names after loading all files
    if len(seq_map) != len(files) and not split:
        raise ValueError(f"Multiple files with the same name detected")

    # Determine reference
    if ref_by_name:
        ref_name = next((n for n in seq_map.keys() if n.lower() == 'reference'), None)
    else:
        ref_name = first_name

    if not ref_name:
        raise ValueError("Reference could not be determined!")

    logging.info(f"Using {ref_name} as reference.")

    ref_seq = seq_map[ref_name]
    # always have reference as first
    names     = [ ref_name ]
    contigs   = [ r['contig'] for r in ref_seq.values()     ]
    sequences = [ [r['sequence'] for r in ref_seq.values()] ]

    for name, v1 in seq_map.items():
        if name == ref_name:
            continue
        if len(v1) != len(sequences[0]):
            raise ValueError(f"{name} contains a different number of contigs than the reference!")
        seq_list = []
        for e, v2 in v1.items():
            if len(v2['sequence']) != len(ref_seq[e]['sequence']):
                raise ValueError(f"Contig {e} in sample {name} is a different length than the reference!")
            seq_list.append(v2['sequence'])

        names.append(name)
        sequences.append(seq_list)

    return sequences, names, contigs

def write_distances(names: List[str], diffs: np.ndarray) -> None:
    """Write distance matrices in wide and long formats."""
    # Wide format
    cols = np.array(['name'] + names)
    rows = np.array(names)
    out = np.column_stack((rows, diffs.astype(str)))
    out = np.vstack((cols, out))
    filename = 'dist_wide.csv'
    np.savetxt(filename, out, delimiter=",", fmt="%s")
    logging.info(f"Saved filed -> {filename}")
    
    # Long format
    n_pairs = 0
    filename = 'dist_long.csv'
    with open('dist_long.csv', 'w') as f:
        f.write("sample1,sample2,diff\n")
        n = len(names)
        for i in range(n):
            for j in range(i+1, n):
                d = int(diffs[i,j])
                f.write(f"{names[i]},{names[j]},{d}\n")
                n_pairs += 1
    logging.info(f"Saved filed -> {filename}")
    logging.info(f"Distance matrices: {n_pairs} pairwise comparisons")

def write_summary(names, stack, gf, cf, variants):
    """
    Write per-sample summary:
      - name
      - sequence length (columns)
      - number of missing sites (Ns)
      - genome fraction
      - core fraction
      - variant count (if provided)
    """
    length = stack.shape[1]
    missing = np.sum(stack == 'N', axis=1)  # per-sample
    lines = ["name,length,missing,genome_fraction,core_fraction,variants"]

    for i, name in enumerate(names):
        row = [
            name,
            str(length),
            str(missing[i]),
            f"{gf[i]:.6f}",
            f"{cf[i]:.6f}",
            str(variants[i])
        ]
        lines.append(",".join(row))

    with open("summary.csv", "w") as f:
        f.write("\n".join(lines) + "\n")
    logging.info(f'Saved file -> summary.csv')

def create_plot(core_fractions: List[float], names: List[str]) -> bool:
    """Create core fraction progression plot."""
    try:
        logging.info("Creating progressive core plot")
        if not core_fractions:
            logging.info("No progression data for plotting")
            return False
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=list(range(1, len(core_fractions)+1)),
            y=core_fractions,
            mode='lines+markers',
            text=names[:len(core_fractions)],
            hovertemplate='<b>%{text}</b><br>Order: %{x}<br>Core Fraction: %{y:.3f}<extra></extra>'
        ))
        fig.update_layout(
            title='Soft-core Genome Fraction vs Sample Addition',
            xaxis_title='Number of sequences included',
            yaxis_title='Soft-core Genome Fraction',
            width=900, height=520, showlegend=False
        )
        pio.write_html(fig, 'core_fraction_plot.html')
        logging.info(f"core_fraction_plot.html: progression from {core_fractions[0]:.3f} to {core_fractions[-1]:.3f}")
        return True
    except ImportError:
        logging.info("Progressive core plot failed, skipping...")
        return False

def stacks_to_fasta(stack_dict, filename):
    with open(filename, 'w') as f:
         for name, stacks in stack_dict.items():
            seq = ''.join([''.join(s.tolist()) for s in stacks])
            f.write(f">{name}\n{seq}\n")

def stacks_to_csv(full_stacks, contig_names, ref_masks, core_masks, const_masks):
    out = np.array(
        ["CHROM", "POS", "FILTER"] + list(full_stacks.keys())
    )
    for i, contig in enumerate(contig_names):
        ref_mask  = ref_masks[i]
        core_mask = core_masks[i]
        const_mask = const_masks[i]

        first = True
        for name, stacks in full_stacks.items():
            stack = stacks[i]
            if first:
                n_sites = stack.shape[0]
                chrom   = np.array([contig]*n_sites)
                pos     = np.array(range(1,n_sites+1))

                ref_mask_filt   = np.where(ref_mask, "ambiguous_reference", "")
                core_mask_filt  = np.where(core_mask, "below_min_cf", "")
                const_mask_filt = np.where(const_mask, "constant", "")
                filt = np.where(~ref_mask | ~core_mask | ~const_mask, ref_mask_filt + core_mask_filt + const_mask_filt, ref_mask_filt + ";" + core_mask_filt + ";" + const_mask_filt)

                block = np.vstack([chrom, pos, filt])
                first = False

            block = np.vstack([block, stack])
        out = np.vstack([out, np.transpose(block)])

    np.savetxt("full.csv", out, delimiter=",", fmt="%s")