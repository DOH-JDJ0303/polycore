from typing import List, Tuple, Dict, Optional
import numpy as np, logging, screed, os
import plotly.graph_objects as go
import plotly.io as pio

def load_sequences(files: List[str]) -> Tuple[List[str], List[str]]:
    sequences, names = [], []
    logging.info(f"Loading {len(files)} FASTA files...")
    for i, filepath in enumerate(files, 1):
        name = 'Reference' if i == 1 else os.path.splitext(os.path.basename(filepath))[0]
        seq_parts = []
        with screed.open(filepath) as handle:
            for record in handle:
                seq_parts.append(record['sequence'].upper())
        seq = ''.join(seq_parts)
        sequences.append(seq)
        names.append(name)
        seq_len = len(seq)
        if i == 1:
            ref_len = seq_len
        if seq_len != ref_len:
            raise ValueError(f"Sample length ({seq_len:,}) differs from the reference {ref_len:,}: {filepath}")
        logging.info(f"  {i}/{len(files)}: {name} ({seq_len:,} bp)")
    logging.info(f"Loaded {len(sequences)} sequences")
    return sequences, names

def write_distances(names: List[str], diffs: np.ndarray) -> None:
    """Write distance matrices in wide and long formats."""
    # Wide format
    filename = 'dist_wide.csv'
    with open('dist_wide.csv', 'w') as f:
        f.write("name," + ",".join(names) + "\n")
        for name, row in zip(names, diffs):
            row_str = ",".join(str(int(x)) if np.isfinite(x) else "" for x in row)
            f.write(f"{name},{row_str}\n")
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

def write_fasta_from_array(array: np.ndarray, names: List[str], filename: str) -> None:
    """
    Write a 2D array of bases (rows = samples, cols = bases) with associated names to a FASTA file.

    Args:
        array: 2D numpy array of shape (n_samples, n_bases), dtype 'U1' or str.
        names: List of sample names (length n_samples).
        filename: Path to output FASTA file.
    """
    if array.shape[0] != len(names):
        raise ValueError("Number of names must match number of rows in array")

    with open(filename, "w") as f:
        for name, row in zip(names, array):
            seq = "".join(row.tolist())
            f.write(f">{name}\n{seq}\n")
    logging.info(f'Saved file -> {filename}')

def write_vcf_from_array(array: np.ndarray, names: List[str], filename: str) -> None:
    """
    Write a VCF file in the same style as `snp-sites -v`.
    
    Parameters
    ----------
    array : np.ndarray
        2D array of shape (n_samples, n_sites).
        Values should be strings: REF base in column 0, ALT base in other entries,
        or '0'/'1' already encoded as genotypes.
    names : List[str]
        Sample names, order matches rows in array.
    filename : str
        Path to output VCF file.
    """
    if array.shape[0] != len(names):
        raise ValueError("Number of names must match number of rows in array")

    n_sites = array.shape[1]

    # Build header
    header = [
        "##fileformat=VCFv4.1",
        f"##contig=<ID=1,length={n_sites}>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    ]
    cols = [
        "#CHROM", "POS", "ID", "REF", "ALT",
        "QUAL", "FILTER", "INFO", "FORMAT"
    ] + names

    with open(filename, "w") as f:
        # Write header lines
        for line in header:
            f.write(line + "\n")
        f.write("\t".join(cols) + "\n")

        # Iterate over sites
        for pos in range(n_sites):
            # Column slice: genotypes at this site
            site = array[:, pos]

            # Define REF as allele in the first sample
            ref = site[0]
            alts = sorted(set(site) - {ref, "0"})

            alt = ",".join(alts) if alts else "."

            # Map alleles to 0 (REF) / 1 (ALT)
            allele_map = {ref: "0"}
            if alts:
                for i, a in enumerate(alts, start=1):
                    allele_map[a] = str(i)

            genotypes = [allele_map.get(x, "0") for x in site]

            row = [
                "1", str(pos + 1), ".", ref, alt,
                ".", ".", ".", "GT"
            ] + genotypes
            f.write("\t".join(row) + "\n")

    logging.info(f"Saved file -> {filename}")

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