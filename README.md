# PolyCore

PolyCore is a Python-based tool for **core genome analysis in polyploid organisms**.  
It loads reference and sample FASTA files, collapses identical sequences, filters by genome/core fractions, and produces:

- Core alignment FASTA files
- Variant alignment FASTA/VCF files
- Pairwise distance matrices (wide and long formats)
- Per-sample summary table
- Progressive core fraction plot (HTML)

---

## Features

- Handles haploid and polyploid genomes (auto-detects ploidy if not specified)
- Soft-core / progressive core fraction calculation
- Collapsing and re-expansion of identical sequences
- Distance matrices with efficient chunking (auto memory-aware)
- Output in CSV, FASTA, and VCF formats
- Interactive visualization with Plotly

---

## Installation

Clone the repository and install:

```bash
git clone https://github.com/WA-DOH/polycore.git
cd polycore
pip install -e .
```
> PolyCore requires Python 3.10+.
Dependencies (numpy, screed, psutil, plotly) are installed automatically.

---
## Usage
Run PolyCore from the command line:
```
polycore --ref reference.fasta --sample sample1.fasta sample2.fasta --progressive
```
### Common options
- `--min-gf` : Minimum genome fraction per sample (default: 0.9)
- `--min-cf` : Minimum fraction of population required per site (default: 0.95)
- `--min-pf` : Minimum fraction of alt allele per site (for SNPs vs SNVs)
- `--min-pn` : Minimum number of samples with alt allele per site
- `--ploidy` : Force ploidy (otherwise auto-detected)
- `--progressive` : Enable soft-core (progressive) calculation

For full options:
```
polycore --help
```

### Outputs
Outputs

- `core.aln` : Core alignment (variants only, FASTA)
- `core.full.aln` : Full core alignment (FASTA)
- `core.vcf` : Variants in VCF format
- `dist_wide.csv` : Pairwise distance matrix (wide)
- `dist_long.csv` : Pairwise distance matrix (long/tidy)
- `summary.csv` : Per-sample statistics
- `core_fraction_plot.html` : Interactive visualization of soft-core genome fraction
