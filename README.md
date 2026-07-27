# PolyCore

PolyCore is a Python-based tool for **core genome analysis in polyploid organisms**.

## Features
Polycore can do the following for haploid, diploid, and triploid organisms:
* Determine the core genome (user defined - "soft" or "hard" core)
* Auto-filter low quality samples by the fraction of sites with missing data (`N`)
* Calculate pairwise nucleotide differences, taking into account ploidy
* Calculate constant sites for use with IQTREE
* Outputs core genome alignments (`FASTA` and `CSV`) and distance matrices (`CSV`)

## Overview

![](docs/flowchart.svg)

---

## Installation

#### Option 1 - Use the container (Preferred Method):

```bash
docker pull public.ecr.aws/o8h2f0o1/polycore:1.1.0
```

#### Option 2 - Clone the repository and install:

```bash
git clone https://github.com/NW-PaGe/polycore.git
cd polycore
pip install -e .
```

---

## Usage

Run PolyCore from the command line:

```bash
polycore \
    --progressive \
    --ref reference.fasta \
    sample1.fasta sample2.fasta
```

### Options

```
PolyCore - Core genome analysis on polyploid organisms

positional arguments:
  input                 Input sequences

options:
  -h, --help            show this help message and exit
  --ref REF             Reference FASTA file
  --mask MASK           Bed file with coordinates of sites to exclude.
  --min-gf MIN_GF       Minimum genome fraction per input
  --min-cf MIN_CF       Minimum fraction with valid data per site
  --progressive
  --ploidy PLOIDY
  --chunk-size CHUNK_SIZE
                        Sites per chunk for pairwise diffs (controls memory)
  --split               Treat each contig in a multi-fasta file as a separate sample.
  --ref-by-name         Treat file/contig labeled 'reference' as the reference (case-insensitive).
  --snippy              Shortcut for Snippy *.full.aln input (sets split and ref-by-name).
  --version             show program's version number and exit
```

---

## Outputs

* `core.aln` : Core alignment (variants only, FASTA)
* `core.full.aln` : Full core alignment (FASTA)
* `full.csv` : Per-site summary for all passing samples
* `dist_wide.csv` : Pairwise distance matrix (wide)
* `dist_long.csv` : Pairwise distance matrix (long/tidy)
* `summary.csv` : Per-sample statistics
* `core_fraction_plot.html` : Interactive visualization of soft-core genome fraction

---

## Inputs
Polycore accepts a multiple sequence alignment as input. This can be supplied as **separate FASTA files** or as **single multi-FASTA**.

#### Separate FASTA files
Each FASTA file is treated as a separate sample. Files with multiple contigs are concatenated into a single sequence.
```
# Separate FASTA files
polycore *.fasta
```

#### A single multi-FASTA file
Each contig is treated as a separate sample.
```
# Single multi-FASTA
polycore --split alignment.fasta
```

### Reference Sequence
Polycore uses a *reference** to calculate each sample's genome fraction. Only sites with valid information (non-`N`) in the reference are considered. The sequence used as the reference is defined in the following manner:

- **Sequence order** - The first sequence in the alignment is treated as the reference. Sequences supplied via the `--ref` parameter are treated as the first sequence.
- **Sequence name** - The contig or FASTA file labelled `reference` is treated as the reference when using the `--ref-by-name` parameter. This supersedes reference selection by sequence order.

> [!NOTE]
> The `--snippy` parameter is the same as running `--split --ref-by-name`.

## Ploidy Determination
Ploidy can be determined automatically or specified by the user using the `--ploidy` parameter. When not specified, Polycore sets the ploidy based on the IUPAC codes observed within the multiple sequence alignment using the rules shown below:

| IUPAC Codes | Ploidy |
| ----------- | ------ |
| `ATCG`      | 1      |
| `RYSWKM`    | 2      |
| `BVDH`      | 3      |

## Masking Invalid Sites
Invalid sites are coded to IUPAC code `N`. Sites that are not valid in the reference are masked to `N` across every sample in the alignment. At all remaining sites, IUPAC codes that exceed the declared ploidy, along with any non-IUPAC codes, are converted to `N` on a per-sample basis.

### Example
```
Input alignment (--ploidy 2):

reference  ACGT A ACGT A ACGT N ACGT
sample_1   ACGT B ACGT R ACGT A ACGT
sample_2   ACGT Y ACGT X ACGT R ACGT

After conversion:

reference  ACGT A ACGT A ACGT N ACGT
sample_1   ACGT N ACGT R ACGT N ACGT
sample_2   ACGT Y ACGT N ACGT N ACGT
```

## Genome Fraction Calculation
The genome fraction of each sample is calculated as the fraction of valid reference sites that are also valid in the sampel:

```
genome fraction = valid sample sites / valid reference sites
```

Samples with genome fractions below the threshold declared by `--min-gf` are excluded from the analysis.

## Core Genome Calculation
The core genome is defined as genome positions with valid IUPAC codes in the minimum fraction of samples (reference inclusive) defined by the core fraction threshold, `--min-cf`. Genome positions that do not meet this minimum threshold are excluded from the core genome analysis. By default, the core genome is calculated in a single step, using all samples that pass the minimum genome fraction threshold. It is also possible to calculate the core genome *progressively* using the `--progressive` parameters. This calculates the core genome following the addition of each sample sequentially, in order of genome fraction size. The final core genome is **indentical** regardless of which method is used. The *progressive* method is provided as a means of measuring how your core genome size changes as your population grows.

## Distance Calculation
Polycore estimates pairwise distances by splitting IUPAC codes into equal **allele fractions** and determing the **total unmatched allele fraction** across comparable sites. This unmatched fraction is then scaled by the ploidy to give the approximate number of unmatched alleles across all chromsome copies. It is critical to understand that IUPAC codes only tell you _which_ alleles were observed, <u> not</u> which chromosome copy the alleles are observed on. This means that **we must assume the location of heterozygous alleles**. Polycore does this by assuming each allele is split across each chromsome copy equally, thus resulting in the possibility of fractional distance values.

The example below demonstrates how IUPAC codes are broken down into allele fractions. Note that IUPAC codes that represent a number of alleles greater than the ploidy are treated as invalid and ignored.
| IUPAC Code  | Haploid (P=1) | Diploid (P=2) | Triploid (P=3) |
| ----------- | ------------: | ------------: | -------------: |
|A            | 1 A         | 1 A         | 1 A          |
|R            | -             | 1/2 A, 1/2 G  | 1/2 A, 1/2 G   |
|B            | -             | -             | 1/3 C, 1/3 G, 1/3 T |

### Example single-site distances

| Sample 1  | Sample 2    | Haploid (P=1) | Diploid (P=2) | Triploid (P=3) |
| --------- | ----------- | ------------: | ------------: | -------------: |
| A         | A           |             0 |             0 |              0 |
| A         | G           |             1 |             2 |              3 |
| A         | R (A/G)     |             - |             1 |            1.5 |
| R (A/G)   | R (A/G)     |             - |             0 |              0 |
| R (A/G)   | W (A/T)     |             - |             1 |            1.5 |
| R (A/G)   | B (C/G/T)   |             - |             - |              2 |
| B (C/G/T) | B (C/G/T)   |             - |             - |              0 |
| B (C/G/T) | D (A/G/T)   |             - |             - |              1 |
| A         | N (unknown) |             - |             - |              - |

`-` indicates that the comparison is not applicable

### Interpretation

For example, at ploidy 3:

```text
A = 3 A
R = 1.5 A + 1.5 G
```

The shared allele contribution is `1.5 A`, resulting in a distance of:

```text
3 - 1.5 = 1.5
```

Thus, an `A` versus `R (A/G)` comparison contributes `1.5` to the distance at a triploid site.

The total pairwise distance is the sum of these per-site distances across all comparable sites.
