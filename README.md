# AMRScan Pipeline

A comprehensive pipeline for identifying antimicrobial resistance (AMR) genes and variants from metagenomic sequencing data.

AMRScan is a robust and user-friendly pipeline designed to process raw sequencing reads and provide a detailed, normalized report of AMR content. It leverages the CARD database and uses a dual-analysis approach:

1.  **HomScan:** Detects the presence and abundance of AMR genes based on homology.
2.  **VarScan:** Detects known resistance-conferring mutations (e.g., SNPs) in target genes.

The pipeline normalizes results against universal single-copy genes (USCGs) and produces summary tables and rich, interactive HTML visualizations for data exploration.

## Features

-   **Dual Detection:** Simultaneously screens for both AMR gene presence (homology) and known resistance variants.
-   **Packaged Databases:** Comes with the CARD and SCG databases pre-packaged for out-of-the-box use, with options to provide custom databases.
-   **Robust Normalization:** Normalizes AMR gene abundance using Reads Per Kilobase per Genome (RPKG) and per single-copy gene (RPKPC/RPKPMC) for accurate sample comparisons.
-   **Intelligent Read Disambiguation:** Implements a "Winner Takes All" (WTA) strategy to resolve reads that map ambiguously to multiple homologous genes.
-   **Interactive Visualizations:** Automatically generates detailed HTML reports for each AMR gene family, showing read coverage, identity, and uniqueness.
-   **Flexible & Controllable:** Offers fine-grained control over analysis parameters, including separate PID cutoffs for homology and variant detection, choice of protein or nucleotide identity, and gene type filtering.
-   **Safe & Smart Reruns:** Protects against accidental overwrites and allows for efficient re-analysis by skipping the time-consuming mapping step.
-   **Easy Installation:** Packaged for simple installation into a dedicated Conda environment.

## Installation

Installation is handled via Conda to ensure all dependencies are managed correctly.

**Prerequisites:**
*   You must have `conda` (or `miniconda`/`mamba`) installed.
*   You must have `git` installed to clone the repository.

Follow these three steps:

**1. Clone the Repository**
First, clone this repository to your local machine.

```bash
git clone https://github.com/hsgweon/amrscan.git
cd amrscan
```

**2. Create the Conda Environment**
Use the provided `environment.yml` file to create a self-contained environment with all necessary software (BWA, Diamond, Samtools, etc.). This command also uses `pip` to install the `amrscan` scripts.

```bash
conda env create -f environment.yml
```
*This may take a few minutes as it downloads and installs all dependencies.*

**3. Activate the Environment**
Activate the newly created environment. You must do this every time you want to use the pipeline.

```bash
conda activate amrscan-env
```

You are now ready to run the pipeline! The `amrscan` command will be available in your PATH.

## Quick Start

To test your installation and run a basic analysis on paired-end reads using the default packaged databases and 16 threads, use the following command:

```bash
amrscan -i test.fastq.gz -o test -t 2
```
This will create a directory named `my_first_run` containing all the results.

## Usage

The main script `amrscan` coordinates the entire workflow.

```bash
amrscan -i <INPUT_FILES> -o <OUTPUT_PREFIX> [OPTIONS]
```

### All Options

#### **Required Arguments**
| Flag | Description |
| :--- | :--- |
| `-i`, `--input-fastqs` | **Required.** Comma-separated list of input FASTQ files. No spaces between paths. (e.g., `read1.fq.gz,read2.fq.gz` or `single_reads.fq.gz`). |
| `-o`, `--output-prefix` | **Required.** A prefix for all output files and the name of the main output directory that will be created. |

#### **Database Arguments**
*By default, the pipeline uses the databases packaged with the installation. Use these flags only if you want to provide your own.*
| Flag | Description |
| :--- | :--- |
| `--db-card` | Path to a custom CARD AMR database FASTA file. |
| `--db-scg` | Path to a custom Single Copy Genes FASTA file. |
| `--db-scg-lengths` | Path to a custom SCG gene lengths TSV file. |
| `--db-card-metadata` | Path to a custom CARD metadata file. |

#### **Analysis Parameters**
| Flag | Description | Default |
| :--- | :--- | :--- |
| `--homscan-pid-cutoff` | Minimum percent identity for HomScan hits (0.0-1.0 scale). | `0.95` |
| `--varscan-pid-cutoff` | Minimum nucleotide percent identity for VarScan hits (0.0-1.0 scale). | `0.95` |
| `--homscan-pid-type` | PID type to use for HomScan filtering and WTA (`protein` or `nucleotide`). | `protein` |
| `--consensus-cutoff` | Minimum fraction of ambiguous hits that must map to the same gene family to reach a consensus. | `0.8` |
| `--homscan-gene-types` | Comma-separated list of gene types for HomScan (e.g., 'H,K'). | `H` |
| `--varscan-gene-types` | Comma-separated list of variant types for VarScan (e.g., 'V,R,O'). | `V,R` |

#### **Performance**
| Flag | Description | Default |
| :--- | :--- | :--- |
| `-t`, `--threads` | Number of threads to use for computationally intensive steps. | `1` |

#### **Workflow Control**
| Flag | Description |
| :--- | :--- |
| `--overwrite` | Overwrite the output directory if it already exists. **Use with caution.** |
| `--skip-mapping` | Skip the initial (and longest) mapping steps. Requires that the pipeline has been run successfully on the same data before. |
| `--debug` | Enable debug mode, which creates detailed log files for each sub-script. |
| `--version` | Show the pipeline version and exit. |
| `-h`, `--help` | Show the help message and exit. |

## Examples

### Example 1: Standard Paired-End Analysis
A standard run on paired-end data, using 32 threads and the default protein PID cutoff of 90% for HomScan.

```bash
amrscan -i sampleA_R1.fastq.gz,sampleA_R2.fastq.gz \
        -o SampleA_results \
        -t 32
```

### Example 2: Advanced Analysis with Custom Parameters
A more advanced run with several custom settings:
- Use a stricter **nucleotide** PID of 95% for HomScan.
- Use a very strict nucleotide PID of 99% for VarScan.
- Analyze both `H` and `K` gene types in HomScan.
- Use a custom CARD metadata file.

```bash
amrscan -i sampleB_single.fastq.gz \
        -o SampleB_custom \
        -t 32 \
        --homscan-pid-type nucleotide \
        --homscan-pid-cutoff 0.95 \
        --varscan-pid-cutoff 0.99 \
        --homscan-gene-types H,K \
        --db-card-metadata /path/to/my/custom_metadata.txt
```

### Example 3: Re-running Analysis with a Different Cutoff
Imagine the run from Example 1 finished, but you want to re-analyze the HomScan results with a stricter protein PID of 98% without re-doing the slow mapping step.

The combination of `--skip-mapping` and `--overwrite` is designed for this. It will safely delete the old analysis results while preserving the mapping data.

```bash
amrscan -i sampleA_R1.fastq.gz,sampleA_R2.fastq.gz \
        -o SampleA_results \
        -t 32 \
        --homscan-pid-cutoff 0.98 \
        --skip-mapping \
        --overwrite
```

## Output Files

All results will be located in the directory specified by `-o`. The key outputs are:

-   **`SampleA_results/`** (Main output directory)
    -   `SampleA_results_homscan.tsv`: The final, normalized summary table for homology-based gene detection.
    -   `SampleA_results_homscan_detailed.tsv`: A more detailed version of the HomScan report with ambiguous hits resolved.
    -   `SampleA_results_varscan.tsv`: The final, normalized summary table for variant detection.
    -   **`SampleA_results_homscan_html/`**: A directory containing interactive HTML coverage plots for each detected AMR gene family.
    -   **`SampleA_results_varscan_html/`**: A directory containing HTML alignment views for each detected AMR variant.
    -   **`logs/`**: A directory containing the main run log and a detailed run summary file.

The `tmp/` directory contains intermediate files and can typically be ignored.

## License
This project is licensed under the MIT License.

## Contact
For questions, bug reports, or suggestions, please open an issue on this GitHub repository.