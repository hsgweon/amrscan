# AMRScan v.1.0.3

A comprehensive pipeline for identifying antimicrobial resistance (AMR) genes and variants from metagenomic sequencing data.

AMRScan is a robust and user-friendly pipeline designed to process raw sequencing reads and provide a detailed, normalised report of AMR content. It leverages the CARD database and uses a dual-analysis approach:

1.  **HomScan:** Detects the presence and abundance of AMR genes based on homology.
2.  **VarScan:** Detects known resistance-conferring mutations (e.g., SNPs) in target genes.

The pipeline normalises results against universal single-copy genes (USCGs) and produces summary tables and rich, interactive HTML visualisations for data exploration.

## Features

-   **Dual Detection:** Simultaneously screens for both AMR gene presence (homology) and known resistance variants.
-   **Packaged Databases:** Comes with the CARD and SCG databases pre-packaged for out-of-the-box use.
-   **Robust Normalisation:** Normalises AMR gene abundance using a suite of metrics (RPK, RPKG, RPKPC, FPK, etc.) for well-informed, multi-level sample comparisons.
-   **Advanced Ambiguity Resolution:** Implements a **Maximum A Posteriori (MAP)** iterative algorithm to statistically resolve ambiguous reads, identifying the most probable gene within a homologous family and reporting the proportion of evidence supporting each candidate.
-   **Prior-Guided Analysis:** Allows the incorporation of external knowledge (e.g., clinical prevalence data) via a user-provided priors file to improve the accuracy of the MAP resolution.
-   **Interactive Visualisations:** Automatically generates detailed HTML reports for each AMR gene family, showing read coverage, identity, and uniqueness.
-   **Flexible & Controllable:** Offers fine-grained control over analysis parameters, including separate PID cutoffs for homology and variant detection, choice of protein or nucleotide identity, and gene type filtering.
-   **Safe & Smart Reruns:** Protects against accidental overwrites and allows for efficient re-analysis by skipping the time-consuming mapping step.
-   **Easy Installation:** Packaged for simple installation into a Conda environment.

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
amrscan -i test.fastq.gz -o test_run -t 2
```
This will create a directory named `test_run` containing all the results.


## Pipeline Workflow

The amrscan command executes a series of scripts in a coordinated workflow:

1.  ***Mapping & QC (BWA & Diamond)***:
-   Input reads are mapped against the CARD AMR database using bwa mem.
-   Input reads are simultaneously mapped against a protein database of Universal Single-Copy Genes (USCGs) using diamond blastx.
2.  ***Normalisation Scaffolding (scgscan_*)***:
-   The total number of bases in the input sample is calculated.
-   The DIAMOND alignments are used to quantify the abundance of USCGs, calculating an overall RPK (Reads Per Kilobase) value that represents the average 'coverage' of a single-copy gene in the sample. This value is crucial for normalisation.
3.  ***VarScan (varscan_*)***:
-   The BWA alignments to the CARD database are processed to identify reads that both map to a variant-type gene and contain the specific, known resistance-conferring mutation (e.g., a SNP).
-   Confirmed variant hits are tabulated, normalised (RPK, RPKG, RPKPC, RPKPMC, FPK, FPKG, FPKPC, FPKPMC), and alignment visualisations are generated.
4.  ***HomScan (homscan_*)***:
-   The BWA alignments are processed again, this time to identify any read that maps to a homology-type AMR gene above a specified identity cutoff.
-   Tabulation & Normalisation: All passing hits are aggregated. This step produces two initial reports: a _homscan.tsv file where ambiguous reads are grouped into 'multiple' categories.
-   MAP Resolution: A final, more sophisticated Maximum A Posteriori (MAP) solver (homscan_resolve_MAP.py) is run. It uses an iterative algorithm to distribute the abundance from ambiguous reads among candidate genes based on the evidence from uniquely mapped reads and optional user-provided priors. This produces the final, most accurate quantitative report.
5.  ***Consolidation***: 
-   The final summary tables from all steps are copied from the temporary directory into the main output directory for easy access.

### Understanding the Normalisation Metrics

Raw read counts are not directly comparable between genes or samples due to variations in gene length and sequencing depth. To address this, **AMRScan** calculates several normalised abundance metrics.

The pipeline uses two fundamental units for these calculations: **Reads** and **Fragments**.
* A **Read** is a single sequence from a FASTQ file.
* A **Fragment** represents the original piece of DNA sequenced. For paired-end data, the two reads (R1 and R2) from the same DNA fragment are counted as a single fragment. 

The primary metric for interpreting AMR abundance in this pipeline is **FPKPMC**, which estimates the copy number of a gene per million bacterial cells. This is the **signature metric** of the pipeline and the default used by the MAP algorithm for the most accurate abundance estimation (configurable).

Here is a breakdown of all the metrics calculated:

| Metric | Calculation | Interpretation |
| :--- | :--- | :--- |
| **RPK / FPK** | Reads / (Gene Length / 1000) | **Reads/Fragments Per Kilobase.** Normalises for gene length. A longer gene will naturally recruit more reads than a shorter one at the same copy number, and this metric corrects for that bias. |
| **RPKG / FPKG** | RPK / (Total Sample Bases / 1,000,000,000) | **Reads/Fragments Per Kilobase per Gigabase.** Normalises RPK/FPK for sequencing depth. This allows for direct comparison of a gene's abundance between two samples that were sequenced to different depths (e.g., one with 2 Gb of data and one with 10 Gb). |
| **RPKPC / FPKPC** | RPK_amr_gene / RPK_uscg_average | **Reads/Fragments Per Kilobase Per single-copy gene Copy.** Normalises the AMR gene's RPK by the average RPK of all Universal Single-Copy Genes (USCGs). USCGs act as a proxy for the number of bacterial genomes in the sample. This metric estimates the average copy number of the AMR gene per bacterium. A value of ~1.0 suggests the gene is present in a single copy in the average genome, while a value of ~2.0 suggests two copies (e.g., on a plasmid and a chromosome). |
| **RPKPMC / FPKPMC** | RPKPC * 1,000,000 | **RPKPC/FPKPC Multiplied by one Million.** This is the RPKPC/FPKPC value scaled up for easier reading, similar to "parts per million". It is useful for comparing the relative prevalence of genes within or between samples, especially when copy numbers are low. |

## Understanding the MAP Resolver
The most challenging aspect of quantifying AMR genes from metagenomes is handling ambiguous reads—reads that map with high identity to multiple different but highly similar reference genes (e.g., different variants of bla_CTX-M or bla_KPC).

AMRScan's ***Maximum A Posteriori (MAP)*** resolver provides a statistical solution to this problem by determining which gene is the most likely source of the ambiguous reads.

### How it Works
The resolver uses an iterative, expectation-maximisation-like algorithm. The core idea is:

1.  ***Initial Evidence***: The algorithm first calculates the abundance for each AMR gene using only the reads that map **uniquely** to it. This forms the initial, most reliable estimate.
2.  ***Iterative Allocation***: It then iterates through all the ambiguous reads. For each ambiguous read, it looks at the current abundance estimates of all the genes it mapped to. It proportionally distributes the ambiguous read's abundance among those candidate genes. For example, if a read maps to Gene A and Gene B, and Gene A currently has 90 "votes" (from unique reads) while Gene B has 10, the ambiguous read's abundance will be split 90% to Gene A and 10% to Gene B.
3.  ***Convergence***: This process is repeated for many iterations. In each round, the abundance estimates are refined based on the newly distributed evidence from the previous round. The algorithm stops when the abundance estimates for all genes stabilise (i.e., they no longer change significantly between iterations).

The final output of this process is not a new set of abundance values, but rather two key pieces of information added to the summary report:

-   Top_ARO: The single ARO within an ambiguous group that received the highest proportion of evidence after the iterative process. If there is a tie, it is reported as 'N/A'.
-   Allocation_Proportions: A user-friendly string showing the percentage of evidence allocated to each candidate gene in the ambiguous group (e.g., ARO_3002312:100.0% or ARO_3002264:33.3%,ARO_3002266:33.3%,ARO_3005493:33.3%).

### Incorporating External Knowledge with Priors
Sometimes, the evidence from unique reads is sparse or non-existent for a group of similar genes. In these cases, the MAP resolver can be guided by external knowledge using a priors file. This allows you to inject information about which genes are more likely to be present based on historical data, clinical relevance, or other studies.

The `--map-priors-file` option takes a simple tab-separated file where you assign a numeric "weight" to specific AROs. These weights act as "pseudo-counts" that are added to the evidence from the sequencing data.
Example Priors File (clinical_priors.tsv):

#### Example Priors File (clinical_priors.tsv):

```tsv
# Priors for common clinically significant AMR genes
# Format: ARO_ID    Weight
# Weights are based on a 5-year study of local hospital isolates.

# --- High Priority / Very Common Carbapenemases ---
ARO_3002312	15.0	# KPC-2, extremely common in our region
ARO_3000589	12.0	# NDM-1
ARO_3001782	9.5	# OXA-48

# --- Common Extended-Spectrum Beta-Lactamases (ESBLs) ---
ARO_3001878	8.0	# CTX-M-15
ARO_3001877	6.0	# CTX-M-14
ARO_3001071	4.0	# SHV-12

# --- Important MCR genes (Colistin Resistance) ---
ARO_3003689	20.0	# mcr-1, high clinical alert status

# --- Tetracycline resistance genes (very common but lower alert) ---
ARO_3000165	2.5	# tet(A)
ARO_3000166	2.5	# tet(B)
```

### MAP Resolver Parameters
You can fine-tune the resolver's behavior with these options:
-   `--map-metric-column`: Specifies which quantitative column should be resolved and distributed. Default: FPKPMC.
-   `--map-priors-file`: Path to your priors file, formatted as shown above.
-   `--map-prior-strength`: A multiplier that controls the influence of the priors. A higher value gives the priors more weight compared to the evidence from the sequencing data. Default: 1.0.
-   `--map-base-prior`: A tiny "pseudo-count" applied to all genes, even those not in the priors file. This prevents division-by-zero errors and ensures that any gene in an ambiguous group has a non-zero chance of being allocated some abundance. Default: 1.0.


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


#### **MAP Resolver Arguments**
| Flag	| Description	| Default |
| :--- | :--- | :--- |
| `--map-priors-file` |	Path to a tab-separated file of priors for the MAP resolver. | None |
| `--map-metric-column` |	The numeric column to use for MAP abundance resolution.	| RPKG |
| `--map-base-prior` |	Baseline prior 'pseudo-count' for genes NOT in the priors file.	| 1.0 |
| `--map-prior-strength` |	Multiplier for the influence of the priors file. | 1.0 |


#### **Performance**
| Flag | Description | Default |
| :--- | :--- | :--- |
| `-t`, `--threads` | Number of threads to use for computationally intensive steps. | `1` |

#### **Workflow Control**
| Flag | Description |
| :--- | :--- |
| `--overwrite` | Overwrite the output directory if it already exists. **Use with caution.** |
| `--skip-mapping` | Skip the initial (and longest) mapping steps. Requires that the pipeline has been run successfully on the same data before. |
| `--skip-to-map-resolve` | Skip all steps and run only the final MAP resolution. |
| `--debug` | Enable debug mode, which creates detailed log files for each sub-script. |
| `--version` | Show the pipeline version and exit. |
| `-h`, `--help` | Show the help message and exit. |

## Examples

### Example 1: Standard Paired-End Analysis
A standard run on paired-end data, using 32 threads and the default protein PID cutoff of 90% for HomScan.

```bash
amrscan -i sampleA_R1.fastq.gz,sampleA_R2.fastq.gz,sampleA_others.fastq.gz \
        -o SampleA_results \
        -t 32
```

### Example 2: Advanced Analysis with Custom Parameters
A more advanced run with several custom settings:
- Use a stricter **nucleotide** PID of 95% for HomScan.
- Use a very strict nucleotide PID of 99% for VarScan.
- Analyse both `H` and `K` gene types in HomScan.
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
Imagine the run from Example 1 finished, but you want to re-analyse the HomScan results with a stricter protein PID of 98% without re-doing the slow mapping step.

The combination of `--skip-mapping` and `--overwrite` is designed for this. It will safely delete the old analysis results while preserving the mapping data.

```bash
amrscan -i sampleA_R1.fastq.gz,sampleA_R2.fastq.gz \
        -o SampleA_results \
        -t 32 \
        --homscan-pid-cutoff 0.98 \
        --skip-mapping \
        --overwrite
```

### Example 4: Re-running only the MAP Resolver with Priors
After a run is complete, you might want to re-calculate the final abundances using a file of clinical priors. The `--skip-to-map-resolve` flag makes this fast and efficient.

```bash
amrscan -i sampleA_R1.fastq.gz,sampleA_R2.fastq.gz \
        -o SampleA_results \
        --map-priors-file /path/to/clinical_priors.tsv \
        --skip-to-map-resolve \
        --overwrite
```

## Output Files

All results will be located in the directory specified by `-o`. The key outputs are:

| File / Directory	| Description |
| :--- | :--- |
| `[prefix]_homscan_MAP.tsv` |	***(Primary Homology Result)*** The final, most accurate quantitative report for homology-based gene detection. It includes the original aggregated metrics for each gene family, plus three new columns: Resolved_[Metric] (e.g., Resolved_RPKG) shows the final abundance after MAP distribution; Top_ARO shows the single ARO within an ambiguous family that was allocated the most abundance; and Allocation_Proportions is a JSON string detailing how the abundance was split among all members of that ambiguous family. |
| `[prefix]_homscan.tsv`	| A "clean" summary where all reads mapping to multiple AROs within the same family are grouped into a single family;multiple line. |
| `[prefix]_varscan.tsv` |	***(Primary Variant Result)*** The final, normalised summary table for confirmed resistance variant detection. |
| `[prefix]_homscan_html/` |	A directory containing interactive HTML coverage plots for each detected AMR gene family. |
| `[prefix]_varscan_html/` |	A directory containing HTML alignment views for each detected AMR variant. |
| `logs/` |	A directory containing the main run log (_run.log) and a detailed run summary file (_run_details.txt). |
| `tmp/` |	Contains all intermediate files from the pipeline. Can typically be ignored but is useful for debugging. |


## License
This project is licensed under the MIT License.

## Contact
For questions, bug reports, or suggestions, please open an issue on this GitHub repository.