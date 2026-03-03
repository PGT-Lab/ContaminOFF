# <img src="contaminoff_figure.png" alt="ContaminOFF Logo" width="100" align="left"/> ContaminOFF
**Ultra-Fast Metagenomic Contamination Filter for Human Sequencing Data**
<br clear="left"/>
<br>

ContaminOFF is a high-performance Python tool designed to rapidly identify and remove bacterial contamination from human whole-genome sequencing (WGS) data, specifically optimized for high-noise sources like saliva samples. 

By coupling Kraken2's taxonomic classification with highly optimized, byte-level I/O processing using `pigz`, ContaminOFF extracts target sequences (e.g., *Homo sapiens*) up to 4x faster than standard parsing tools. It minimizes RAM usage and saves days of computational time when processing large cohorts.


## The ContaminOFF Saliva Database

We provide the **ContaminOFF Saliva DB** — a strictly curated, highly specific database tailored for the oral/salivary microbiome. 
* **Zero Human Cross-Talk:** Masked against the human reference genome to prevent false-positive bacterial assignments.
* **Lightweight & Fast:** Contains only the essential clinically relevant taxa found in oral environments, drastically reducing RAM requirements during classification.

**Download the ContaminOFF Saliva Database:**
```bash
wget [https://figshare.com/articles/dataset/ContaminOFF_saliva_microbiome_and_Homo_sapiens_database_for_kraken2_/31460278](https://figshare.com/articles/dataset/ContaminOFF_saliva_microbiome_and_Homo_sapiens_database_for_kraken2_/31460278)
tar -xzvf contaminoff_db.tar.gz
```

## Key Features

* **Ultra-Fast I/O:** Bypasses standard Python text-parsing bottlenecks by using multithreaded `pigz` and byte-level matching ($O(1)$ complexity).
* **Dual Processing Modes:** Process a single sample pair or automatically batch-process an entire directory of FASTQ files.
* **Automated QC Reporting:** Generates high-resolution composition pie charts, taxonomic abundance bar plots, and TSV summaries.

## Dependencies

ContaminOFF requires the following tools to be accessible in your system's PATH:
* `kraken2`
* `pigz`
* Python 3.6+ with libraries: `pandas`, `matplotlib`, `seaborn`

**Quick Conda Installation:**
```bash
conda create -n contaminoff_env -c conda-forge -c bioconda kraken2 pigz pandas matplotlib seaborn python
conda activate contaminoff_env
```

## Quick Start & Usage

ContaminOFF offers two flexible ways to process your data, depending on your pipeline needs. Both modes support all taxonomic tuning and plotting flags.

**Note: All modes accept both `.fq.gz` and `.fastq.gz` file extensions, and strictly require the files to be gzip compressed (`.gz`).**

### 1. Single-Sample Mode
Use this mode to process a specific pair of paired-end FASTQ files.
```bash
python contaminoff.py \
    -r1 WC-001_R1.fq.gz \
    -r2 WC-001_R2.fq.gz \
    -o ./clean_results \
    -d /path/to/ContaminOFF_DB \
    -t 16 \
    -p WC-001_report \
    --tax-level G \
    --top-taxa 20
```

### 2. Multi-Sample (Batch) Mode
Point ContaminOFF to a directory, and it will automatically find all paired `_R1`/`_R2` (or `_1`/`_2`) files and process them sequentially without overflowing your RAM.
```bash
python contaminoff.py \
    -i /path/to/raw_fastq_folder/ \
    -o ./clean_cohort_results \
    -d /path/to/ContaminOFF_DB \
    -t 16 \
    --tax-level G
    --top-taxa 20
```

## Expected Output
For every sample processed, ContaminOFF generates four optimized files in your output directory:
1. **`*_clean_R1.fq.gz` & `*_clean_R2.fq.gz`**: The filtered, purely human FASTQ files, highly compressed via `pigz`.
2. **`*.report`**: The standard Kraken2 taxonomic report.
3. **`*_summary.tsv`**: A data table containing the exact percentage of Human, Bacteria, and Unclassified reads, plus the absolute number of extracted reads.
4. **`*_QC_plots.png`**: A image containing a composition pie chart and a bar plot of the top contaminating taxa.

## Command-Line Arguments

| Argument | Description | Default |
| :--- | :--- | :--- |
| **Input Modes** | *(Must choose one)* | |
| `-r1`, `--read1` | Input FASTQ R1 (.gz) [Single Mode] | - |
| `-i`, `--input-dir` | Directory containing FASTQ pairs (.gz) [Multi Mode] | - |
| **Core Arguments** | | |
| `-r2`, `--read2` | Input FASTQ R2 (.gz) [Required if using `-r1`] | - |
| `-o`, `--out-dir` | Output directory for clean FASTQs and reports | **Required** |
| `-d`, `--db` | Path to your Kraken2 database | **Required** |
| `-t`, `--threads` | Number of CPU threads to use | `8` |
| `-p`, `--prefix` | Output prefix (Ignored in Multi Mode) | `contaminoff_out` |
| **Tuning & QC** | | |
| `--taxid` | Target TaxID to keep (default is *Homo sapiens*) | `9606` |
| `--tax-level` | Taxonomic level for the QC bar plot (D, P, C, O, F, G, S) | `G` (Genus) |
| `--top-taxa` | Number of top contaminants to show in the QC plot | `10` |
| `--keep-kraken` | Flag to prevent deletion of the massive `.kraken` file | False |

## Citation

If you use ContaminOFF in your research, please cite our upcoming paper:
*Centurion VB et al., ContaminOFF: Ultra-Fast Metagenomic Contamination Filter for Human Sequencing Data

## License and Third-Party Software

ContaminOFF is distributed under a GPL-3 license. Additionally, ContaminOFF relies on the following third-party software and libraries to function:
* [Kraken2](https://ccb.jhu.edu/software/kraken2/)
* [pigz](https://zlib.net/pigz/)
* [Pandas](https://pandas.pydata.org/)
* [Matplotlib](https://matplotlib.org/)
* [Seaborn](https://seaborn.pydata.org/)

## Appreciation

This tool was developed by the Psychiatric Genetics Team (PGT) at the Universidade Federal de São Paulo (UNIFESP) as part of the research for the **Brazilian High-Risk Cohort (BHRC)** project.

This work was supported by the **São Paulo Research Foundation (FAPESP)** under grant numbers: **2021/05332-8** and **2025/02111-1**.
