<p align="center">
  <img src="slogan_contaminoff_transparent.jpg" alt="ContaminOFF Logo" width="600"/>
</p>

# ContaminOFF
**Ultra-Fast Metagenomic Contamination Filter for Human Sequencing Data**

ContaminOFF is a high-performance Python tool designed to rapidly identify and remove bacterial contamination from human whole-genome sequencing (WGS) data, specifically optimized for high-noise sources like saliva samples. 

By coupling Kraken2's taxonomic classification with highly optimized, byte-level I/O processing using `pigz`, ContaminOFF extracts target sequences (e.g., *Homo sapiens*) up to 4x faster than standard parsing tools. It minimizes RAM usage and saves days of computational time when processing large cohorts.

## ✨ Key Features
* ⚡ **Ultra-Fast I/O:** Bypasses standard Python text-parsing bottlenecks by using multithreaded `pigz` and byte-level matching ($O(1)$ complexity).
* 🔄 **Dual Processing Modes:** Process a single sample pair or automatically batch-process an entire directory of FASTQ files.
* 📊 **Automated QC Reporting:** Generates high-resolution composition pie charts, taxonomic abundance bar plots, and TSV summaries on the fly.
* 🔒 **Strict Validation:** Enforces `.gz` compression to protect users from accidental I/O system crashes and storage bloat.

## 🛠️ Dependencies
ContaminOFF requires the following tools to be accessible in your system's PATH:
* `kraken2`
* `pigz`
* Python 3.6+ with libraries: `pandas`, `matplotlib`, `seaborn`

**Quick Conda Installation:**
```bash
conda create -n contaminoff_env -c conda-forge -c bioconda kraken2 pigz pandas matplotlib seaborn python
conda activate contaminoff_env




conda create -n contaminoff_env -c conda-forge -c bioconda kraken2 pigz pandas matplotlib seaborn python
conda activate contaminoff_env
