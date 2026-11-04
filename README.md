# QuicBac

A streamlined Python pipeline designed to identify, extract, and taxonomically classify 16S rRNA genes from metagenomic assemblies. 

The tool uses `barrnap` to locate rRNA genes, `bedtools` to extract the 16S sequences, and `blastn` against the NCBI 16S RefSeq database to assign taxonomy. Finally, it uses the NCBI Entrez API to resolve the full lineage (Phylum, Class, Order, Family) and generates publication-ready stacked bar charts of the taxonomic composition.
<img width="2752" height="1536" alt="Gemini_Generated_Image_hkvpfjhkvpfjhkvp" src="https://github.com/user-attachments/assets/4ab1cc2c-aee0-49e0-8962-9d0e7ea724d9" />

## Features
* **Automated Extraction:** Finds and extracts full and partial 16S sequences directly from your assembly FASTA.
* **High-Confidence Taxonomy:** Uses BLASTn and the Yarza et al. (2014) identity thresholds for robust taxonomic classification.
* **Lineage Resolution:** Automatically fetches high-level taxonomy (Phylum through Family) via NCBI E-utilities.
* **Visualizations:** Automatically generates stacked bar plots of the community composition.

## Prerequisites & Installation

It is highly recommended to use `conda` or `mamba` to install the required dependencies.

1. **Create and activate a new environment:**
```bash
conda create -n 16s_pipeline -c bioconda -c conda-forge barrnap bedtools "blast>=2.14.0" pandas matplotlib requests
conda activate 16s_pipeline
```
*(Note: We specify BLAST >= 2.14.0 to avoid integer overflow issues with modern NCBI databases).*
Clone the repository.
```bash
git clone https://github.com/BenMinch/quicbac
```

## Database Setup

Before running the pipeline for the first time, you must download the NCBI 16S ribosomal RNA database.

1. **Download and decompress the database:**
```bash
mkdir -p ~/blastdb
cd ~/blastdb
update_blastdb.pl --decompress 16S_ribosomal_RNA taxdb
```

## Usage

Run the pipeline by providing your input metagenomic assembly FASTA file and an output prefix.

```bash
python quicbac.py -i assembly.fasta -o my_results -t 16 --db blastdb/16S_ribosomal_RNA
```

### Optional Arguments
* `-t, --threads`: Number of threads to use for Barrnap and BLAST (default: 8).
* `--db`: The name or absolute path of the BLAST database to use (default: `16S_ribosomal_RNA`).
* `--min-qcov`: Minimum query coverage percentage for a valid hit (default: 80.0).

## Output Files

The pipeline generates the following files based on your provided `-o` prefix (e.g., `my_results`):

1. **`my_results.16S.fasta`**: The extracted 16S rRNA nucleotide sequences.
2. **`my_results.taxonomy.tsv`**: A tab-separated table containing the BLAST metrics, taxonomic assignments, NCBI TaxIDs, and full lineage information.
3. **`my_results_tax_barplot.pdf`**: A figure containing four stacked bar charts representing the community composition at the Phylum, Class, Order, and Family levels.

## Taxonomic Confidence Thresholds
Taxonomic ranks are assigned dynamically based on percent identity (Yarza et al., 2014):
* `species`: $\ge$ 98.65%
* `genus`: $\ge$ 94.5%
* `family`: $\ge$ 86.5%
* Sequences falling below these thresholds are marked as `unclassified_below_family`.
* Sequences with query coverage below the `--min-qcov` threshold are marked as `low_coverage`.
