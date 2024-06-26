# zebrafish-tRNA-seq

This repository provides code to reproduce the results presented in [https://doi.org/10.1101/2024.01.30.575011](https://doi.org/10.1101/2024.01.30.575011) on OSX/Linux.

## Install Dependencies

Please install the following software:

- [Snakemake](https://snakemake.readthedocs.io/en/stable/#getting-started)
- fastp
- cutadapt
- umi_tools
- cmalign
- segemehl
- python

and the following Python 3 modules:

- pandas
- sklearn
- numpy
- matplotlib
- seaborn
- yaml
- biopython
- scipy
- BCBio

We recommend installation via Conda, as explained below.

### Install Conda

The installation instructions are linked below:

[Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

### Setting Up the Conda Environment

1. **Create a new Conda environment with all necessary packages:**

    ```bash
    conda create -n zebrafish python=3.8 \
        bioconda::snakemake \
        bioconda::fastp bioconda::cutadapt bioconda::umi_tools bioconda::cmalign bioconda::segemehl \
        conda-forge::numpy conda-forge::pandas conda-forge::matplotlib-base conda-forge::seaborn-base conda-forge::pyyaml \
        conda-forge::biopython conda-forge::scipy bioconda::bcbio-nextgen conda-forge::scikit-learn
    ```

2. **Activate the Conda environment:**

    ```bash
    conda activate zebrafish
    ```

## Download this Workflow

Download this Snakemake workflow, e.g., with `git clone`:

```bash
git clone git@github.com:mwaldl/zebrafish-tRNA-seq.git
```

## Set Up Input Files

Required inputs are the demultiplexed fastq sequencing files, a file with meta data for each sample and a config file that among other things specifies genome sources. To set up this input files:

- Create `raw/rawr_reads` folder and copy your (unzipped) FASTQ files there. The FASTQ files for reproducing the results presented in Rappol et al. 2024 can be downloaded here: https://dataview.ncbi.nlm.nih.gov/object/PRJNA1061456.
- Set up the c`config/samples.tsv` file with metadata for each sample. An example for reproducing the results presented in our paper is included with the workflow. The columns `fastq`, `treatment`, `timepoint`, `timepoint_name` and `replicate` are required. The name listed in the column `fastq` should correspond to the FASTQ file name (e.g., `EV01001` and `raw/rawr_reads/EV01001.R1.fastq`).
- Edit the `config/config.yaml` to reflect your desired parameters (see comments in config file for details).
- Canonical tRNA positions are annotated based on the tRNA Rfam alignment within the `per_ref_nt_count` rule in `workflow/rules/coverage.smk`: The Stockholm file `raw/canonical_tRNA/RF00005.stockholm_withcanonical annotation.txt` includes manually annotated canonical positions after the #=GC SS_cons line. Zebrafish tRNAs are aligned to the Rfam alignment and canonical positions are inferred by the corresponding positions in X14835.1_6927-7002, which represents all canonical positions. The assignment can be changed by editing the `alignment_to_canonical_positions_mapper` function within the previously mentioned rule. No editing is required if you want to use the annotated canonical positions as in our publication.
- The tRNA Rfam alignment and the corresponding covariance model need to be provided within the path specified in the `condig/config.yaml` file. By default:
  - `rfam_alignment: 'raw/canonical_tRNA/RF00005.stockholm.txt'`
  - `rfam_cm: 'raw/canonical_tRNA/RF00005.cm'`
    (The version used in our paper is provided with the workflow).
- The genomic tRNA reference genome is generated based on tRNAscanSE scans and the zebrafish genome. Fasta files for extracting the genomic tRNAs are downloaded automatically. The assembly ID of the zebrafish genome needs to be provided in the config file (`assembly_id: 'GCF_000002035.6_GRCz11'`). The required annotations from tRNAscanSE can, for example, be downloaded from [gtrnadb.ucsc.edu](https://gtrnadb.ucsc.edu/genomes/eukaryota/Dreri11/) ("Download tRNAscan-SE Results"). The following three files need to be copied to the paths specified in the config file:
  - `hc_genomic_tRNAs_from_GtRNAdb: 'raw/references/danRer11-tRNAs/danRer11-mature-tRNAs.fa'`
  - `downloaded_tRNAscanSE_summary: 'raw/references/danRer11-tRNAs/danRer11-tRNAs-detailed.out'`
  - `downloaded_tRNAscanSE_name_mapping: 'raw/references/danRer11-tRNAs/danRer11-tRNAs_name_map.txt'`
- Mitochondrial tRNAs need to be provided at the path specified in the config file (e.g., `mitochondrial_tRNAs_from_mt_tRNADB: 'raw/references/tRNAdb/mit_from_tRNAdb.fst'`, downloaded from tRNADB, included with the workflow).
- A manually curated reference genome, including SNPs found during our analysis, is provided at `raw/references/manual.fa`. This can be edited to fit your needs. Make sure you follow the same nameing scheme.

## How to Run the Pipeline

After you downloaded the worklflow, edited inoutfiles and parameters as needed and installed the dependencies via Conda, you can 

1. move to the worklflow directory

    ```bash
    cd zebrafish-tRNA-seq

    ```

2. activate the conda enviornment (if not already active)
    ```bash
    conda activate zebrafish

    ```

3. execute the pipeline
    ```bash
    snakemake 'results/all.txt' -c2

    ```
    Replace '2' with the number of cores you wish to use, at least 2.
