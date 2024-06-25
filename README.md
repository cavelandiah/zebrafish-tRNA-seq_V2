# zebrafish-tRNA-seq

This repository provides code to reproduce the results presented in [https://doi.org/10.1101/2024.01.30.575011 ] on OSX/Linux.

## download this workflow

Download this Snakemake wokflow, e.g. with git clone:
`git clone git@github.com:mwaldl/zebrafish-tRNA-seq.git`

## set up input files

Required inputs are the demultiplexed fastq sequencing files, a file with meta data for each sample and a config file that among other things specifies genome sources. To set up this input files:

- Create `raw/rawr_reads` folder and copy your (unziped) fastq files there. The fasq files for reproducing the results presented in Rappol et al. 2024 can be downloaded here: https://dataview.ncbi.nlm.nih.gov/object/PRJNA1061456.
- Set up the `config/samples.tsv` file with metadata for each sample. An example for reproducing the results presented in our paper is included with the worklflow. The columns `fastq`, `treatment`, `timepoint`, `timepoint_name` and `replicate`  are required. The name listed in the column `fastq` should correspond to the fastq file name (e.g. `EV01001` and `raw/rawr_reads/EV01001.R1.fastq`)
- Edit the `config/config.yaml` to reflect your desired parameters. (TODO: explain parameters)
- Canonical tRNA positions are annotated based on the tRNA Rfam alignment within `per_ref_nt_count` rule in `workflow/rules/coverage.smk`: The stockholm file `raw/canonical_tRNA/RF00005.stockholm_withcanonical annotation.txt` includes manually annotated canonical positions after the #=GC SS_cons line. Zebrafish tRNAs are aligned to the Rfam alignment and canonical positions are infered by the corresponding positions in X14835.1_6927-7002, which represents all canonical positions. The assignment can be changed by editing the `alignment_to_canonical_positions_mapper` function within the previously mentioned rule. No editing is requird if you want to use the annotated canonical positions as in our publication.
- The tRNA Rfam alignment and the corresponding covariance model need to be provided within the path specified in the `condig/config.yaml` file, by default `rfam_alignment: 'raw/canonical_tRNA/RF00005.stockholm.txt'` and `rfam_alignment: 'raw/canonical_tRNA/RF00005.cm'` (the version used in our paper is provided with the worklflow).
- The genomic tRNA refernce genome is generated based on tRNAscanSE scans and the zebrafish genome. Fasta files for extracting the genomic tRNAs are downloaded automatically. The assembly ID of the zebrafish genome needs to be provided in the config file (`assembly_id: 'GCF_000002035.6_GRCz11'`). The required annotations from tRNAscanSE can for example be downloaded from https://gtrnadb.ucsc.edu/genomes/eukaryota/Dreri11/ ("Download tRNAscan-SE Results"). The following three files need to be copied to the paths specified in the config file:
  - `hc_genomic_tRNAs_from_GtRNAdb: 'raw/references/danRer11-tRNAs/danRer11-mature-tRNAs.fa'`
  - `downloaded_tRNAscanSE_summary: 'raw/references/danRer11-tRNAs/danRer11-tRNAs-detailed.out'`
  - `downloaded_tRNAscanSE_name_mapping: 'raw/references/danRer11-tRNAs/danRer11-tRNAs_name_map.txt'`
- Mitochondrial tRNAs need to be provided at the path specified in the config file (eg. `mitochondrial_tRNAs_from_mt_tRNADB: 'raw/references/tRNAdb/mit_from_tRNAdb.fst'`, downloaded from tRNADB, included with the worklflow).
- A manually currated reference genome, including snps found during our analysis, is provided at `raw/references/manual.fa`.

## install dependencies

### dependencies are the following tools:

- Snakemake (https://snakemake.readthedocs.io/en/stable/#getting-started)
- fastp
- cutadapt
- umi_tools
- cmalign
- segemehl


### as well as python3 and the following python modules
- pandas
- sklearn
- numpy
- matplotlib
- seaborn
- yaml
- shutil
- biopython
- scipy
- BCBio

# How to run pipelin

call snakemake:
`snakemake 'results/all.txt' -c2`
