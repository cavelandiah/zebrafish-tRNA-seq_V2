# zebrafish-tRNA-seq

This repository provides code to reproduce the results on OSX/Linux presented in [https://doi.org/10.1101/2024.01.30.575011 ].

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

create `raw/rawr_reads` folder and copy your fastq files there. The fasq files for reproducing the results presented in Rappol et al. 2024 can be downloaded here: https://dataview.ncbi.nlm.nih.gov/object/PRJNA1061456.

set up the samples.csv file with metadata for each sample

call snakemake:
`snakemake 'results/all.txt' -c2`
