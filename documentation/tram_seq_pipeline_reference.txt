## TRAM-seq Pipeline Reference

A step-by-step guide to the rules, inputs, outputs, and key parameters in the TRAM-seq Snakemake workflow.

---

### 1. Configuration

```yaml
reads_filter: "umi"                # Strategy for handling UMIs (“umi” to collapse duplicates, “none” to keep all)
ref_set: "selected"                # Which reference set to use: “all” or “selected”
inosine_no_filter:                 # Anticodon families exempt from coverage filtering
  - Ala-AGC
  - Arg-ACG
  - Ile-AAT
  - Leu-AAG
  - Pro-AGG
  - Ser-AGA
  - Thr-AGT
  - Val-AAC

min_coverage_per_ref:              # Threshold for selecting high-coverage references
  - ["max all fraction", 0.0005]

poly_T_processing: "ONtttt"        # Trim poly-T tails of length ≥4 from reads
selected_alignment:                # Path to final alignments of selected tRNAs
  "resources/references/alignment/selected_tRNAs.fa"
min_raw_abundance_alignment:       # Path to alignment input for low-abundance filter (no CCA)
  "resources/references/alignment/min_raw_abundance_tRNAs.fa"
min_raw_abundance_refs:            # Path to raw low-abundance tRNA FASTA
  "resources/references/min_raw_abundance_tRNAs.fa"
rfam_cm:                           # Rfam covariance model for tRNA
  "raw/canonical_tRNA/RF00005.cm"

min_raw_abundance_count_for_alignment:
  - ["max all count", 10]          # Require ≥10 reads per tRNA for inclusion

raw_reads_folder:                  # Location of input FASTQ files
  "raw/raw_reads"
```

---

### 2. Reference Fasta Generation

#### 2.1 `all_refs.fa`
- **Rule:** `get_merged_tRNA_refs` (in `ref_genome.smk`)  
- **Input:** Annotated zebrafish tRNAs (genomic + mitochondrial)  
- **Output:**  
  ```
  resources/references/all_refs.fa
  ```

#### 2.2 `selected_refs.fa`
- **Rule:** `get_selected_refs`  
- **Inputs:**
  - `resources/references/all_refs.fa`
  - Coverage filter summary:
    ```
    resources/min_coverage_refs/pre-filter_${reads_filter}/min_cov_refs.yaml
    ```
- **Coverage filtering (`get_min_cov_refs` in `filter_mappings.smk`):**
  1. Reads from:
     ```
     resources/coverage/pre-filter_${reads_filter}/raw/min_coverage_summary_DM-with-max.tsv
     ```
  2. Apply `min_coverage_per_ref` (0.0005 fraction).
  3. Exempt `inosine_no_filter` anticodons from this filter.
- **Output:**
  ```
  resources/references/selected_refs.fa
  ```

---

### 3. Coverage Summary Computation

#### 3.1 Per-sample stats (`get_coverage_per_ref`)
- **Inputs:** Three SAM files per sample under:
  ```
  resources/filtered-mappings/pre-filter_${reads_filter}/raw/{unique,random,mapped}/${sample}.sam
  ```
- **Processing:**
  1. **Random**: Shuffle, keep last per QNAME.
  2. **Unique**: Drop all reads sharing a QNAME, keep only unique.
  3. **Mapped (all)**: Remove PCR duplicates by (UMI, RNAME), keep first.
  4. Concatenate counts into DataFrame and compute:
     - `unique count`, `random count`, `all count`
     - Fractions (each / its total across refs)
     - Ratios: `unique/random`, `random/all`
     - Extract `anticodon` from `RNAME`
- **Output:**
  ```
  resources/coverage/pre-filter_${reads_filter}/raw/per-ref/${sample}.tsv
  ```

#### 3.2 Merge per-sample (`get_coverage_summary`)
- **Inputs:** All per-ref TSVs for DM samples  
- **Output:**
  ```
  resources/coverage/pre-filter_${reads_filter}/${ref_set}/coverage_summary_DM.tsv
  ```

#### 3.3 Add min/max (`get_min_coverage_summary_with_max`)
- **Input:** `coverage_summary_DM.tsv`  
- **Output:**
  ```
  resources/coverage/pre-filter_${reads_filter}/raw/min_coverage_summary_DM-with-max.tsv
  ```
- **Adds:** Minimum and maximum per metric across samples

---

### 4. Read Mapping

#### 4.1 Pre-mapping & Trimming (`mapping_segemehl`)
- **Rule:** `mapping_segemehl` (in `mapping.smk`)  
- **Reference index:**  
  ```
  resources/index/${ref_set}.{idx,ctidx,gaidx}
  ```
- **Input FASTQ:**
  ```
  resources/filtered-reads/a_u_t_q_${sample}.fastq
  ```
  - Poly-T trimming (`poly_T_processing = "ONtttt"`)
- **Segemehl params:**  
  - If `ref_set == "all"`, use `config['pre_mapping']`; else `config['mapping']`.
- **Outputs:**
  - Mapped SAM:
    ```
    resources/mapping/${ref_set}/polyTtrimming_${poly_T}/${sample}.sam
    ```
  - Unmapped FASTQ & stats TSV

#### 4.2 UMI Deduplication (`filter_umi`)
- **Input:**
  ```
  resources/filtered-mappings/pre-filter_none/raw/mapped/${sample}.sam
  ```
- **Action:** Drop duplicates by `(umi, RNAME)`, keep first.
- **Output:**
  ```
  resources/filtered-mappings/pre-filter_umi/raw/mapped/${sample}.sam
  ```

---

### 5. Final tRNA Alignments

#### 5.1 Sort Selected Alignments (`get_sorted_selected_alignment`)
- **Rule:** `get_sorted_selected_alignment` (in `alignments.smk`)
- **Inputs:**
  - `config['selected_alignment']`  
  - `min_cov_refs.yaml` from step 2
- **Output:**
  ```
  resources/references/alignment/selected_tRNAs.fa_sorted.fa
  ```

#### 5.2 Add CCA Tail (`add_CCA_to_alignment`)
- **Inputs:**
  - `min_cov_refs.yaml`
  - `min_raw_abundance_tRNAs.fa` (no-CCA alignment via `cmalign`)
- **Output:**
  ```
  resources/references/alignment/selected_tRNAs.fa
  ```

##### Building `min_raw_abundance_tRNAs.fa`
- **Rule:** `get_min_raw_coverage_fasta`
- **Inputs:**
  - `resources/min_coverage_refs/pre-filter_umi/min_cov_refs_alignment.yaml`
  - `resources/references/all_refs.fa`
- **Filter:** Keep tRNAs with `max all count ≥ 10`

---

### 6. Quality Control & Reports

| Report                     | Rule                                 | Output                                                                                 |
|----------------------------|--------------------------------------|----------------------------------------------------------------------------------------|
| Edit-distance histogram    | `plot_dist_hist` (in `clustering.smk`)| `qc/cluster/${ref_set}/editdist_hist.pdf`                                             |
| FastQC (raw reads)         | `raw_fastqc` (`preprocessing_reads.smk`)| `qc/fastqc/raw_${sample}.json`                                                         |
| Adapter & quality trimming | `trim_adapter_only`, `quality_filter_after_t_trimming` | `resources/filtered-reads/a_${sample}.fastq` and `qc/fastqc/a_u_t_q_${sample}.html` |
| Preprocessing summary      | `preprocessing_summary` (`qc.smk`)   | `qc/preprocessing_summary/preprocessing.pdf`                                           |
| Mapping summary            | `mapping_summary_all_samples`        | `qc/mapping-summary/mapping-${ref_set}.pdf`                                           |
| Duplication rate           | `duplication_rate` (`duplication_rate.smk`) | `qc/umi/${ref_set}_duplication_rate.pdf`                                               |
| Multimapper vs. dist       | `multimapper_vs_dist` (`clustering.smk`)| `qc/cluster/pre-filter_${reads_filter}/${ref_set}/dist_multimapper.pdf`                |

---

### 7. Downstream Analyses

- **Abundance PCA**  
  - **Rule:** `get_all_abundance_pca_plots`  
  - **Output:**  
    ```
    results/abundance/${reads_filter}_${ref_set}_pca-done.txt
    ```

- **Modification & mismatch plots**  
  - **Rule:** `get_all_mismatch_plots_CAVH` (`modifications.smk`)  
  - **Outputs:**  
    - Heatmaps and line-plots under  
      `results/modifications/pre-filter_${reads_filter}/${ref_set}/…`

- **Coverage plots**  
  - **Rule:** `get_all_coverage_plots` (`coverage_plots.smk`)  
  - **Output:**  
    ```
    results/coverage_plots/done.txt
    ```

---

## Notes & Clarifications

- **`reads_filter`**:  
  - `"umi"` collapses by UMI + reference  
  - `"none"` keeps all reads  

- **`ref_set`**:  
  - `"all"` uses every annotated tRNA  
  - `"selected"` uses only high-coverage subset  

- **`poly_T_processing = "ONtttt"`**:  
  Trim reads ending in ≥4 consecutive T’s (common RT artifact).

- **Raw vs. Selected folders**:  
  Early filters operate on the complete (“raw”) data; later steps restrict to the “selected” references.

- **Hard-coded thresholds**:  
  - `min_coverage_per_ref = 0.0005` (fraction)  
  - `min_raw_abundance_count_for_alignment = 10` (reads)  

Consider centralizing these parameters in `config.yaml` and documenting their biological rationale.
