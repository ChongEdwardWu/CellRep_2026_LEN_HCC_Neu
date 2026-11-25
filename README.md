# CellRep_2026_LEN_HCC_Neu

Analysis code for mouse single-cell and bulk RNA sequencing, human nCounter immune signatures, spatial cell distances, and TCGA-LIHC survival analyses.

The repository contains path-sanitized code for preprocessing, quality control, Seurat integration, immune-cell annotation, neutrophil reference mapping, pySCENIC analysis, differential expression, enrichment analysis, and figure generation.

## Repository structure

```text
LEN_HCC_Neu/
├── R/
│   ├── QC/
│   │   ├── scRNA_QC_loop_CB_step1.R
│   │   └── scRNA_QC_loop_CB_step2.R
│   ├── BulkRNA/
│   │   └── Step01_Neu_RNA_DEG.R
│   ├── Human/
│   │   ├── Step1_Human_nCounter_signature_scores.R
│   │   ├── Step2_Human_PD_delta_heatmaps.R
│   │   ├── Step2b_Human_PD_signature_statistics.R
│   │   ├── Step3_Human_PD_signature_ranking.R
│   │   ├── Paired_Human_nCounter_panels.R
│   │   └── Human_Fig1C_PFS.R
│   ├── Spatial/
│   │   └── Distance_and_phenotype_analysis.R
│   ├── Survival/
│   │   └── KM_survival.R
│   └── Analysis/
│       ├── Step1_Len_Data_Integration.R
│       ├── Step1_5_Len_pySCENIC.R
│       ├── Step2_Len_Clustering.R
│       ├── Step3_Len_Neu_sub_clustering.R
│       ├── Step4_Len_Neu_DEG.R
│       ├── Step5_Len_LENvsCTR_Figures.R
│       ├── Step6_Len_Neu_query_native_display.R
│       └── Step7_Len_CD45_lineage_remodeling.R
├── python/
│   └── run_SCENIC_mouse.py
├── shell/
│   └── run_cellbender_batch.sh
├── env/
│   └── .gitkeep
├── run_clinical.R
└── README.md
```

## Analysis workflow

### Core workflow

1. `shell/run_cellbender_batch.sh`
   - Runs `cellbender remove-background` on Cell Ranger raw-feature matrices.

2. `R/QC/scRNA_QC_loop_CB_step1.R`
   - Reads CellBender-filtered matrices, calculates sample-level quality-control metrics, performs exploratory clustering, and calls doublets with `scDblFinder`.

3. `R/QC/scRNA_QC_loop_CB_step2.R`
   - Applies final filtering, harmonizes gene symbols, calculates cell-cycle scores, and creates sample-level Seurat objects with `SCTransform`.

4. `R/Analysis/Step1_Len_Data_Integration.R`
   - Performs multi-sample SCT-RPCA integration, dimensionality reduction, and graph-based clustering.

5. `R/Analysis/Step1_5_Len_pySCENIC.R` and `python/run_SCENIC_mouse.py`
   - Export merged RNA counts, run pySCENIC, and add regulon activity to the Seurat object.

6. `R/Analysis/Step2_Len_Clustering.R`
   - Assigns broad immune-cell annotations and summarizes cluster composition.

7. `R/Analysis/Step3_Len_Neu_sub_clustering.R`
   - Subsets neutrophils, maps discrete Ng et al. maturation states, and prepares the neutrophil object used by downstream analyses.

8. `R/Analysis/Step4_Len_Neu_DEG.R`
   - Performs neutrophil differential-expression, regulon-activity, score, and gene-set enrichment analyses.
   - For binary regulon activity, the two-sided Wilcoxon rank-sum test accounts for ties at both 0 and 1 and uses a continuity correction. Bonferroni adjustment uses the full number of regulons in the assay (214 for Data S3).

9. `R/Analysis/Step5_Len_LENvsCTR_Figures.R`
   - Generates neutrophil-focused figure panels.

10. `R/Analysis/Step6_Len_Neu_query_native_display.R`
    - Recalculates a query-native neutrophil UMAP, performs de novo five-cluster analysis, retains upstream Ng state labels, transfers tumour-derived Xue et al. mouse neutrophil states, and calculates the Xue core-state maturation score.

11. `R/Analysis/Step7_Len_CD45_lineage_remodeling.R`
    - Performs matched-background differential-expression analyses across eight immune lineages, identifies genes with concordant attenuation under PERK inhibition and myeloid Eif2ak3 deletion, runs GO Biological Process enrichment, and generates the myeloid and lymphoid Venn/GO displays.

### Human nCounter workflow

Run the clinical workflow from the repository root with R and the packages `readxl`, `ggplot2`, `pheatmap`, and `survival`. The command below uses `nCounter_Normalization.xlsx` from the [Yamauchi et al. data deposit](https://doi.org/10.6084/m9.figshare.22215760.v3) and Supplementary Data 3 (`43856_2023_390_MOESM5_ESM.xlsx`) from [the original article](https://doi.org/10.1038/s43856-023-00390-x).

```bash
Rscript run_clinical.R "data/human/nCounter_Normalization.xlsx" "data/human/43856_2023_390_MOESM5_ESM.xlsx" "results/human" "data/human/43856_2023_390_MOESM5_ESM.xlsx"
```

The first argument supplies normalized expression (sheet `NC_9785_N1_Normalization`). The second supplies clinical annotations (sheet `Fig1b_annotation` or `Patient annotation`). The third is a new output directory. The fourth identifies biopsy timepoints: the `Fig3b` sheet in Supplementary Data 3 identifies the eight week-4 biopsies among the 20 follow-up samples. The remaining 12 are progression biopsies. Alternatively, supply a CSV with one row per expression sample and the columns `sample_id`, `patient_id`, and `timepoint` (`baseline`, `week4`, or `PD`).

| Panels | Calculation |
| --- | --- |
| S1A | Marker color shows log2(mean normalized expression + 1). |
| 1A, S1B | Steps 1 and 2 average log2-transformed marker expression within each signature, standardize the signature against 39 baseline biopsies using the sample SD, and calculate within-patient changes from baseline to progression. |
| S1B | Step 2b compares changes between the PD and non-PD response groups. It computes the rank-biserial effect from all cross-group patient pairs and exact within-group percentile bootstrap intervals for the difference of medians. Step 3 orders signatures by the composite score and writes the myeloid panel to `step3/SFig_rank_lollipop_MYeloid.pdf`. |
| 1B, S1C–F | Each log-transformed marker is standardized across 40 biopsies from 20 paired patients using the population SD. Marker Z scores are averaged within signatures. Baseline and progression scores are compared by a paired t test within each best-response group (3 PD and 9 non-PD patients). |
| 1C | The Kaplan–Meier curves and log-rank test use all 20 paired patients, grouped by whether their baseline-referenced neutrophil signature increased. |

### Spatial distance analysis

`R/Spatial/Distance_and_phenotype_analysis.R` calculates neutrophil-to-nearest-CD8Tc distances within each `QupathID` and scales them to the within-sample maximum (0–100). The input RDS contains `QupathID`, `SampleName`, `Cell ID` (or `CellID`), `Cell X Position`, `Cell Y Position`, and `Phenotype`, with labels `Neutrophil` and `CD8Tc`. It exports a distance CSV and session information.

```bash
Rscript R/Spatial/Distance_and_phenotype_analysis.R "data/spatial/cells.rds" "results/spatial"
```

### TCGA-LIHC survival analysis

`R/Survival/KM_survival.R` combines CIBERSORT neutrophil fractions and PERK-UPR ssGSEA scores for Kaplan–Meier overall-survival analysis, using optimized cutpoints (`minprop = 0.3`) and 72-month follow-up. Inputs are a `log2(TPM + 1)` expression TSV (Ensembl IDs in the first column and TCGA sample barcodes as column names) and a survival TSV with `sample`, `OS` (0/1), and `OS.time` (days). It exports Kaplan–Meier PDFs, scores, groups, cutpoints, and session information.

```bash
Rscript R/Survival/KM_survival.R "data/TCGA-LIHC.star_tpm.tsv" "data/TCGA-LIHC.survival.tsv" "results/TCGA_LIHC"
```

Both scripts write to a new or empty output directory.

## Main software dependencies

### R

- `Seurat`, `SeuratObject`, `SeuratWrappers`
- `SingleCellExperiment`, `scater`, `scran`, `scuttle`, `DropletUtils`
- `scDblFinder`, `bluster`, `BiocSingular`, `BiocParallel`
- `SingleR`, `UCell`, `destiny`, `presto`
- `SCENIC`, `SCopeLoomR`
- `msigdbr`, `fgsea`, `clusterProfiler`, `org.Mm.eg.db`, `AnnotationDbi`, `GOSemSim`
- `tidyverse`, `Matrix`, `openxlsx2`, `patchwork`, `ggridges`, `ggVennDiagram`, `ragg`
- `readxl`, `ggplot2`, `pheatmap`, `survival` (human nCounter workflow)
- `phenoptr`, `rtree` (spatial distances)
- `data.table`, `IOBR`, `GSVA`, `msigdbr`, `org.Hs.eg.db`, `AnnotationDbi`, `survival`, `survminer`, `maxstat`, `ggplot2` (TCGA-LIHC survival)

### Python and command-line tools

- `cellbender`
- `pyscenic`
- `loompy`
- `pandas`

Step 6 requires Seurat 5.1 or later.

## External inputs

- Cell Ranger and CellBender outputs
- pySCENIC motif annotations, transcription-factor lists, and ranking databases
- a processed Ng et al. neutrophil reference object
- the Ng et al. supplementary gene list used by the original Step 3 workflow
- a processed Xue et al. mouse liver-tumour reference object
- the annotated CD45-positive Seurat object required by Step 7
- any additional gene-set resources referenced by the core workflow

Replace placeholder paths in the core scripts and supply Step 6 and Step 7 paths through their command-line arguments.

