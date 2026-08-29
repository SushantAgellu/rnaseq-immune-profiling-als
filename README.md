# RNA-seq Immune Profiling in ALS

**Semester 3 — Data Analysis Course Project (SLU)**  
Differential gene expression re-analysis of peripheral immune cells in **Amyotrophic Lateral Sclerosis (ALS)** using the publicly available **GSE60424** RNA-seq dataset.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**[View the full rendered analysis report](https://sushantagellu.github.io/rnaseq-immune-profiling-als/03_DataAnalysis/02_Visualization.html)** — no cloning or knitting required.

## TL;DR

Reprocessed public RNA-seq data (GSE60424) to compare **7 sorted immune cell types** from ALS patients vs. healthy controls using **DESeq2**. Neutrophils and monocytes carried 2–3× more differentially expressed genes than the adaptive immune cell types, with a consistent **down-regulation of interferon-response genes** (IFI44, ISG15, OAS2/3, etc.) shared across both — pointing to a blunted innate antiviral response in ALS. 17 genes replicate as "consensus DEGs" across both cell types.

## Key Results

| Global PCA | Neutrophil Volcano | Consensus Heatmap |
|---|---|---|
| ![Global PCA](03_DataAnalysis/02_Visualization/02_PCA_Global.png) | ![Neutrophil Volcano](03_DataAnalysis/02_Visualization/06_Volcano_Neutrophils.png) | ![Consensus Heatmap](03_DataAnalysis/02_Visualization/10_Consensus_Heatmap_Neutrophils_Monocytes.png) |
| Cell-type identity dominates global variance (PC1 = 56%), motivating a per-cell-type analysis | Neutrophils show the largest DEG signature (38 genes) of any cell type | 17 genes replicate across Neutrophils and Monocytes — the highest-confidence ALS signature in this study |

See [Section 9](#9-key-findings) for the full DEG table across all 7 cell types.

## Table of Contents

1. [What is your Question?](#1-what-is-your-question)
2. [What is your Response Variable?](#2-what-is-your-response-variable)
3. [What are your Predictor Variables?](#3-what-are-your-predictor-variables)
4. [Unit of Replication](#4-unit-of-replication)
5. [Statistical Approach](#5-statistical-approach)
6. [Repository Structure](#6-repository-structure)
7. [How to Reproduce](#7-how-to-reproduce)
8. [Key Design Decisions](#8-key-design-decisions)
9. [Key Findings](#9-key-findings)
10. [Data Source](#10-data-source)
11. [Key Versions](#11-key-versions)
12. [References](#12-references)
13. [Code of Ethics & Data Use Compliance](#13-code-of-ethics--data-use-compliance)
14. [License](#14-license)
15. [Author](#15-author)


## 1. What is your Question?

**Which genes are differentially expressed in circulating immune cells of ALS patients compared to healthy controls, and which cell types carry the strongest transcriptional signature?**

Secondary question: *Are any differentially expressed genes shared across multiple immune cell types (i.e., consensus DEGs), making them more robust candidate markers?*



## 2. What is your Response Variable?

**Dependent response variable — Gene expression level:**
- Measured as RNA-seq read counts
- One response value per gene per sample
- Discrete, non-negative integers
- Modeled with a **negative binomial distribution** via DESeq2 (not Gaussian — counts violate normality and are overdispersed)



## 3. What are your Predictor Variables?

**Primary predictor (independent variable of interest):**
- `disease` — two-level factor: `ALS` vs `Healthy Control`

**Stratification variable** (analysis run separately within each level):
- `celltype` — 7 levels: Neutrophils, Monocytes, CD4 T cells, CD8 T cells, B cells, NK cells, Whole Blood

**Available covariates** (not used in the primary model due to small sample size):
- `age` (years)
- `gender` (Male / Female)
- `race`

**Design formula used:** `~ disease`



## 4. Unit of Replication

- **One individual RNA-seq library = one unit of replication**
- Total libraries in GSE60424: **134**
- Libraries used in this ALS-vs-HC analysis: **49** (3 ALS + 4 HC per cell type × 7 cell types)
- Each library contributes one expression profile of ~50,045 genes (Ensembl IDs)



## 5. Statistical Approach

**Test:** Wald test on a negative binomial generalized linear model (DESeq2 v1.48.1)

**Why DESeq2:**
- RNA-seq counts are integer-valued and overdispersed → t-test, Wilcoxon, and Poisson are poor fits
- Dispersion shrinkage stabilizes variance estimates at small sample sizes (n = 3 vs 4)
- Field-standard method for small-sample RNA-seq differential expression (Love, Huber & Anders, 2014)

**Pipeline per cell type:**
1. Filter low-count genes (≥ 10 counts in ≥ 2 samples)
2. Median-of-ratios normalization
3. Gene-wise dispersion estimation + empirical-Bayes shrinkage
4. Negative binomial GLM fit with design `~ disease`
5. Wald test on the disease coefficient
6. Benjamini–Hochberg FDR correction
7. DEG classification: `padj < 0.05` AND `|log2 fold change| > 1`



## 6. Repository Structure

The repository is organized into numbered folders matching the order of execution:

rnaseq-immune-profiling-als/
├── README.md                           # This file
├── LICENSE                             # MIT license
├── .gitignore
├── 001_Data_analysis.Rproj             # RStudio project file
├── ALS_Final_Report.pdf                # Written project report
├── ALS_Final_Presentation.pptx         # Slide deck
│
├── 00_DataFiles/                       # Raw inputs from GEO
│   ├── GSE60424_metadata.csv
│   └── GSE60424_expression_matrix_from_supp.csv
│
├── 01_DataCleaning/                    # STEP 1: clean raw data
│   ├── 01_DataCleaning.Rmd
│   ├── 01_DataCleaning.html
│   ├── GSE60424_counts_clean.csv       # 50,045 genes × 134 samples
│   ├── GSE60424_metadata_clean.csv
│   └── GSE60424_sample_map.csv         # links GSM IDs ↔ lib IDs
│
├── 02_DataExploration/                 # STEP 2: QC and exploration
│   ├── 02_DataExploration.Rmd
│   ├── PCA_plot_ALS.png
│   └── PCA_plot_T1D.png
│
└── 03_DataAnalysis/                    # STEP 3: analysis + visualization
       ├── 01_StatisticalAnalysis.Rmd      # DESeq2 pipeline
       ├── 01_StatisticalAnalysis.html
       ├── 01_StatisticalAnalysis/         # auto-created outputs
       ├── 02_Visualization.Rmd            # Plotting pipeline
       ├── 02_Visualization.html
       └── 02_Visualization/               # auto-created outputs



## 7. How to Reproduce

**Prerequisites:** R ≥ 4.5, RStudio, and the following packages:

```r
# CRAN packages
install.packages(c("data.table", "dplyr", "ggplot2", "ggrepel",
                   "pheatmap", "RColorBrewer", "tibble"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "org.Hs.eg.db", "AnnotationDbi"))
```

**Before running:** update the hard-coded paths at the top of each `.Rmd` (`clean_dir`, `output_dir`, etc.) to match your local machine.

**Execution order:**

1. Open `001_Data_analysis.Rproj` in RStudio
2. Knit `01_DataCleaning/01_DataCleaning.Rmd` → produces cleaned CSVs
3. Knit `03_DataAnalysis/01_StatisticalAnalysis.Rmd` → produces DEG tables and VST counts
4. Knit `03_DataAnalysis/02_Visualization.Rmd` → produces all figures

Each script consumes the outputs of the previous one. Output subfolders are created automatically.


## 8. Key Design Decisions

**Why ALS vs Healthy Control only?** GSE60424 contains several disease cohorts (MS, T1D, Sepsis, ALS). I restricted to ALS vs HC to keep the comparison focused. The same pipeline generalizes to any of the other comparisons by changing two parameters.

**Why stratify by cell type?** Global PCA showed that cell-type identity dominates transcriptional variation (PC1 = 56% of global variance). A pooled analysis across cell types would have drowned any disease signal in cell-type differences. Running DESeq2 separately within each cell type isolates the disease effect.

**Why these thresholds (padj < 0.05, |LFC| > 1)?** padj < 0.05 is the standard FDR cutoff; |LFC| > 1 requires a 2-fold change, which is a biologically meaningful effect size. Both are DESeq2 documentation defaults.

**Why focus on Neutrophils and Monocytes for detail plots?** They produced the highest DEG counts (38 and 34 respectively) — roughly 2–3× more than the adaptive immune compartments — consistent with published ALS immunology implicating innate immune dysregulation.

**Why a consensus heatmap?** Genes significant in both Neutrophils AND Monocytes (17 total) are more likely to reflect real ALS biology rather than cell-specific noise or single-donor artifacts. This cross-cell-type intersection is the highest-confidence output of the analysis.

**Why DESeq2?**

Raw RNA-seq counts are overdispersed integers — variance increases with the mean, violating the constant-variance assumption of linear models. DESeq2 uses a **negative binomial distribution** with a log link function, which correctly handles overdispersed count data. Gene-wise dispersion estimates are shrunk toward a fitted mean–dispersion trend, which stabilizes variance estimates at small sample sizes (n = 3 vs 4 per stratum in this study).

**Why a Wald test with Benjamini–Hochberg correction?**
The Wald test directly evaluates whether the disease coefficient in the negative binomial GLM differs significantly from zero — a straightforward test for a two-level factor. Benjamini–Hochberg controls the false discovery rate across the ~8,000 genes tested per cell type, which is appropriate for genome-scale data where many real discoveries are expected and Bonferroni would be overly conservative.

**Why a consensus (set intersection) across cell types?**
DESeq2 identifies significant genes within one cell type, but does not tell you which genes reflect a shared biological response versus cell-specific artifacts. By intersecting the DEG sets from Neutrophils and Monocytes, the consensus approach filters for genes whose ALS signal replicates in a second independent cell population — a stronger standard than single-cell-type significance alone.
DESeq2 is the field-standard method for small-sample RNA-seq differential expression and is more statistically appropriate than parametric alternatives.


## 9. Key Findings

| Cell Type   | N Samples | Up | Down | Total DEGs |
|-|--|-|||
| Neutrophils | 7         | 10 | 28   | **38**     |
| Monocytes   | 7         | 18 | 16   | **34**     |
| NK          | 7         | 16 | 12   | 28         |
| Whole Blood | 7         | 9  | 16   | 25         |
| CD4         | 7         | 12 | 2    | 14         |
| B-cells     | 7         | 11 | 3    | 14         |
| CD8         | 7         | 12 | 1    | 13         |

*Thresholds: padj < 0.05 and |log2FC| > 1.*

**Biological takeaways:**
- **Innate > adaptive:** Neutrophils and monocytes carry ~2–3× more DEGs than T and B cells.
- **17 consensus DEGs** are shared between Neutrophils and Monocytes.
- **Interferon-response genes** (IFI44, IFI44L, IFIT1, IFI6, OAS2, OAS3, OASL, CMPK2, XAF1, ISG15) are consistently **down-regulated in ALS** in both cell types — suggesting a blunted type-I interferon response in ALS innate immunity. This immune imbalance may have diagnostic potential as a biomarker signature.



## 10. Data Source

**Dataset:** GSE60424 — Gene Expression Omnibus  
**Original publication:** Linsley PS, Speake C, Whalen E, Chaussabel D (2014). *Copy number loss of the interferon gene cluster in melanomas is linked to reduced T cell infiltrate and poor patient prognosis.* PLoS One 9(10): e109760.  
**Platform:** Illumina HiSeq 2000  
**URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60424



## 11. Key Versions
- R 4.5.0
- DESeq2 1.48.1
- ggplot2 4.0.2
- pheatmap 1.0.13
- org.Hs.eg.db 3.21.0



## 12. References 

1. **Love MI, Huber W, Anders S.** (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15(12): 550.
2. **Linsley PS, Speake C, Whalen E, Chaussabel D.** (2014). Copy number loss of the interferon gene cluster in melanomas. *PLoS One* 9(10): e109760. [GSE60424 deposit]
3. **Beers DR, Appel SH.** (2019). Immune dysregulation in amyotrophic lateral sclerosis. *Lancet Neurology* 18(2): 211–220.
4. **Murdock BJ, et al.** (2017). Correlation of peripheral immunity with rapid ALS progression. *JAMA Neurology* 74(12): 1446–1454.

## 13. Code of Ethics & Data Use Compliance

In this project, only publicly archived human RNA-seq data was utilized, and the work is therefore governed by codes of ethics covering
(a) the data source
(b) human genomic data sharing
(c) responsible conduct in human-genetics research.

### 1. NCBI GEO Data Use Policy

**Source:** [https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html)

>**GEO Disclaimer, paragraph 1, sentence 2:** *"NCBI places no restrictions on the use or distribution of the GEO data. However, some submitters may claim patent, copyright, or other intellectual property rights in all or a portion of the data they have submitted."*

GEO: GSE60424 Dataset (open-access archive) - GEO is the open-access archive, which we use to host the dataset in this study. Its policy reads that there are no limitations on the use or Distribution of the GEO data by NCBI, although some submitters may assert patent, copyright or other intellectual property rights in all or part of any submitted text. Along with all applicable federal, Tribal, state and local laws, regulations, statutes, guidance and institutional policies at the time of deposit—human-derived submissions to unrestricted-access GEO must already comply with the original consent.

**For this Repository and the Project:**
- Data downloaded only from the official GEO repository — no scraping or unauthorized redistribution was performed .
- The original data submitters (Linsley et al., 2014) are credited in this README and in the project report.
- Only the obtained analytical outputs (filtered counts, DEG tables, normalized expression for visualization) are committed to this repository — not the raw GEO deposit.
- Re-identification of the donors was not at all attempted respecting IPR and Donor Privacy.

### 2. NIH Genomic Data Sharing (GDS) Policy

**Source:** [https://grants.nih.gov/grants/guide/notice-files/not-od-14-124.html](https://grants.nih.gov/grants/guide/notice-files/not-od-14-124.html)

> **GEO FAQ, "Human Subject Guidelines" section:** *"If you plan to submit genomic data from human specimens... it is your responsibility to ensure that the submitted information does not compromise participant privacy, and is in accord with the original consent, in addition to all applicable laws, regulations, and institutional policies."*

This is the framework GEO operates under.The GDS Policy covers all NIH-funded research that generates large-scale human or non-human genomic data and the use of such data for subsequent research; and large-scale data explicitly includes transcriptomic, metagenomic, epigenomic and gene expression.

**For this Repository and the Project:**
- The data are used only for the educational/research purpose stated in this README — characterizing transcriptional differences between ALS and Healthy Control donors.
- Data are not re-distributed; the cleaned files committed to this repository are derived analytical outputs (filtered count matrices, DEG tables, normalized expression for visualization), not the raw deposit.


## 14. License

This project's code is released under the [MIT License](LICENSE). The underlying RNA-seq data belongs to GEO accession GSE60424 (Linsley et al., 2014) and is subject to the data-use terms in [Section 13](#13-code-of-ethics--data-use-compliance).

## 15. Author

**Sushant Agellu** — [GitHub](https://github.com/SushantAgellu)

