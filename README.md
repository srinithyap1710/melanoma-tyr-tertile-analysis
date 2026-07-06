# melanoma-tyr-tertile-analysis
Differential gene expression analysis of TCGA Skin Cutaneous Melanoma (SKCM) on Making TYR expression tertiles using R, TCGAbiolinks, and DESeq2.

Biological Question-

Do melanoma tumors with high TYR expression exhibit distinct transcriptomic profiles compared to tumors with low TYR expression?**

This project explores how differences in TYR expression are associated with global gene expression changes across melanoma samples.

---

Project Overview-

This project investigates whether melanoma tumors with high TYR expression exhibit different global gene expression patterns than tumors with low TYR expression. of the **TYR (Tyrosinase)** gene using RNA-seq data from **The Cancer Genome Atlas (TCGA)**.

TYR is a melanocyte-specific gene involved in melanin synthesis and is widely used as a marker of melanoma differentiation. Samples were divided into **Low**, **Medium**, and **High** TYR expression groups using tertiles, followed by differential expression analysis between the High and Low groups.

---

 Objectives of the project-

- Download RNA-seq data from TCGA-SKCM
- Explore the SummarizedExperiment object
- Extract TYR gene expression
- Divide samples into TYR expression tertiles
- Compare High vs Low TYR expression groups
- Perform differential expression analysis using DESeq2
- Visualize results using volcano plots and heatmaps

---

Dataset which is used-

Project: TCGA-SKCM (Skin Cutaneous Melanoma)

**Data Type:**
- Transcriptome Profiling
- Gene Expression Quantification
- STAR - Counts workflow

Data were downloaded using the **TCGAbiolinks** package.

---

## Workflow

1. Download TCGA melanoma RNA-seq data
2. Prepare data using `GDCprepare()`
3. Extract raw count matrix
4. Retrieve TYR gene expression
5. Create Low, Medium, and High TYR expression groups using tertiles
6. Remove the Medium group
7. Perform differential expression analysis with DESeq2
8. Generate volcano plot
9. Generate heatmap of top differentially expressed genes

---

## Tools and Packages

- R
- Bioconductor
- TCGAbiolinks
- DESeq2
- SummarizedExperiment
- ggplot2
- dplyr
- EnhancedVolcano
- pheatmap

---

## Repository Structure

```
melanoma-tyr-tertile-analysis/

├── data/
│   ├── raw/
│   └── processed/
│
├── figures/
│
├── results/
│
├── scripts/
│   └── 01_download_data.R
│
├── README.md

```

---

## Output

The project generates:

- TYR expression histogram
- TYR density plot
- TYR boxplot
- TYR tertile distribution
- Volcano plot
- Heatmap of top differentially expressed genes
- Complete DESeq2 differential expression results
- Significant differentially expressed genes

---

## Learning Outcomes

Through this project I learned how to:

- Access TCGA data using TCGAbiolinks
- Work with SummarizedExperiment objects
- Process RNA-seq count data
- Perform differential expression analysis using DESeq2
- Create publication-style visualizations
- Organize a reproducible bioinformatics workflow using GitHub

---

## Author

**Srinithya Reddy**

B.Tech Biotechnology

Mahindra University
