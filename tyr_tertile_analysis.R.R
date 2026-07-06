


library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)

library(dplyr)
library(ggplot2)
library(tidyr)

library(EnhancedVolcano)
library(pheatmap)



dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)



query <- GDCquery(
  project = "TCGA-SKCM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)

melanoma <- GDCprepare(query)



saveRDS(
  melanoma,
  "data/raw/tcga_skcm.rds"
)



melanoma <- readRDS(
  "data/raw/tcga_skcm.rds"
)



class(melanoma)

dim(melanoma)

melanoma

assayNames(melanoma)


counts <- assay(melanoma)

dim(counts)

head(counts)



gene_info <- as.data.frame(
  rowData(melanoma)
)

head(gene_info)

write.csv(
  gene_info,
  "data/processed/gene_information.csv",
  row.names = FALSE
)



sample_info <- as.data.frame(
  colData(melanoma)
)

head(sample_info)

write.csv(
  sample_info,
  "data/processed/sample_information.csv",
  row.names = FALSE
)


gene_info %>%
  filter(gene_name == "TYR")



tyr_counts <- counts["ENSG00000077498.9", ]

tyr_df <- data.frame(
  Sample = colnames(counts),
  TYR = as.numeric(tyr_counts)
)

head(tyr_df)

summary(tyr_df$TYR)

write.csv(
  tyr_df,
  "data/processed/TYR_expression.csv",
  row.names = FALSE
)


ggplot(
  tyr_df,
  aes(TYR)
) +
  geom_histogram(
    bins = 30,
    fill = "steelblue"
  )

ggsave(
  "figures/TYR_histogram.png",
  width = 6,
  height = 4
)




ggplot(
  tyr_df,
  aes(TYR)
) +
  geom_density(
    fill = "orange"
  )

ggsave(
  "figures/TYR_density.png",
  width = 6,
  height = 4
)



ggplot(
  tyr_df,
  aes(y = TYR)
) +
  geom_boxplot(
    fill = "pink"
  )

ggsave(
  "figures/TYR_boxplot.png",
  width = 5,
  height = 5
)



cutpoints <- quantile(
  tyr_df$TYR,
  probs = c(0, 0.33, 0.67, 1)
)

cutpoints

tyr_df$Group <- cut(
  tyr_df$TYR,
  breaks = cutpoints,
  labels = c(
    "Low",
    "Medium",
    "High"
  ),
  include.lowest = TRUE
)

table(tyr_df$Group)



write.csv(
  tyr_df,
  "data/processed/TYR_groups.csv",
  row.names = FALSE
)



ggplot(
  tyr_df,
  aes(Group,
      fill = Group)
) +
  geom_bar()

ggsave(
  "figures/TYR_groups.png",
  width = 5,
  height = 4
)



selected <- tyr_df %>%
  filter(Group != "Medium")



counts_subset <- counts[
  ,
  selected$Sample
]


metadata <- data.frame(
  Group = selected$Group
)

rownames(metadata) <- selected$Sample



dds <- DESeqDataSetFromMatrix(
  countData = counts_subset,
  colData = metadata,
  design = ~ Group
)

dds <- DESeq(dds)



res <- results(dds)

res_df <- as.data.frame(res)



write.csv(
  res_df,
  "results/DESeq2_results.csv"
)



sig <- res_df %>%
  filter(
    padj < 0.05,
    abs(log2FoldChange) > 1
  )

write.csv(
  sig,
  "results/significant_genes.csv"
)



EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj"
)



png(
  "figures/volcano_plot.png",
  width = 1800,
  height = 1400,
  res = 300
)

EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj"
)

dev.off()



top_genes <- rownames(sig)[1:30]

norm_counts <- counts(
  dds,
  normalized = TRUE
)

heatmap_data <- norm_counts[
  top_genes,
]

png(
  "figures/heatmap.png",
  width = 1800,
  height = 1500,
  res = 300
)

pheatmap(
  heatmap_data,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Top 30 Differentially Expressed Genes"
)

dev.off()
cat("Analysis completed successfully!\n")
