
 Parkinson’s Microarray Analysis (limma)
 Dataset: GSE99039


# Load libraries
library(GEOquery)
library(limma)
library(tidyverse)

# Setting working directory
setwd("~/parkinsons_project")

# Create results folder
dir.create("results", showWarnings = FALSE)


# Step 1 — Load the local GEO file


gset <- getGEO(
  filename = "GSE99039_series_matrix.txt.gz",
  GSEMatrix = TRUE
)[[1]]


# Step 2 — Extract phenotype data


pdata <- pData(gset)

# Keep only CONTROL and IPD samples
keep <- pdata$`disease label:ch1` %in% c("CONTROL", "IPD")
pdata_clean <- pdata[keep, ]


# Step 3 — Extract expression matrix


expr <- exprs(gset)
expr_clean <- expr[, rownames(pdata_clean)]


# Step 4 —(if needed)Reduce to small subset


control_samples <- rownames(
  pdata_clean[pdata_clean$`disease label:ch1` == "CONTROL", ]
)[1:20]

ipd_samples <- rownames(
  pdata_clean[pdata_clean$`disease label:ch1` == "IPD", ]
)[1:20]

final_samples <- c(control_samples, ipd_samples)

expr_final  <- expr_clean[, final_samples]
pdata_final <- pdata_clean[final_samples, ]


# Step 5 — Create design matrix


group <- factor(pdata_final$`disease label:ch1`)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)


# Step 6 — limma differential expression


fit <- lmFit(expr_final, design)

contrast.matrix <- makeContrasts(IPD - CONTROL, levels = design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


# Step 7 — Get results


res <- topTable(fit2, number = Inf)

head(res)


# Step 8 — Save results


write.csv(res, "results/limma_results_GSE99039.csv")



# Step 6 — Prepare ranked gene list for GSEA


# Load annotation package for GPL570
BiocManager::install("hgu133plus2.db", ask = FALSE, update = FALSE)
library(hgu133plus2.db)
library(dplyr)

# Add probe IDs as a column
res$probe_id <- rownames(res)

# Map probe IDs to gene symbols
res$symbol <- mapIds(
  hgu133plus2.db,
  keys = res$probe_id,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

# Remove probes without gene symbol
res <- res[!is.na(res$symbol), ]

# If multiple probes map to same gene, keep the one with highest absolute logFC
res <- res %>%
  arrange(desc(abs(logFC))) %>%
  distinct(symbol, .keep_all = TRUE)

# Create ranked gene list (named numeric vector)
gene_list <- res$logFC
names(gene_list) <- res$symbol

# Sort for GSEA
gene_list <- sort(gene_list, decreasing = TRUE)

# Save for future use
write.csv(gene_list, "results/gene_list_for_gsea.csv")

head(gene_list)

# Step 7 — GSEA Pathway Enrichment


library(fgsea)
library(msigdbr)
library(dplyr)

# Get Hallmark gene sets for humans
msig <- msigdbr(species = "Homo sapiens", category = "H")

# Convert to list format required by fgsea
pathways <- msig %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Run GSEA using ranked gene list
set.seed(1)
gene_list <- gene_list + rnorm(length(gene_list), 0, 1e-6)

fgseaRes <- fgsea(
  pathways = pathways,
  stats    = gene_list
)

# Sort by adjusted p-value (FDR)
fgseaRes <- fgseaRes %>% arrange(padj)

# View top pathways
head(fgseaRes)
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, paste, collapse = ",")
# Save results
write.csv(fgseaRes, "results/gsea_results.csv", row.names = FALSE)


# Plot enrichment of top 5 pathway
top_pathways <- fgseaRes %>%
  dplyr::arrange(padj) %>%
  dplyr::slice(1:5)

write.csv(top_pathways, "results/top_5_pathways.csv", row.names = FALSE)

# Plot enrichment of top pathway
plotEnrichment(pathways[[fgseaRes$pathway[1]]], gene_list)

# Extract core genes from the top 5 significant pathway
core_genes <- unique(unlist(strsplit(top_pathways$leadingEdge, ",")))

write.csv(core_genes, "results/core_genes_top_pathways.csv", row.names = FALSE)
