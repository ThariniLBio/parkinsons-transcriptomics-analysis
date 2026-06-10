
# Step 1 — Restore correct sample IDs as column names


colnames(expr_final) <- pdata_final$geo_accession



# Step 2 — Create matrix of only core genes


mat <- expr_final[
  rownames(expr_final) %in% res$probe_id[res$symbol %in% core_genes],
]

write.csv(mat, "results/heatmap_gene_matrix.csv")

mat_scaled <- t(scale(t(mat)))



# Step 3 — Create proper annotation (Condition per sample)


annotation_col <- data.frame(
  Condition = pdata_final$`disease label:ch1`
)

rownames(annotation_col) <- pdata_final$geo_accession



# Step 4 — Plot heatmap (only core genes)


library(pheatmap)

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean"
)