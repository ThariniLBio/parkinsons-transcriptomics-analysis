# To Find Top 20 Target Genes from Top 5 Significant Pathway

# Step 1 Keep only pathway genes

target_table <- res[res$symbol %in% core_genes, ]

# Sort by significance

target_table <- target_table[order(target_table$adj.P.Val), ]

# Take top 20

top20_targets <- target_table[1:20, ]

# View important columns

top20_targets[, c("symbol","logFC","adj.P.Val")]


# Save

write.csv(top20_targets, "results/top20_targets_genes.csv", row.names = FALSE)
