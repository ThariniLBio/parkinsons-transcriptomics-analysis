# Step - 1 Install once if not present

install.packages("ggplot2")

library(ggplot2)

# Step - 2 Add columns needed for volcano plot

res$negLogP <- -log10(res$adj.P.Val)

# Step - 3 Define significance

res$threshold <- "Not Significant"
res$threshold[res$adj.P.Val < 0.05 & res$logFC > 1]  <- "Upregulated"
res$threshold[res$adj.P.Val < 0.05 & res$logFC < -1] <- "Downregulated"

# Step - 4 Volcano plot

ggplot(res, aes(x = logFC, y = negLogP, color = threshold)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot - Parkinson’s (IPD vs Control)",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )

# Step - 5 To save

ggsave("volcano_plot.png", width = 8, height = 6)
