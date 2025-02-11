
BiocManager::install('Maaslin2')
library(Maaslin2)


ps <- ps.glom.norm.top
samples.to.keep <- grep('FR|human|canis', sample_names(ps), value = T)

ps <- prune_samples(samples.to.keep, ps)

# Prepare data for MaAsLin2
otutab   <- as.data.frame(t(as.data.frame(otu_table(ps))))
otutab   <- data.frame(ID = rownames(otutab), otutab)

metadata <- as.data.frame(sample_data(ps))
metadata <- data.frame(ID = rownames(metadata), metadata)

# Run MaAsLin2
maaslin_res <- Maaslin2(
  input_data = otutab,
  input_metadata = metadata,
  min_abundance = 0.01,
  min_prevalence = 0.0,
  min_variance = 0.0,
  normalization = 'NONE',
  transform = 'NONE',
  output = "maaslin_results",
  fixed_effects = c("source"),
  reference = c("source,FR")
)


# View results
maaslin_res$results

## coef barplot
ggplot(maaslin_res$results[maaslin_res$results$qval < 0.05, ],
       aes(x = feature, y = coef, fill = value)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Significant Associations", x = "Source", y = "Effect Size")
  # + facet_grid(rows=vars(feature, value), scales='free')

## Heatmap Significant associations (−log(qval)*sign(coeff))
# Prepare the data
results <- maaslin_res$results

# Transform results table to create the heatmap values
results$heatmap_value <- -log10(results$qval) * sign(results$coef)

setDT(results)
results[,is.significant := fifelse(qval < 0.05, T, F)]
results[,sign           := fifelse(coef < 0, '-', '+')]

# Number of significant interactions
res_sum <- results[,.N,by=.(value, is.significant, sign)]

ggplot(res_sum, aes(x=value, y=N, fill=sign)) +
  geom_col(position = 'dodge') +
  theme_minimal() +
  scale_fill_manual(
    values = c("blue", "red")) +
  facet_grid(cols=vars(is.significant))


# Filter for significant results
significant_results <- results[qval < 0.05, ]

# Reshape data for heatmap
heatmap_data <- dcast(significant_results, feature ~ value, value.var = "heatmap_value", fill = 0)

# Convert to matrix
rownames(heatmap_data) <- heatmap_data$feature
heatmap_data$feature <- NULL
heatmap_matrix <- as.matrix(heatmap_data)

# Melt for ggplot
heatmap_df <- melt(heatmap_matrix)
colnames(heatmap_df) <- c("Feature", "Group", "Value")

# Create a column for "+" and "-" signs
heatmap_df$Sign <- ifelse(heatmap_df$Value > 0, "+", ifelse(heatmap_df$Value < 0, "-", ""))

# Create the heatmap
ggplot(heatmap_df, aes(x = Group, y = Feature, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, name = "-log(qval) * sign(coef)"
  ) +
  geom_text(aes(label = Sign), color = "black", size = 6) + # Add "+" and "-" signs
  theme_minimal() +
  theme(
    axis.text.x = element_text(face='italic', angle = 45, hjust = 1),
    axis.text.y = element_text(face='italic'),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(
    title = "Significant Associations Heatmap",
    fill = "Effect Size"
  )










































