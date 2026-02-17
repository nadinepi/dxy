# Final Analysis Script
# Script to process Dxy data from Pixy (50kb windows)
# Generates two plots for the report

library(tidyverse)
library(scales)

# Input files
input_file <- "results/pixy/ProjTaxa_dxy.txt"
gff_file   <- "data/ProjTaxa.gff"
# Output location
output_dir <- "results/figures"

# Create output folder if missing
if(!dir.exists(output_dir)) dir.create(output_dir)

pop_colors <- c(
  "Pop8N vs PopK" = "#C77CFF",
  "Pop8N vs PopLesina" = "#56B4E9", 
  "PopK vs PopLesina" = "#009E73"
)

# Loading the data from Pixy
dxy_data <- read_delim(input_file, delim = "\t", col_types = cols()) %>%
  filter(!is.na(avg_dxy)) %>%
  mutate(
    comparison = paste(pop1, "vs", pop2),
    # Calculate midpoint for plotting
    window_mid = (window_pos_1 + window_pos_2) / 2
  )

# Calculating Z-scores to identify outliers
dxy_data <- dxy_data %>%
  group_by(comparison) %>%
  mutate(
    z_score = (avg_dxy - mean(avg_dxy)) / sd(avg_dxy),
    is_outlier = z_score > 3 
  )

# Statistical test for Chromosome 5 vs Z Chromosome
stats_file <- file.path(output_dir, "statistics_summary.txt")

stats_data <- dxy_data %>%
  mutate(chrom_type = ifelse(chromosome == "chrZ", "Z Chromosome", "Autosome"))

means <- stats_data %>%
  group_by(chrom_type) %>%
  summarise(
    mean_dxy = mean(avg_dxy, na.rm = TRUE),
    n = n()
  )

# Wilcoxon rank sum test
wilcox_res <- wilcox.test(avg_dxy ~ chrom_type, data = stats_data)

# Save results to file
sink(stats_file)
cat("Statistical Analysis: Chromosome 5 vs Z Chromosome\n\n")
cat("1. Mean DXY values by chromosome type:\n")
print(means)
cat("\n")

cat("2. Wilcoxon Rank Sum Test Result:\n")
print(wilcox_res)
cat("\n")

# Interpretation
if(wilcox_res$p.value < 2.2e-16) {
    cat("Highly Significant (p < 2.2e-16)\n")
} else {
    cat("Result: p-value =", wilcox_res$p.value, "\n")
}
sink()

cat("Statistical summary saved to:", stats_file, "\n")

# Figure 1: Genome-wide Absolute Divergence
p1 <- ggplot(dxy_data, aes(x = window_mid/1e6, y = avg_dxy)) +
  geom_point(aes(color = comparison), alpha = 0.6, size = 0.5) +
  
  # Highlight the outliers
  geom_point(data = filter(dxy_data, is_outlier == TRUE),
             color = "firebrick", size = 0.8, alpha = 0.8) +
  
  # Add smoothed trend line
  geom_smooth(aes(color = comparison), method = "loess", se = FALSE, span = 0.1, color = "black", linewidth = 0.5) +
  scale_color_manual(values = pop_colors) +
  facet_grid(comparison ~ chromosome, scales = "free_x", space = "free_x") +
  labs(
    title = "Genome-wide Absolute Divergence (Red Points = Z-score > 3)",
    x = "Position (Mb)",
    y = expression(D[XY])
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none", 
    strip.text = element_text(face="bold", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 16)
  )

ggsave(file.path(output_dir, "Figure1_Genome_Dxy.pdf"), p1, width = 12, height = 8)


# Figure 2: Zoom on Chromosome Z inversion
# Filter data for the specific region
zoom_data <- dxy_data %>%
  filter(chromosome == "chrZ", window_pos_1 >= 11500000, window_pos_2 <= 20000000)

p2 <- ggplot(zoom_data, aes(x = window_mid/1e6, y = avg_dxy, color = comparison)) +
  geom_point(alpha = 0.6, size = 1.5) +
  # trend line
  geom_smooth(method = "loess", se = FALSE, span = 0.2, linewidth = 1.5) +
  scale_color_manual(values = pop_colors) +
  labs(
    title = "Divergence on Chromosome Z (11.5-20Mb)",
    x = "Position (Mb)",
    y = expression(D[XY]),
    color = ""
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "bottom", 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.key.size = unit(1.0, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 3, size = 3)))

ggsave(file.path(output_dir, "Figure2_ChrZ_Zoom_50kb.pdf"), p2, width = 10, height = 7)

cat("Done! Files are in:", output_dir, "\n")
