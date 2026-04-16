# Fig.5e

library(ggplot2)
library(ggpubr)
library(dplyr)

# load
# df <- read.csv("../cohort/copy_number_of_HER2_amp_with_ctDNA_for_R.csv")
df <- read.csv("/path/to/your/copy_number_of_HER2_amp_with_ctDNA_for_R.csv")

# RFS_Event
df$ctDNA_Result_4W <- factor(df$ctDNA_Result_4W, levels = c("Negative","Positive"))

# pvalue
pval <- wilcox.test(copy_number ~ ctDNA_Result_4W, data = df)$p.value

# median, IQR
stats <- df %>%
  group_by(ctDNA_Result_4W) %>%
  summarise(
    median_raw = median(copy_number, na.rm = TRUE),
    Q1_raw = quantile(copy_number, 0.25, na.rm = TRUE),
    Q3_raw = quantile(copy_number, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    median = sprintf("%.2f", median_raw),
    Q1     = sprintf("%.2f", Q1_raw),
    Q3     = sprintf("%.2f", Q3_raw),
    IQR_range = paste0(Q1, "–", Q3)
  ) %>%
  select(ctDNA_Result_4W, median, Q1, Q3, IQR_range)

print(stats)

# n
n_counts <- df %>%
  group_by(ctDNA_Result_4W) %>%
  summarise(n = n()) %>%
  arrange(factor(ctDNA_Result_4W, levels = c("Negative","Positive")))

# label
new_labels <- paste0(n_counts$ctDNA_Result_4W, " (N = ", n_counts$n, ")")

df$ctDNA_Result_4W <- factor(df$ctDNA_Result_4W,
                             levels = c("Negative","Positive"),
                             labels = new_labels)

# boxplot
g=ggplot(df, aes(x = ctDNA_Result_4W, y = copy_number, fill = ctDNA_Result_4W)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 60)) +
  scale_y_continuous(breaks = seq(0, 60, 10)) +
  labs(title = "",
       x = "4-week ctDNA status",
       y = "Copy number") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12, family = "Arial"),
        axis.title.y = element_text(size = 12, family = "Arial"),
        axis.text.x  = element_text(size = 12, family = "Arial"),
        axis.text.y  = element_text(size = 12, family = "Arial")
  ) +
  
  annotate("text",
           x = 1.5,
           y = max(df$copy_number) * 1.0,
           label = sprintf("Mann-Whitney U : p = %.4f", pval),
           size = 3.4,
           family = "Arial")




ggplot2::ggsave(
  "Fig5e.pdf",
  plot   = g,
  device = grDevices::cairo_pdf,   
  width  = 8, height = 6, units = "in",
  bg     = "white"
)