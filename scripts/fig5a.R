# Fig.5a

library(ggplot2)
library(ggpubr)
library(dplyr)

# load
# df <- read.csv("../cohort/DFS_of_HER2_amp_with_copy_number.csv")
df <- read.csv("/path/to/your/DFS_of_HER2_amp_with_copy_number.csv")

# RFS_Event
df$RFS_Event <- factor(df$RFS_Event, levels = c(0,1),
                       labels = c("No recurrence (N = 38)", "Recurrence (N = 10)"))

# pvalue
pval <- wilcox.test(copy_number ~ RFS_Event, data = df)$p.value

# median, IQR
stats <- df %>%
  group_by(RFS_Event) %>%
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
  select(RFS_Event, median, Q1, Q3, IQR_range)

print(stats)

# boxplot
g=ggplot(df, aes(x = RFS_Event, y = copy_number, fill = RFS_Event)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 60)) +
  scale_y_continuous(breaks = seq(0, 60, 10)) +
  labs(title = "pStage I-III",
       x = "",
       y = "Copy number") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12, family = "Arial"),
        axis.title.y = element_text(size = 12, family = "Arial"),
        axis.text.x  = element_text(size = 12, family = "Arial"),
        axis.text.y  = element_text(size = 12, family = "Arial"),
        plot.title = element_text(
          size = 14,
          family = "Arial",
          hjust = 0.5,
          vjust = 0.5
        )
  ) +
  
  annotate("text",
           x = 1.5,
           y = max(df$copy_number) * 1.0,
           label = sprintf("Mann-Whitney U : p = %.4f", pval),
           size = 3.4,
           family = "Arial")

print(g)

ggplot2::ggsave(
  "Fig5a.pdf",
  plot   = g,
  device = grDevices::cairo_pdf,   
  width  = 8, height = 6, units = "in",
  bg     = "white"
)