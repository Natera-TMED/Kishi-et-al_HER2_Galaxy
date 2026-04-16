# Fig.5f

library(ggplot2)

# load
# df <- read.csv("../cohort/copy_number_of_HER2_amp_with_MTM_for_R.csv")
df <- read.csv("/path/to/your/copy_number_of_HER2_amp_with_MTM_for_R.csv")

#  correlation coefficient, pvalue
pearson_test  <- cor.test(df$copy_number, df$MTM_initial, method = "pearson")
spearman_test <- cor.test(df$copy_number, df$MTM_initial, method = "spearman")

cor_pearson  <- pearson_test$estimate
p_pearson    <- pearson_test$p.value

cor_spearman <- spearman_test$estimate
p_spearman   <- spearman_test$p.value

# scatter plot
g=ggplot(df, aes(x = copy_number, y = MTM_initial)) +
  geom_point(color = "steelblue", size = 2.5) +
  annotate("text", 
           x = Inf, y = Inf, 
           label = paste0(
             "Pearson r = ", round(cor_pearson, 4),
             " (p = ", round(p_pearson, 4), ")\n",
             "Spearman ρ = ", round(cor_spearman, 4),
             " (p = ", round(p_spearman, 4), ")"
           ),
           hjust = 1.05, vjust = 1.2, size = 4, color = "black") +
  labs(
    x = "copy number",
    y = "preoperative MTM"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )+
  
  scale_x_continuous(
    limits = c(0, NA),
    breaks = pretty(c(0, max(df$copy_number, na.rm = TRUE)))
  )


ggplot2::ggsave(
  "Fig5f.pdf",
  plot   = g,
  device = grDevices::cairo_pdf,   
  width  = 8, height = 6, units = "in",
  bg     = "white"
)