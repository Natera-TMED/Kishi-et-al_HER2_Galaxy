# Fig.2b

library(survival)
library(survminer)
library(dplyr)

# setting
axis_title_size  <- 12
axis_text_size   <- 12
legend_text_size <- 12
title_size       <- 12
annot_size       <- 3.4
table_axis_title <- 12
table_axis_text  <- 12
font_family      <- "Arial"

# load
# df <- read.csv("../cohort/Prognostic_differences_associated_with_HER2amp_for_R.csv")
df <- read.csv("/path/to/your/Prognostic_differences_associated_with_HER2amp_for_R.csv")
df <- df %>% mutate(DFS_month_3year = DFS_3year / 30.4)

# pStage I
df_stage <- df %>%
  filter(pStage == "I") %>%
  droplevels()

# survival object
surv_obj <- with(df_stage, Surv(DFS_month_3year, DFS_event_3year))
fit <- survfit(surv_obj ~ `HER2_amp`, data = df_stage)

# log-rank & Cox
lg  <- survdiff(surv_obj ~ `HER2_amp`, data=df_stage)
p_logrank <- pchisq(lg$chisq, df = length(lg$n) - 1, lower.tail = FALSE)

cx <- coxph(surv_obj ~ `HER2_amp`, data=df_stage)
hr <- exp(coef(cx))[1]
ci <- exp(confint(cx))[1, ]
p_cox <- summary(cx)$coefficients[1, "Pr(>|z|)"] # 今回は使わないが、Coxのp値の計算

# 3year-DFS
sf36 <- summary(fit, times = 36)

km36 <- data.frame(
  strata = sf36$strata,
  surv   = sf36$surv,
  lower  = sf36$lower,
  upper  = sf36$upper,
  n_risk = sf36$n.risk
)

# label
labs_map <- c("negative" = "HER2 amp(-)", "positive" = "HER2 amp(+)")
km36$label <- labs_map[sub(".*=", "", km36$strata)]

# Round DFS to one decimal place
fmt_pct <- function(x) sprintf("%.1f", x * 100)
km36$surv_pct  <- fmt_pct(km36$surv)
km36$lower_pct <- fmt_pct(km36$lower)
km36$upper_pct <- fmt_pct(km36$upper)

lines_3yr <- sprintf("%s: %s%% (95%%CI %s–%s)",
                     km36$label, km36$surv_pct, km36$lower_pct, km36$upper_pct)

# pvalue, HR
label_text <- paste0(
  sprintf("log-rank p=%.4f\nCox HR=%.2f (95%%CI %.2f–%.2f)",
          p_logrank, hr, ci[1], ci[2]),
  "\n3-year DFS\n  ", paste(lines_3yr, collapse = "\n  ")
)

# Kaplan-Meier
g <- ggsurvplot(
  fit,
  data = df_stage,
  palette = "lancet",
  risk.table = TRUE,
  title = "pStage I",
  xlab = "Time (Months)",
  ylab = "Disease-free survival",
  legend.title = "",
  legend.labs = c("HER2 amp(-)", "HER2 amp(+)"),
  linetype = c("dotdash", "solid"),
  xlim = c(0, 36),
  break.time.by = 6,
  risk.table.height = 0.30
)

g$plot <- g$plot +
  theme(
    legend.position = c(0.8, 0.5),
    text         = element_text(family = font_family),
    axis.title.x = element_text(size = axis_title_size),
    axis.title.y = element_text(size = axis_title_size),
    axis.text.x  = element_text(size = axis_text_size),
    axis.text.y  = element_text(size = axis_text_size),
    legend.text  = element_text(size = legend_text_size),
    plot.title   = element_text(size = title_size, hjust = 0.5)
  )

# risk table
g$table <- g$table +
  theme(
    text         = element_text(family = font_family),
    axis.title.x = element_text(size = table_axis_title),
    axis.text.x  = element_text(size = table_axis_text)
  )

# pvalue, HR
g$plot <- g$plot +
  annotate("text",
           x = 0, y = 0.18,
           hjust = 0,
           size = annot_size,
           label = label_text
  )

print(g)

library(gridExtra)
g_comb <- arrangeGrob(g$plot, g$table, ncol = 1, heights = c(0.70, 0.30))

ggplot2::ggsave(
  "Fig2b.pdf",
  plot   = g_comb,
  device = grDevices::cairo_pdf,   
  width  = 8, height = 6, units = "in",
  bg     = "white"
)