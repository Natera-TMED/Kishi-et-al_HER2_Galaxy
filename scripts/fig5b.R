# Fig.5b

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
# df <- read.csv("../cohort/ctDNA_of_HER2_amplification_all_Stage_for_R.csv")
df <- read.csv("/path/to/your/ctDNA_of_HER2_amplification_all_Stage_for_R.csv")

df <- df %>% mutate(DFS_month_3year = DFS_3year / 30.4)

# survival object
surv_obj <- with(df, Surv(DFS_month_3year, DFS_event_3year))
fit <- survfit(surv_obj ~ `ctDNA_Result_4W`, data = df)

# log-rank & Cox
lg  <- survdiff(surv_obj ~ `ctDNA_Result_4W`, data=df)
p_logrank <- pchisq(lg$chisq, df = length(lg$n) - 1, lower.tail = FALSE)

cx <- coxph(surv_obj ~ `ctDNA_Result_4W`, data=df)
hr <- exp(coef(cx))[1]
ci <- exp(confint(cx))[1, ]
p_cox <- summary(cx)$coefficients[1, "Pr(>|z|)"]

# 3year-DFS
sf36 <- summary(fit, times = 36)

km36_raw <- data.frame(
  strata = sf36$strata,
  surv   = sf36$surv,
  lower  = sf36$lower,
  upper  = sf36$upper,
  n_risk = sf36$n.risk
)

all_strata <- names(fit$strata)
km36 <- merge(
  data.frame(strata = all_strata),
  km36_raw,
  by = "strata",
  all.x = TRUE,
  sort = FALSE
)

# label
group_val <- sub(".*=", "", km36$strata)  # "Negative" / "Positive"
labs_map  <- c("Negative" = "4W ctDNA(-)", "Positive" = "4W ctDNA(+)")
km36$label <- labs_map[group_val]

# # Round DFS to one decimal place
fmt_pct <- function(x) ifelse(is.na(x), "NA", sprintf("%.1f", x * 100))
km36$surv_pct  <- fmt_pct(km36$surv)
km36$lower_pct <- fmt_pct(km36$lower)
km36$upper_pct <- fmt_pct(km36$upper)

# 3year-DFS
lines_3yr <- sprintf("%s: %s%% (95%%CI %s–%s)",
                     km36$label, km36$surv_pct, km36$lower_pct, km36$upper_pct)

# pvalue, HR
label_text <- paste0(
  sprintf("log-rank p=%.2e\nCox HR=%.2f (95%%CI %.2f–%.2f)",
          p_logrank, hr, ci[1], ci[2]),
  "\n3-year DFS\n  ", paste(lines_3yr, collapse = "\n  ")
)

# Kaplan-Meier
g <- ggsurvplot(
  fit,
  data = df,
  palette = "lancet",
  risk.table = TRUE,
  title = "pStage I-III",
  xlab = "Time (Months)",
  ylab = "Disease-free survival",
  legend.title = "",
  legend.labs = c("4W ctDNA(-)", "4W ctDNA(+)"),
  linetype = c("dotdash", "solid"),
  xlim = c(0, 36),
  break.time.by = 6,
  risk.table.height = 0.30
)

g$plot <- g$plot +
  theme(
    legend.position = c(0.83, 0.3),
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
           x = -1.5, y = 0.16,   # 目盛りの数字で位置を指定
           hjust = 0,
           size = annot_size,
           label = label_text
  )

print(g)


library(gridExtra)
g_comb <- arrangeGrob(g$plot, g$table, ncol = 1, heights = c(0.70, 0.30))

ggplot2::ggsave(
  "Fig5b.pdf",
  plot   = g_comb,
  device = grDevices::cairo_pdf,   
  width  = 8, height = 6, units = "in",
  bg     = "white"
)