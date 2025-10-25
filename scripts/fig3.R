rm(list=ls())
library(data.table)
library(survival)
library(broom)
library(forestplot)
library(grid)
library(data.table)
library(dplyr)      
library(tibble)

forest_plot_cox <- function(model = NULL, formula = NULL, data = NULL) {
  # Helper for R‐style significance codes (only ***, **, *)
  sig_code <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.0001)     "****"
    else if (p < 0.001) "***"
    else if (p < 0.01)  "**"
    else if (p < 0.05)  "*"
    else                ""
  }
  
  format_p <- function(p) {
    if (is.na(p)) return("") # Handle missing values
    
    sig <- sig_code(p) # Get significance stars
    
    if (p < 0.0001) {
      val_str <- "<0.0001"
    } else {
      val_str <- sprintf("%.4f", p)
    }
    
    return(paste0(val_str, sig)) # Combine value and stars
  }
  
  fit <- if (!is.null(model)) model else survival::coxph(formula, data = data)
  
  tr     <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  tr$fdr <- p.adjust(tr$p.value, method = "fdr")
  
  df       <- if (!is.null(data)) data else stats::model.frame(fit)
  terms_in <- attr(terms(fit), "term.labels")
  
  lab_list <- list()
  meanv    <- lowerv <- upperv <- numeric()
  
  for (v_raw in terms_in) {
    v_clean <- gsub("`", "", v_raw)
    
    if (v_clean %in% names(df) && is.factor(df[[v_clean]])) {
      lvls <- levels(df[[v_clean]])
      for (i in seq_along(lvls)) {
        lvl <- lvls[i]
        if (i == 1) {
          n0 <- sum(df[[v_clean]] == lvl, na.rm = TRUE)
          lab_list[[length(lab_list)+1]] <-
            c(v_clean,
              paste0(lvl, " (N=", n0, ")"),
              "reference", "", "")
          meanv  <- c(meanv, 1); lowerv <- c(lowerv, 1); upperv <- c(upperv, 1)
        } else {
          term_nm <- paste0(v_raw, lvl)
          rx      <- tr[tr$term == term_nm, ]
          n1      <- sum(df[[v_clean]] == lvl, na.rm = TRUE)
          
          p_str   <- format_p(rx$p.value)
          fdr_str <- format_p(rx$fdr)

          lab_list[[length(lab_list)+1]] <-
            c("",
              paste0(lvl, " (N=", n1, ")"),
              sprintf("%.2f (%.2f–%.2f)", rx$estimate, rx$conf.low, rx$conf.high),
              p_str, fdr_str)
          meanv  <- c(meanv,  rx$estimate)
          lowerv <- c(lowerv, rx$conf.low)
          upperv <- c(upperv, rx$conf.high)
        }
      }
      
    } else if (grepl(":", v_raw, fixed = TRUE)) {
      parts_raw   <- strsplit(v_raw, ":", fixed = TRUE)[[1]]
      parts_clean <- gsub("`", "", parts_raw)
      
      lv1   <- levels(df[[parts_clean[1]]])
      lv2   <- levels(df[[parts_clean[2]]])
      ref1  <- lv1[1]; ref2 <- lv2[1]
      n_ref <- sum(df[[parts_clean[1]]] == ref1 &
                     df[[parts_clean[2]]] == ref2, na.rm = TRUE)
      lab_list[[length(lab_list)+1]] <-
        c(v_clean,
          paste0(ref1, " & ", ref2, " (N=", n_ref, ")"),
          "reference", "", "")
      meanv  <- c(meanv, 1); lowerv <- c(lowerv, 1); upperv <- c(upperv, 1)
      
      pattern <- paste0("^", parts_raw[1], ".*:", parts_raw[2])
      rx_all  <- tr[grepl(pattern, tr$term), ]
      for (k in seq_len(nrow(rx_all))) {
        rx         <- rx_all[k, ]
        tsplit     <- strsplit(rx$term, ":", fixed = TRUE)[[1]]
        combo1_raw <- tsplit[1]
        combo2_raw <- tsplit[2]
        lvl1 <- sub(paste0("^", parts_raw[1]), "", combo1_raw)
        lvl1 <- gsub("`", "", lvl1)
        lvl2 <- sub(paste0("^", parts_raw[2]), "", combo2_raw)
        lvl2 <- gsub("`", "", lvl2)
        n_cmp <- sum(df[[parts_clean[1]]] == lvl1 &
                       df[[parts_clean[2]]] == lvl2, na.rm = TRUE)
        
        p_str   <- format_p(rx$p.value)
        fdr_str <- format_p(rx$fdr)

        lab_list[[length(lab_list)+1]] <-
          c("",
            paste0(lvl1, " & ", lvl2, " (N=", n_cmp, ")"),
            sprintf("%.2f (%.2f–%.2f)", rx$estimate, rx$conf.low, rx$conf.high),
            p_str, fdr_str)
        meanv  <- c(meanv,  rx$estimate)
        lowerv <- c(lowerv, rx$conf.low)
        upperv <- c(upperv, rx$conf.high)
      }
      
    } else {
      n0 <- if (v_clean %in% names(df)) sum(!is.na(df[[v_clean]])) else nrow(df)
      rx     <- tr[tr$term == v_raw, ]
      
      p_str   <- format_p(rx$p.value)
      fdr_str <- format_p(rx$fdr)

      lab_list[[length(lab_list)+1]] <-
        c(v_clean,
          paste0("N=", n0),
          sprintf("%.2f (%.2f–%.2f)", rx$estimate, rx$conf.low, rx$conf.high),
          p_str, fdr_str)
      meanv  <- c(meanv,  rx$estimate)
      lowerv <- c(lowerv, rx$conf.low)
      upperv <- c(upperv, rx$conf.high)
    }
  }
  
  lab_mat <- do.call(rbind, lab_list)
  header  <- c("Predictor", "N", "HR (95% CI)", "p-value", "FDR")
  lab_mat <- rbind(header, lab_mat)
  meanv   <- c(NA, meanv)
  lowerv  <- c(NA, lowerv)
  upperv  <- c(NA, upperv)
  
  rows <- which(lab_mat[,1] != "" & seq_len(nrow(lab_mat)) != 1)
  hrzl_sep <- setNames(
    lapply(rows, function(i) grid::gpar(col = "#DDDDDD", lwd = 1)),
    as.character(rows)
  )
  
  is_summary <- c(TRUE, rep(FALSE, nrow(lab_mat) - 1))
  
  forestplot::forestplot(
    labeltext  = lab_mat,
    mean       = meanv,
    lower      = lowerv,
    upper      = upperv,
    zero       = 1,
    xlog       = TRUE,
    graph.pos  = 4,
    boxsize    = 0.2,
    hrzl_lines = hrzl_sep,
    is.summary = is_summary,
    txt_gp     = forestplot::fpTxtGp(
      label = grid::gpar(cex = 0.8),
      ticks = grid::gpar(cex = 0.7),
      xlab  = grid::gpar(cex = 0.8)
    ),
    col        = forestplot::fpColors(box = "black", line = "black"),
    xlab       = "Hazard Ratio"
  )
}

dataset=fread('../Demographics for Natera.csv')
# dataset=fread('/path/to/your/Demographics for Natera.csv')

dataset[,DFS.event:=`DFS event=1 free=0`]
dataset[,histological_type:=`Degree of Histological Differentiation por/sig or not`]
dataset[,`Histological type`:=histological_type]
dataset[`Histological type`=='',`Histological type`:=NA]
dataset[`MSI category`=='',`MSI category`:=NA]
dataset[RAS=='',RAS:=NA]
dataset[RAS=='Not Done',RAS:=NA]
dataset[RAS=='Unknown',RAS:=NA]
dataset[`BRAF category`=='',`BRAF category`:=NA]
dataset[`ctDNA Result 4W`=='TNP',`ctDNA Result 4W`:=NA]
dataset[`ctDNA Result 4W`=='NC',`ctDNA Result 4W`:=NA]
dataset[`ctDNA Result 4W`=='',`ctDNA Result 4W`:=NA]
dataset[,MSI:=factor(`MSI category`,levels=c('MSI-HIGH','MSS'))]
dataset[,BRAF:=`BRAF category`]
dataset[,ACT:=`Perioperative treatmentChemotherapy`]
dataset[,`ctDNA status at 4 weeks post surgery`:=`ctDNA Result 4W`]
dataset[,location:=`Location 解析4 colon vs rectum`]
dataset[location=='COLON',location:='Colon']
dataset[location=='RECTUM',location:='Rectum']
dataset[BRAF=='BRAF neg',BRAF:='Negative (wild type)']
dataset[BRAF=='BRAF posi',BRAF:='Positive (mutant)']
dataset[`HER2 amplification`=='negative',`HER2 amplification`:='Negative']
dataset[`HER2 amplification`=='positive',`HER2 amplification`:='Positive']
dataset[`Histological type`=='Por/Sig',`Histological type`:='Poorly differentiated /Signet ring']

cox_fit <- coxph(Surv(DFS, DFS.event) ~ Sex + Age + location +`Stage I,II,III`+RAS+BRAF+MSI + `HER2 amplification` + `ctDNA status at 4 weeks post surgery` + `Histological type`, data = dataset)

summary(cox_fit)
forest_plot_cox(cox_fit)

