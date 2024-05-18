#!/usr/bin/env Rscript

packs <- c("ggplot2", "data.table", "R.utils", "optparse", "rjson", "stringi", "ggpubr")

for (p in packs) {
  if (suppressPackageStartupMessages(!require(p, character.only = T))) {
    install.packages(p,  repos = c(CRAN = "http://cran.r-project.org"))
  }
}

option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output file name [default=%default]", metavar="character"),
  make_option(c("-m","--method"), type="character", default="inv_var",
              help="meta-analysis method [default=%default]", metavar="character"),
  make_option(c("-l","--loo"), action="store_true", default=FALSE,
              help="use leave-one-out results"),
  make_option("--conf", type="character", default=NULL,
              help="meta-analysis config json", metavar="character"),
  make_option("--pval_thresh", type="character", default="5e-8, 5e-6",
              help="comma separated list of p-value thresholds used to filter the data for qc", metavar="numeric"),
  make_option("--region", type="numeric", default=1.0,
              help="region size in megabases used when counting unique loci hits", metavar="numeric"),
  make_option(c("-c","--chr_col"), type="character", default="#CHR",
              help="chromosome column [default=%default]", metavar="character"),
  make_option(c("-b","--bp_col"), type="character", default="POS",
              help="bp column [default=%default]", metavar="character"),
  make_option(c("-r","--ref_col"), type="character", default="REF",
              help="ref column [default=%default]", metavar="character"),
  make_option(c("-a","--alt_col"), type="character", default="ALT",
              help="alt column [default=%default]", metavar="character"),
  make_option(c("--af_alt_col_suffix"), type="character", default="_af_alt",
              help="af alt column suffix [default=%default]", metavar="character"),
  make_option(c("--pheno"), type="character", default="pheno",
              help="phenotype name [default=%default]", metavar="character"),
  make_option(c("-w","--weighted"), action="store_true", default=FALSE,
              help="do inverse variance weighted linear regression"),
  make_option(c("--keep-hla"), action="store_true", default=FALSE,
              help="do not remove HLA region variants from QC"),
);

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser, positional_arguments = 0)

options(bitmapType = "cairo")

print(str(opt))
chr_col <- opt$options$chr_col
bp_col <- opt$options$bp_col
ref_col <- opt$options$ref_col
alt_col <- opt$options$alt_col
af_alt_col_suffix <- opt$options$af_alt_col_suffix
pheno <- opt$options$pheno
method <- opt$options$method
pval_thresh <- sort(unique(as.numeric(unlist(strsplit(gsub("[[:blank:]]", "", opt$options$pval_thresh), ",")))), decreasing = T)
leave <- opt$options$loo
region <- round(opt$options$region * 10^6 / 2, 0)
weighted <- opt$options$weighted
keep_hla <- opt$options$keep_hla

file <- opt$options$file
conf <- rjson::fromJSON(file = opt$options$conf)

qc_dt <- data.table(pheno = pheno)
for (s in conf$meta) {
  qc_dt[, (paste0(s$name, "_n_cases")) := s$n_cases]
  qc_dt[, (paste0(s$name, "_n_controls")) := s$n_controls]
}

studies <- sapply(conf$meta, function(x) x$name)

study_pval_cols <- sapply(conf$meta, function(x) paste(x$name, x$pval, sep = "_"))
meta_pval_col <- paste("all", method, "meta_p", sep = "_")

study_beta_cols <- sapply(conf$meta, function(x) paste(x$name, x$effect, sep = "_"))
meta_beta_col <- paste("all", method, "meta_beta", sep = "_")

het_p_col <- paste("all", method, "het_p", sep = "_")

af_cols <- paste0(studies, af_alt_col_suffix)

pval_cols <- c(study_pval_cols, meta_pval_col)
meta_pval_cols <- meta_pval_col

keep_cols <- c(chr_col, bp_col, ref_col, alt_col, af_cols,
               study_pval_cols, meta_pval_col,
               study_beta_cols, meta_beta_col,
               het_p_col)

if (weighted) {
  study_sebeta_cols <- sapply(conf$meta, function(x) paste(x$name, x$se, sep = "_"))
  meta_sebeta_col <- paste("all", method, "meta_sebeta", sep = "_")
  keep_cols <- c(keep_cols, study_sebeta_cols, meta_sebeta_col)
}

if (leave) {
  leave_pval_cols <- paste("leave", studies, method, "meta_p", sep = "_")
  leave_beta_cols <- paste("leave", studies, method, "meta_beta", sep = "_")
  #leave_het_p_cols <- paste("leave", studies, method, "meta_het_p", sep = "_")
  keep_cols <- c(keep_cols, leave_pval_cols, leave_beta_cols)
  pval_cols <- c(pval_cols, leave_pval_cols)
  meta_pval_cols <- c(meta_pval_cols, leave_pval_cols)
}

message("Reading file ", file, " ...")
data <- fread(file, header = T, select = keep_cols)

# Remove variants not in reference
data <- data[! is.na(get(study_pval_cols[[1]]))]

# Remove variants with weak pvals that will be never used
pass <- apply(data[, ..pval_cols], 1, function(x) any(x < max(pval_thresh), na.rm = T))
data <- data[pass]

# Remove HLA region variants
if (! keep_hla) {
  data <- data[! (data[[chr_col]] == 6 & data[[bp_col]] >= 25e6 & data[[bp_col]] <= 35e6)]
}

for (pval_thresh_i in pval_thresh) {
  
  message("Calculating qc metrics with p-value threshold ", pval_thresh_i, " ...")
  
  output_prefix <- paste0(ifelse(test = is.null(opt$options$out), yes = file, no = opt$options$out), ".qc.", pval_thresh_i)
  
  qc_dt_i <- copy(qc_dt)
  
  pass <- apply(data[, ..pval_cols], 1, function(x) any(x < pval_thresh_i, na.rm = T))
  DATA <- data[pass]
  
  # Get gw signif hits
  sig_loc_list <- list()
  for (pval_col in pval_cols) {
    tempdata <- DATA[DATA[[pval_col]] < pval_thresh_i]
    n_sig_loci <- 0
    message("Finding significant loci according to \"", pval_col, "\"")
    while (nrow(tempdata) > 0) {
      maxrow <- tempdata[which.min(tempdata[[pval_col]])]
      maxrow_u_bp <- maxrow[[bp_col]] + region
      maxrow_l_bp <- maxrow[[bp_col]] - region
      maxrow_chr <- maxrow[[chr_col]]
      in_region <- tempdata[[chr_col]] == maxrow_chr & tempdata[[bp_col]] >= maxrow_l_bp & tempdata[[bp_col]] <= maxrow_u_bp
      if (anyNA(maxrow[, ..study_pval_cols])) {
        message("Top hit in a region is not found from all of the studies: ",
                paste(maxrow[, 1:4], collapse = "-"))
        #temp_sig_loci <- na.omit(tempdata[in_region])
        temp_sig_loci <- tempdata[in_region]
        tempdt <- data.table(N = apply(is.na(temp_sig_loci[, ..study_pval_cols]), 1, sum),
                             pval = temp_sig_loci[[pval_col]])
        temp_sig_loci <- temp_sig_loci[order(tempdt$N, tempdt$pval)]
        if (! equals(temp_sig_loci[1,], maxrow)) {
          maxrow <- temp_sig_loci[1,]
          message("Replaced top hit in region with the next best SNP found in more studies: ",
                  paste(maxrow[, 1:4], collapse = "-"))
        } else {
          message("Could not find better replacement for the top hit in the region. This may effect the heterogeneity calculations. Number of significant SNPs in region: ", nrow(temp_sig_loci))
        }
      }
      tempdata <- tempdata[! in_region]
      n_sig_loci <- n_sig_loci + 1
      if (is.null(sig_loc_list[[pval_col]])) {
        sig_loc_list[[pval_col]] <- maxrow
      } else {
        sig_loc_list[[pval_col]] <- rbind(sig_loc_list[[pval_col]], maxrow)
      }
    }
    message("Found ", n_sig_loci, " significant loci.")
    qc_dt_i[, (sapply(strsplit(pval_col, "_"), function(x) paste(c(x[1:(length(x)-1)], "N_hits"), collapse = "_"))) := n_sig_loci]
  }
  rm(tempdata)
  
  ref_hits <- sig_loc_list[[pval_cols[1]]]
  
  beta_cols <- c(study_beta_cols[-1], meta_beta_col)
  sebeta_cols <- c(study_sebeta_cols[-1], meta_sebeta_col)
  for (i in 1:length(beta_cols)) {
    tryCatch(
      expr = {
        if (weighted) {
          weights <- 1 / ref_hits[[sebeta_cols[i]]]^2
          m <- lm(as.formula(paste(beta_cols[i], "~", study_beta_cols[1], "+ 0")), data = ref_hits, weights = weights)
        } else {
          m <- lm(as.formula(paste(beta_cols[i], "~", study_beta_cols[1], "+ 0")), data = ref_hits)
        }
        ms <- summary(m)
        qc_dt_i[, (paste(study_beta_cols[1], "vs", beta_cols[i], "slope", sep = "_")) := round(ms$coefficients[1,1], 3)]
        qc_dt_i[, (paste(study_beta_cols[1], "vs", beta_cols[i], "r2", sep = "_")) := round(ms$r.squared, 3)]
        qc_dt_i[, (paste(study_beta_cols[1], "vs", beta_cols[i], "r2adj", sep = "_")) := round(ms$adj.r.squared, 3)]
      },
      error = function(e) {
        qc_dt_i[, (paste(study_beta_cols[1], "vs", beta_cols[i], "slope", sep = "_")) := NA]
        qc_dt_i[, (paste(study_beta_cols[1], "vs", beta_cols[i], "r2", sep = "_")) := NA]
        qc_dt_i[, (paste(study_beta_cols[1], "vs", beta_cols[i], "r2adj", sep = "_")) := NA]
      }
    )
  }
  
  for (pval_col in meta_pval_cols) {
    tryCatch(
      expr = {
        qc_dt_i[, (paste("pct_pval_stronger_in", pval_col, "vs", studies[1], sep = "_")) := round(sum(ref_hits[[pval_cols[1]]] > ref_hits[[pval_col]], na.rm = T) / nrow(ref_hits), 3)]
      },
      error = function(e) {
        qc_dt_i[, (paste("pct_pval_stronger_in", pval_col, "vs", studies[1], sep = "_")) := NA]
      }
    )
  }
  
  tryCatch(
    expr = {
      ref_hits[, het_p_fdr := p.adjust(ref_hits[[het_p_col]], method = "fdr")]
      qc_dt_i[, het_p_fdr_signif_in_meta_pct := round(sum(ref_hits[["het_p_fdr"]] < 0.05, na.rm = T) / sum(!is.na(ref_hits[["het_p_fdr"]])), 3)]
    },
    error = function(e) {
      qc_dt_i[, het_p_fdr_signif_in_meta_pct := NA]
    }
  )
  
  fwrite(qc_dt_i, paste0(output_prefix, ".tsv"), col.names = T, row.names = F, quote = F, sep = "\t", na = "NA")
  
  # Don't try plotting if less than two significant SNPs
  if (is.null(ref_hits)) {
    message("Skipping plotting due to no significant hits.")
    next
  }
  
  for (h in names(sig_loc_list)) {
    fwrite(sig_loc_list[[h]], paste(output_prefix, h, "hits", sep = "."), col.names = T, row.names = F, quote = F, sep = "\t", na = "NA")
  }
  
  if (nrow(ref_hits) < 2) {
    message("Skipping plotting due to less than 2 significant hits.")
    next
  }
  
  # Change afs to mafs
  for (af_col in af_cols) {
    af <- ref_hits[[af_col]]
    af[af > 0.5 & ! is.na(af)] <- 1 - af[af > 0.5 & ! is.na(af)]
    ref_hits[, (af_col) := af]
  }
  rm(af)
  
  # P-value comparisons scatter plot
  pval_plots <- lapply(2:length(study_pval_cols), function(i) {
    tempcols <- study_pval_cols[c(1,i)]
    af_col <- af_cols[i]
    tempdata <- na.omit(cbind(-log10(ref_hits[, ..tempcols]), ref_hits[, ..af_col]))
    tempdata[[af_col]] = cut(tempdata[[af_col]], breaks = c(0.5, 0.2, 0.05, 0.01, 0))
    color_vector <- c("#ca0020", "#f4a582", "#92c5de", "#0571b0")
    names(color_vector) <- levels(tempdata[[af_col]])
    p <- ggplot(tempdata, aes_string(x = tempcols[1], y = tempcols[2], color = af_col)) +
      geom_point() +
      #geom_smooth(method = "lm", se = F) +
      geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
      expand_limits(x = 0, y = 0) +
      labs(x = paste(studies[1], "mlogp"),
           y = paste(studies[i], "mlogp"),
           color = "Ext maf") +
      scale_color_manual(labels = names(color_vector), values = color_vector, drop = F) +
      theme_bw() +
      theme(axis.title = element_text(size = 7),
            axis.text = element_text(size = 7))
    
    min_limit <- min(ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[1], ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[1])
    max_limit <- max(ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[2], ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2])
    p <- p + xlim(min_limit, max_limit) + ylim(min_limit, max_limit)
    
    p
  })
  
  
  # Effect size comparisons scatter plot
  beta_plots <- lapply(2:length(study_beta_cols), function(i) {
    af_col <- af_cols[i]
    tempcols <- c(study_beta_cols[c(1,i)], af_col)
    tempdata <- na.omit(ref_hits[, ..tempcols])
    neg_betas <- tempdata[[tempcols[1]]] < 0
    new_beta1 <- tempdata[[tempcols[1]]]
    new_beta1[neg_betas] <- -new_beta1[neg_betas]
    new_beta2 <- tempdata[[tempcols[2]]]
    new_beta2[neg_betas] <- -new_beta2[neg_betas]
    tempdata[, (tempcols[1]) := new_beta1]
    tempdata[, (tempcols[2]) := new_beta2]
    tempdata[[af_col]] = cut(tempdata[[af_col]], breaks = c(0.5, 0.2, 0.05, 0.01, 0))
    color_vector <- c("#ca0020", "#f4a582", "#92c5de", "#0571b0")
    names(color_vector) <- levels(tempdata[[af_col]])
    p <- ggplot(tempdata, aes_string(x = tempcols[1], y = tempcols[2], color = af_col)) +
      geom_point() +
      stat_smooth(method = "lm", se = F, formula = "y ~ x + 0", fullrange = T, size = .8) +
      geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
      labs(x = paste(studies[1], "beta"),
           y = paste(studies[i], "beta"),
           color = "Ext maf") +
      scale_color_manual(labels = names(color_vector), values = color_vector, drop = F) +
      theme_bw() +
      theme(axis.title = element_text(size = 7),
            axis.text = element_text(size = 7))
    
    min_limit <- min(ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[1], ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[1])
    max_limit <- max(ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[2], ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2])
    p <- p + xlim(min_limit, max_limit) + ylim(min_limit, max_limit)
    
    p
  })
  
  # Het p histogram
  tempdata <- -log10(na.omit(ref_hits[, ..het_p_col]))
  het_hist <- ggplot(tempdata, aes_string(x = het_p_col)) +
    geom_histogram() +
    xlab(paste("all", method, "het_mlogp", sep = "_")) +
    theme_bw() +
    theme(axis.title = element_text(size = 7),
          axis.text = element_text(size = 7))
  
  # P-value comparison histogram
  p_cols <- c(meta_pval_col, study_pval_cols[1])
  tempdata <- -log10(na.omit(ref_hits[, ..p_cols]))
  diff <- data.table(tempdata[[1]] - tempdata[[2]])
  pval_hist <- ggplot(diff, aes(x = V1)) + 
    geom_histogram() +
    xlab(paste0("all_", method, "_meta_mlogp - ", studies[1], "_mlogp")) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    theme_bw() +
    theme(axis.title = element_text(size = 7),
          axis.text = element_text(size = 7))
  
  # Leave-one-out qc
  if (leave) {
    p_cols <- c(study_pval_cols[1], leave_pval_cols[-1])
    tempdata <- -log10(ref_hits[, ..p_cols])
    loo_hists <- lapply(2:length(p_cols), function(i) {
      diff <- na.omit(data.table(tempdata[[i]] - tempdata[[1]]))
      p <- ggplot(diff, aes(x = V1)) +
        geom_histogram() +
        xlab(paste0("leave_", studies[i], "_", method, "_meta_mlogp - ", studies[1], "_mlogp")) +
        geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
        theme_bw() +
        theme(axis.title = element_text(size = 7),
              axis.text = element_text(size = 7))
      p
    })
  }
  
  pval_plots_arranged <- ggarrange(plotlist = pval_plots, common.legend = T, legend = "bottom", nrow = 1)
  beta_plots_arranged <- ggarrange(plotlist = beta_plots, legend = "none", nrow = 1)
  hists_arranged <- ggarrange(het_hist, pval_hist, legend = "none", nrow = 1)
  
  plotlist <- list(pval_plots_arranged, beta_plots_arranged, hists_arranged)
  
  if (leave) {
    loo_hists_arranged <- ggarrange(plotlist = loo_hists, legend = "none", nrow = 1)
    plotlist <- append(plotlist, list(loo_hists_arranged))
  }
  
  plots_arranged <- ggarrange(plotlist = plotlist, ncol = 1) +
    theme(axis.title = element_text(size = 7),
          axis.text = element_text(size = 7))

  ggsave(filename = paste0(output_prefix, ".pdf"),
         plot = plots_arranged,
         device = "pdf",
         dpi = 300,
         width = 7,
         height = 9)
  
}
