#!/usr/bin/env Rscript

packs <- c("ggplot2", "data.table", "R.utils", "optparse", "rjson")

for (p in packs) {
  if (!require(p, character.only = T)) {
    install.packages(p,  repos = c(CRAN = "http://cran.r-project.org"))
  }
}

option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-m","--method"), type="character", default="inv_var",
              help="meta-analysis method [default= %default]", metavar="character"),
  make_option(c("-l","--loo"), action="store_true",
              help="use leave-one-out results [default= %default]"),
  make_option("--conf", type="character", default=NULL,
              help="meta-analysis config json", metavar="character"),
  make_option("--pval_thresh", type="character", default=5e-8,
              help="meta-analysis config json", metavar="numeric"),
  make_option(c("-c","--chr_col"), type="character", default="#CHR",
              help="chromosome column [default= %default]", metavar="character"),
  make_option(c("-b","--bp_col"), type="character", default="POS",
              help="bp column [default= %default]", metavar="character"),
  make_option(c("-r","--ref_col"), type="character", default="REF",
              help="ref column [default= %default]", metavar="character"),
  make_option(c("-a","--alt_col"), type="character", default="ALT",
              help="alt column [default= %default]", metavar="character")
);

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser, positional_arguments = 0)

options(bitmapType = "cairo")

print(str(opt))
chr_col <- opt$options$chr_col
bp_col <- opt$options$bp_col
ref_col <- opt$options$ref_col
alt_col <- opt$options$alt_col
method <- opt$options$method
pval_thresh <- opt$options$pval_thresh
leave <- opt$options$loo

file <- opt$options$file
conf <- rjson::fromJSON(file = opt$options$conf)

studies <- sapply(conf$meta, function(x) x$name)

study_pval_cols <- sapply(conf$meta, function(x) paste(x$name, x$pval, sep = "_"))
meta_pval_col <- paste("all", method, "meta_p", sep = "_")

study_beta_cols <- sapply(conf$meta, function(x) paste(x$name, x$effect, sep = "_"))
meta_beta_col <- paste("all", method, "meta_beta", sep = "_")

het_p_col <- paste("all", method, "het_p", sep = "_")

keep_cols <- c(chr_col, bp_col, ref_col, alt_col, study_pval_cols, meta_pval_col, study_beta_cols, meta_beta_col, het_p_col)

if (leave) {
  leave_pval_cols <- paste("leave", studies, method, "meta_p", sep = "_")
  leave_beta_cols <- paste("leave", studies, method, "meta_beta", sep = "_")
  #leave_het_p_cols <- paste("leave", studies, method, "meta_het_p", sep = "_")
  keep_cols <- c(keep_cols, leave_pval_cols, leave_beta_cols)
}

data <- fread(file, header = T, select = keep_cols)

# Get gw signif hits, 1MB regions, each study
# TODO
tempdata <- data[data[[meta_pval_col]] < 5e-8]

pass <- data[[study_pval_cols[1]]] < pval_thresh
while (sum(pass) < 5) {
  pval_thresh <- pval_thresh + 1e-8
  pass <- data[[study_pval_cols[1]]] < pval_thresh
}

data <- data[pass]

qc_dt <- data.table()

m <- lm(as.formula(paste(study_beta_cols[1], "~", meta_beta_col)), data = data)
ms <- summary(m)
qc_dt[, (paste(studies[1], "vs_meta_beta_slope", sep = "_")) := round(ms$coefficients[2,1], 3)]
qc_dt[, (paste(studies[1], "vs_meta_beta_r2", sep = "_")) := round(ms$r.squared, 3)]
qc_dt[, (paste(studies[1], "vs_meta_beta_r2adj", sep = "_")) := round(ms$adj.r.squared, 3)]

output_prefix <- ifelse(test = is.null(opt$options$out), yes = file, no = opt$options$out)
output_prefix <- paste0(output_prefix, "_qc.pval_lt_", pval_thresh)
pdf(paste0(output_prefix, ".pdf"))

# P-value comparisons scatter plot
for (i in 2:length(study_pval_cols)) {
  tempcols <- study_pval_cols[c(1,i)]
  tempdata <- -log10(na.omit(data[, ..tempcols]))
  p <- ggplot(tempdata, aes_string(x = tempcols[1], y = tempcols[2])) +
    geom_point() +
    #geom_smooth(method = "lm", se = F) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    expand_limits(x = 0, y = 0) +
    xlab(paste0(studies[1], "_mlogp")) +
    ylab(paste0(studies[i], "_mlogp")) +
    theme_bw()
  plot(p)
}

# Effect size comparisons scatter plot
for (i in 2:length(study_beta_cols)) {
  tempcols <- study_beta_cols[c(1,i)]
  tempdata <- na.omit(data[, ..tempcols])
  p <- ggplot(tempdata, aes_string(x = tempcols[1], y = tempcols[2])) +
    geom_point() +
    #geom_smooth(method = "lm", se = F) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    #expand_limits(x = 0, y = 0) +
    #xlab(paste0("abs(", studies[1], "_beta)")) +
    #ylab(paste0("abs(", studies[i], "_beta)")) +
    xlab(paste0(studies[1], "_beta")) +
    ylab(paste0(studies[i], "_beta")) +
    theme_bw()
  plot(p)

}

# Het p histogram
tempdata <- -log10(na.omit(data[, ..het_p_col]))
p <- ggplot(tempdata, aes_string(x = het_p_col)) + 
  geom_histogram() +
  xlab(paste("all", method, "het_mlogp", sep = "_")) +
  ggtitle("Heterogeneity -log10(p-value) distribution") +
  theme_bw()
plot(p)

# P-value comparison histogram
pval_cols <- c(meta_pval_col, study_pval_cols[1])
tempdata <- -log10(na.omit(data[, ..pval_cols]))
diff <- data.table(tempdata[[1]] - tempdata[[2]])
p <- ggplot(diff, aes(x = V1)) + 
  geom_histogram(binwidth = .5) +
  xlab(paste0("all_", method, "_meta_mlogp - ", studies[1], "_mlogp")) +
  ggtitle(paste("Is", studies[1], "p-value getting stronger in meta-analysis (all studies)")) +
  theme_bw()
plot(p)

qc_dt[, (paste("pct_pval_stronger_in", studies[1], "vs_meta", sep = "_")) := round(sum(diff > 0) / length(diff[[1]]), 3)]


# Leave-one-out qc
if (leave) {
  pval_cols <- c(study_pval_cols[1], leave_pval_cols[-1])
  tempdata <- -log10(data[, ..pval_cols])
  for (i in 2:length(pval_cols)) {
    diff <- na.omit(data.table(tempdata[[i]] - tempdata[[1]]))
    p <- ggplot(diff, aes(x = V1)) +
      geom_histogram(binwidth = .5) +
      xlab(paste0("leave_", studies[i], "_", method, "_meta_mlogp - ", studies[1], "_mlogp")) +
      ggtitle(paste0("Is ", studies[1], " p-value getting stronger in meta-analysis (drop ", studies[i], ")")) +
      theme_bw()
    plot(p)
    
    qc_dt[, (paste("pct_pval_stronger_in", studies[1], "vs_meta_leave", studies[i], sep = "_")) := round(sum(diff > 0) / length(diff[[1]]), 3)]
    
  }
}

graphics.off()

fwrite(qc_dt, paste0(output_prefix, ".tsv"), col.names = T, row.names = F, quote = F, sep = "\t")
