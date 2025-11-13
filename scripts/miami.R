#!/usr/bin/env Rscript

# Script to generate Miami plots from GWAS summary statistics
# Tailored for FinnGen GWAS meta-analysis summary statistics

packs <- c("qqman", "optparse", "data.table", "R.utils")

for (p in packs) {
  if (!require(p, character.only = TRUE)) {
    print(p)
    install.packages(p,  repos = c(CRAN = "http://cran.r-project.org"))
  }
}

option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = NULL,
              help = "Summary statistics file",
              metavar = "character"),
  make_option(c("-o", "--out"), type = "character",
              help = "output file name [default= %default]",
              metavar = "character"),
  make_option(c("--chr_col"), type = "character", default = "#CHR",
              help = "chromosome column [default= %default]",
              metavar = "character"),
  make_option(c("--pos_col"), type = "character", default = "POS",
              help = "pos column [default= %default]",
              metavar = "character"),
  make_option(c("-p", "--pval_cols"), type = "character", default = "P",
              help = "Two p-value columns, comma-separated",
              metavar = "character"),
  make_option("--pvalue_type", type = "character", default = "detect",
              help = "Type of p-values: 'p' for p-values, 'mlogp' for -log10 p-values, 'detect' for automatically detecting based on data [default= %default]",
              metavar = "character"),
  make_option("--highlight", action = "store_true", default = FALSE,
              help = "Highlight variants not genome-wide significant in the other p-value column")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser, positional_arguments = 0)

options(bitmapType = "cairo")

print(str(opt))
pos_col <- opt$options$pos_col
chr_col <- opt$options$chr_col

pval_cols <- unlist(strsplit(opt$options$pval_cols, ","))

if (length(pval_cols) != 2) {
  stop("Please provide exactly two p-value columns, comma-separated")
}

pvalue_type <- opt$options$pvalue_type
if (! pvalue_type %in% c("p", "mlogp", "detect")) {
  stop("Invalid pvalue-type: ", pvalue_type)
}

file <- opt$options$file
print(paste("reading file:", file))
data <- fread(file, header = TRUE, select = c(pval_cols, c(pos_col, chr_col)))

if (!is.null(opt$options$out)) {
  output_prefix <- opt$options$out
} else {
  output_prefix <- file
}

if (any(! append(pval_cols, c(pos_col, chr_col)) %in% colnames(data))) {
  stop("All required columns do not exist in the data: ", paste(c(pval_cols, pos_col, chr_col), collapse = ", "))
}

if (pvalue_type == "detect") {
  # Detect p-value type based on data
  pvalue_type <- "p"
  for (pcol in pval_cols) {
    if (any(data[[pcol]] > 1)) {
      pvalue_type <- "mlogp"
      break
    }
  }
}

if (pvalue_type == "p") {
  # Convert p-values to -log10
  # If any p-value is zero, set it to smallest non-zero value allowed by R
  data[, (pval_cols) := lapply(.SD, function(x) ifelse(x==0, 5e-324, x)), .SDcols = pval_cols]
  data[, (pval_cols) := lapply(.SD, function(x) -log10(x)), .SDcols = pval_cols]
}
data <- na.omit(data)[data[, Reduce(`|`, lapply(.SD, `>`, -log10(0.01))), .SDcols = pval_cols]]

data[, (chr_col):= gsub("^chr", "", get(chr_col))]
data[, (chr_col):= gsub("X", "23", get(chr_col))]
data[, (chr_col):= gsub("Y", "24", get(chr_col))]
data[, (chr_col):= gsub("MT|M", "25", get(chr_col))]
data[, (chr_col):= as.numeric(get(chr_col))]

data[, SNP := as.character(1:.N)]  # Dummy SNP column for manhattan function

maxp <- max(data[, max(.SD), .SDcols = pval_cols])  # maximum -log10 p-value across both columns

highlight <- opt$options$highlight
if (highlight) {
  hilight_snps1 <- data[get(pval_cols[1]) > -log10(5e-8) & get(pval_cols[1]) > get(pval_cols[2]), SNP]
  hilight_snps2 <- data[get(pval_cols[2]) > -log10(5e-8) & get(pval_cols[2]) > get(pval_cols[1]), SNP]
  if (length(hilight_snps1) == 0) {
    hilight_snps1 <- NULL
  }
  if (length(hilight_snps2) == 0) {
    hilight_snps2 <- NULL
  }
} else {
  hilight_snps1 <- NULL
  hilight_snps2 <- NULL
}

png(paste(output_prefix, pcol, "miami.png", sep = "_"), 1300, 800)
par(mfrow = c(2, 1))
par(mar = c(1, 6, 3, 3))
manhattan(data,
          chr = chr_col,
          bp = pos_col,
          p = pval_cols[1],
          snp = "SNP",
          ylim = c(0, maxp + 0.2 * maxp),
          cex.lab = 2.5,
          font.lab = 2,
          font.axis = 2,
          cex.axis = 1.6,
          suggestiveline = F,
          chrlabs = c(1:22, "X"),
          logp = F,
          highlight = hilight_snps1,
          main = pval_cols[1],
          cex.main = 1.5,
          las = 1,
          col = c(rgb(53, 0, 212, maxColorValue = 255), rgb(40, 40, 40, maxColorValue = 255)),
          font = 4,
          xlab = "")
par(mar = c(5, 6, 1.6, 3))
manhattan(data,
          chr = chr_col,
          bp = pos_col,
          p = pval_cols[2],
          snp = "SNP",
          ylim = c(maxp + 0.2 * maxp,0),
          cex.lab = 2.5,
          font.lab = 2,
          font.axis = 2,
          cex.axis = 1.6,
          suggestiveline = F,
          logp = F,
          highlight = hilight_snps2,
          sub = pval_cols[2],
          cex.sub = 1.5,
          font.sub = 2,
          las = 1,
          col = c("grey20", "grey80"),
          font = 4,
          xlab = "",
          xaxt = "n")
dev.off()