# META ANALYSIS tools

Tools for doing x-way meta-analysis of GWAS summary statistics including data munging, harmonization to reference, meta-analysis, visualization and QC.

## Table of Contents

* [Variant matching across studies](#variant-matching-across-studies)
* [Running single trait meta-analysis](#running-single-trait-meta-analysis)
* [Data Preparation and Munging](#data-preparation-and-munging)
  * [Munging Summary Statistics](#munging-summary-statistics)
  * [Harmonizing Summary Statistics](#harmonizing-summary-statistics)
* [Visualization](#visualization)
  * [QQ and Manhattan Plots](#qq-and-manhattan-plots)
  * [Miami Plots](#miami-plots)
* [Meta-analysis QC](#meta-analysis-qc)
* [Workflow Execution (WDL)](#workflow-execution-wdl)
  * [Meta-Analysis Workflow](#meta-analysis-workflow)
  * [Munging Workflow](#munging-workflow)
  * [Liftover Workflow](#liftover-workflow)
* [Utility Scripts](#utility-scripts)
  * [Creating Meta-Analysis Configuration Files](#creating-meta-analysis-configuration-files)
  * [Copy Cromwell Outputs](#copy-cromwell-outputs)
  * [Remove Leave-One-Out Columns](#remove-leave-one-out-columns)
* [Docker Image](#docker-image)

## Variant matching across studies

Variants are matched using chr pos ref and alt. With FinnGen, all 37 build results need to first be lifted over to 38. [lift.wdl](wdl/lift.wdl) can be used to liftover results first if needed.

**IMPORTANT**: Variants in the summary statistics files need to be ordered by chr (1-22, x, y, mt) and position. Chromosome can be indicated with numbers 1-25 or 1-22, X, Y, MT, with or without 'chr' prefix and they will be internally coded to numerical values.

## Running single trait meta-analysis

[meta_analysis.py](scripts/meta_analysis.py) is the main script for running meta-analysis for a single trait. Meta-analyses to be performed are specified with json configuration file. Example configuration with three studies in [data/conf.json](data/conf.json). Script will try to align using both strands as well as by flipping ref vs. alt.

```bash
usage: meta_analysis.py [-h] [--not_quiet] [--leave_one_out] [--is_het_test]
                        [--pairwise_with_first] [--sep SEP] [--chrom CHROM]
                        [--flip_indels]
                        config_file path_to_res methods [methods ...]

Run x-way meta-analysis

positional arguments:
  config_file           Configuration file
  path_to_res           Result file
  methods               Methods to use in calculating meta-analysis statistics.
                        Allowed values: n, inv_var, variance.

optional arguments:
  -h, --help            show this help message and exit
  --not_quiet           Print matching variants to stdout
  --leave_one_out       Do leave-one-out meta-analysis
  --is_het_test         Do heterogeneity tests based on Cochrans Q and output het_p
  --pairwise_with_first
                        Do pairwise meta-analysis with the first given study
  --sep SEP             Input file field separator (default: tab)
  --chrom CHROM         Restrict to given chromosome
  --flip_indels         Try variant aligning by flipping indels also. By default
                        indels are not flipped
```

The configuration file should be a json file with two or more studies with these elements in each:

```json
            "name": "FINNGEN",
            "file": "/path/to/GWAS/summary/statistics.gz",
            "n_cases": 6570 , # number of cases. Used only if sample size weighted meta-analysis is used
            "n_controls": 48378, # number of controls. Used only if sample size weighted meta-analysis is used
            "chr": "CHR", # chromosome column name in the file
            "pos": "POS", # position column name in the file
            "ref": "Allele1", # reference allele column name in the file
            "alt": "Allele2", # alternate allele column name in the file
            "effect": "BETA", # effect size column name in the file
            "effect_type": "beta", # effect type. Allowed values: beta, OR. In case of OR the value will be log transformed to beta.
            "pval": "p.value", # effect size column name in the file
            "se": "SE" # effect standard error. This column is optional if not using inv_var method. Otherwise additional p-value will be added using this as a weight for z-score.
```

[meta_analysis.py](scripts/meta_analysis.py) supports 3 different meta-analysis methods

* `N`: purely weight by sample size and use z-score from p-value,
* `variance`: weight z-score from p-value by variance,
* `inv_var`: regular inverse variance weighted betas meta-analysis.

`inv_var` is recommended if betas and variances are comparable. In case of combining data from different models (e.g.) linear vs. logistic you should use sample size weighted meta.

## Data Preparation and Munging

### Munging Summary Statistics

[munge.py](scripts/munge.py) standardizes and cleans GWAS summary statistics files to prepare them for meta-analysis.

```bash
usage: munge.py [-h] --chr-col CHR_COL --pos-col POS_COL --ref-col REF_COL
                --alt-col ALT_COL --af-col AF_COL --effect-col EFFECT_COL
                --se-col SE_COL --pval-col PVAL_COL
                [--effect-type {beta,or}] [--se-type {se,ci}]
                [--pval-type {p,mlog10p}] [--af-allele {alt,ref}]
                [--flip-alleles] [--recalculate-se]
                [--rounding-precision ROUNDING_PRECISION]
                [--filt-col FILT_COL] [--filt-threshold FILT_THRESHOLD]
                [--delim DELIM] [--verbose]
                sumstat_file

Munge a summary statistics file.

positional arguments:
  sumstat_file          Input summary statistics file. Can be plain text or
                        gzipped (.gz).

Column Mapping (Required):
  --chr-col CHR_COL     Name of the chromosome column.
  --pos-col POS_COL     Name of the position column.
  --ref-col REF_COL     Name of the reference/non-effect allele column.
  --alt-col ALT_COL     Name of the alternate/effect allele column.
  --af-col AF_COL       Name of the allele frequency column.
  --effect-col EFFECT_COL
                        Name of the effect size column (beta or OR).
  --se-col SE_COL       Name of the standard error or confidence interval
                        column.
  --pval-col PVAL_COL   Name of the p-value column.

Transformation and Filtering:
  --effect-type {beta,or}
                        Type of effect size reported. (default: beta)
  --se-type {se,ci}     Type of error reported (standard error or confidence
                        interval). (default: se)
  --pval-type {p,mlog10p}
                        Type of p-value reported (raw p-value or -log10(p)).
                        (default: p)
  --af-allele {alt,ref}
                        Allele described by the frequency in --af-col.
                        (default: alt)
  --flip-alleles        Flip REF/ALT alleles and invert beta/AF. (default:
                        False)
  --recalculate-se      Recalculate SE from beta and p-value. (default: False)
  --rounding-precision ROUNDING_PRECISION
                        Decimal places to round numeric outputs to. (default:
                        4)
  --filt-col FILT_COL   Name of the column to use for filtering. (default:
                        None)
  --filt-threshold FILT_THRESHOLD
                        Rows with a value in --filt-col LESS than or EQUAL to
                        this are removed. (default: None)

optional arguments:
  --delim DELIM         Delimiter used in the summary statistics file.
                        (default: tab)
  --verbose             Enable verbose output for debugging. (default: False)
```

The script performs the following transformations:

* Standardizes column names to `#CHR`, `POS`, `REF`, `ALT`, `af_alt`, `beta`, `sebeta`, `pval`
* Converts OR to log(OR) if `--effect-type or` is specified
* Converts CI to SE if `--se-type ci` is specified
* Calculates `-log10(p)` and adds it as `mlogp` column
* Applies QC filters (valid chromosomes, p-values in [0,1], non-zero betas, valid allele frequencies)
* Optionally filters based on custom column thresholds

Example usage:

```bash
python3 munge.py biomarkers.tsv.bgz \
  --chr-col chr --pos-col pos --ref-col ref --alt-col alt \
  --af-col af_EUR --effect-col beta_EUR --se-col se_EUR \
  --pval-col neglog10_pval_EUR --effect-type beta --pval-type mlog10p \
  | sort -k 1,1g -k2,2g | bgzip > munged.gz
```

### Harmonizing Summary Statistics

[harmonize.py](scripts/harmonize.py) harmonizes GWAS summary statistics to a reference (typically gnomAD) by matching variants and optionally filtering based on allele frequencies.

```bash
usage: harmonize.py [-h] [--chr_col CHR_COL] [--pos_col POS_COL]
                    [--ref_col REF_COL] [--alt_col ALT_COL]
                    [--af_col AF_COL] [--beta_col BETA_COL]
                    [--require_gnomad] [--passing_only]
                    [--gnomad_min_an GNOMAD_MIN_AN]
                    [--gnomad_max_abs_diff GNOMAD_MAX_ABS_DIFF]
                    [--pre_aligned] [--keep_best_duplicate]
                    file_in file_ref

Harmonize GWAS summary stats to reference

positional arguments:
  file_in               GWAS summary stats
  file_ref              GnomAD reference file

optional arguments:
  --chr_col CHR_COL     Chromosome column (default: #CHR)
  --pos_col POS_COL     Position column (default: POS)
  --ref_col REF_COL     Reference allele column (default: REF)
  --alt_col ALT_COL     Alternative allele column (default: ALT)
  --af_col AF_COL       Allele frequency allele column (default: af_alt)
  --beta_col BETA_COL   Beta column (default: beta)
  --require_gnomad      Filter out variants not in gnomAD
  --passing_only        Filter out non-passing variants in gnomAD
  --gnomad_min_an GNOMAD_MIN_AN
                        Minimum AN in gnomAD (default: 0)
  --gnomad_max_abs_diff GNOMAD_MAX_ABS_DIFF
                        Maximum absolute difference between variant and gnomAD
                        AF (default: 1.0)
  --pre_aligned         Input summary stats are already aligned to reference
                        (disables flipping of alleles to try find best match)
  --keep_best_duplicate
                        If duplicate variants (by chr:pos:ref:alt) keep
                        variant with smallest AF difference to reference.
                        Otherwise discard both
```

The script aligns variants across both strands and by flipping ref/alt, matches them to gnomAD reference, and can filter based on:

* Presence in gnomAD (`--require_gnomad`)
* gnomAD filter status (`--passing_only`)
* Minimum allele number (`--gnomad_min_an`)
* Maximum AF difference from gnomAD (`--gnomad_max_abs_diff`)

## Visualization

### QQ and Manhattan Plots

[qqplot.R](scripts/qqplot.R) generates quantile-quantile (QQ) plots and Manhattan plots from GWAS summary statistics.

```text
Usage: qqplot.R [options]

Options:
  -f CHARACTER, --file=CHARACTER
    dataset file name

  -o CHARACTER, --out=CHARACTER
    output file name [default=NULL]

  -c CHARACTER, --chrcol=CHARACTER
    chromosome column [default=CHR]

  -p CHARACTER, --pval_col=CHARACTER
    pvalue column [default=P]. This can be a comma separated list and
    plots will be generated for each of these

  -b CHARACTER, --bp_col=CHARACTER
    bp column [default=BP]

  -l INTEGER, --loglog_pval=INTEGER
    -log10 p-val threshold for using log-log scale in manhattan plot
    [default=10]

  -y INTEGER, --loglog_ylim=INTEGER
    -log10 p-val limit for y-axis of log-log manhattan [default=324]

  -m CHARACTER, --minrep_col=CHARACTER
    if given then chr:bp:ref:alt identifier assumed and chr and bp are
    read from there [default=NULL]
```

The script generates both standard and log-log scale Manhattan plots along with QQ plots to visualize GWAS results.

### Miami Plots

[miami.R](scripts/miami.R) generates Miami plots (mirrored Manhattan plots) to compare two sets of p-values from the same dataset.

```text
Usage: miami.R [options]

Options:
  -f CHARACTER, --file=CHARACTER
    Summary statistics file

  -o CHARACTER, --out=CHARACTER
    output file name [default=NULL]

  --chr_col=CHARACTER
    chromosome column [default=#CHR]

  --pos_col=CHARACTER
    pos column [default=POS]

  -p CHARACTER, --pval_cols=CHARACTER
    Two p-value columns, comma-separated

  --pvalue_type=CHARACTER
    Type of p-values: 'p' for p-values, 'mlogp' for -log10 p-values,
    'detect' for automatically detecting based on data [default=detect]

  --highlight
    Highlight variants not genome-wide significant in the other p-value
    column
```

Miami plots are particularly useful for comparing:

* Meta-analysis results vs. single study results
* Leave-one-out meta-analysis results
* Different ancestry groups
* Different phenotype definitions

## Meta-analysis QC

The meta-analysis workflow produces QC statistics and plots with [qc.R](scripts/qc.R).

```text
Usage: qc.R [options]

Options:
  -f CHARACTER, --file=CHARACTER
    dataset file name

  -o CHARACTER, --out=CHARACTER
    output file name [default=NULL]

  -m CHARACTER, --method=CHARACTER
    meta-analysis method [default=inv_var]

  -l, --loo
    use leave-one-out results

  --conf=CHARACTER
    meta-analysis config json

  --pval_thresh=NUMERIC
    comma separated list of p-value thresholds used to filter the data for qc

  --region=NUMERIC
    region size in megabases used when counting unique loci hits

  -c CHARACTER, --chr_col=CHARACTER
    chromosome column [default=#CHR]

  -b CHARACTER, --bp_col=CHARACTER
    bp column [default=POS]

  -r CHARACTER, --ref_col=CHARACTER
    ref column [default=REF]

  -a CHARACTER, --alt_col=CHARACTER
    alt column [default=ALT]

  --af_alt_col_suffix=CHARACTER
    af alt column suffix [default=_af_alt]

  --pheno=CHARACTER
    phenotype name [default=pheno]

  -h, --help
    Show this help message and exit
```

[qc.R](scripts/qc.R) first estimates genome-wide significant hits from input with specified region size (default: 1MB) and then calculates linear regression for those hits against using meta-analysis summary statistics.

By default, the workflow compares each study and the calculated meta-analysis summary statistics against each other and produces metrics and plots which can be used for checking whether different studies are agreeing with each other.

The QC report contains the following fields (example of a 3-way meta-analysis with studies FINNGEN, UKBB and ESTBB):

Column | Description
------ | -----------
pheno | phenotype name
FINNGEN_n_cases | FINNGEN number of cases
FINNGEN_n_controls | FINNGEN number of controls
UKBB_n_cases | UKBB number of cases
UKBB_n_controls | UKBB number of controls
ESTBB_n_cases | ESTBB number of cases
ESTBB_n_controls | ESTBB number of controls
FINNGEN_N_hits | Estimated number of genome-wide significant hits in FINNGEN
UKBB_N_hits | Estimated number of genome-wide significant hits in UKBB
ESTBB_N_hits | Estimated number of genome-wide significant hits in ESTBB
all_inv_var_meta_N_hits | Estimated number of genome-wide significant hits in the 3-way meta-analysis
leave_FINNGEN_inv_var_meta_N_hits | Estimated number of genome-wide significant hits in 2-way meta-analysis (FINNGEN left out)
leave_UKBB_inv_var_meta_N_hits | Estimated number of genome-wide significant hits in 2-way meta-analysis (UKBB left out)
leave_ESTBB_inv_var_meta_N_hits | Estimated number of genome-wide significant hits in 2-way meta-analysis (ESTBB left out)
FINNGEN_beta_vs_UKBB_beta_slope | Linear regression slope of FINNGEN vs UKBB betas
FINNGEN_beta_vs_UKBB_beta_r2 | Linear regression r2 of FINNGEN vs UKBB betas
FINNGEN_beta_vs_UKBB_beta_r2adj | Linear regression adjusted r2 of FINNGEN vs UKBB betas
FINNGEN_beta_vs_ESTBB_beta_slope | Linear regression slope of FINNGEN vs ESTBB betas
FINNGEN_beta_vs_ESTBB_beta_r2 | Linear regression r2 of FINNGEN vs ESTBB betas
FINNGEN_beta_vs_ESTBB_beta_r2adj | Linear regression adjusted r2 of FINNGEN vs ESTBB betas
FINNGEN_beta_vs_all_inv_var_meta_beta_slope | Linear regression slope of FINNGEN vs 3-way meta-analysis betas
FINNGEN_beta_vs_all_inv_var_meta_beta_r2 | Linear regression r2 of FINNGEN vs 3-way meta-analysis betas
FINNGEN_beta_vs_all_inv_var_meta_beta_r2adj | Linear regression adjusted r2 of FINNGEN vs 3-way meta-analysis betas
pct_pval_stronger_in_all_inv_var_meta_p_vs_FINNGEN | Percentage of p-values that are stronger in 3-way meta-analysis vs FINNGEN
pct_pval_stronger_in_leave_FINNGEN_inv_var_meta_p_vs_FINNGEN | Percentage of p-values that are stronger in 2-way meta-analysis (FINNGEN left-out) vs FINNGEN
pct_pval_stronger_in_leave_UKBB_inv_var_meta_p_vs_FINNGEN | Percentage of p-values that are stronger in 2-way meta-analysis (UKBB left-out) vs FINNGEN
pct_pval_stronger_in_leave_ESTBB_inv_var_meta_p_vs_FINNGEN | Percentage of p-values that are stronger in 2-way meta-analysis (FINNGEN left-out) vs FINNGEN
het_p_fdr_signif_in_meta_pct | Percentage of heterogeneity p-values (FDR corrected) that are siginificant with threshold 0.05.

## Workflow Execution (WDL)

### Meta-Analysis Workflow

[meta.wdl](wdl/meta.wdl) is the main WDL workflow for running multi-study meta-analysis at scale. It:

* Runs meta-analysis per chromosome in parallel
* Combines chromosome-level results
* Adds rsIDs from reference database
* Applies post-filtering
* Generates QC metrics and plots

### Munging Workflow

[munge.wdl](wdl/munge.wdl) automates the data preparation for meta-analysis:

* Cleans and filters summary statistics
* (Optionally) lifts over variants to GRCh38
* Harmonizes to gnomAD reference
* Applies QC filters
* Generates diagnostic plots

### Liftover Workflow

[lift.wdl](wdl/lift.wdl) lifts over summary statistics from GRCh37 to GRCh38:

* Converts summary statistics to VCF format
* Performs liftover using Picard LiftoverVcf
* Converts back to summary statistics format
* If beta/AF columns are provided, flips them as needed based on strand changes during liftover

## Utility Scripts

### Creating Meta-Analysis Configuration Files

[create_meta_confs.py](scripts/create_meta_confs.py) generates JSON configuration files for meta-analysis from a mapping file containing phenotype and study information.

```bash
usage: create_meta_confs.py [-h] [--studies STUDIES [STUDIES ...]]
                            [--min_studies MIN_STUDIES]
                            [--required_studies REQUIRED_STUDIES [REQUIRED_STUDIES ...]]
                            [--complete] [--continuous]
                            [--bucket BUCKET]
                            [--study_cols_json STUDY_COLS_JSON]
                            in_mapping_file sumstat_filelist_name
                            json_filelist_name

Generate meta-analysis configuration files from a mapping file

positional arguments:
  in_mapping_file       Input mapping file (tab-separated, can be gzipped)
  sumstat_filelist_name Output file listing summary statistic files
  json_filelist_name    Output file listing generated JSON configs

optional arguments:
  --studies STUDIES [STUDIES ...]
                        List of study names to include
  --min_studies MIN_STUDIES
                        Minimum number of studies required per phenotype
                        (default: 2)
  --required_studies REQUIRED_STUDIES [REQUIRED_STUDIES ...]
                        Studies that must be present for a phenotype
  --complete            Only include phenotypes with all studies present
  --continuous          Treat phenotypes as continuous (n_controls = 0)
  --bucket BUCKET       GCS bucket to upload files to
  --study_cols_json STUDY_COLS_JSON
                        JSON file specifying study column mappings
```

The script reads a mapping file where each row represents a phenotype and columns contain study-specific information (file paths, sample sizes, column names). It generates:

* Individual JSON configuration files for each phenotype in `jsons/` directory
* A file listing all summary statistic file paths
* A file listing all generated JSON configuration paths

### Copy Cromwell Outputs

[copy_cromwell_outputs.sh](scripts/copy_cromwell_outputs.sh) is a utility script for copying cromwell workflow output files to a GCS bucket with appropriate directory structure.

```bash
Usage: bash copy_cromwell_outputs.sh <cromwell_meta_output_prefix> <bucket>
Example: bash copy_cromwell_outputs.sh /cromwell/output/path/6223ac01.meta_analysis. gs://my-bucket
```

This script organizes meta-analysis outputs from Cromwell into a structured directory layout in the destination bucket, making it easier to find and access results.

### Remove Leave-One-Out Columns

[remove_loo_cols.sh](scripts/remove_loo_cols.sh) removes leave-one-out columns from summary statistics files.

```bash
Usage: bash remove_loo_cols.sh gs://path/to/summary/stats/file
```

This script identifies and removes columns starting with `leave_` prefix (generated by `--leave_one_out` option in meta_analysis.py). Can be useful to reduce file size when LOO results are not needed.

## Docker Image

The project includes a [Dockerfile](docker/Dockerfile) that builds a container image with all required dependencies for running the meta-analysis tools.

```dockerfile
FROM eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8.2

ADD scripts/*.py /usr/local/bin/
ADD scripts/*R /usr/local/bin/

RUN chmod a+x /usr/local/bin/*.R && chmod a+x /usr/local/bin/*.py

RUN R -e "install.packages(c('openxlsx'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

WORKDIR /
```

The Docker image:

* Based on the FinnGen bioinformatics base image (includes R, Python, common bioinformatics tools), but can be replaced with almost any image with R and Python 3 installed
* Includes all Python scripts (munge.py, meta_analysis.py, harmonize.py, create_meta_confs.py)
* Includes all R scripts (qc.R, qqplot.R, miami.R)
* Installs additional R packages (openxlsx for Excel report generation)
* Makes all scripts executable and available in PATH

To build the Docker image:

```bash
docker build -t meta-analysis:latest -f docker/Dockerfile .
```

The Docker image is used by the WDL workflows to ensure consistent execution environment across all tasks.
