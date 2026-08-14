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
usage: meta_analysis.py [-h] [--not_quiet] [--leave_one_out]
                        [--leave_one_cohort_out] [--leave_one_population_out]
                        [--is_het_test] [--pairwise_with_first] [--sep SEP]
                        [--chrom CHROM] [--flip_indels]
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
  --leave_one_cohort_out
                        Do leave-one-cohort-out meta-analysis using the optional
                        "cohort" config field. Studies without a cohort are
                        always retained.
  --leave_one_population_out
                        Do leave-one-population-out meta-analysis using the
                        optional "population" config field. Studies without a
                        population are always retained.
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
            "se": "SE", # effect standard error. This column is optional if not using inv_var method. Otherwise additional p-value will be added using this as a weight for z-score.
            "extra_cols": ["AF.Cases","AF.Controls"], # optional. Additional columns to carry over to the output, prefixed with the study name
            "cohort": "FINNGEN", # optional. Cohort label used by --leave_one_cohort_out
            "population": "EUR" # optional. Population label used by --leave_one_population_out
```

With `--leave_one_out` the output gets a `leave_<STUDY>_*` column block per study.
`--leave_one_cohort_out` and `--leave_one_population_out` add one block per distinct
value of the `cohort` / `population` config field, named `leave_cohort_<VALUE>_*` and
`leave_population_<VALUE>_*`. Studies that do not define the field are always kept in
the meta-analysis, and the flag is skipped with a message to stderr if no study defines it.

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
                        6)
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
* Converts CI to SE if `--se-type ci` is specified. Note that CI is assumed to be 95% and formatted as `[lower],[upper]`
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

  -w, --weighted
    do inverse variance weighted linear regression

  --keep_hla
    do not remove HLA region variants from QC

  --hla_region=NUMERIC
    HLA region boundaries

  -h, --help
    Show this help message and exit
```

[qc.R](scripts/qc.R) first estimates genome-wide significant hits from input with specified region size (default: 1MB) and then calculates linear regression for those hits against using meta-analysis summary statistics. With `--weighted` the regression is inverse variance weighted using the standard errors; the [meta.wdl](wdl/meta.wdl) workflow always runs it this way.

Variants in the HLA region (chromosome 6, `--hla_region` boundaries, default 20-40Mb) are excluded from the QC by default because their long-range LD inflates hit counts and regression fits. Use `--keep_hla` to include them.

Besides the report table and QC plots, [qc.R](scripts/qc.R) writes `<out>.<pval_thresh>.forest_plots.pdf` with one forest plot per unique significant hit, showing beta and 95% CI with the p-value annotated for each study, the meta-analysis (labelled `META`, drawn at the bottom) and, with `--loo`, each per-study leave-one-out meta-analysis.

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

The `leave_*` fields cover per-study leave-one-out only. The `leave_cohort_<VALUE>_*` and
`leave_population_<VALUE>_*` columns produced by `--leave_one_cohort_out` /
`--leave_one_population_out` are written to the summary statistics but are not included in
the QC report or plots. `--loo` is ignored if the file has no per-study `leave_*_meta_p`
columns, or if fewer than 3 studies are present.

## Workflow Execution (WDL)

### Meta-Analysis Workflow

[meta.wdl](wdl/meta.wdl) is the main WDL workflow for running multi-study meta-analysis at scale. It:

* Runs meta-analysis per chromosome in parallel
* Combines chromosome-level results
* Adds rsIDs from reference database
* Applies post-filtering
* Generates QC metrics and plots: the QC report and plots from [qc.R](scripts/qc.R) (run with `--weighted`), per-hit forest plots, QQ/Manhattan plots from [qqplot.R](scripts/qqplot.R) and Miami plots from [miami.R](scripts/miami.R)

Which meta-analysis options are used is controlled by `meta_analysis.run_range.opts` in
[meta.json](wdl/meta.json), e.g. `--is_het_test --leave_one_out --leave_one_cohort_out
--leave_one_population_out`. `--loo` is passed to [qc.R](scripts/qc.R) automatically when
`meta_analysis.plots.pvals_to_plot` contains a `leave_` column.

### Munging Workflow

[munge.wdl](wdl/munge.wdl) automates the data preparation for meta-analysis:

* Cleans and filters summary statistics
* (Optionally) lifts over variants to GRCh38
* Harmonizes to gnomAD reference
* Applies QC filters
* Generates diagnostic plots

### Liftover Workflow

[lift.wdl](wdl/lift.wdl) lifts over summary statistics from GRCh37 to GRCh38:

* Validates that all configured input columns exist in the summary statistics header and fails early listing the missing and available columns
* Converts summary statistics to VCF format
* Performs liftover using Picard LiftoverVcf
* Converts back to summary statistics format
* If beta/AF columns are provided, flips them as needed based on strand changes during liftover

The `af_col` and `beta_col` inputs are optional and each accepts a comma-separated list of
column names, so files carrying per-population or per-model effect and frequency columns
can be lifted in one pass:

```json
"liftover.af_col": "af_controls_EUR,af_controls_AFR",
"liftover.beta_col": "beta_EUR,beta_AFR"
```

## Utility Scripts

### Creating Meta-Analysis Configuration Files

[create_meta_confs.py](scripts/create_meta_confs.py) generates JSON configuration files for meta-analysis from a mapping file containing phenotype and study information.

```bash
usage: create_meta_confs.py [-h] [--studies STUDIES [STUDIES ...]]
                            [--bucket BUCKET]
                            [--sumstat_filelist_name SUMSTAT_FILELIST_NAME]
                            [--json_filelist_name JSON_FILELIST_NAME]
                            [--phenotype_col PHENOTYPE_COL]
                            [--link_col_suffix LINK_COL_SUFFIX]
                            [--n_cases_col_suffix N_CASES_COL_SUFFIX]
                            [--n_controls_col_suffix N_CONTROLS_COL_SUFFIX]
                            [--study_meta_json STUDY_META_JSON] [--continuous]
                            [--min_studies MIN_STUDIES] [--complete]
                            [--required_studies REQUIRED_STUDIES [REQUIRED_STUDIES ...]]
                            in_mapping_file

Create a json config file for meta-analysis

positional arguments:
  in_mapping_file       A tab-delimited text file with phenotype mapping
                        information between studies

optional arguments:
  -h, --help            show this help message and exit
  --studies STUDIES [STUDIES ...]
                        List of studies to include in the meta-analysis
                        configs. (Default: all studies in the mapping file)
  --bucket BUCKET       GCS bucket path
  --sumstat_filelist_name SUMSTAT_FILELIST_NAME
                        Name of file listing the sumstats (default:
                        sumstat_files.txt)
  --json_filelist_name JSON_FILELIST_NAME
                        Name of file listing the jsons (default:
                        conf_jsons.txt)
  --phenotype_col PHENOTYPE_COL
                        Column name containing the phenotype name. (Default:
                        'phenotype')
  --link_col_suffix LINK_COL_SUFFIX
                        Suffix for the column name containing the sumstat URI.
                        (Default: '_link')
  --n_cases_col_suffix N_CASES_COL_SUFFIX
                        Suffix for the column name containing the number of
                        cases. (Default: '_n_cases')
  --n_controls_col_suffix N_CONTROLS_COL_SUFFIX
                        Suffix for the column name containing the number of
                        controls. (Default: '_n_controls')
  --study_meta_json STUDY_META_JSON
                        JSON file (keyed by study name) with per-study sumstat
                        column names and optionally 'cohort'/'population'
                        grouping metadata passed through to the configs
  --continuous          Phenotypes are continuous (assume number of controls
                        equals 0)
  --min_studies MIN_STUDIES
                        Minimum number of studies required for a phenotype to
                        be included in the meta-analysis. (Default: 2)
  --complete            Only include phenotypes with all studies present
  --required_studies REQUIRED_STUDIES [REQUIRED_STUDIES ...]
                        List of studies that are required for a phenotype to
                        be included in the meta-analysis
```

The script reads a mapping file where each row represents a phenotype and columns contain study-specific information (file paths, sample sizes, column names). It generates:

* Individual JSON configuration files for each phenotype in `jsons/` directory
* A file listing all summary statistic file paths (`--sumstat_filelist_name`, default `sumstat_files.txt`)
* A file listing all generated JSON configuration paths (`--json_filelist_name`, default `conf_jsons.txt`)

`--bucket` copies the two list files and the whole `jsons/` directory to the given GCS path.
Note that the JSON list is written with the bucket prefix, so it is only populated when
`--bucket` is given.

The `--study_meta_json` file is keyed by study name. Besides the sumstat column
names (`chr`, `pos`, `ref`, `alt`, `effect`, `effect_type`, `pval`, optional
`se`/`extra_cols`), each study entry may include optional `"cohort"` and
`"population"` fields. When present they are passed through verbatim into the
generated study configs, enabling `--leave_one_cohort_out` /
`--leave_one_population_out` in [meta_analysis.py](scripts/meta_analysis.py).

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

This script identifies and removes columns starting with `leave_` prefix (generated by the `--leave_one_out`, `--leave_one_cohort_out` and `--leave_one_population_out` options in meta_analysis.py). Can be useful to reduce file size when LOO results are not needed.

## Docker Image

The project includes a [Dockerfile](docker/Dockerfile) that builds a container image with all required dependencies for running the meta-analysis tools.

The Docker image:

* Based on the FinnGen bioinformatics base image (includes R, Python, common bioinformatics tools), but can be replaced with almost any image with R and Python 3 installed
* Includes all Python scripts (munge.py, meta_analysis.py, harmonize.py, create_meta_confs.py)
* Includes all R scripts (qc.R, qqplot.R, miami.R)
* Installs the Python dependencies listed in [docker/requirements.txt](docker/requirements.txt) (numpy, scipy) with pip
* Installs the R dependencies listed in [docker/install_packages.R](docker/install_packages.R) (data.table, ggplot2, ggpubr, optparse, qqman, R.utils, rjson, stringi, openxlsx). Already present packages are skipped and the build fails if any package is still missing afterwards
* Makes all scripts executable and available in PATH

To build the Docker image:

```bash
docker build -t meta-analysis:latest -f docker/Dockerfile .
```

The Docker image is used by the WDL workflows to ensure consistent execution environment across all tasks.
