# META ANALYSIS tools

Tools for doing x way meta-analysis

## Variant matching across studies

Variants are matched using chr pos ref and alt. For this reason all 37 build results need to first be lifted over to 38 to match FinnGen. [lift.wdl](wdl/lift.wdl) can be used to liftover results first if needed.

IMPORTANT: Studies need to be ordered by chr (1-22, x,y,mt) and position. Chromosome can be indicated with numbers 1-25 or 1-22, X, Y, MT, with or without 'chr' prefix and they will be internally coded to numerical values.

## Running single trait meta-analysis

[meta_analysis.py](scripts/meta_analysis.py) is the main script for running meta-analysis for a single trait. Meta-analyses to be performed are specified with json configuration file. Example configuration with three studies in (data/conf.json). Script will try to align using both strands as well as by flipping ref vs. alt.

```bash
usage: meta_analysis.py [-h] [--not_quiet] [--leave_one_out] [--is_het_test]
                        [--pairwise_with_first] [--dont_allow_space]
                        [--chrom CHROM]
                        config_file path_to_res methods

Run x-way meta-analysis

positional arguments:
  config_file           Configuration file
  path_to_res           Result file
  methods               List of meta-analysis methods to compute separated by commas. Allowed values [n,inv_var,variance]

optional arguments:
  -h, --help            show this help message and exit
  --not_quiet           Print matching variants to stdout
  --leave_one_out       Do leave-one-out meta-analysis
  --is_het_test         Do heterogeneity tests based on Cochrans Q and output
                        het_p
  --pairwise_with_first
                        Do pairwise meta-analysis with the first given study
  --dont_allow_space    Do not allow space as field delimiter
  --chrom CHROM         Restrict to given chromosome
  --flip_indels         Try variant aligning by flipping indels also. By default indels are not flipped
```

The configuration file should be a json file with two or more studies with these elements in each:

```json
            "name":"FINNGEN",
            "file":"/Users/mitja/projects/finngen/META_ANALYSIS/I9_AF.gz",
            "n_cases": 6570 , # number of cases. Used only if sample size weighted meta-analysis is used
            "n_controls": 48378, # number of controls. Used only if sample size weighted meta-analysis is used
            "chr":"CHR", # chromosome column name in the file
            "pos":"POS", # position column name in the file
            "ref":"Allele1", # reference allele column name in the file
            "alt":"Allele2", # alternate allele column name in the file
            "effect":"BETA", # effect size column name in the file
            "effect_type":"beta", # effect type. Allowed values: beta, OR. In case of OR the value will be log transformed to beta.
            "pval":"p.value" # effect size column name in the file
            "se":"SE" # effect standard error. This column is optional if not using inv_var method. Otherwise additional p-value will be added using this as a weight for z-score.
```

[meta_analysis.py](scripts/meta_analysis.py) supports 3 different meta-analysis methods

* `N`: purely weight by sample size and use z-score from p-value,
* `variance`: weight z-score from p-value by variance,
* `inv_var`: regular inverse variance weighted betas meta-analysis.

`inv_var` is recommended if betas and variances are comparable. In case of combining data from different models (e.g.) linear vs. logistic you should use sample size weighted meta.

## Meta-analysis QC

The meta-analysis workflow produces QC statistics and plots with [qc.R](scripts/qc.R).

```
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
-------| -----------
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
