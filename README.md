# META ANALYSIS tools

Tools for doing x way meta-analysis

## Variant matching across studies

Variants are matched using chr pos ref and alt. For this reason all 37 build results need to first be lifted over to 38.
<https://github.com/FINNGEN/commons/tree/master/liftover> can be used to liftover results first if needed.

IMPORTANT: Studies need to be ordered by chr (1-22, x,y,mt) and position. Chromosome can be indicated with numbers 1-25 or 1-22, X, Y, MT, with or without 'chr' prefix and they will be internally coded to numerical values.

## Running single trait meta-analysis

[meta_analysis.py](scripts/meta_analysis.py) is the main script for running meta-analysis for a single trait. Meta-analyses to be performed are specified with json configuration file. Example configuration in (data/conf.json). Script will try to align using both strands as well as by flipping ref vs. alt.

```bash
usage: meta_analysis.py [-h] [--not_quiet] [--leave_one_out] [--is_het_test]
                        [--pairwise_with_first] [--dont_allow_space]
                        [--chrom CHROM]
                        config_file path_to_res methods

Run x-way meta-analysis

positional arguments:
  config_file           Configuration file
  path_to_res           Result file
  methods               List of meta-analysis methods to compute separated by
                        commas.Allowed values [n,inv_var,variance]

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
```

The configuration file should be a json file with these elements:

```json
            "name":"FINNGEN",
            "file":"/Users/mitja/projects/finngen/META_ANALYSIS/I9_AF.gz",
            "n_cases": 6570 , # number of cases. Used only if sample size weighted meta-analysis is used
            "n_controls": 48378, # number of controls. Used only if sample size weighted meta-analysis is used
            "chr":"CHR", #chromosome column name in the file
            "pos":"POS", #position column name in the file
            "ref":"Allele1", #reference allele column name in the file
            "alt":"Allele2", #alternate allele column name in the file
            "effect":"BETA", #effect size column name in the file
            "effect_type":"beta", #is the effect column beta or or. In case of OR the value will be log transformed to beta.
            "pval":"p.value" # effect size column name in the file
            "se":"SE" <- this parameter is optional. If given for compared studies additional p-value will be added using this as a weight for z-score.
```

[meta_analysis.py](scripts/meta_analysis.py) supports 3 different meta-analysis methods

* `N`: purely weight by sample size and use z-score from p-value,
* `variance`: weight z-score from p-value by variance,
* `inv_var`: regular inverse variance weighted betas meta-analysis.

`inv_var` is recommeneded if betas and variances are comparable. In case of combining data frmo different models (e.g.) linear vs. logistic you should use sample size weighted meta.
