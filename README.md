Epistasis for antibiotic resistance
===================================


This code runs the analyses described in Guerrero _et al_. (2019) "Proteostasis environment shapes higher-order epistasis operating on antibiotic resistance". The goal is to infer interactions acting on three phenotypes measured for nine strains of __E. coli__ (including transgenic strains and evolved mutants).

The analyses use functions from the `glm.net` and `tidyverse` packages. The plots call functions from `treemapify` and `cowplot`.

Load data
---------
The script `ms1_load_data.R` reads two data files ('GGG Whole Raw Set Clean.xlsx' and 'zero_growth_rate_reps.csv') carries out basic transformations. The `allphenos` data frame contains all the data for downstream analysis, coded in the following columns: 
'context'-- one of three proteostasis environments (WT, GroEL, lon). 
'species'-- one of three DHFR backbones, from _E.coli_, _C. muridarum_, or _L. grayi_.
'binary'-- the haplotype at the three sites of interest (coded from 000 to 111)
'pheno'-- the type of phenotype (IC50, abundance, or growth)
'W'-- the value of the phenotype
'logW'--the log10 of variable 'W'. 

Linear models
-------------
The script `ms1_linear_models.R` "sources" the previous step (loads data) and carries out the bulk of the regularized regressions. The key functions are:
`run_long_lm` -- carry out Elastic Net on the model  ~ species * context * P21L * A26T * L28R. It assumes a gaussian error distribution and intercept of zero. It does n-fold cross-validation (n = number of observations). The dependent variable can be set through parameter `yname`. If `yname=="normy"`, then it assumes the dependent variable is already normalized. Otherwise it normalizes it (recommended for regularized regressions). The function returns the coefficients at lambda_min + 1 std.error (calculated through the cross-validation).

`one_cs_lm` -- carries out Elastic Net on the model ~ P21L * A26T * L28R (meant to be fit to subsets of the data, by context-species). It always normalizes the dependent variable first (which can be set through parameter `yvar`). It also does n-fold cross-validation and returns coefficients at lambda_min + 1 std.error.

Plots and supplement
--------------------
The scripts `ms1_plots.R` and `ms1_supplemental.R` source the linear models script and produce unpolished PDF versions of all the plots in the manuscript. 