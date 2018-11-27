
<!-- README.md is generated from README.Rmd. Please edit that file -->
Scripts
=======

Below are descriptions of each script and the intended order of running the scripts.

Phenotypic Analysis
-------------------

1.  `phenotypic_analysis.R` - trial analysis, variance components, and heritability.
2.  `phenotypic_variance_analysis.R` - calculation of empirical family mean, genetic variance, and superior progeny mean.

Predictions
-----------

1.  `popvar_predictions.R` - generate predictions using the package `PopVar`.
2.  `run_predictions.sh` - Shell script to run the `popvar_predictions.R` script. Formatted for execution on the Minnesota Supercomputing Institute. You may need to modify this script for your own system.
3.  `prediction_analysis.R` - analysis of the predictions from `PopVar`.
4.  `prediction_accuracy_analysis.R` - calculation and analysis of prediction accuracy and bias.
