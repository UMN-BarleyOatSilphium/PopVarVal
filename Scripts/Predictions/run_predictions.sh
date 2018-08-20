#!/bin/bash

#PBS -l walltime=80:00:00,mem=62gb,nodes=1:ppn=24
#PBS -N PVV_popvar_predictions
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/PopVarVal

module load R/3.4.0

# Prediction of all BP families
Rscript Predictions/PVV_popvar_predictions.R
