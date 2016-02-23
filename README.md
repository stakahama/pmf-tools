PMF Tools
===

[![DOI](https://zenodo.org/badge/19334/stakahama/pmf-tools.svg)](https://zenodo.org/badge/latestdoi/19334/stakahama/pmf-tools)

## Introduction

These set of scripts have been used for PMF analysis of FTIR, STXM-NEXAFS, and ACSM spectra. Provided are scripts for running PMF and visualizing outputs (factor profiles, explained variation, Q-values, correlations in factor strengths/contributions). There are additional scripts provided for cluster analysis and visualization on the same data matrix using the Ward algorithm (Ward, 1963) implemented in R. Additional scripts for fixed-size moving window analysis and visualization of residual element correlations for FTIR data sets are included.

The primary requirements are matrix.dat, std_dev.dat, samples.txt, and variables.txt; helper functions to generate these matrices are provided for some measurement types.

## Instructions

---Instructions for setup of OS---

From the PMF_CD.zip provided by P. Paatero,

copy PMF_CD/PMF2key.key to C:/PMF/
copy PMF_CD/PMFx/pmf2wopt.exe to C:/PMF/
copy PMF_CD/PMFx/Exampl_2/pmf2def.ini to C:/PMF/

In Windows (on PC or as virtual machine):
- install R, add C:\Program Files\R\R-x.xx\bin to Environment Path
- requires installation of R package RJSONIO

In either Windows (guest) or Mac/Linux OS (host), create FOLDER/ which includes:
- a .json file containing variables {nFactors, FPEAK, Seed}. Example:
        {
		  "nFactors":[2,3,4,5]
		  "FPEAK":[-1.8,0.2,1.8]
		  "Seed":[1,10,100]
	    }
- input matrices (matrix.dat and std_dev.dat) and vectors (samples.txt and variables.txt)

---Instructions for running PMF---

*Must be run in Windows*

Steps:
- open DOS terminal in Windows
- type: `Rscript /path/to/pmf-run.r /path/to/userinputs.json`


---Instructions for post-processing---

*Postprocessing functions can be run in Windows or Mac/Linux.*

All scripts require the syntax `Rscript /path/to/script/scriptname.r /path/to/userinputs.json`.

Required R packages:
- fields, reshape2, gridBase

General scripts:
- pmf-qvalues.r (Q-values)
- pmf-g-correlations.r (Inter-correlations among G-factors)
- pmf-explvar.r (Explained Variations)
- pmf-explvar-plots.r (Explained Variations plots)

Measurement-specific postprocessing functions:
- ftir-factorplots.r
- ftir-residualmatrix.r
- ftir-factor-group.r
- ftir-cluster-analysis.r

---Instructions for pre-processing---

*FSMWA can be used to guide the number of factors considered for FTIR*

- ftir-fsmwa.r (Fixed-Size Moving Window Analysis)

## Other scripts

Older versions of the code read user inputs from "userinputs.r". Current version is updated to read from a .json file passed as the first command-line argument to scripts. The scripts below have not been updated.

- ftir-makematrix.r
- ftir-generate-F-file.r (for separate FTIR analysis program)
- ftir-factor-classification.r
- ftir-residualmatrix_runno.r
- ftir-residualmatrix2_runno.r
- stxm-cluster-analysis.r (a separate classification program is maintained by the author)
- stxm-factorplots.r
- acsm-extractmatrix.r
- acsm-factor-classification.r (based on .db file created from J. Jimenez's Wiki data)
- acsm-residualmatrix_runno.r
- ...

Scripts for classification analyses require separate .db files.
