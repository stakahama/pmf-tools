PMF Tools
===

[![DOI](https://zenodo.org/badge/19334/stakahama/pmf-tools.svg)](https://zenodo.org/badge/latestdoi/19334/stakahama/pmf-tools)

These set of scripts have been used for PMF analysis of FTIR, STXM-NEXAFS, and ACSM spectra. Provided are scripts for running PMF and visualizing outputs (factor profiles, explained variation, Q-values, correlations in factor strengths/contributions). There are additional scripts provided for cluster analysis and visualization on the same data matrix using the Ward algorithm (Ward, 1963) implemented in R. Additional scripts for fixed-size moving window analysis and visualization of residual element correlations for FTIR data sets are included.

The primary requirements are matrix.dat, std_dev.dat, samples.txt, and variables.txt; helper functions to generate these matrices are provided for some measurement types.

---Instructions---

From the PMF_CD.zip provided by P. Paatero,

copy PMF_CD/PMF2key.key to C:/PMF/
copy PMF_CD/PMFx/pmf2wopt.exe to C:/PMF/
copy PMF_CD/PMFx/Exampl_2/pmf2def.ini to C:/PMF/

In Windows (on PC or as virtual machine):
- install R, add C:\Program Files\R\R-x.xx\bin to Environment Path

In either Windows (guest) or Mac/Linux OS (host), create FOLDER/ which includes:
- userinputs.json containing variables: {nFactors, FPEAK, Seeds}
- input matrices (matrix.dat and std_dev.dat) and vectors (samples.txt and variables.txt)

In Windows:
- open DOS terminal in Windows
- type: Rscript --no-save pmf-run.r

Postprocessing functions can be run in Windows or Mac/Linux (can be set up to be called from batch.r):
- pmf-qvalues.r (Q-values)
- pmf-g-correlations.r (Inter-correlations among G-factors)
- pmf-explvar.r (Explained Variations)
- pmf-explvar-plots.r (Explained Variations plots)

Measurement-specific postprocessing functions:
- ftir-factorplots.r
- ftir-residualmatrix.r
- ftir-factor-classification.r (please request .db file)
- ftir-generate-F-file.r (for separate FTIR analysis program)
- stxm-factorplots.r

Measurement-specific preprocessing functions:
- ftir-fsmwa.r (Fixed-Size Moving Window Analysis)
- ftir-makematrix.r
- acsm-extractmatrix.r

Cluster analysis (exports clusters and figures; requires gridBase package in R):
- ftir-cluster-analysis.r
- stxm-cluster-analysis.r

Classification analysis
- ftir-factor-classification.r (requires .db file provided separately)
- acsm-factor-classification.r (to add; based on .db file created from J. Jimenez's Wiki data)
