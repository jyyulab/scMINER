## ---------------------------
##
## Script name: scMINER v1.0
##
## Purpose of script: designed for preprocessing, QC, clustering,
##    and hidden driver analysis of single-cell RNA-seq data
##
## Author: The scMINER software is developed and maintained by the Yu Laboratory @ St. Jude
##
## Date Created: 2023-03-24
##
## ---------------------------
##
## Notes: scMINER enables mutual information-based cell clustering,
##    cell-type-specific gene regulatory network (GRN) reverse engineering and
##    protein activity inference, to identify hidden transcriptional factors (TFs)
##    and signaling factors (SIGs) driving cellular lineage differentiation and
##    tissue specific specification.
##
##    scMINER software consists of three components:
##
##      (1) MICA: Mutual Information based Clustering analysis (https://github.com/jyyulab/MICA).
##      (2) SJARACNe:a scalable solution of ARACNe that improves the computational performance to reconstruct the regulatory network (https://github.com/jyyulab/SJARACNe).
##      (3) MINIE: Mutual Information-based Network Inference Engine (MINIE.R).
##
##    scMINER provides supporting functions to prepare datasets for running MICA and SJARACNe and import results into R.
##
## ---------------------------



