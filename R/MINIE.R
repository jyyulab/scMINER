## ---------------------------
##
## Script name: MINIE (Mutual Information-based Network Inference Engine)
##
## Purpose of script: to reconstruct cell-type-specific GRNs for driver activity inference and target network rewiring analysis
##
## Author: The scMINER software is developed and maintained by the Yu Laboratory @ St. Jude
##
## Date Created: 2023-03-24
##
## ---------------------------
##
## Notes:
##    (1) MINIE takes inputs of a gene expression profile and cell cluster labels.
##    (2) Then MINIE invokes SJARACNe to reconstruct cell-type-specific transcriptional factor and signaling networks.
##    (3) Finally, with the predicted targets of a driver for a cell cluster, MINIE calculates the driver activity by performing a column-wise normalization to ensure each cell is on a similar expression level, followed by averaging the expression of the driver's target genes.
##
## ---------------------------


