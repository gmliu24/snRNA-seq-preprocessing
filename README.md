# snRNA-seq-preprocessing
Run QC on drprnull42d and w111842d cells separately 
QC using median absolute deviations (MAD) method, Ambient RNA correction using SoupX package, Normalization using log10 transformation
Outputs AnnData .h5ad files for both samples separately
Code adpated from https://biostatsquid.com/scrnaseq-preprocessing-workflow-seurat/#step3

Didn't remove doublets
