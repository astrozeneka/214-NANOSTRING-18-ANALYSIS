## Molecular subtype identification from NanoString nCounter data

This pipeline is intended to predict the molecular subtype of colorectal cancer samples using a targeted ~770-gene panel from **NanoString nCounter data**.

Raw probe counts are first cleaned of QC controls and converted to length-normalized TPM (then log2) using transcript lengths from Ensembl biomaRt. 

The log2-TPM matrix is turned into MSigDBv7 functional spectra and passed to a DeepCC model trained on CRC TCGA, which assigns each sample a CMS label (CMS1–CMS4, or NA below the 0.5 probability cutoff) together with its class probabilities.
