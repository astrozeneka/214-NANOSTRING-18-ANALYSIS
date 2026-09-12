
import pandas as pd

if __name__ == '__main__':
    rnaseq_log2tpm_df = pd.read_csv("../211-ONCOGENES/tpm_finished_data/RNA-seq -18-samples-tpm.csv", index_col=0)
    nanostring_raw_df = pd.read_csv("data/Nanostring rawdata for 18 patients.csv", index_col=0)

    rnaseq_genes = set(rnaseq_log2tpm_df.index)
    nanostring_genes = set(nanostring_raw_df.index)

    print(f"Number of genes in RNA-seq data: {len(rnaseq_genes)}")
    print(f"Number of genes in Nanostring data: {len(nanostring_genes)}")
    print(f"Number of genes in the intersection: {len(rnaseq_genes & nanostring_genes)}")
    print(f"Number of genes in the union: {len(rnaseq_genes | nanostring_genes)}")
    print(f"Jaccard similarity: {len(rnaseq_genes & nanostring_genes) / len(rnaseq_genes | nanostring_genes):.4f}")
