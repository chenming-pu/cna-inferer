import numpy as np
import scipy.sparse as sp
import pandas as pd

def sliding_window_aggregate(adata, window_size=100):
    """
    Efficient sliding window aggregation for gene expression data.

    This function bins genes along the genome (assumed to be sorted by chromosome and position)
    and computes average expression per bin for each cell using sparse matrix operations.

    Parameters:
    ----------
    adata : AnnData
        AnnData object with raw gene expression in .X and genomic information in .var
        (must contain 'chromosome' and 'start' columns).
    window_size : int
        Number of genes per bin.

    Returns:
    -------
    adata : AnnData
        The same AnnData object with the following additions:
        - .obsm["X_binned"] : [n_cells × n_bins] matrix of binned expression values
        - .uns["bin_info"]  : DataFrame with bin-level metadata (chromosome, start position)

    Example:
    --------
    >>> adata = sliding_window_aggregate(adata, window_size=100)
    >>> Xb = adata.obsm["X_binned"]
    >>> bin_info = adata.uns["bin_info"]
    """

    # 1) Number of genes and bins
    n_genes = adata.n_vars
    n_bins = int(np.ceil(n_genes / window_size))

    # 2) Assign each gene to a bin index
    bins = np.repeat(np.arange(n_bins), window_size)[:n_genes]

    # 3) Construct a sparse gene × bin matrix
    bin_matrix = sp.csr_matrix(
        (np.ones(n_genes), (np.arange(n_genes), bins)),
        shape=(n_genes, n_bins)
    )
    bin_matrix = bin_matrix.multiply(1 / window_size)  # Normalize for average

    # 4) Multiply: [cells × genes] @ [genes × bins] → [cells × bins]
    Xbinned = adata.X @ bin_matrix

    # Convert to dense if needed and store in .obsm
    if sp.issparse(Xbinned):
        adata.obsm["X_binned"] = Xbinned.toarray()
    else:
        adata.obsm["X_binned"] = Xbinned

    # 5) Construct bin-level metadata from gene annotations
    bin_info = pd.DataFrame({
        "chromosome": adata.var["chromosome"].values,
        "start_gene": adata.var["start"].values
    }, index=bins).groupby(level=0).first()
    
    adata.uns["bin_info"] = bin_info

    return adata
