import scanpy as sc
import pandas as pd
from cna_inferer.annotation import annotate_genes_from_gtf
from cna_inferer.aggregation import sliding_window_aggregate
from cna_inferer.segmentation import call_cnas

def process_and_call_cnas(
    h5_file,
    gtf_file="Homo_sapiens.GRCh38.104.gtf",
    window_size=150,
    z_thresh=1.5,
    min_bins=3,
    n_bkps=4,
    model="l2",
    n_jobs=1
):
    """
    Load scRNA-seq data, annotate genes, preprocess, bin, and call CNAs.

    Parameters
    ----------
    h5_file : str
        Path to .h5 or .h5ad file.
    gtf_file : str
        Path to GTF annotation file.
    window_size : int
        Number of genes per bin.
    z_thresh : float
        Z-score threshold for calling gain/loss.
    min_bins : int
        Minimum length (in bins) of CNA to retain.
    n_bkps : int
        Number of breakpoints for segmentation.
    model : str
        Cost model to use in segmentation.
    n_jobs : int
        Number of parallel jobs.

    Returns
    -------
    adata : AnnData
        Annotated object with CNA calls.
    """
    # Load data depending on file type
    if h5_file.endswith(".h5ad"):
        adata = sc.read_h5ad(h5_file)
    elif h5_file.endswith(".h5"):
        adata = sc.read_10x_h5(h5_file)
    else:
        raise ValueError(f"Unsupported file type: {h5_file}")

    adata.var_names_make_unique()

    # Annotate genes with chromosome and start position
    adata = annotate_genes_from_gtf(adata, gtf_file)
    adata.var["chromosome"] = adata.var["chromosome"].astype(str)

    # Quality control: filter cells and genes
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    sc.pp.filter_cells(adata, min_counts=500)
    sc.pp.filter_cells(adata, max_counts=50000)
    adata = adata[adata.obs['pct_counts_mt'] < 5, :]

    sc.pp.filter_genes(adata, min_cells=3)

    # Normalize and log transform
    adata.raw = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Sort genes by chromosome and position
    chrom_order = [str(i) for i in range(1, 23)] + ["X", "Y"]
    adata.var["chromosome"] = pd.Categorical(adata.var["chromosome"], categories=chrom_order, ordered=True)
    adata.var = adata.var.sort_values(["chromosome", "start"])
    adata._inplace_subset_var(adata.var.index)

    # Sliding window aggregation
    adata = sliding_window_aggregate(adata, window_size=window_size)

    # Call CNAs
    adata = call_cnas(
        adata,
        window_size=window_size,
        z_thresh=z_thresh,
        min_bins=min_bins,
        n_bkps=n_bkps,
        model=model,
        n_jobs=n_jobs
    )

    # Global z-score (optional: for plotting or scoring)
    Xbinned = adata.obsm["X_binned"]
    z = (Xbinned - Xbinned.mean(0)) / Xbinned.std(0)
    zs = z.ravel()  # not returned, but available for downstream usage

    # Record number of CNA events per cell
    event_counts = [len(segments) for segments in adata.obs["cna_calls"]]
    adata.obs["n_cna_events"] = event_counts

    print(f"📊 Total CNA events in {h5_file}:", sum(event_counts))
    return adata
