# cna_inferer/segmentation.py

import pandas as pd
import numpy as np
import ruptures as rpt
from joblib import Parallel, delayed
from .aggregation import sliding_window_aggregate
import multiprocessing

def segment_binseg(zsig: np.ndarray,
                   n_bkps: int = 8,
                   model: str = "l2",
                   z_thresh: float = 1.5,
                   min_length: int = 3):
    """
    Apply binary segmentation (Binseg) to a z-score signal and identify gain/loss events.

    Parameters:
    -----------
    zsig : np.ndarray
        1D array of z-score values for one cell.
    n_bkps : int
        Maximum number of breakpoints.
    model : str
        Cost model to use in segmentation ('l2', 'rbf', etc.).
    z_thresh : float
        Absolute z-score threshold to call a gain or loss.
    min_length : int
        Minimum segment length to consider as a CNA event.

    Returns:
    --------
    events : list of dicts
        Each dict has keys: start, end, type ('gain' or 'loss').
    """
    algo = rpt.Binseg(model=model).fit(zsig)
    bkps = algo.predict(n_bkps=n_bkps)
    events = []
    start = 0
    for end in bkps:
        seg = zsig[start:end]
        mean_z = seg.mean()
        length = end - start
        if mean_z >= z_thresh and length >= min_length:
            events.append({"start": start, "end": end - 1, "type": "gain"})
        elif mean_z <= -z_thresh and length >= min_length:
            events.append({"start": start, "end": end - 1, "type": "loss"})
        start = end
    return events


def call_cnas(adata,
              window_size: int = 100,
              z_thresh: float = 1.5,
              min_bins: int = 3,
              n_bkps: int = 8,
              model: str = "l2",
              n_jobs: int = 1):
    """
    Parallelized unsupervised CNA calling using z-score segmentation.

    Parameters:
    -----------
    adata : AnnData
        Input single-cell expression dataset (with .var containing 'chromosome' and 'start').
    window_size : int
        Number of genes per bin in sliding window aggregation.
    z_thresh : float
        Z-score threshold for calling CNAs.
    min_bins : int
        Minimum number of bins per CNA segment.
    n_bkps : int
        Maximum number of breakpoints per cell.
    model : str
        Cost model to use in segmentation ('l2', 'rbf', etc.).
    n_jobs : int
        Number of parallel jobs. If None, auto-detected from CPU cores.

    Returns:
    --------
    adata : AnnData
        Modified AnnData with:
        - .obs["cna_calls"]: list of CNA events per cell
        - .uns["cna_events"]: combined CNA calls across all cells in a DataFrame
    """
    # Auto-determine number of cores if not set
    if n_jobs is None:
        n_jobs = max(multiprocessing.cpu_count() - 1, 1)
        print(f"🔧 Auto-selected n_jobs = {n_jobs}")

    # 1) Aggregate if not already done
    if "X_binned" not in adata.obsm:
        adata = sliding_window_aggregate(adata, window_size)

    # 2) Compute z-score matrix
    Xb = adata.obsm["X_binned"]
    Xz = (Xb - Xb.mean(axis=1, keepdims=True)) / (Xb.std(axis=1, keepdims=True) + 1e-6)

    # 3) Segment each cell in parallel
    def _process_cell(idx):
        if idx % 500 == 0:
            print(f"✅ Processed {idx} cells...")
        cell = adata.obs_names[idx]
        segs = segment_binseg(
            Xz[idx],
            n_bkps=n_bkps,
            model=model,
            z_thresh=z_thresh,
            min_length=min_bins
        )
        return cell, segs

    results = Parallel(n_jobs=n_jobs)(
        delayed(_process_cell)(i) for i in range(adata.n_obs)
    )

    # 4) Collect results
    all_events = []
    adata.obs["cna_calls"] = [[] for _ in adata.obs_names]
    for cell, segs in results:
        adata.obs.at[cell, "cna_calls"] = segs
        for ev in segs:
            all_events.append({
                "cell":      cell,
                "start_bin": ev["start"],
                "end_bin":   ev["end"],
                "type":      ev["type"]
            })

    # 5) Create summary DataFrame
    if all_events:
        df = pd.DataFrame(all_events)
        bin_info = adata.uns["bin_info"]
        df["chromosome"] = df["start_bin"].map(bin_info["chromosome"])
        df["start_gene"] = df["start_bin"].map(bin_info["start_gene"])
    else:
        df = pd.DataFrame(columns=[
            "cell", "start_bin", "end_bin", "type", "chromosome", "start_gene"
        ])

    adata.uns["cna_events"] = df
    return adata
