# cna_inferer/annotation.py

import pandas as pd
from gtfparse import read_gtf

def annotate_genes_from_gtf(adata, gtf_path):
    """
    Annotate genes in an AnnData object using coordinates from a GTF file.

    Parameters
    ----------
    adata : AnnData
        The input AnnData object. Must contain either 'gene_ids' column or use .var_names as gene symbols.
    gtf_path : str
        Path to the GTF file.

    Returns
    -------
    adata : AnnData
        Annotated AnnData object with added columns: chromosome, start, end.
    """
    # 1) Read GTF (as a Polars or Pandas DataFrame)
    gtf = read_gtf(gtf_path)

    # If using Polars, convert to Pandas
    try:
        import polars as pl
        if isinstance(gtf, pl.DataFrame):
            gtf = gtf.to_pandas()
    except ImportError:
        pass

    # 2) Filter for 'gene' features and keep relevant columns
    genes = gtf[gtf['feature'] == 'gene'][['gene_id', 'gene_name', 'seqname', 'start', 'end']]
    genes = genes.rename(columns={'seqname': 'chromosome'}).set_index('gene_id')

    # 3) Join with adata.var by gene_id or gene_name
    if 'gene_ids' in adata.var.columns:
        adata.var = adata.var.join(genes, on='gene_ids', how='left')
    else:
        adata.var['gene_name'] = adata.var_names
        adata.var = adata.var.join(genes, on='gene_name', how='left', rsuffix='_gtf')

    return adata
