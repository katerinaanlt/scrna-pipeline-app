"""Reading, merging and exporting datasets (.h5ad and .csv <-> AnnData)."""
from __future__ import annotations

import io
import os
import tempfile
from typing import Dict, List

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def read_h5ad_bytes(data: bytes) -> ad.AnnData:
    """Read an .h5ad upload (Streamlit gives us bytes, anndata wants a path)."""
    fd, path = tempfile.mkstemp(suffix=".h5ad")
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
        adata = ad.read_h5ad(path)
    finally:
        os.unlink(path)
    return adata


def read_csv_bytes(data: bytes, transpose: bool = False) -> ad.AnnData:
    """Convert a numeric CSV/TSV expression matrix to AnnData.

    Expected layout: first column = row names, header = column names, rows = cells and
    columns = genes. Tick *transpose* if your file has genes as rows instead.
    """
    try:
        frame = pd.read_csv(io.BytesIO(data), index_col=0, sep=None, engine="python")
    except (pd.errors.ParserError, pd.errors.EmptyDataError, UnicodeDecodeError) as exc:
        raise ValueError(f"Could not parse the file as CSV/TSV: {exc}") from exc
    if frame.empty:
        raise ValueError("The file contains no data.")
    numeric = frame.select_dtypes(include="number")
    if numeric.shape[1] != frame.shape[1]:
        raise ValueError("All columns must be numeric (apart from the first, which holds row names).")
    if transpose:
        frame = frame.T
    adata = ad.AnnData(
        X=sparse.csr_matrix(frame.to_numpy(dtype=np.float32)),
        obs=pd.DataFrame(index=frame.index.astype(str)),
        var=pd.DataFrame(index=frame.columns.astype(str)),
    )
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    return adata


def load_file(filename: str, data: bytes, transpose: bool = False) -> ad.AnnData:
    """Dispatch on the file extension."""
    lower = filename.lower()
    if lower.endswith(".h5ad"):
        return read_h5ad_bytes(data)
    if lower.endswith((".csv", ".tsv", ".txt")):
        return read_csv_bytes(data, transpose=transpose)
    raise ValueError("Unsupported file type. Please upload a .h5ad or .csv file.")


def merge_samples(samples: Dict[str, ad.AnnData], join: str = "inner") -> ad.AnnData:
    """Concatenate samples into one AnnData with a ``sample`` column in ``obs``.

    ``join='inner'`` keeps only genes present in every sample; ``'outer'`` keeps all genes
    (missing values filled with 0). Cell names are made unique by appending the sample name.
    """
    if len(samples) < 2:
        raise ValueError("Select at least two samples to merge.")
    if join not in ("inner", "outer"):
        raise ValueError("join must be 'inner' or 'outer'.")
    adatas: List[ad.AnnData] = []
    for adata in samples.values():
        adata = adata.copy()
        adata.var_names_make_unique()
        adatas.append(adata)
    merged = ad.concat(
        adatas,
        join=join,
        label="sample",
        keys=list(samples),
        index_unique="-",
        fill_value=0,
    )
    if merged.n_vars == 0:
        raise ValueError("The samples share no genes; try the 'outer' join.")
    return merged


def adata_to_h5ad_bytes(adata: ad.AnnData) -> bytes:
    """Serialise an AnnData object to .h5ad bytes (for a download button)."""
    fd, path = tempfile.mkstemp(suffix=".h5ad")
    os.close(fd)
    try:
        adata.write_h5ad(path)
        with open(path, "rb") as handle:
            return handle.read()
    finally:
        os.unlink(path)
