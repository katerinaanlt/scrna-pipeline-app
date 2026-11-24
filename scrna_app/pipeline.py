"""Analysis backend: the ``AnnDataHandler`` from the design report's class diagram.

Every step works on a copy first and only replaces the stored dataset when it succeeds, so a
failed step (bad parameters, too few cells...) never leaves the session half-processed.
"""
from __future__ import annotations

from typing import Dict, List, Optional

import pandas as pd
import scanpy as sc
from anndata import AnnData

from . import data_io


class PipelineError(ValueError):
    """A user-facing problem (missing prerequisite step, invalid parameters...)."""


class AnnDataHandler:
    def __init__(self, adata: Optional[AnnData] = None, description: str = ""):
        self.adata: Optional[AnnData] = adata
        self.description = description
        self.deg_skipped: List[str] = []

    # ---- dataset management -------------------------------------------------
    @property
    def has_data(self) -> bool:
        return self.adata is not None

    def _require(self) -> AnnData:
        if self.adata is None:
            raise PipelineError("No dataset loaded. Load data first.")
        return self.adata

    def set_adata(self, adata: AnnData, description: str = "") -> None:
        self.adata = adata
        self.description = description
        self.deg_skipped = []

    def load_file(self, filename: str, data: bytes, transpose: bool = False) -> None:
        self.set_adata(data_io.load_file(filename, data, transpose=transpose), filename)

    # ---- preprocessing ------------------------------------------------------
    def filter_data(self, min_genes: int = 200, min_cells: int = 3) -> AnnData:
        """Drop cells with < min_genes detected genes and genes seen in < min_cells cells."""
        adata = self._require().copy()
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        if adata.n_obs == 0 or adata.n_vars == 0:
            raise PipelineError("Filtering removed every cell or gene. Lower the thresholds.")
        self.adata = adata
        return adata

    def normalize(
        self,
        target_sum: float = 1e4,
        n_top_genes: int = 2000,
        max_value: float = 10.0,
        already_normalized: bool = False,
    ) -> AnnData:
        """normalize_total -> log1p -> highly variable genes -> scale.

        The full log-normalised matrix is stored in ``adata.raw`` *before* subsetting to the
        highly variable genes, so any gene can still be plotted and tested for DEG later.
        """
        adata = self._require().copy()
        if not already_normalized:
            sc.pp.normalize_total(adata, target_sum=target_sum)
            sc.pp.log1p(adata)
        adata.raw = adata
        sc.pp.highly_variable_genes(adata, n_top_genes=min(n_top_genes, adata.n_vars), flavor="seurat")
        adata = adata[:, adata.var["highly_variable"]].copy()
        sc.pp.scale(adata, max_value=max_value)
        self.adata = adata
        return adata

    def run_pca(self, n_comps: int = 50) -> AnnData:
        adata = self._require()
        n = min(n_comps, adata.n_obs - 1, adata.n_vars - 1)
        if n < 2:
            raise PipelineError("Not enough cells/genes left for PCA.")
        sc.tl.pca(adata, n_comps=n)
        return adata

    def _require_pca(self) -> AnnData:
        adata = self._require()
        if "X_pca" not in adata.obsm:
            raise PipelineError("Run normalisation + PCA first.")
        return adata

    def run_umap(self, n_neighbors: int = 15, n_pcs: int = 40) -> AnnData:
        adata = self._require_pca()
        sc.pp.neighbors(
            adata,
            n_neighbors=min(n_neighbors, adata.n_obs - 1),
            n_pcs=min(n_pcs, adata.obsm["X_pca"].shape[1]),
        )
        sc.tl.umap(adata)
        return adata

    def run_tsne(self, n_pcs: int = 40, perplexity: float = 30.0) -> AnnData:
        adata = self._require_pca()
        perplexity = min(perplexity, max(5.0, (adata.n_obs - 1) / 3))
        sc.tl.tsne(adata, n_pcs=min(n_pcs, adata.obsm["X_pca"].shape[1]), perplexity=perplexity)
        return adata

    def cluster_data(self, resolution: float = 1.0, key: str = "leiden") -> AnnData:
        """Leiden clustering on the neighbour graph; higher resolution = more clusters."""
        adata = self._require()
        if "neighbors" not in adata.uns:
            raise PipelineError("Compute the neighbour graph (UMAP step) before clustering.")
        sc.tl.leiden(adata, resolution=resolution, key_added=key,
                     flavor="igraph", n_iterations=2, directed=False)
        return adata

    # ---- differential expression -------------------------------------------
    def run_deg(self, groupby: str = "leiden", method: str = "t-test", n_genes: int = 10) -> pd.DataFrame:
        """Rank genes per group (vs. all other cells) and return the top ``n_genes`` table."""
        adata = self._require()
        if groupby not in adata.obs:
            raise PipelineError(f"'{groupby}' is not a column of the cell metadata.")
        adata.obs[groupby] = adata.obs[groupby].astype("category")
        counts = adata.obs[groupby].value_counts()
        usable = [str(g) for g, n in counts.items() if n >= 2]
        self.deg_skipped = [str(g) for g, n in counts.items() if n < 2]
        if len(usable) < 2:
            raise PipelineError("Need at least two groups with 2+ cells each.")
        n_total = adata.raw.n_vars if adata.raw is not None else adata.n_vars
        sc.tl.rank_genes_groups(adata, groupby, groups=usable, method=method,
                                n_genes=min(100, n_total))
        return self.deg_table(n_genes)

    def deg_table(self, n_genes: int = 10) -> Optional[pd.DataFrame]:
        """Top genes per group as a table (columns = groups, rows = rank), or None if not computed."""
        result = None if self.adata is None else self.adata.uns.get("rank_genes_groups")
        if result is None:
            return None
        return pd.DataFrame(result["names"]).head(n_genes)

    def suggest_genes(self, n: int = 4) -> List[str]:
        """Best-ranked gene of the first groups, as a sensible default for the plot pages."""
        table = self.deg_table(1)
        if table is None:
            return []
        return list(dict.fromkeys(table.iloc[0].tolist()))[:n]

    # ---- introspection helpers for the UI ------------------------------------
    def gene_names(self) -> List[str]:
        adata = self._require()
        return list((adata.raw.var_names if adata.raw is not None else adata.var_names))

    def categorical_columns(self) -> List[str]:
        obs = self._require().obs
        cols = []
        for col in obs.columns:
            series = obs[col]
            if isinstance(series.dtype, pd.CategoricalDtype) or pd.api.types.is_string_dtype(series):
                if series.nunique() <= 100:
                    cols.append(str(col))
        return cols

    def embeddings(self) -> List[str]:
        adata = self._require()
        return [name for name in ("umap", "tsne") if f"X_{name}" in adata.obsm]

    def status(self) -> Dict[str, bool]:
        adata = self.adata
        if adata is None:
            return {}
        return {
            "Normalised": adata.raw is not None,
            "PCA": "X_pca" in adata.obsm,
            "UMAP": "X_umap" in adata.obsm,
            "t-SNE": "X_tsne" in adata.obsm,
            "Clusters": "leiden" in adata.obs,
            "DEG": "rank_genes_groups" in adata.uns,
        }
