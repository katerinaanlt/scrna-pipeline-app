"""``VisualizationManager`` from the design report: Scanpy plots returned as matplotlib figures."""
from __future__ import annotations

import threading
from typing import List, Optional

import matplotlib

matplotlib.use("Agg")  # headless: figures are rendered by Streamlit, never by a GUI window
import matplotlib.pyplot as plt  # noqa: E402
import scanpy as sc  # noqa: E402
from anndata import AnnData  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

# pyplot keeps global state and Streamlit serves each session in its own thread.
_LOCK = threading.Lock()


class VisualizationManager:
    def __init__(self, adata: AnnData):
        self.adata = adata

    def _render(self, func, *args, **kwargs) -> Figure:
        with _LOCK:
            plt.close("all")
            func(*args, show=False, **kwargs)
            fig = plt.gcf()
            plt.close(fig)  # detach from pyplot; the Figure object stays usable
        return fig

    def plot_embedding(self, basis: str, color: Optional[List[str] | str] = None,
                       legend_loc: str = "right margin") -> Figure:
        plotter = {"umap": sc.pl.umap, "tsne": sc.pl.tsne}[basis]
        return self._render(plotter, self.adata, color=color, legend_loc=legend_loc)

    def plot_umap(self, color=None, legend_loc: str = "right margin") -> Figure:
        return self.plot_embedding("umap", color, legend_loc)

    def plot_tsne(self, color=None, legend_loc: str = "right margin") -> Figure:
        return self.plot_embedding("tsne", color, legend_loc)

    def plot_clusters(self, key: str = "leiden", basis: str = "umap") -> Figure:
        return self.plot_embedding(basis, key)

    def plot_violin(self, genes: List[str], groupby: str) -> Figure:
        return self._render(sc.pl.violin, self.adata, keys=genes, groupby=groupby)

    def plot_dotplot(self, genes: List[str], groupby: str) -> Figure:
        return self._render(sc.pl.dotplot, self.adata, var_names=genes, groupby=groupby)

    def plot_heatmap(self, genes: List[str], groupby: str) -> Figure:
        return self._render(sc.pl.heatmap, self.adata, var_names=genes, groupby=groupby)

    def plot_rank_genes(self, n_genes: int = 10) -> Figure:
        return self._render(sc.pl.rank_genes_groups, self.adata, n_genes=n_genes, sharey=False)
