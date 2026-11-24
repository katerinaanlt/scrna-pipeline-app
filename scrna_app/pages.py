"""``UserInterface`` from the design report: one render function per sidebar menu entry.

Streamlit re-runs the whole script on every interaction, so the analysis state lives in
``st.session_state`` (as an :class:`AnnDataHandler` that owns the AnnData object).
"""
from __future__ import annotations

from typing import Callable, Optional

import streamlit as st

from . import data_io
from .demo import make_synthetic_adata
from .pipeline import AnnDataHandler, PipelineError
from .plotting import VisualizationManager
from .utils import parse_gene_list, resolve_genes, sample_name, unique_names


# ---- shared helpers ----------------------------------------------------------------
def get_handler() -> AnnDataHandler:
    if "handler" not in st.session_state:
        st.session_state["handler"] = AnnDataHandler()
    return st.session_state["handler"]


def require_data() -> Optional[AnnDataHandler]:
    handler = get_handler()
    if not handler.has_data:
        st.info("No dataset loaded yet. Go to **Data** (or **Merge**) first.")
        return None
    return handler


def run_step(label: str, func: Callable, *args, **kwargs) -> bool:
    """Run an analysis step with a spinner and friendly error reporting."""
    try:
        with st.spinner(f"{label}…"):
            func(*args, **kwargs)
    except PipelineError as exc:
        st.error(str(exc))
        return False
    except Exception as exc:  # UI boundary: never show a raw traceback to a biologist
        st.error(f"{label} failed: {exc}")
        return False
    st.success(f"{label}: done")
    return True


def show_figure(fig) -> None:
    st.pyplot(fig)


def show_summary(handler: AnnDataHandler) -> None:
    adata = handler.adata
    c1, c2, c3 = st.columns(3)
    c1.metric("Cells", f"{adata.n_obs:,}")
    c2.metric("Genes", f"{adata.n_vars:,}")
    c3.metric("Source", handler.description or "-")
    status = handler.status()
    st.caption("  ·  ".join(f"{'✅' if ok else '⬜'} {name}" for name, ok in status.items()))


def sidebar_summary() -> None:
    handler = get_handler()
    if handler.has_data:
        st.sidebar.success(f"{handler.adata.n_obs:,} cells × {handler.adata.n_vars:,} genes")
    else:
        st.sidebar.warning("No dataset loaded")


# ---- pages -------------------------------------------------------------------------
def page_data() -> None:
    st.header("Data")
    handler = get_handler()
    uploaded = st.file_uploader("Upload a dataset (.h5ad or .csv)", type=["h5ad", "csv", "tsv"])
    transpose = st.checkbox("CSV has genes as rows (transpose)", value=False,
                            help="Default layout: rows = cells, columns = genes.")
    left, right = st.columns(2)
    if left.button("Load file", disabled=uploaded is None):
        run_step("Loading file", handler.load_file, uploaded.name, uploaded.getvalue(), transpose)
    if right.button("Load demo dataset"):
        run_step("Generating demo data", handler.set_adata, make_synthetic_adata(), "synthetic demo")

    if handler.has_data:
        st.divider()
        show_summary(handler)
        tab_obs, tab_var = st.tabs(["Cell metadata (obs)", "Gene metadata (var)"])
        tab_obs.dataframe(handler.adata.obs.head(50))
        tab_var.dataframe(handler.adata.var.head(50))


def page_merge() -> None:
    st.header("Merge samples")
    st.write("Upload two or more `.h5ad` files. A `sample` column is added to the cell metadata "
             "so you can colour plots by origin.")
    files = st.file_uploader("Samples", type=["h5ad"], accept_multiple_files=True)
    join = st.radio("Genes to keep", ["inner", "outer"], horizontal=True,
                    help="inner = genes present in every sample; outer = all genes (filled with 0).")
    if st.button("Merge", disabled=len(files) < 2):
        names = unique_names(sample_name(f.name) for f in files)

        def _merge():
            samples = {n: data_io.read_h5ad_bytes(f.getvalue()) for n, f in zip(names, files)}
            get_handler().set_adata(data_io.merge_samples(samples, join=join),
                                    f"merge of {len(files)} samples")

        run_step("Merging samples", _merge)
    handler = get_handler()
    if handler.has_data and "sample" in handler.adata.obs:
        show_summary(handler)
        st.write(handler.adata.obs["sample"].value_counts().rename("cells per sample"))


def page_preprocessing() -> None:
    st.header("Preprocessing")
    handler = require_data()
    if handler is None:
        return

    st.subheader("1 · Filtering")
    c1, c2 = st.columns(2)
    min_genes = c1.number_input("Min genes per cell", 0, 10000, 200, step=50)
    min_cells = c2.number_input("Min cells per gene", 0, 1000, 3)
    if st.button("Apply filter"):
        run_step("Filtering", handler.filter_data, int(min_genes), int(min_cells))

    st.subheader("2 · Normalisation and PCA")
    c1, c2, c3, c4 = st.columns(4)
    target = c1.number_input("Target sum", 100, 1_000_000, 10_000, step=1000)
    n_hvg = c2.number_input("Highly variable genes", 50, 20000, 2000, step=100)
    max_value = c3.number_input("Clip scaled values at", 1.0, 100.0, 10.0)
    n_comps = c4.number_input("PCA components", 2, 200, 50)
    already = st.checkbox("Data is already log-normalised (skip normalise + log1p)")
    if st.button("Normalise + PCA"):
        def _norm_pca():
            handler.normalize(float(target), int(n_hvg), float(max_value), already)
            handler.run_pca(int(n_comps))
        run_step("Normalisation and PCA", _norm_pca)

    st.subheader("3 · Embedding")
    c1, c2, c3 = st.columns(3)
    n_neighbors = c1.number_input("Neighbours", 2, 200, 15)
    n_pcs = c2.number_input("PCs to use", 2, 100, 40)
    perplexity = c3.number_input("t-SNE perplexity", 5.0, 100.0, 30.0)
    left, right = st.columns(2)
    if left.button("Compute UMAP"):
        run_step("UMAP", handler.run_umap, int(n_neighbors), int(n_pcs))
    if right.button("Compute t-SNE"):
        run_step("t-SNE", handler.run_tsne, int(n_pcs), float(perplexity))

    st.subheader("4 · Clustering (Leiden)")
    resolution = st.slider("Resolution", 0.1, 3.0, 1.0, 0.1,
                           help="Higher values produce more clusters.")
    if st.button("Run Leiden clustering"):
        run_step("Clustering", handler.cluster_data, float(resolution))

    st.divider()
    show_summary(handler)


def page_deg() -> None:
    st.header("Differential expression")
    handler = require_data()
    if handler is None:
        return
    groups = handler.categorical_columns()
    if not groups:
        st.info("No categorical column to group by. Run clustering first.")
        return
    groupby = st.selectbox("Group cells by", groups,
                           index=groups.index("leiden") if "leiden" in groups else 0)
    method = st.selectbox("Method", ["t-test", "wilcoxon"])
    n_genes = st.slider("Top genes per group", 5, 50, 10)
    if st.button("Run differential expression"):
        run_step("Differential expression", handler.run_deg, groupby, method, n_genes)
        if handler.deg_skipped:
            st.warning("Skipped groups with fewer than 2 cells: " + ", ".join(handler.deg_skipped))

    table = handler.deg_table(n_genes)
    if table is not None:
        st.subheader(f"Top {n_genes} genes per group")
        st.dataframe(table)
        st.download_button("Download table (CSV)", table.to_csv(index=False).encode("utf-8"),
                           file_name="deg_top_genes.csv", mime="text/csv")
        show_figure(VisualizationManager(handler.adata).plot_rank_genes(min(n_genes, 25)))


def _gene_input(handler: AnnDataHandler, key: str):
    default = ", ".join(handler.suggest_genes())
    text = st.text_input("Genes (comma separated)", value=default, key=key)
    found, missing = resolve_genes(parse_gene_list(text), handler.gene_names())
    if missing:
        st.warning("Not found in this dataset: " + ", ".join(missing))
    return found


def page_plots() -> None:
    st.header("Plots")
    handler = require_data()
    if handler is None:
        return
    viz = VisualizationManager(handler.adata)
    kind = st.selectbox("Plot type", ["UMAP / t-SNE", "Violin", "Dot plot", "Heatmap"])

    if kind == "UMAP / t-SNE":
        bases = handler.embeddings()
        if not bases:
            st.info("Compute UMAP or t-SNE first (Preprocessing).")
            return
        basis = st.radio("Embedding", bases, horizontal=True)
        mode = st.radio("Colour by", ["Cluster / metadata", "Gene expression"], horizontal=True)
        if mode == "Cluster / metadata":
            cols = handler.categorical_columns()
            if not cols:
                st.info("No categorical metadata available yet.")
                return
            key = st.selectbox("Column", cols, index=cols.index("leiden") if "leiden" in cols else 0)
            show_figure(viz.plot_embedding(basis, key))
        else:
            genes = _gene_input(handler, "umap_genes")
            if genes:
                show_figure(viz.plot_embedding(basis, genes))
        return

    cols = handler.categorical_columns()
    if not cols:
        st.info("Run clustering first so there are groups to compare.")
        return
    groupby = st.selectbox("Group by", cols, index=cols.index("leiden") if "leiden" in cols else 0)
    genes = _gene_input(handler, f"genes_{kind}")
    if not genes:
        st.info("Enter at least one valid gene name.")
        return
    handler_fn = {"Violin": viz.plot_violin, "Dot plot": viz.plot_dotplot, "Heatmap": viz.plot_heatmap}[kind]
    show_figure(handler_fn(genes, groupby))


def page_export() -> None:
    st.header("Export")
    handler = require_data()
    if handler is None:
        return
    st.write("Download the processed dataset (including embeddings, clusters and DEG results) "
             "to continue in Scanpy or share it.")
    if st.button("Prepare .h5ad"):
        with st.spinner("Writing file…"):
            st.session_state["export_bytes"] = data_io.adata_to_h5ad_bytes(handler.adata)
    if "export_bytes" in st.session_state:
        st.download_button("Download processed.h5ad", st.session_state["export_bytes"],
                           file_name="processed.h5ad", mime="application/octet-stream")


PAGES = {
    "Data": page_data,
    "Merge": page_merge,
    "Preprocessing": page_preprocessing,
    "Differential expression": page_deg,
    "Plots": page_plots,
    "Export": page_export,
}
