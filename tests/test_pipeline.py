import pytest

pytest.importorskip("scanpy")
pytest.importorskip("igraph")

from scrna_app.demo import make_synthetic_adata  # noqa: E402
from scrna_app.pipeline import AnnDataHandler, PipelineError  # noqa: E402
from scrna_app.plotting import VisualizationManager  # noqa: E402


@pytest.fixture(scope="module")
def processed():
    handler = AnnDataHandler(make_synthetic_adata(n_cells=400, n_genes=800, markers_per_group=25), "demo")
    handler.filter_data(min_genes=50, min_cells=3)
    handler.normalize(n_top_genes=300)
    handler.run_pca(n_comps=20)
    handler.run_umap(n_neighbors=15, n_pcs=20)
    handler.cluster_data(resolution=0.5)
    return handler


def test_steps_require_their_prerequisites():
    empty = AnnDataHandler()
    with pytest.raises(PipelineError):
        empty.filter_data()
    fresh = AnnDataHandler(make_synthetic_adata(n_cells=100, n_genes=200, markers_per_group=10))
    with pytest.raises(PipelineError):
        fresh.run_umap()
    with pytest.raises(PipelineError):
        fresh.cluster_data()


def test_filter_that_removes_everything_leaves_data_untouched():
    handler = AnnDataHandler(make_synthetic_adata(n_cells=100, n_genes=200, markers_per_group=10))
    before = handler.adata.shape
    with pytest.raises(PipelineError):
        handler.filter_data(min_genes=10_000)
    assert handler.adata.shape == before


def test_normalisation_keeps_all_genes_in_raw(processed):
    assert processed.adata.raw is not None
    assert processed.adata.n_vars <= 300 < processed.adata.raw.n_vars
    assert "X_pca" in processed.adata.obsm and "X_umap" in processed.adata.obsm


def test_clusters_recover_the_planted_groups(processed):
    obs = processed.adata.obs
    assert 2 <= obs["leiden"].nunique() <= 12
    purity = obs.groupby("leiden", observed=True)["true_group"].agg(lambda s: s.value_counts(normalize=True).iloc[0])
    assert purity.mean() > 0.8


def test_deg_finds_marker_genes(processed):
    table = processed.run_deg("leiden", "t-test", n_genes=10)
    assert table.shape[0] == 10 and table.shape[1] == processed.adata.obs["leiden"].nunique()
    for column in table.columns:
        assert sum(name.startswith("MARKER") for name in table[column]) >= 3
    assert processed.deg_table(5).shape[0] == 5


def test_deg_rejects_unknown_column(processed):
    with pytest.raises(PipelineError):
        processed.run_deg("does_not_exist")


def test_plots_return_figures(processed):
    from matplotlib.figure import Figure

    viz = VisualizationManager(processed.adata)
    genes = processed.suggest_genes(3)
    assert genes
    assert isinstance(viz.plot_clusters("leiden"), Figure)
    assert isinstance(viz.plot_umap(genes[0]), Figure)
    assert isinstance(viz.plot_violin(genes, "leiden"), Figure)
    assert isinstance(viz.plot_dotplot(genes, "leiden"), Figure)
    assert isinstance(viz.plot_heatmap(genes, "leiden"), Figure)
    assert isinstance(viz.plot_rank_genes(5), Figure)


def test_tsne_runs(processed):
    processed.run_tsne(n_pcs=20, perplexity=20)
    assert "tsne" in processed.embeddings()
