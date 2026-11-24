import numpy as np
import pytest

pytest.importorskip("anndata")

from scrna_app import data_io  # noqa: E402
from scrna_app.demo import make_synthetic_adata  # noqa: E402


def test_h5ad_roundtrip():
    adata = make_synthetic_adata(n_cells=60, n_genes=120, markers_per_group=10)
    restored = data_io.read_h5ad_bytes(data_io.adata_to_h5ad_bytes(adata))
    assert restored.shape == adata.shape
    assert list(restored.var_names) == list(adata.var_names)


def test_csv_cells_by_genes_and_transpose():
    csv = b"cell,G1,G2,G3\nc1,1,0,3\nc2,0,2,0\n"
    adata = data_io.load_file("x.csv", csv)
    assert adata.shape == (2, 3) and list(adata.var_names) == ["G1", "G2", "G3"]
    flipped = data_io.load_file("x.csv", b"gene,c1,c2\nG1,1,0\nG2,0,2\n", transpose=True)
    assert flipped.shape == (2, 2) and list(flipped.obs_names) == ["c1", "c2"]


def test_csv_rejects_non_numeric_and_unknown_extension():
    with pytest.raises(ValueError):
        data_io.load_file("x.csv", b"cell,G1,label\nc1,1,a\nc2,2,b\n")
    with pytest.raises(ValueError):
        data_io.load_file("x.xlsx", b"")


def test_merge_adds_sample_column_and_unique_cells():
    a = make_synthetic_adata(n_cells=40, n_genes=100, markers_per_group=10, seed=1)
    b = make_synthetic_adata(n_cells=30, n_genes=100, markers_per_group=10, seed=2)
    merged = data_io.merge_samples({"s1": a, "s2": b})
    assert merged.n_obs == 70
    assert set(merged.obs["sample"]) == {"s1", "s2"}
    assert merged.obs_names.is_unique


def test_merge_needs_two_samples_and_valid_join():
    a = make_synthetic_adata(n_cells=20, n_genes=100, markers_per_group=10)
    with pytest.raises(ValueError):
        data_io.merge_samples({"only": a})
    with pytest.raises(ValueError):
        data_io.merge_samples({"a": a, "b": a}, join="weird")
