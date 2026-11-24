import numpy as np
import pytest

from scrna_app.demo import make_synthetic_counts


def test_shapes_and_names():
    counts, obs, genes = make_synthetic_counts(n_cells=200, n_genes=300, n_groups=3, markers_per_group=10)
    assert counts.shape == (200, 300)
    assert len(genes) == len(set(genes)) == 300
    assert list(obs.columns) == ["true_group"] and obs.index.is_unique
    assert obs["true_group"].nunique() == 3


def test_markers_are_enriched_in_their_own_group():
    counts, obs, genes = make_synthetic_counts(seed=1)
    for g in range(4):
        cols = [i for i, name in enumerate(genes) if name.startswith(f"MARKER{g}_")]
        own = (obs["true_group"] == f"group_{g}").to_numpy()
        assert counts[own][:, cols].mean() > 2 * counts[~own][:, cols].mean()


def test_deterministic_for_a_seed():
    a, _, _ = make_synthetic_counts(seed=3)
    b, _, _ = make_synthetic_counts(seed=3)
    assert np.array_equal(a, b)


def test_rejects_impossible_requests():
    with pytest.raises(ValueError):
        make_synthetic_counts(n_genes=50, n_groups=4, markers_per_group=25)
    with pytest.raises(ValueError):
        make_synthetic_counts(n_groups=1)
