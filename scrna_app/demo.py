"""Synthetic scRNA-seq data so the app (and the tests) work without any real dataset or network."""
from __future__ import annotations

from typing import List, Tuple

import numpy as np
import pandas as pd


def make_synthetic_counts(
    n_cells: int = 600,
    n_genes: int = 1000,
    n_groups: int = 4,
    markers_per_group: int = 25,
    seed: int = 0,
) -> Tuple[np.ndarray, pd.DataFrame, List[str]]:
    """Poisson counts with ``n_groups`` cell populations, each with its own marker genes.

    Marker genes are named ``MARKER<group>_<k>`` and are ~8x more expressed in their own group,
    so clustering and differential expression have a known ground truth to recover.
    Returns ``(counts[cells x genes], obs DataFrame with 'true_group', gene_names)``.
    """
    if n_groups < 2:
        raise ValueError("n_groups must be at least 2")
    if n_groups * markers_per_group > n_genes:
        raise ValueError("n_genes is too small for the requested number of marker genes")

    rng = np.random.default_rng(seed)
    groups = np.arange(n_cells) % n_groups          # balanced groups
    rng.shuffle(groups)

    base = rng.gamma(shape=0.6, scale=1.0, size=n_genes)
    means = np.tile(base, (n_groups, 1))
    gene_names = [f"Gene{i:04d}" for i in range(n_genes)]
    for g in range(n_groups):
        lo, hi = g * markers_per_group, (g + 1) * markers_per_group
        means[g, lo:hi] = np.maximum(base[lo:hi], 0.3) * 8
        gene_names[lo:hi] = [f"MARKER{g}_{k:02d}" for k in range(markers_per_group)]

    depth = rng.lognormal(mean=0.0, sigma=0.3, size=n_cells)
    counts = rng.poisson(means[groups] * depth[:, None]).astype(np.float32)

    obs = pd.DataFrame(
        {"true_group": [f"group_{g}" for g in groups]},
        index=[f"cell_{i:04d}" for i in range(n_cells)],
    )
    return counts, obs, gene_names


def make_synthetic_adata(**kwargs):
    """Same as :func:`make_synthetic_counts` but returned as an AnnData object (sparse X)."""
    import anndata as ad
    from scipy import sparse

    counts, obs, genes = make_synthetic_counts(**kwargs)
    return ad.AnnData(X=sparse.csr_matrix(counts), obs=obs, var=pd.DataFrame(index=pd.Index(genes)))
