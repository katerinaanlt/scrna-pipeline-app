"""Small dependency-free helpers (no Scanpy/Streamlit imports, so they are trivially testable)."""
from __future__ import annotations

import os
import re
from typing import Iterable, List, Tuple


def parse_gene_list(text: str) -> List[str]:
    """Split free text ("PPY, INS  GCG;SST") into unique gene names, preserving order."""
    tokens = [t for t in re.split(r"[,\s;]+", text.strip()) if t]
    return list(dict.fromkeys(tokens))


def resolve_genes(requested: Iterable[str], available: Iterable[str]) -> Tuple[List[str], List[str]]:
    """Match requested gene names against the dataset.

    Exact matches win; otherwise a case-insensitive match is accepted (so ``ppy`` finds ``PPY``,
    and ``Ins1`` finds ``INS1`` in mouse/human data). Returns ``(found, missing)`` where *found*
    uses the dataset's canonical spelling.
    """
    available = list(available)
    exact = set(available)
    folded = {}
    for name in available:
        folded.setdefault(name.casefold(), name)
    found: List[str] = []
    missing: List[str] = []
    for gene in requested:
        if gene in exact:
            match = gene
        else:
            match = folded.get(gene.casefold())
        if match is None:
            missing.append(gene)
        elif match not in found:
            found.append(match)
    return found, missing


def sample_name(filename: str) -> str:
    """'data/Sample_1.h5ad' -> 'Sample_1'."""
    return os.path.splitext(os.path.basename(filename))[0]


def unique_names(names: Iterable[str]) -> List[str]:
    """Make names unique by suffixing duplicates: ['a', 'a'] -> ['a', 'a_2']."""
    seen = {}
    result = []
    for name in names:
        seen[name] = seen.get(name, 0) + 1
        result.append(name if seen[name] == 1 else f"{name}_{seen[name]}")
    return result
