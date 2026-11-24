from scrna_app.utils import parse_gene_list, resolve_genes, sample_name, unique_names


def test_parse_gene_list_handles_mixed_separators_and_duplicates():
    assert parse_gene_list(" PPY, INS  GCG;SST,,PPY ") == ["PPY", "INS", "GCG", "SST"]
    assert parse_gene_list("   ") == []


def test_resolve_genes_exact_and_case_insensitive():
    found, missing = resolve_genes(["PPY", "ins", "NOPE"], ["PPY", "INS", "GCG"])
    assert found == ["PPY", "INS"]
    assert missing == ["NOPE"]


def test_resolve_genes_does_not_duplicate():
    found, _ = resolve_genes(["ppy", "PPY"], ["PPY"])
    assert found == ["PPY"]


def test_sample_names():
    assert sample_name("data/Sample_1.h5ad") == "Sample_1"
    assert unique_names(["a", "b", "a", "a"]) == ["a", "b", "a_2", "a_3"]
