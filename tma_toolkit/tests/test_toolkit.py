"""End-to-end test on a synthetic 200-core TMA:  pytest tma_toolkit/tests -q"""

import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, PKG)
sys.path.insert(0, HERE)

from make_synthetic_tma import make  # noqa: E402
from tma_toolkit import TMACohort, validate_core_metadata  # noqa: E402


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = str(tmp_path_factory.mktemp("tma"))
    meta, truth, missing = make(root)
    out = os.path.join(root, "tma.h5ad")
    subprocess.run([sys.executable, os.path.join(PKG, "scripts", "build_tma_object.py"),
                    "--slides", os.path.join(root, "slides.csv"),
                    "--core-metadata", os.path.join(root, "core_metadata.csv"),
                    "--project-dir", root, "--out", out, "--min-cells-per-core", "100"], check=True)
    return root, out, meta, truth, missing


def test_metadata_validation_catches_spelling_variants():
    df = pd.DataFrame({"core_id": ["a", "b"], "slide_id": ["s", "s"], "core_row": ["A", "A"],
                       "core_col": [1, 2], "organ": ["colon", "Colon "]})
    with pytest.warns(UserWarning, match="spelling variants"):
        out = validate_core_metadata(df)
    assert out.organ.nunique() == 1


def test_every_present_core_gets_the_right_label(built):
    root, out, meta, truth, missing = built
    c = TMACohort.read(out)
    adata = c.adata
    # cells keep their true core (row/col recovered despite rotation and a missing column)
    t = truth.reindex(adata.obs_names)
    agree = (adata.obs["core_id"].astype(str).to_numpy() == t.core_id.to_numpy()).mean()
    assert agree > 0.99, agree
    found = set(c.cores.index[c.cores.found_on_slide])
    assert found == set(meta.core_id) - {meta.core_id[i] for i in missing}


def test_selection_and_core_sets(built):
    root, out, *_ = built
    c = TMACohort.read(out)
    colon = c.select(organ="colon")
    assert colon.n_cores > 0 and set(colon.cores.organ) == {"colon"}      # "Colon " typo unified
    tum = c.select(tissue_type="tumor")
    both = colon & tum
    assert set(both.core_ids) == set(colon.core_ids) & set(tum.core_ids)
    old = c.select(query="age >= 60")
    assert (old.cores.age >= 60).all()
    rng = c.select(age=(50, 55))
    assert rng.cores.age.between(50, 55).all()
    c.save_core_set("colon_tumor", both, "test")
    c.write(out)
    c2 = TMACohort.read(out)
    assert c2.core_set("colon_tumor").core_ids == both.core_ids
    sub = both.adata()
    assert set(sub.obs.core_id.astype(str)) == set(both.core_ids)
    per = c.select(organ=["colon", "lung"]).sample(2, by="organ")
    assert per.cores.groupby("organ", observed=True).size().max() <= 2


def test_comparisons(built, tmp_path):
    root, out, meta, truth, missing = built
    c = TMACohort.read(out)
    c.adata.obs["cell_type"] = pd.Categorical(truth.reindex(c.adata.obs_names).cell_type)
    sel = c.select(organ=["colon", "lung", "liver"])
    comp = sel.compare_composition("cell_type", group_col="organ")
    assert comp.padj.min() < 0.05                                       # organs differ by design
    comp_pat = sel.compare_composition("cell_type", group_col="organ", unit="patient_id")
    assert len(comp_pat) > 0
    de = c.select().differential("tissue_type", "tumor", "normal")
    assert de.index[0] == "GENE_A" and de.loc["GENE_A", "log2FC"] > 0      # the planted tumour gene
    de_ct = c.select().differential("tissue_type", "tumor", "normal", label="cell_type", label_value="T_cell")
    assert "GENE_A" in de_ct.index[:3]
    corr = sel.similarity()
    rel = sel.similarity_by_relation(corr)
    assert {"same_organ", "same_patient_id"} <= set(rel.columns)
    g = sel.gene_by_core(["CD3E", "EPCAM"], stat="frac")
    assert g.shape[1] == 2
    for name, fig in [("gallery", sel.plot_cores("cell_type", max_cores=10)),
                      ("composition", sel.plot_composition("cell_type", group_col="organ")),
                      ("similarity", sel.plot_similarity(corr))]:
        fig.savefig(tmp_path / f"{name}.png")
    assert os.path.exists(os.path.join(root, "figures", "dearray_TMA2.png"))
