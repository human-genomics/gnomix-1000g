"""Pedigree labels used by the research browser."""

import pandas as pd

from gnomix1000g.common import ANCESTRIES
from gnomix1000g.site import sample_rows


def test_family_roles_include_duos_and_multiple_generations():
    samples = pd.DataFrame([
        ("father", "0", "0"), ("mother", "0", "0"),
        ("child", "father", "mother"), ("grandchild", "child", "outside-panel"),
        ("unrelated", "outside-panel", "0"),
    ], columns=["IID", "PAT", "MAT"])
    glob = samples[["IID"]].assign(Population="ASW", SuperPop="AFR", gnofix_in_final=False)
    for ancestry in ANCESTRIES:
        glob[ancestry] = 0.0
    rows = {r["id"]: r for r in sample_rows(glob, samples, set())}

    assert rows["father"]["family"] == ["Trio parent"]
    assert rows["mother"]["family"] == ["Trio parent"]
    assert rows["child"]["family"] == ["Trio child", "Duo parent"]
    assert rows["grandchild"]["family"] == ["Duo child"]
    assert rows["unrelated"]["family"] == ["Unrelated"]
    assert rows["child"]["kids"] == ["grandchild"]
    assert rows["grandchild"]["par"] == ["child"]
