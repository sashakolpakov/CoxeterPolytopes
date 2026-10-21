import json
from pathlib import Path
import shutil

import pytest

from scripts.certify_rational_p17_canonical_cusps import build_certificate


@pytest.mark.skipif(
    shutil.which("cddexec_gmp") is None,
    reason="the exact P17 dual-cone replay requires GMP cddlib",
)
def test_p17_canonical_negative_branches_are_bounded_or_one_cusp_orbit():
    certificate = build_certificate()
    frozen = json.loads(
        Path("data/enumerations/rational_p17_canonical_cusps.json").read_text()
    )
    assert certificate == frozen
    assert certificate["nef_extreme_ray_count"] == 807
    assert certificate["negative_nef_branch_census"] == {
        "infeasible": 735,
        "bounded": 64,
        "unbounded": 8,
        "largest_exact_bounded_a": "123",
    }
    assert certificate["weighted_diagram_automorphism_order"] == 8
    assert all(certificate["checks"].values())
    assert certificate["standard_cusp_lemma"]["conclusion"] == "K=(-3,1^17)"
