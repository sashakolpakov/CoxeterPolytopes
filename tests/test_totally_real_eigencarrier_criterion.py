import json
from pathlib import Path

from scripts.certify_totally_real_eigencarrier_criterion import (
    carrier_dimension_ceiling,
    period_dimension,
    trace_rank,
    verify,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "data/candidates/totally_real_eigencarrier_criterion.json"


def test_general_rank_and_period_formulas() -> None:
    assert trace_rank(2, 9) == 18
    assert trace_rank(2, 10) == 20
    assert trace_rank(2, 11) == 22
    assert carrier_dimension_ceiling(22, 2) == 9
    assert carrier_dimension_ceiling(24, 2) == 10
    assert carrier_dimension_ceiling(22, 3) == 5
    assert period_dimension(22, 18) == 2
    assert period_dimension(22, 20) == 0
    assert period_dimension(24, 22) == 0


def test_frozen_totally_real_eigencarrier_audit() -> None:
    archived = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    summary = verify(archived)
    assert summary["rank_ladder"] == [18, 20, 22]
    assert summary["dimension_ceilings"] == {
        "K3_degree_2": 9,
        "OG10_degree_2": 10,
        "formula": "floor((b2-2)/[F:Q])-1",
    }
    assert summary["baseline_statuses"] == [
        "complete compact 11-facet eigen-nef chamber",
        "complete norm-two chamber has infinite volume",
        "complete norm-two movable chamber has infinite volume",
    ]
