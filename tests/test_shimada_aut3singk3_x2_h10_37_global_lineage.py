from __future__ import annotations

import hashlib
import json
from pathlib import Path

from scripts.certify_shimada_aut3singk3_x2_h10_37_global_lineage import (
    _candidate_relations,
    _relation_histogram,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT
    / "data/enumerations/shimada_aut3singk3_x2_h10_37_global_lineage.json"
)
CANDIDATE = (
    ROOT
    / "data/enumerations/shimada_aut3singk3_x2_rational_involution_h10_37.json"
)


def test_x2_h10_37_bounded_global_lineage_certificate() -> None:
    raw = CERTIFICATE.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == (
        "52b4800b07608804b7834aad65a4671b3255ef918c07779497742f11091a984d"
    )
    record = json.loads(raw)
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))

    assert record["all_checks_pass"] is True
    assert record["claim"]["worldwide_priority"] == "open"
    assert record["corpus"]["ordinary_record_count"] == 1214
    assert record["corpus"]["doubling_record_count"] == 2959
    assert len(record["corpus"]["ordinary_37_node_records"]) == 2
    assert sum(
        row["record_count"]
        for row in record["corpus"]["doubling_37_node_histogram_classes"]
    ) == 43
    assert record["corpus"]["all_37_node_record_count"] == 45
    assert record["corpus"]["candidate_histogram_match_count"] == 0

    relations = _candidate_relations(
        candidate["rational_involution_fold"]["facet_gram"]
    )
    assert _relation_histogram(relations) == {
        "0": 121,
        "1": 55,
        "2": 405,
        "3": 15,
        "4": 70,
    }
    assert record["candidate"]["relation_histogram"] == _relation_histogram(
        relations
    )

    competing_histograms = [
        row["relation_histogram"]
        for row in record["corpus"]["ordinary_37_node_records"]
    ] + [
        row["relation_histogram"]
        for row in record["corpus"]["doubling_37_node_histogram_classes"]
    ]
    assert record["candidate"]["relation_histogram"] not in competing_histograms
