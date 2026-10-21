from __future__ import annotations

from collections import Counter
import json

from scripts.certify_shimada_three_discriminant_h18_minus4_reflection_family import (
    EXPECTED_COUNTS,
    EXPECTED_SMITH,
    OUTPUT,
    SCHEMA,
    _sha256,
)


def _load() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_all_934_fresh_h18_classes_are_primitive_minus4_reflections() -> None:
    payload = _load()
    assert payload["schema"] == SCHEMA
    assert payload["summary"]["all_checks_pass"] is True
    assert payload["summary"]["case_count"] == 934
    assert payload["summary"]["case_count_by_discriminant"] == {
        str(key): value for key, value in EXPECTED_COUNTS.items()
    }
    assert payload["summary"]["reflection_denominator_histogram"] == {"2": 934}
    assert payload["summary"]["vector_square_histogram"] == {"-4": 934}
    assert payload["summary"]["vector_divisibility_histogram"] == {"1": 934}
    assert payload["summary"]["fixed_lattice_determinant_histogram"] == {
        "44": 216,
        "60": 144,
        "64": 574,
    }


def test_every_record_replays_the_formula_and_has_a_coxeter_obstruction() -> None:
    payload = _load()
    records = payload["records"]
    assert Counter(record["discriminant"] for record in records) == EXPECTED_COUNTS
    for record in records:
        discriminant = record["discriminant"]
        vector = record["primitive_minus4_vector_primal_row"]
        assert len(vector) == 20
        assert next(value for value in vector if value) > 0
        assert record["vector_square"] == -4
        assert record["vector_divisibility"] == 1
        assert record["reflection_action_denominator"] == 2
        assert record["fixed_lattice_determinant_absolute"] == 4 * discriminant
        assert record["fixed_lattice_smith_diagonal_absolute"] == EXPECTED_SMITH[
            discriminant
        ]
        assert record["certified_non_Coxeter_pair_count"] >= 1
        assert all(record["checks"].values())
    assert min(record["certified_non_Coxeter_pair_count"] for record in records) == 1
    assert sum(
        record["certified_non_Coxeter_pair_count"] == 1 for record in records
    ) == 38


def test_pinned_inputs_have_not_changed() -> None:
    payload = _load()
    for source in payload["sources"].values():
        assert _sha256(__import__("pathlib").Path(source["primary_atlas"])) == source[
            "primary_atlas_sha256"
        ]
        for key in ("inventory", "manifest"):
            path = OUTPUT.parents[2] / source[key]
            assert _sha256(path) == source[f"{key}_sha256"]

