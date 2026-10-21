from __future__ import annotations

import json

from scripts.certify_jacobian_kummer_h12_fixed_lattice import (
    OUTPUT,
    build_certificate,
)


def test_archived_fixed_lattice_certificate_replays_exactly() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    rebuilt = build_certificate()
    assert rebuilt == archived
    assert archived["all_checks_pass"] is True


def test_fixed_lattice_identity_and_reflection_types() -> None:
    certificate = json.loads(OUTPUT.read_text(encoding="utf-8"))
    fixed = certificate["fixed_lattice"]
    target = certificate["target_identification"]
    reflections = certificate["facet_reflections"]

    assert fixed["rank"] == 13
    assert fixed["signature"] == [1, 12, 0]
    assert fixed["determinant"] == 256
    assert fixed["smith_diagonal"] == [1] * 8 + [2] * 2 + [4] * 3
    assert fixed["canonical_discriminant_form"] == (
        "u(2) + <3/4> + <3/4> + <7/4>"
    )
    assert target["lattice"] == "U(2) + D5(-1)^2 + <-4>"
    assert reflections["primitive_square_divisibility_histogram"] == {
        "-2,1": 16,
        "-4,2": 36,
        "-4,4": 96,
        "-8,4": 24,
    }
    assert reflections["all_reflections_integral_involutive_isometries"] is True
