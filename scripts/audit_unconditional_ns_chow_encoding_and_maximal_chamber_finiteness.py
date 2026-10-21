#!/usr/bin/env python3
"""Independent hostile audit of the unconditional NS/Chow synthesis.

This replay does not import the synthesis certifier or either branch owner/
audit producer.  It pins their terminal bytes, rebuilds the sensitive exact
matrix calculations, and enforces the provenance boundary as a directed
dependency graph.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from math import comb, factorial
from pathlib import Path
import re
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
AUDIT_REPORT = (
    ROOT
    / "results/unconditional_ns_chow_encoding_and_maximal_chamber_finiteness_"
    "independent_hostile_audit.md"
)
OUTPUT = (
    ROOT
    / "data/enumerations/unconditional_ns_chow_encoding_and_maximal_chamber_"
    "finiteness_independent_hostile_audit.json"
)
SCHEMA = (
    "coxeter.unconditional-ns-chow-and-maximal-chamber."
    "independent-hostile-audit.v1"
)


SYNTHESIS_OWNER_PINS = {
    "results/unconditional_ns_chow_encoding_and_maximal_chamber_finiteness.md": (
        "ab17efccf6961be32a45d11ae6683f7d6637c42162626653b145ce0529f2336b"
    ),
    "scripts/certify_unconditional_ns_chow_encoding_and_maximal_chamber_finiteness.py": (
        "6711f6c587c6f20c0df596ed5eaf5d385caff4cc360c0087392a920835141438"
    ),
    "data/enumerations/unconditional_ns_chow_encoding_and_maximal_chamber_finiteness.json": (
        "c61c6b16d9268455ec82134c7c012b9690f40bfe4761bd9d1a9d987f67353565"
    ),
    "tests/test_unconditional_ns_chow_encoding_and_maximal_chamber_finiteness.py": (
        "3f71fcb06b34a5ca52cf4c13297aaecee19ed6c31012bfc8b77aa3b97fa90814"
    ),
}

CLASSICAL_OWNER_PINS = {
    "results/unconditional_maximal_arithmetic_reflection_chamber_finiteness.md": (
        "6eee91512840287ce47a195a317fa4bbfd19678a95b52869ec0a45775d61e67b"
    ),
    "scripts/certify_unconditional_maximal_arithmetic_reflection_chamber_finiteness.py": (
        "c32c7969f900093db388718f3ac769f669d42e43b6c2a95c76297eb34899f541"
    ),
    "data/enumerations/unconditional_maximal_arithmetic_reflection_chamber_finiteness.json": (
        "2bce48c0a7156caea3f4b4890411d770450912b4d4b2d3827ce93cc7e9730acc"
    ),
    "tests/test_unconditional_maximal_arithmetic_reflection_chamber_finiteness.py": (
        "4105f4f3d08dfdcbc1d98f272aae6168a99221f6e426f1bcf1fb5b50303d6702"
    ),
}

CLASSICAL_AUDIT_PINS = {
    "results/unconditional_maximal_arithmetic_reflection_chamber_finiteness_independent_hostile_audit.md": (
        "1f1592ef207580b765d56bdca9461429456ef9fa46eb180b0ae65d69dc68ef3c"
    ),
    "scripts/audit_unconditional_maximal_arithmetic_reflection_chamber_finiteness.py": (
        "e6fc5f7543093b6167e6d46612b64c7247007cc16b06f081e14260a283193400"
    ),
    "data/enumerations/unconditional_maximal_arithmetic_reflection_chamber_finiteness_independent_hostile_audit.json": (
        "87ec814ae05d29b04ac068296c4f995144264fed353fa02f607eb22841cff734"
    ),
    "tests/test_unconditional_maximal_arithmetic_reflection_chamber_finiteness_independent_hostile_audit.py": (
        "2ab582e72677922fa772ff57008cb17e307556f5bd595ca767f6130fe7b251ea"
    ),
}

WEAK_OWNER_PINS = {
    "results/universal_weak_chow_action_carrier_and_leverage_no_go.md": (
        "df38b26e0cffd8a391834be922869153af5badf41a632e74a48b472b332e101a"
    ),
    "scripts/certify_universal_weak_chow_action_carrier.py": (
        "d3a19cfe7069b811721a77a624952d92ab5bc7c5574d9e332d2ff8c7f0205eba"
    ),
    "data/enumerations/universal_weak_chow_action_carrier.json": (
        "d75897bbf8d80a662a7b9d3aa3ae729052494e20633d835566db70870bb0fc9e"
    ),
    "tests/test_universal_weak_chow_action_carrier.py": (
        "6da7b9e4590d7e531a7f3502ba2c971e806cbd3f259cb4722169553b7532c802"
    ),
}

WEAK_AUDIT_PINS = {
    "results/universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.md": (
        "91345590b506ef7c7d4287be5a4d94c7563ba7613167ffbe34272e6b7fc422c1"
    ),
    "scripts/audit_universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile.py": (
        "98765eaae927d9265cf69591f02c31d1d0adc71653936d2199f53eb75b16a347"
    ),
    "data/enumerations/universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.json": (
        "a9c19a9f2dce89c03182c9dc7ad8248a0ab75e32e5e58cd1bca6fac7683764f4"
    ),
    "tests/test_universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.py": (
        "46be65edc68679ed6f43da75f59393c69e83630810f8c2f7404c3c40381cdc42"
    ),
}

SCOPE_PINS = {
    "results/unconditional_ag_absw_universality_leverage_no_go.md": (
        "fcc3acb8bee4751bef083cd00a179f9b542602561ba1bf81d89f374dd5523150"
    ),
}

AUDIT_REPORT_SHA256 = (
    "75e93b1938eaf0d89ba4d382f61ca089eeb86ea08102918845d30d6d582d5a4f"
)

PIN_GROUPS = {
    "synthesis_owner": SYNTHESIS_OWNER_PINS,
    "classical_owner": CLASSICAL_OWNER_PINS,
    "classical_independent_audit": CLASSICAL_AUDIT_PINS,
    "weak_chow_owner": WEAK_OWNER_PINS,
    "weak_chow_independent_audit": WEAK_AUDIT_PINS,
    "scope_no_go": SCOPE_PINS,
}
ALL_INPUT_PINS = {
    relative: expected
    for pins in PIN_GROUPS.values()
    for relative, expected in pins.items()
}

PRIMARY_SOURCE_PINS = {
    "absw_2008": {
        "path": "/private/tmp/absw.pdf",
        "sha256": "974ac7c5141ca1c45b515f145215ad115aba867e42551a9b103233ea52930ee4",
        "locus": "Introduction and Theorem 6.1",
    },
    "vinberg_1971": {
        "path": "/private/tmp/vinberg_1971_definition_fields.pdf",
        "sha256": "16a4f497688e6f968b4c63d81a6b049b3f096f1cd4f43e458b2882a5c44ded38",
        "locus": "Theorems 1, 3, and 5",
    },
    "vinberg_dimension": {
        "path": "/private/tmp/vinberg_absence_complete.pdf",
        "sha256": "fd25f8f892eed27e4303a8f3252fd384d8954587a6c1c67f1efed2c38aabe44f",
        "locus": "compact and arithmetic noncompact dimension exclusions",
    },
    "nikulin_2007": {
        "path": "/private/tmp/nikulin_finiteness.pdf",
        "sha256": "9a820fef0b88158b283057f7d97acc89130610a9a4aacc6c53a8215e23ce874a",
        "locus": "Theorem 2.2",
    },
    "nikulin_1982": {
        "path": "/private/tmp/nikulin_1981_eng.pdf",
        "sha256": "696652b6a67c6ed420de6c1cb61d9422e3619e28e081a78131a5323ce2edee9d",
        "locus": "Theorems 2.1.3(2) and 5.2.1; frame constant 62",
    },
    "nikulin_2011": {
        "path": "/private/tmp/nikulin_transition_journal.pdf",
        "sha256": "716e6e7333efbfaa22060b60dc97e8303363a4914409ccdd56ef7a90c923c7af",
        "locus": "Theorem 3.2 and transition bounds",
    },
    "belolipetsky_linowitz_v2": {
        "path": "/tmp/belolipetsky_linowitz_2013.pdf",
        "sha256": "026e0f635f71e214807e0b505bf7297673eb9e0d47a6a91b7bd6dec134eab97f",
        "locus": "Corollary 5.1",
    },
}

DIMENSION_MAX = 29
RANK_MAX = DIMENSION_MAX + 1
FIELD_DEGREE_MAX = 25
TRACE_RANK_MAX = FIELD_DEGREE_MAX * RANK_MAX
BLOWUP_MAX = 4 * TRACE_RANK_MAX - 3
FRAME_CONSTANT = 62
CONJUGATE_HEIGHT = FRAME_CONSTANT**RANK_MAX


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def boolean_checks(checks: dict[str, Any]) -> dict[str, bool]:
    return {key: bool(value) for key, value in checks.items()}


def validate_input_pins() -> dict[str, dict[str, str]]:
    observed: dict[str, dict[str, str]] = {}
    for group, pins in PIN_GROUPS.items():
        group_observed = {
            relative: sha256_file(ROOT / relative) for relative in pins
        }
        if group_observed != pins:
            raise AssertionError(
                f"terminal {group} pins changed: {group_observed} != {pins}"
            )
        observed[group] = group_observed
    if sha256_file(AUDIT_REPORT) != AUDIT_REPORT_SHA256:
        raise AssertionError("independent synthesis audit report changed")
    return observed


def validate_optional_primary_sources() -> None:
    for record in PRIMARY_SOURCE_PINS.values():
        path = Path(record["path"])
        if path.exists() and sha256_file(path) != record["sha256"]:
            raise AssertionError(f"audit-time primary source changed: {path}")


def read_json(relative: str) -> dict[str, Any]:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def local_link_audit(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    pattern = re.compile(r"!?\[[^\]]*\]\((<?[^)\s>]+>?)")
    broken: list[str] = []
    escaped: list[str] = []
    checked: list[str] = []
    root = ROOT.resolve()
    for match in pattern.finditer(text):
        raw = match.group(1).strip("<>")
        if raw.startswith(("http://", "https://", "mailto:", "#")):
            continue
        target = raw.split("#", 1)[0]
        if not target:
            continue
        # A TeX interval such as ``[k:\mathbf Q](n+1)`` has Markdown-like
        # punctuation but is not a link.  All repository artifacts audited
        # here have a filename suffix.
        if not Path(target).suffix:
            continue
        checked.append(target)
        resolved = (path.parent / target).resolve()
        try:
            resolved.relative_to(root)
        except ValueError:
            escaped.append(target)
            continue
        if not resolved.is_file():
            broken.append(target)
    return {
        "local_links_checked": len(checked),
        "broken_local_links": sorted(set(broken)),
        "links_escaping_repository": sorted(set(escaped)),
        "all_local_links_resolve": not broken and not escaped,
    }


def synthesis_scope_audit() -> dict[str, Any]:
    report_path = (
        ROOT
        / "results/unconditional_ns_chow_encoding_and_maximal_chamber_finiteness.md"
    )
    report = report_path.read_text(encoding="utf-8")
    flattened_report = re.sub(r"\s+", " ", report)
    required_markers = (
        "standard admissible Vinberg lattice-form model",
        "No maximality is needed for Theorem 1.1",
        "projected wall normals",
        "raw ambient trace",
        "The first does **not** prove the second",
        "does not invoke the ABSW finiteness conclusion",
        "finite Vieta boxes",
        "not an AG-only proof of finiteness",
    )
    marker_checks = {
        marker: marker in flattened_report for marker in required_markers
    }

    wrapper = read_json(
        "data/enumerations/"
        "unconditional_ns_chow_encoding_and_maximal_chamber_finiteness.json"
    )
    payload = wrapper["payload"]
    weak = payload["unconditional_weak_carrier"]
    maximal = payload["unconditional_maximal_chambers"]
    provenance = payload["provenance_boundary"]
    dependencies = payload["dependency_sha256"]

    machine_checks = {
        "owner_payload_digest_valid": wrapper["payload_sha256"] == digest(payload),
        "machine_input_keeps_lattice_form": (
            "admissible Vinberg lattice-form" in weak["input"]
        ),
        "weak_carrier_requires_no_maximality": (
            weak["requires_reflection_maximality"] is False
        ),
        "weak_carrier_is_exact_only_as_eigenslice": (
            weak["exact_semialgebraic_eigenslice_return"] is True
            and weak["effective_or_nef_return"] is False
        ),
        "machine_bounds_are_29_25_750_2997": (
            weak["maximum_dimension_imported"] == DIMENSION_MAX
            and weak["maximum_field_degree_imported"] == FIELD_DEGREE_MAX
            and weak["maximum_trace_rank"] == TRACE_RANK_MAX
            and weak["maximum_blowups"] == BLOWUP_MAX
        ),
        "maximal_input_is_literal": (
            "reflection-maximal" in maximal["input"]
            and maximal["finite_unmarked_normalized_chambers"] is True
        ),
        "classical_route_has_no_carrier_hypothesis": (
            maximal["AG_carrier_hypothesis"] is False
        ),
        "provenance_denies_weak_to_finite": (
            provenance["finiteness_derived_from_weak_carrier"] is False
        ),
        "provenance_marks_classical_imports": (
            provenance["classical_dimension_and_field_inputs_imported"] is True
            and provenance["global_AG_only_finiteness_claimed"] is False
        ),
        "classical_audit_report_pinned_by_owner": (
            dependencies[
                "results/unconditional_maximal_arithmetic_reflection_chamber_finiteness_independent_hostile_audit.md"
            ]
            == CLASSICAL_AUDIT_PINS[
                "results/unconditional_maximal_arithmetic_reflection_chamber_finiteness_independent_hostile_audit.md"
            ]
        ),
        "weak_audit_report_pinned_by_owner": (
            dependencies[
                "results/universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.md"
            ]
            == WEAK_AUDIT_PINS[
                "results/universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.md"
            ]
        ),
    }

    links = local_link_audit(report_path)
    checks = boolean_checks({**marker_checks, **machine_checks})
    return {
        "input": weak["input"],
        "marker_checks": marker_checks,
        "machine_checks": machine_checks,
        "link_audit": links,
        "all_checks_pass": all(checks.values()) and links["all_local_links_resolve"],
    }


def inertia(matrix: sp.Matrix) -> list[int]:
    values = matrix.eigenvals()
    positive = sum(m for value, m in values.items() if value.is_positive is True)
    negative = sum(m for value, m in values.items() if value.is_negative is True)
    zero = matrix.rows - positive - negative
    return [positive, negative, zero]


def reflection(gram: sp.Matrix, index: int) -> sp.Matrix:
    root = sp.eye(gram.rows).col(index)
    square = (root.T * gram * root)[0]
    return sp.eye(gram.rows) - 2 * root * (root.T * gram) / square


def cycle(operator: sp.Matrix, intersection: sp.Matrix) -> sp.Matrix:
    return intersection.inv() * operator.T


def chow_composition_replay() -> dict[str, Any]:
    gram = sp.Matrix([[-2, 0, 1], [0, -2, 2], [1, 2, -2]])
    source = [reflection(gram, index) for index in range(3)]
    injection = sp.Matrix(
        [
            [0, 0, 3],
            [0, 0, -1],
            [4, 0, -2],
            [4, 0, -2],
            [0, 4, -4],
            [0, 4, -4],
        ]
    )
    intersection = sp.diag(1, -1, -1, -1, -1, -1)
    embedded = injection.T * intersection * injection
    projector = injection * embedded.inv() * injection.T * intersection

    def lift(operator: sp.Matrix) -> sp.Matrix:
        return (
            injection
            * operator
            * embedded.inv()
            * injection.T
            * intersection
        )

    lifted = [lift(operator) for operator in source]
    cycles = [cycle(operator, intersection) for operator in lifted]
    projector_cycle = cycle(projector, intersection)
    g = lifted[0]
    h = lifted[2]
    c_g = cycles[0]
    c_h = cycles[2]

    checks = {
        "source_determinant_two": gram.det() == 2,
        "embedding_identity": embedded == 16 * gram,
        "projector_idempotent": projector**2 == projector,
        "projector_is_summand_unit": projector * injection == injection,
        "lifts_restrict_to_source": all(
            lifted[i] * injection == injection * source[i] for i in range(3)
        ),
        "noncommuting_pair": g * h != h * g,
        "cycles_act_as_declared_operators": all(
            cycles[i].T * intersection == lifted[i] for i in range(3)
        ),
        "standard_composition_order": (
            c_h * intersection * c_g == cycle(g * h, intersection)
        ),
        "reversed_order_is_detectably_wrong": (
            c_g * intersection * c_h != cycle(g * h, intersection)
        ),
        "reflection_cycles_self_transpose": all(
            coefficient == coefficient.T for coefficient in cycles
        ),
        "reflection_cycles_square_to_projector_cycle": all(
            coefficient * intersection * coefficient == projector_cycle
            for coefficient in cycles
        ),
        "m2_relation": (lifted[0] * lifted[1]) ** 2 == projector,
        "m3_relation": (lifted[0] * lifted[2]) ** 3 == projector,
        "parabolic_word_nontrivial": lifted[1] * lifted[2] != projector,
        "parabolic_nilpotency_index_three": (
            (lifted[1] * lifted[2] - projector) ** 2 != sp.zeros(6)
            and (lifted[1] * lifted[2] - projector) ** 3 == sp.zeros(6)
        ),
    }
    checks = boolean_checks(checks)
    return {
        "source": "arithmetic (2,3,infinity) triangle",
        "ambient": "I_(1,5)",
        "composition_rule": "C_A J C_B = C_(BA)",
        "standard_action_rule": "Gamma_g o Gamma_h = Gamma_(gh)",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def quadratic_eigenslice_replay() -> dict[str, Any]:
    root_three = sp.sqrt(3)
    evaluation = sp.Matrix(
        [
            [1, root_three, 0, 0, 0, 0],
            [0, 0, 1, root_three, 0, 0],
            [0, 0, 0, 0, 1, root_three],
            [1, -root_three, 0, 0, 0, 0],
            [0, 0, 1, -root_three, 0, 0],
            [0, 0, 0, 0, 1, -root_three],
        ]
    )
    archimedean_gram = sp.diag(root_three, -1, -1, -root_three, -1, -1)
    trace_gram = sp.simplify(evaluation.T * archimedean_gram * evaluation)
    block = sp.Matrix([[0, 3], [1, 0]])
    multiplication = sp.diag(block, block, block)
    diagonal_action = sp.diag(
        root_three,
        root_three,
        root_three,
        -root_three,
        -root_three,
        -root_three,
    )
    plus_projector_arch = sp.diag(1, 1, 1, 0, 0, 0)
    plus_projector = evaluation.inv() * plus_projector_arch * evaluation

    plus_point = sp.Matrix([2, -1, -1, 0, 0, 0])
    trace_point = evaluation.inv() * plus_point
    walls = sp.Matrix.hstack(
        sp.Matrix([0, 0, 1, 0, 0, 0]),
        sp.Matrix([1, 0, 2, 0, 0, 0]),
    )
    arch_walls = evaluation * walls
    trace_pairings = sp.simplify(trace_point.T * trace_gram * walls)
    distinguished_pairings = sp.simplify(
        plus_point[:3, :].T
        * sp.diag(root_three, -1, -1)
        * arch_walls[:3, :]
    )
    raw_wall_gram = sp.simplify(walls.T * trace_gram * walls)
    projected_walls = plus_projector * walls
    projected_wall_gram = sp.simplify(
        projected_walls.T * trace_gram * projected_walls
    )
    distinguished_wall_gram = sp.simplify(
        arch_walls[:3, :].T
        * sp.diag(root_three, -1, -1)
        * arch_walls[:3, :]
    )
    raw_normalized_square = sp.simplify(
        4 * raw_wall_gram[0, 1] ** 2
        / (raw_wall_gram[0, 0] * raw_wall_gram[1, 1])
    )
    projected_normalized_square = sp.simplify(
        4 * projected_wall_gram[0, 1] ** 2
        / (projected_wall_gram[0, 0] * projected_wall_gram[1, 1])
    )

    checks = {
        "trace_gram_integral": all(entry.is_Integer for entry in trace_gram),
        "trace_signature_1_5": inertia(trace_gram) == [1, 5, 0],
        "field_action_splits_by_embeddings": (
            evaluation * multiplication == diagonal_action * evaluation
        ),
        "field_action_trace_self_adjoint": (
            multiplication.T * trace_gram == trace_gram * multiplication
        ),
        "plus_projector_idempotent": plus_projector**2 == plus_projector,
        "plus_projector_is_root_three_eigenspace": (
            multiplication * plus_projector == root_three * plus_projector
        ),
        "distinguished_signature_1_2": (
            inertia(sp.diag(root_three, -1, -1)) == [1, 2, 0]
        ),
        "conjugate_negative_definite": (
            inertia(sp.diag(-root_three, -1, -1)) == [0, 3, 0]
        ),
        "positive_interior_point": (
            plus_point[:3, :].T
            * sp.diag(root_three, -1, -1)
            * plus_point[:3, :]
        )[0]
        > 0,
        "trace_pairings_equal_distinguished_pairings": (
            trace_pairings == distinguished_pairings
        ),
        "oriented_wall_pairings_positive": all(
            value > 0 for value in trace_pairings
        ),
        "projected_wall_gram_is_distinguished_gram": (
            projected_wall_gram == distinguished_wall_gram
        ),
        "raw_trace_gram_differs_from_projected_gram": (
            raw_wall_gram != projected_wall_gram
        ),
        "raw_normalized_gram_differs_from_projected": (
            raw_normalized_square != projected_normalized_square
        ),
    }
    checks = boolean_checks(checks)
    return {
        "field": "Q(sqrt(3))",
        "native_rank": 3,
        "trace_rank": 6,
        "trace_signature": inertia(trace_gram),
        "raw_wall_gram": [
            [str(raw_wall_gram[i, j]) for j in range(2)] for i in range(2)
        ],
        "projected_wall_gram": [
            [str(projected_wall_gram[i, j]) for j in range(2)]
            for i in range(2)
        ],
        "raw_normalized_pair_square": str(raw_normalized_square),
        "projected_normalized_pair_square": str(projected_normalized_square),
        "wall_pairings": [str(value) for value in trace_pairings],
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def simple_cycle_slot_bound(vertices: int) -> int:
    return sum(
        factorial(vertices) // factorial(vertices - length) // length
        for length in range(2, vertices + 1)
    )


def arithmetic_bound_and_box_replay() -> dict[str, Any]:
    cycle_slots = simple_cycle_slot_bound(RANK_MAX)
    vieta = [
        comb(FIELD_DEGREE_MAX, index) * CONJUGATE_HEIGHT**index
        for index in range(1, FIELD_DEGREE_MAX + 1)
    ]
    classical_owner = read_json(
        "data/enumerations/"
        "unconditional_maximal_arithmetic_reflection_chamber_finiteness.json"
    )["payload"]
    classical_audit = read_json(
        "data/enumerations/"
        "unconditional_maximal_arithmetic_reflection_chamber_finiteness_"
        "independent_hostile_audit.json"
    )["payload"]
    owner_box = classical_owner["finite_field_box"]
    audit_box = classical_audit["finite_algebraic_integer_box"]

    checks = {
        "dimension_to_rank": DIMENSION_MAX + 1 == RANK_MAX == 30,
        "degree_rank_product": FIELD_DEGREE_MAX * RANK_MAX == TRACE_RANK_MAX == 750,
        "blowup_arithmetic": 4 * TRACE_RANK_MAX - 3 == BLOWUP_MAX == 2997,
        "cycle_height": CONJUGATE_HEIGHT == 62**30,
        "simple_cycle_slots": (
            cycle_slots == 24893912605687593731774060679329
        ),
        "vieta_vector_has_25_positive_entries": (
            len(vieta) == 25 and all(value > 0 for value in vieta)
        ),
        "owner_box_matches_independent_recomputation": (
            owner_box["uniform_conjugate_height"] == str(CONJUGATE_HEIGHT)
            and owner_box["simple_directed_cycle_slot_bound"] == cycle_slots
            and owner_box["Vieta_coefficient_bounds"]
            == [str(value) for value in vieta]
        ),
        "classical_audit_box_matches_independent_recomputation": (
            audit_box["conjugate_height"] == str(CONJUGATE_HEIGHT)
            and audit_box["uniform_Vieta_bounds"]
            == [str(value) for value in vieta]
            and audit_box["uniform_Vieta_bounds_sha256"]
            == digest([str(value) for value in vieta])
        ),
        "factorization_is_ABSW_conclusion_free": (
            classical_owner["two_terminal_routes"][
                "uses_ABSW_finiteness_conclusion_in_factorization"
            ]
            is False
            and classical_audit["finite_pair_factorization"][
                "uses_ABSW_finiteness_conclusion"
            ]
            is False
        ),
    }
    checks = boolean_checks(checks)
    return {
        "maximum_dimension": DIMENSION_MAX,
        "maximum_native_rank": RANK_MAX,
        "maximum_field_degree": FIELD_DEGREE_MAX,
        "maximum_trace_rank": TRACE_RANK_MAX,
        "maximum_blowups": BLOWUP_MAX,
        "conjugate_height": str(CONJUGATE_HEIGHT),
        "simple_cycle_slots": cycle_slots,
        "Vieta_bounds_sha256": digest([str(value) for value in vieta]),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def finite_gallery_calibration() -> dict[str, Any]:
    gram = sp.Matrix([[2, -1], [-1, 2]])
    generators = [reflection(gram, index) for index in range(2)]
    identity = sp.eye(2)

    def key(matrix: sp.Matrix) -> tuple[str, ...]:
        return tuple(str(entry) for entry in matrix)

    group: dict[tuple[str, ...], sp.Matrix] = {key(identity): identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = current * generator
            candidate_key = key(candidate)
            if candidate_key not in group:
                group[candidate_key] = candidate
                frontier.append(candidate)

    interior = gram.inv() * sp.ones(2, 1)
    positive_roots = (
        sp.Matrix([1, 0]),
        sp.Matrix([0, 1]),
        sp.Matrix([1, 1]),
    )
    sign_patterns = {
        tuple(
            int(sp.sign((matrix * interior).T.dot(gram * root)))
            for root in positive_roots
        )
        for matrix in group.values()
    }
    checks = {
        "two_reflections_are_involutions": all(
            generator**2 == identity for generator in generators
        ),
        "coxeter_m3_relation": (generators[0] * generators[1]) ** 3 == identity,
        "A2_group_has_six_elements": len(group) == 6,
        "regular_interior_has_no_zero_pairing": all(
            all(sign != 0 for sign in pattern) for pattern in sign_patterns
        ),
        "one_orbit_meets_all_six_chambers": len(sign_patterns) == 6,
    }
    checks = boolean_checks(checks)
    return {
        "calibration": "finite A2 complete mirror arrangement",
        "group_order": len(group),
        "chamber_sign_patterns": [list(pattern) for pattern in sorted(sign_patterns)],
        "general_bridge": [
            "finite-covolume discrete reflection overgroups have finite index",
            "finite-index commensurability preserves arithmeticity",
            "reflection maximality forces W(P)=R(L)",
            "local finiteness plus a generic compact path gives a finite gallery",
            "complete chambers form one reflection-group orbit",
        ],
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def reachable(edges: list[tuple[str, str]], start: str, target: str) -> bool:
    adjacency: dict[str, set[str]] = {}
    for left, right in edges:
        adjacency.setdefault(left, set()).add(right)
    seen = {start}
    frontier = [start]
    while frontier:
        node = frontier.pop()
        for successor in adjacency.get(node, set()):
            if successor == target:
                return True
            if successor not in seen:
                seen.add(successor)
                frontier.append(successor)
    return start == target


def provenance_replay() -> dict[str, Any]:
    edges = [
        ("admissible_lattice_form_input", "weak_chow_carrier"),
        ("weak_chow_carrier", "four_R_minus_three_formula"),
        ("four_R_minus_three_formula", "uniform_2997_carrier"),
        ("classical_dimension_bound", "uniform_2997_carrier"),
        ("classical_degree_bound", "uniform_2997_carrier"),
        ("classical_dimension_bound", "finite_native_pairs"),
        ("classical_degree_bound", "finite_native_pairs"),
        ("bounded_integral_frame_cycles", "finite_native_pairs"),
        ("finite_native_pairs", "finite_native_similarity_hosts"),
        ("Nikulin_fixed_data", "finite_native_similarity_hosts"),
        ("finite_native_similarity_hosts", "finite_maximal_chambers"),
        ("reflection_maximality", "finite_maximal_chambers"),
        ("complete_mirror_gallery", "finite_maximal_chambers"),
        ("ABSW_or_Nikulin_2007", "finite_maximal_chambers"),
    ]
    requirements = {
        "uniform_2997_carrier": [
            "four_R_minus_three_formula",
            "classical_dimension_bound",
            "classical_degree_bound",
        ],
        "factorized_finite_maximal_chambers": [
            "finite_native_pairs",
            "Nikulin_fixed_data",
            "reflection_maximality",
            "complete_mirror_gallery",
        ],
        "direct_finite_maximal_chambers": ["ABSW_or_Nikulin_2007"],
    }
    checks = {
        "weak_carrier_reaches_its_embedding_formula": reachable(
            edges, "weak_chow_carrier", "four_R_minus_three_formula"
        ),
        "classical_bounds_feed_uniform_2997": reachable(
            edges, "classical_dimension_bound", "uniform_2997_carrier"
        )
        and reachable(edges, "classical_degree_bound", "uniform_2997_carrier"),
        "weak_carrier_does_not_reach_native_pairs": not reachable(
            edges, "weak_chow_carrier", "finite_native_pairs"
        ),
        "weak_carrier_does_not_reach_similarity_hosts": not reachable(
            edges, "weak_chow_carrier", "finite_native_similarity_hosts"
        ),
        "weak_carrier_does_not_reach_chamber_finiteness": not reachable(
            edges, "weak_chow_carrier", "finite_maximal_chambers"
        ),
        "factorized_finish_requires_maximality": (
            "reflection_maximality"
            in requirements["factorized_finite_maximal_chambers"]
        ),
        "direct_classical_finish_has_no_carrier": (
            "weak_chow_carrier"
            not in requirements["direct_finite_maximal_chambers"]
        ),
    }
    checks = boolean_checks(checks)
    return {
        "edges": [list(edge) for edge in edges],
        "requirements": requirements,
        "forbidden_edge": ["weak_chow_carrier", "finite_maximal_chambers"],
        "finiteness_derived_from_weak_carrier": False,
        "uniform_2997_imports_classical_bounds": True,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def branch_archive_cross_checks() -> dict[str, bool]:
    weak_owner = read_json(
        "data/enumerations/universal_weak_chow_action_carrier.json"
    )
    weak_audit = read_json(
        "data/enumerations/"
        "universal_weak_chow_action_carrier_and_leverage_no_go_"
        "independent_hostile_audit.json"
    )
    classical_audit = read_json(
        "data/enumerations/"
        "unconditional_maximal_arithmetic_reflection_chamber_finiteness_"
        "independent_hostile_audit.json"
    )
    checks = {
        "weak_owner_passes": weak_owner["all_checks_pass"] is True,
        "weak_audit_digest_valid": (
            weak_audit["payload_sha256"] == digest(weak_audit["payload"])
        ),
        "weak_audit_passes": weak_audit["payload"]["all_checks_pass"] is True,
        "classical_audit_digest_valid": (
            classical_audit["payload_sha256"]
            == digest(classical_audit["payload"])
        ),
        "classical_audit_verdict_is_classical_not_AG_only": (
            classical_audit["payload"]["verdict"]
            == "PASS_UNCONDITIONAL_CLASSICAL_NOT_AG_ONLY"
        ),
    }
    return boolean_checks(checks)


def build_payload() -> dict[str, Any]:
    validate_optional_primary_sources()
    observed_pins = validate_input_pins()
    scope = synthesis_scope_audit()
    chow = chow_composition_replay()
    eigenslice = quadratic_eigenslice_replay()
    arithmetic = arithmetic_bound_and_box_replay()
    gallery = finite_gallery_calibration()
    provenance = provenance_replay()
    branch_checks = branch_archive_cross_checks()
    report_links = local_link_audit(AUDIT_REPORT)

    checks = {
        "all_21_terminal_input_bytes_match": len(ALL_INPUT_PINS) == 21,
        "synthesis_scope_passes": scope["all_checks_pass"],
        "chow_composition_passes": chow["all_checks_pass"],
        "quadratic_eigenslice_passes": eigenslice["all_checks_pass"],
        "arithmetic_bound_and_box_passes": arithmetic["all_checks_pass"],
        "finite_gallery_calibration_passes": gallery["all_checks_pass"],
        "provenance_separation_passes": provenance["all_checks_pass"],
        "branch_archives_cross_check": all(branch_checks.values()),
        "audit_report_links_resolve": report_links["all_local_links_resolve"],
    }
    checks = boolean_checks(checks)

    source_payload = {
        name: {
            "sha256": record["sha256"],
            "locus": record["locus"],
            "temporary_path": record["path"],
            "required_for_clean_replay": False,
        }
        for name, record in PRIMARY_SOURCE_PINS.items()
    }
    return {
        "schema": SCHEMA,
        "verdict": "PASS_AFTER_SCOPE_REPAIRS",
        "terminal_input_sha256": observed_pins,
        "audit_report_sha256": AUDIT_REPORT_SHA256,
        "repairs_verified": [
            "weak all-input prose retains the admissible Vinberg lattice-form model",
            "projected wall-normal Gram is distinguished from raw ambient trace Gram",
            "machine archive input retains the lattice-form qualifier",
        ],
        "literal_scope_audit": scope,
        "independent_chow_composition_replay": chow,
        "independent_quadratic_eigenslice_replay": eigenslice,
        "independent_arithmetic_bound_and_box_replay": arithmetic,
        "independent_group_to_chamber_calibration": gallery,
        "provenance_separation_replay": provenance,
        "branch_archive_cross_checks": branch_checks,
        "primary_source_pins": source_payload,
        "replay_boundary": {
            "imports_synthesis_owner_producer": False,
            "imports_classical_owner_or_audit_producer": False,
            "imports_weak_owner_or_audit_producer": False,
            "published_literature_machine_proved": False,
            "temporary_primary_sources_required_for_clean_replay": False,
        },
        "audit_report_link_check": report_links,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def encoded_audit() -> str:
    payload = build_payload()
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-output", action="store_true")
    parser.add_argument("--check-output", action="store_true")
    args = parser.parse_args()

    encoded = encoded_audit()
    if args.write_output:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(encoded, encoding="utf-8")
        print(OUTPUT)
        return
    if args.check_output:
        if not OUTPUT.is_file() or OUTPUT.read_text(encoding="utf-8") != encoded:
            raise SystemExit("frozen independent synthesis audit differs from replay")
        print("unconditional NS/Chow synthesis independent hostile audit: PASS")
        return
    print(encoded, end="")


if __name__ == "__main__":
    main()
