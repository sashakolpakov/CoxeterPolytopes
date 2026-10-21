#!/usr/bin/env python3
"""Independent hostile audit of unconditional arithmetic chamber finiteness."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import comb
from pathlib import Path
import re
import sys
from typing import Any, Iterable


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


ROOT = Path(__file__).resolve().parents[1]
OWNER_REPORT = (
    ROOT
    / "results/unconditional_maximal_arithmetic_reflection_chamber_finiteness.md"
)
OWNER_SCRIPT = (
    ROOT
    / "scripts/certify_unconditional_maximal_arithmetic_reflection_chamber_finiteness.py"
)
OWNER_ARCHIVE = (
    ROOT
    / "data/enumerations/unconditional_maximal_arithmetic_reflection_chamber_finiteness.json"
)
OWNER_TEST = (
    ROOT
    / "tests/test_unconditional_maximal_arithmetic_reflection_chamber_finiteness.py"
)
SCOPE_NO_GO_REPORT = (
    ROOT / "results/unconditional_ag_absw_universality_leverage_no_go.md"
)
AUDIT_REPORT = (
    ROOT
    / "results/unconditional_maximal_arithmetic_reflection_chamber_finiteness_independent_hostile_audit.md"
)
AUDIT_ARCHIVE = (
    ROOT
    / "data/enumerations/unconditional_maximal_arithmetic_reflection_chamber_finiteness_independent_hostile_audit.json"
)


OWNER_SHA256 = {
    "report": "6eee91512840287ce47a195a317fa4bbfd19678a95b52869ec0a45775d61e67b",
    "script": "c32c7969f900093db388718f3ac769f669d42e43b6c2a95c76297eb34899f541",
    "archive": "2bce48c0a7156caea3f4b4890411d770450912b4d4b2d3827ce93cc7e9730acc",
    "test": "4105f4f3d08dfdcbc1d98f272aae6168a99221f6e426f1bcf1fb5b50303d6702",
}
SCOPE_NO_GO_SHA256 = (
    "fcc3acb8bee4751bef083cd00a179f9b542602561ba1bf81d89f374dd5523150"
)
AUDIT_REPORT_SHA256 = (
    "1f1592ef207580b765d56bdca9461429456ef9fa46eb180b0ae65d69dc68ef3c"
)


DIMENSION_MAX = 29
RANK_MAX = 30
FIELD_DEGREE_MAX = 25
FRAME_CONSTANT = 62
CONJUGATE_HEIGHT = FRAME_CONSTANT**RANK_MAX
BLOWUP_MAX = 4 * FIELD_DEGREE_MAX * RANK_MAX - 3


SOURCE_PINS = {
    "absw_2008": {
        "path": "/private/tmp/absw.pdf",
        "sha256": "974ac7c5141ca1c45b515f145215ad115aba867e42551a9b103233ea52930ee4",
        "locator": "Introduction definition and Theorem 6.1",
    },
    "vinberg_1971": {
        "path": "/private/tmp/vinberg_1971_definition_fields.pdf",
        "sha256": "16a4f497688e6f968b4c63d81a6b049b3f096f1cd4f43e458b2882a5c44ded38",
        "locator": "Theorems 1, 3, and 5, printed pp. 45-46 and 55",
    },
    "vinberg_dimension": {
        "path": "/private/tmp/vinberg_absence_complete.pdf",
        "sha256": "fd25f8f892eed27e4303a8f3252fd384d8954587a6c1c67f1efed2c38aabe44f",
        "locator": "Theorems 1 and 4, printed p. 69",
    },
    "nikulin_2007": {
        "path": "/private/tmp/nikulin_finiteness.pdf",
        "sha256": "9a820fef0b88158b283057f7d97acc89130610a9a4aacc6c53a8215e23ce874a",
        "locator": "Theorem 2.2, printed p. 56",
    },
    "nikulin_1982": {
        "path": "/private/tmp/nikulin_1981_eng.pdf",
        "sha256": "696652b6a67c6ed420de6c1cb61d9422e3619e28e081a78131a5323ce2edee9d",
        "locator": "Theorem 2.1.3(2), p. 104; Theorem 5.2.1, pp. 119-122",
    },
    "nikulin_2011": {
        "path": "/private/tmp/nikulin_transition_journal.pdf",
        "sha256": "716e6e7333efbfaa22060b60dc97e8303363a4914409ccdd56ef7a90c923c7af",
        "locator": "Theorem 3.2 and Theorem 5.5",
    },
    "belolipetsky_linowitz_v2": {
        "path": "/tmp/belolipetsky_linowitz_2013.pdf",
        "sha256": "026e0f635f71e214807e0b505bf7297673eb9e0d47a6a91b7bd6dec134eab97f",
        "locator": "Corollary 5.1, PDF p. 8",
    },
    "dolgachev_esselmann_secondary": {
        "path": "/tmp/dolgachev_reflection_groups_ag.pdf",
        "sha256": "9f29e1e87eee55f1d66433ccbabb5be0833d3c691f9e67a1643e63f7bf44fc26",
        "locator": "Theorem 4.7",
    },
}


NONINTEGRAL_CARTAN = (
    (Fraction(2), Fraction(1, 2), Fraction(1, 6)),
    (Fraction(2), Fraction(2), Fraction(1, 3)),
    (Fraction(6), Fraction(3), Fraction(2)),
)
SYMMETRIZER = (Fraction(36), Fraction(9), Fraction(1))


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def simple_cycle_slot_bound(vertices: int) -> int:
    """Count ordered simple directed cycles modulo cyclic rotation."""

    falling = 1
    total = 0
    for length in range(1, vertices + 1):
        falling *= vertices - length + 1
        if length >= 2:
            total += falling // length
    return total


def coefficient_bounds(max_degree: int, height: int) -> list[int]:
    """Uniform Vieta bounds valid for every degree at most max_degree."""

    return [
        comb(max_degree, index) * height**index
        for index in range(1, max_degree + 1)
    ]


def coefficient_box_union_upper_bound(max_degree: int, height: int) -> int:
    """Upper bound for monic integral coefficient tuples of degrees <= d."""

    total = 0
    for degree in range(1, max_degree + 1):
        possibilities = 1
        for index in range(1, degree + 1):
            bound = comb(degree, index) * height**index
            possibilities *= 2 * bound + 1
        total += possibilities
    return total


def cartan_cycle(indices: Iterable[int]) -> Fraction:
    cycle = tuple(indices)
    if len(cycle) < 2:
        raise ValueError("cycle needs at least two vertices")
    return_value = Fraction(1)
    for index, vertex in enumerate(cycle):
        return_value *= NONINTEGRAL_CARTAN[vertex][cycle[(index + 1) % len(cycle)]]
    return return_value


def cyclic_integrality_countercontrol() -> dict[str, Any]:
    size = len(NONINTEGRAL_CARTAN)
    for i in range(size):
        for j in range(size):
            if SYMMETRIZER[i] * NONINTEGRAL_CARTAN[i][j] != (
                SYMMETRIZER[j] * NONINTEGRAL_CARTAN[j][i]
            ):
                raise AssertionError("Cartan sample is not symmetrizable")

    nonintegral_entries = sorted(
        {
            str(NONINTEGRAL_CARTAN[i][j])
            for i in range(size)
            for j in range(size)
            if i != j and NONINTEGRAL_CARTAN[i][j].denominator != 1
        }
    )
    two_cycles = [
        NONINTEGRAL_CARTAN[i][j] * NONINTEGRAL_CARTAN[j][i]
        for i in range(size)
        for j in range(i + 1, size)
    ]
    forward = cartan_cycle((0, 1, 2))
    reverse = cartan_cycle((0, 2, 1))
    edge_product = two_cycles[0] * two_cycles[1] * two_cycles[2]
    if any(value.denominator != 1 for value in two_cycles):
        raise AssertionError("two-cycle product lost integrality")
    if forward.denominator != 1 or reverse.denominator != 1:
        raise AssertionError("three-cycle product lost integrality")
    if forward != reverse or forward * forward != edge_product:
        raise AssertionError("b_C^2 identity failed")

    return {
        "symmetrizer": [str(value) for value in SYMMETRIZER],
        "nonintegral_individual_entries": nonintegral_entries,
        "integral_two_cycles": [str(value) for value in two_cycles],
        "forward_three_cycle": str(forward),
        "reverse_three_cycle": str(reverse),
        "cycle_square": str(forward * forward),
        "product_of_two_cycles": str(edge_product),
    }


def directed_edges(walk: tuple[int, ...]) -> list[tuple[int, int]]:
    return list(zip(walk, walk[1:]))


def decompose_closed_walk(walk: tuple[int, ...]) -> list[tuple[int, ...]]:
    """Split a loop-free closed walk into simple directed cycles."""

    if len(walk) < 3 or walk[0] != walk[-1]:
        raise ValueError("expected a nontrivial closed walk")
    if any(left == right for left, right in directed_edges(walk)):
        raise ValueError("diagonal steps are excluded")

    pending = walk
    cycles: list[tuple[int, ...]] = []
    while len(pending) > 1:
        first_position: dict[int, int] = {}
        repeated: tuple[int, int] | None = None
        for right, vertex in enumerate(pending):
            if vertex in first_position:
                repeated = (first_position[vertex], right)
                break
            first_position[vertex] = right
        if repeated is None:
            raise AssertionError("closed walk unexpectedly has no repeated vertex")
        left, right = repeated
        cycle = pending[left : right + 1]
        if len(cycle) < 3 or len(set(cycle[:-1])) != len(cycle) - 1:
            raise AssertionError("decomposition emitted a nonsimple cycle")
        cycles.append(cycle)
        pending = pending[: left + 1] + pending[right + 1 :]
    return cycles


def closed_walk_decomposition_replay(
    vertices: int = 4, maximum_edges: int = 7
) -> dict[str, int]:
    walks = 0
    cycles = 0
    maximum_components = 0
    for start in range(vertices):
        for edge_count in range(2, maximum_edges + 1):
            for interior in product(range(vertices), repeat=edge_count - 1):
                walk = (start,) + interior + (start,)
                if any(left == right for left, right in directed_edges(walk)):
                    continue
                pieces = decompose_closed_walk(walk)
                original_edges = Counter(directed_edges(walk))
                decomposed_edges: Counter[tuple[int, int]] = Counter()
                for piece in pieces:
                    decomposed_edges.update(directed_edges(piece))
                if decomposed_edges != original_edges:
                    raise AssertionError("closed-walk edge product did not factor")
                walks += 1
                cycles += len(pieces)
                maximum_components = max(maximum_components, len(pieces))
    return {
        "vertices": vertices,
        "maximum_edges": maximum_edges,
        "closed_walks_checked": walks,
        "simple_cycle_factors_checked": cycles,
        "maximum_factors_in_one_walk": maximum_components,
    }


def validate_optional_primary_sources() -> None:
    for record in SOURCE_PINS.values():
        path = Path(record["path"])
        if path.exists() and sha256_file(path) != record["sha256"]:
            raise AssertionError(f"audit-time source hash changed: {path}")


def validate_owner_and_scope() -> dict[str, Any]:
    owner_paths = {
        "report": OWNER_REPORT,
        "script": OWNER_SCRIPT,
        "archive": OWNER_ARCHIVE,
        "test": OWNER_TEST,
    }
    for key, path in owner_paths.items():
        actual = sha256_file(path)
        if actual != OWNER_SHA256[key]:
            raise AssertionError(f"owner {key} hash changed: {actual}")
    if sha256_file(SCOPE_NO_GO_REPORT) != SCOPE_NO_GO_SHA256:
        raise AssertionError("AG universality-leverage no-go bytes changed")
    if sha256_file(AUDIT_REPORT) != AUDIT_REPORT_SHA256:
        raise AssertionError("independent audit report bytes changed")

    owner_text = OWNER_REPORT.read_text(encoding="utf-8")
    required_owner = (
        "Theorem 1.1 (unconditional maximal-chamber finiteness)",
        "Individual Cartan coefficients need not be elements",
        "every cyclic product",
        "Theorem 2.1.3(2)",
        "finite-index commensurability invariant",
        "k=\\mathbf Q(b_C:C\\text{ a simple frame cycle})",
        "Theorem 5.2.1",
        "W(P)=R(L)",
        "not an AG-relative result",
    )
    for needle in required_owner:
        if needle not in owner_text:
            raise AssertionError(f"owner scope marker missing: {needle}")
    if re.search(r"(?<!\\)qquad", owner_text):
        raise AssertionError("owner report has a stripped TeX command")

    no_go_text = SCOPE_NO_GO_REPORT.read_text(encoding="utf-8")
    required_no_go = (
        "weak carriers are universal and therefore give no",
        "A prime wall cannot itself be its regular reflection action",
        "no unconditional AG-only ABSW theorem may",
        "unconditional classical chamber theorem",
    )
    for needle in required_no_go:
        if needle not in no_go_text:
            raise AssertionError(f"AG scope marker missing: {needle}")

    archived = json.loads(OWNER_ARCHIVE.read_text(encoding="utf-8"))
    if archived["payload_sha256"] != digest(archived["payload"]):
        raise AssertionError("owner archive payload digest mismatch")
    return archived["payload"]


def build_certificate() -> dict[str, Any]:
    validate_optional_primary_sources()
    owner = validate_owner_and_scope()

    if DIMENSION_MAX + 1 != RANK_MAX:
        raise AssertionError("dimension/rank conversion changed")
    if CONJUGATE_HEIGHT != 62**30:
        raise AssertionError("cycle-height ceiling changed")
    if BLOWUP_MAX != 2997:
        raise AssertionError("weak rational-surface bound changed")

    cycle_slots = simple_cycle_slot_bound(RANK_MAX)
    if cycle_slots != 24893912605687593731774060679329:
        raise AssertionError("rank-30 simple-cycle slot count changed")

    bounds = coefficient_bounds(FIELD_DEGREE_MAX, CONJUGATE_HEIGHT)
    owner_box = owner["finite_field_box"]
    if owner_box["frame_constant"] != FRAME_CONSTANT:
        raise AssertionError("owner frame constant changed")
    if owner_box["maximum_cycle_length"] != RANK_MAX:
        raise AssertionError("owner cycle-length bound changed")
    if owner_box["uniform_conjugate_height"] != str(CONJUGATE_HEIGHT):
        raise AssertionError("owner conjugate-height bound changed")
    if owner_box["maximum_element_degree"] != FIELD_DEGREE_MAX:
        raise AssertionError("owner degree bound changed")
    if owner_box["simple_directed_cycle_slot_bound"] != cycle_slots:
        raise AssertionError("owner cycle-slot count differs")
    if owner_box["Vieta_coefficient_bounds"] != [str(value) for value in bounds]:
        raise AssertionError("owner Vieta vector differs")

    union_bound = coefficient_box_union_upper_bound(
        FIELD_DEGREE_MAX, CONJUGATE_HEIGHT
    )
    union_bound_text = str(union_bound)
    walk_replay = closed_walk_decomposition_replay()
    cyclic_control = cyclic_integrality_countercontrol()

    source_payload = {
        key: {
            "sha256": record["sha256"],
            "locator": record["locator"],
            "temporary_path": record["path"],
            "required_for_clean_replay": False,
        }
        for key, record in SOURCE_PINS.items()
    }

    return {
        "schema": (
            "coxeter.unconditional-maximal-arithmetic-reflection-chamber-"
            "finiteness.independent-hostile-audit.v1"
        ),
        "verdict": "PASS_UNCONDITIONAL_CLASSICAL_NOT_AG_ONLY",
        "owner_sha256": OWNER_SHA256,
        "scope_no_go_sha256": SCOPE_NO_GO_SHA256,
        "audit_report_sha256": AUDIT_REPORT_SHA256,
        "terminal_scope": {
            "all_finite_volume_arithmetic_inputs_dimension_at_least_two": True,
            "maximum_dimension": DIMENSION_MAX,
            "maximum_ground_field_degree": FIELD_DEGREE_MAX,
            "finite_ground_field_set": True,
            "finite_native_pair_set": True,
            "finite_unmarked_normalized_chambers_for_reflection_maximal_inputs": True,
            "direct_ABSW_or_Nikulin_2007_route": True,
            "AG_carrier_existence_hypothesis": False,
        },
        "noncompact_refinement": {
            "ground_field": "Q",
            "dimensions": "n<=19 or n=21",
            "cusp_to_k_isotropy": True,
            "nondistinguished_definiteness_forces_degree_one": True,
            "Esselmann_integral_reflective_lattice_theorem_imported": True,
            "owner_machine_certificate_claims_to_replay_this_import": False,
        },
        "finite_pair_factorization": {
            "uses_ABSW_finiteness_conclusion": False,
            "reflective_host_chain": "W(P) <= R(L) <= O+(L) <= O(L)",
            "bounded_frame_taken_from_complete_RL_chamber": True,
            "individual_Cartan_entry_integrality_required": False,
            "cyclic_Cartan_product_integrality_required": True,
            "cyclic_integrality_countercontrol": cyclic_control,
            "cycle_square_identity": "b_C^2=product(c_ij along C)",
            "all_conjugates_bound": "|sigma(b_C)|<=62^m<=62^30",
            "full_rank_frame_ground_field_inheritance": True,
            "finite_index_minimal_Vinberg_field_invariance": True,
            "field_generation": "k=Q(simple cyclic frame products)",
            "closed_walk_decomposition_replay": walk_replay,
            "maximum_simple_cycle_slots": cycle_slots,
            "fixed_data_invoked_only_after_finite_native_pairs": True,
        },
        "finite_algebraic_integer_box": {
            "maximum_degree": FIELD_DEGREE_MAX,
            "conjugate_height": str(CONJUGATE_HEIGHT),
            "uniform_Vieta_bounds": [str(value) for value in bounds],
            "uniform_Vieta_bounds_sha256": digest([str(value) for value in bounds]),
            "degree_1_through_25_coefficient_box_union_upper_bound": union_bound_text,
            "coefficient_box_union_upper_bound_decimal_digits": len(union_bound_text),
            "coefficient_box_union_upper_bound_sha256": hashlib.sha256(
                union_bound_text.encode("ascii")
            ).hexdigest(),
            "finite": True,
        },
        "fixed_data_and_chambers": {
            "Nikulin_5_2_1_fixed_field_fixed_rank_similarity_classes": True,
            "projective_modules_in_scope": True,
            "similarity_preserves_reflections_and_normalized_Gram": True,
            "reflection_maximality_forces_W_equals_R": True,
            "arithmetic_only_and_absolute_reflection_maximality_agree": True,
            "complete_mirror_chambers_one_group_orbit": True,
        },
        "weak_rational_surface_corollary": {
            "maximum_blowups": BLOWUP_MAX,
            "bound_source": "imported arithmetic n and field-degree bounds",
            "effective_or_nef_return_claimed": False,
            "genuine_reflection_automorphisms_claimed": False,
            "independent_AG_finiteness_leverage": False,
        },
        "AG_scope": {
            "unconditional_weak_numerical_carrier_exists": True,
            "weak_carrier_nonemptiness_restricts_inputs": False,
            "universal_strong_intrinsic_carrier_proved": False,
            "positive_AG_only_ABSW_proof_claimed": False,
            "only_presently_closed_positive_global_route": "classical arithmetic",
        },
        "source_pins": source_payload,
        "replay_boundary": {
            "owner_and_scope_report_bytes_verified": True,
            "owner_archive_payload_digest_verified": True,
            "finite_field_box_recomputed_independently": True,
            "temporary_source_hashes_checked_when_present": True,
            "temporary_source_files_required": False,
            "published_literature_machine_proved": False,
        },
        "logic": [
            "Vinberg n<30 and Belolipetsky-Linowitz degree<=25 are imported",
            "complete reflective chamber -> connected full-rank frame with A=62",
            "integral cyclic products plus conjugate bounds -> finite alphabet",
            "field inheritance and finite-index invariance -> alphabet generates k",
            "bounded simple-cycle tuple length -> finitely many fields and native pairs",
            "Nikulin fixed-(k,m) theorem -> finitely many native similarity hosts",
            "reflection maximality and chamber transitivity -> finite chamber classes",
        ],
    }


def encoded_certificate() -> str:
    payload = build_certificate()
    wrapped = {"payload": payload, "payload_sha256": digest(payload)}
    return json.dumps(wrapped, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-output", action="store_true")
    parser.add_argument("--write-output", action="store_true")
    args = parser.parse_args()

    encoded = encoded_certificate()
    if args.check_output:
        if not AUDIT_ARCHIVE.exists() or AUDIT_ARCHIVE.read_text(
            encoding="utf-8"
        ) != encoded:
            raise SystemExit("frozen independent certificate differs from exact replay")
        print("unconditional arithmetic chamber independent hostile audit: PASS")
        return
    if args.write_output:
        AUDIT_ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        AUDIT_ARCHIVE.write_text(encoded, encoding="utf-8")
        print(AUDIT_ARCHIVE)
        return
    print(encoded, end="")


if __name__ == "__main__":
    main()
