#!/usr/bin/env python3
"""Certify three single doublings of the Jacobian--Kummer H12 chamber.

The frozen JK-H12 chamber has 172 facets and is already certified to have
finite volume and cusps.  A facet is a doubling wall when every finite angle
at that facet has even Coxeter order.  This script enumerates the complete
doubling-wall stock, constructs three representatives by exact reflection,
and proves every displayed boundary inequality is an irredundant integral
Coxeter facet on the same full Picard lattice.

Finite volume and noncompactness are inherited from the exact two-copy
convex gluing.  No direct enumeration of the enormous descendant vertex
sets is needed for that implication.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "data/enumerations/jacobian_kummer_h12_full_picard_k3_carrier.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/jacobian_kummer_h12_single_doublings.json"
)

SOURCE_SHA256 = (
    "f2c57b9fb0423560c9e238541b0ac1a81ca141917c864929501ebdfa7f294ff2"
)
SCHEMA = "coxeter.jacobian_kummer.h12_single_doublings.v1"
REPRESENTATIVES = (30, 60, 76)  # zero-based: walls 31, 61, 77


Vector = tuple[Fraction, ...]
IntegralVector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fraction(value: int | Fraction) -> Fraction:
    return value if isinstance(value, Fraction) else Fraction(value)


def _payload_number(value: int | Fraction) -> int | str:
    rational = _fraction(value)
    return rational.numerator if rational.denominator == 1 else str(rational)


def _vector_payload(vector: Sequence[int | Fraction]) -> list[int | str]:
    return [_payload_number(value) for value in vector]


def _matrix(payload: Sequence[Sequence[int]]) -> Matrix:
    return tuple(tuple(int(value) for value in row) for row in payload)


def _vector(payload: Sequence[int]) -> IntegralVector:
    return tuple(int(value) for value in payload)


def _covector(gram: Matrix, vector: Sequence[int | Fraction]) -> Vector:
    return tuple(
        sum((_fraction(gram[row][column]) * _fraction(vector[column])
             for column in range(len(vector))), Fraction(0))
        for row in range(len(vector))
    )


def _pair_from_covector(
    left: Sequence[int | Fraction], right_covector: Sequence[int | Fraction]
) -> Fraction:
    return sum(
        (_fraction(value) * _fraction(covalue)
         for value, covalue in zip(left, right_covector, strict=True)),
        Fraction(0),
    )


def _pair(
    gram: Matrix,
    left: Sequence[int | Fraction],
    right: Sequence[int | Fraction],
) -> Fraction:
    return _pair_from_covector(left, _covector(gram, right))


def _content(vector: Iterable[int]) -> int:
    return math.gcd(*(abs(value) for value in vector if value))


def _primitive(vector: Sequence[int | Fraction]) -> IntegralVector:
    values = tuple(_fraction(value) for value in vector)
    denominator = math.lcm(*(value.denominator for value in values))
    integral = tuple(int(value * denominator) for value in values)
    content = _content(integral)
    return tuple(value // content for value in integral)


def _unoriented_key(vector: Sequence[int | Fraction]) -> IntegralVector:
    primitive = _primitive(vector)
    first = next(value for value in primitive if value)
    return primitive if first > 0 else tuple(-value for value in primitive)


def _reflection_apply(
    gram: Matrix,
    root: Sequence[int | Fraction],
    vector: Sequence[int | Fraction],
) -> Vector:
    root_square = _pair(gram, root, root)
    coefficient = 2 * _pair(gram, vector, root) / root_square
    return tuple(
        _fraction(value) - coefficient * _fraction(root_value)
        for value, root_value in zip(vector, root, strict=True)
    )


def _reflection_integral(gram: Matrix, root: IntegralVector) -> bool:
    square = _pair(gram, root, root)
    covector = _covector(gram, root)
    return square < 0 and all(
        (2 * value / square).denominator == 1 for value in covector
    )


def _relation(
    left: Sequence[int | Fraction],
    left_covector: Sequence[int | Fraction],
    right: Sequence[int | Fraction],
    right_covector: Sequence[int | Fraction],
) -> tuple[str, Fraction]:
    product = _pair_from_covector(left, right_covector)
    if product < 0:
        raise ArithmeticError("the doubled chamber is not acute-angled")
    ratio = product * product / (
        _pair_from_covector(left, left_covector)
        * _pair_from_covector(right, right_covector)
    )
    finite = {
        Fraction(0): "m2",
        Fraction(1, 4): "m3",
        Fraction(1, 2): "m4",
        Fraction(3, 4): "m6",
        Fraction(1): "parallel",
    }
    if ratio in finite:
        return finite[ratio], ratio
    if ratio > 1:
        return "ultraparallel", ratio
    raise ArithmeticError(f"non-Coxeter normalized product {ratio}")


def _load_source() -> tuple[dict[str, Any], Matrix, tuple[IntegralVector, ...], IntegralVector]:
    actual = _sha256(SOURCE)
    if actual != SOURCE_SHA256:
        raise ArithmeticError(f"pinned JK-H12 carrier changed: {actual}")
    payload = json.loads(SOURCE.read_text(encoding="utf-8"))
    chamber = payload["chamber"]
    gram = _matrix(chamber["carrier"]["gram"])
    wall_records = chamber["oriented_wall_system"]["wall_records"]
    roots = tuple(_vector(record["oriented_coordinates"]) for record in wall_records)
    interior = _vector(
        chamber["oriented_wall_system"]["strict_interior_coordinates"]
    )
    if not (
        payload["all_checks_pass"]
        and chamber["dimension"] == 12
        and chamber["facet_count"] == 172
        and chamber["finite_volume"]
        and gram == tuple(tuple(row) for row in gram)
        and exact_inertia(sp.Matrix(gram)) == (1, 12, 0)
        and len(roots) == 172
        and all(_content(root) == 1 for root in roots)
        and all(_pair(gram, root, root) < 0 for root in roots)
        and all(_pair(gram, interior, root) > 0 for root in roots)
        and all(_reflection_integral(gram, root) for root in roots)
    ):
        raise ArithmeticError("the frozen JK-H12 source failed its exact replay")
    return payload, gram, roots, interior


def _seed_relations(
    gram: Matrix, roots: Sequence[IntegralVector]
) -> tuple[tuple[str, ...], ...]:
    covectors = tuple(_covector(gram, root) for root in roots)
    rows: list[tuple[str, ...]] = []
    for left, root in enumerate(roots):
        row = []
        for right, other in enumerate(roots):
            if left == right:
                row.append("m1")
            else:
                row.append(
                    _relation(
                        root,
                        covectors[left],
                        other,
                        covectors[right],
                    )[0]
                )
        rows.append(tuple(row))
    return tuple(rows)


def _doubling_walls(relations: Sequence[Sequence[str]]) -> tuple[int, ...]:
    odd = {"m3", "m5", "m7", "m9"}
    return tuple(
        index
        for index, row in enumerate(relations)
        if not any(token in odd for token in row)
    )


def _seed_facet_witness(
    gram: Matrix, interior: IntegralVector, root: IntegralVector
) -> Vector:
    coefficient = _pair(gram, interior, root) / _pair(gram, root, root)
    witness = tuple(
        _fraction(value) - coefficient * root_value
        for value, root_value in zip(interior, root, strict=True)
    )
    if not (
        _pair(gram, witness, root) == 0
        and _pair(gram, witness, witness) > 0
    ):
        raise ArithmeticError("orthogonal projection did not give a facet witness")
    return witness


def _boundary_count(
    gram: Matrix,
    roots: Sequence[IntegralVector],
    selected: int,
) -> int:
    doubling_root = roots[selected]
    keys: set[IntegralVector] = set()
    for copy in range(2):
        for index, root in enumerate(roots):
            if index == selected:
                continue
            image: Sequence[int | Fraction] = root
            if copy:
                image = _reflection_apply(gram, doubling_root, root)
            keys.add(_unoriented_key(image))
    return len(keys)


def _construct_double(
    gram: Matrix,
    seed_roots: Sequence[IntegralVector],
    seed_interior: IntegralVector,
    selected: int,
    relations: Sequence[Sequence[str]],
) -> dict[str, Any]:
    doubling_root = seed_roots[selected]
    if "m3" in relations[selected]:
        raise ArithmeticError("the selected wall has an odd incident angle")
    reflected_interior = _reflection_apply(gram, doubling_root, seed_interior)
    interior: Vector = tuple(
        _fraction(left) + _fraction(right)
        for left, right in zip(seed_interior, reflected_interior, strict=True)
    )
    seed_witnesses = tuple(
        _seed_facet_witness(gram, seed_interior, root) for root in seed_roots
    )

    boundary: dict[IntegralVector, dict[str, Any]] = {}
    for copy in range(2):
        for index, seed_root in enumerate(seed_roots):
            if index == selected:
                continue
            root_image: Sequence[int | Fraction] = seed_root
            witness: Sequence[int | Fraction] = seed_witnesses[index]
            if copy:
                root_image = _reflection_apply(gram, doubling_root, seed_root)
                witness = _reflection_apply(gram, doubling_root, witness)
            root = _primitive(root_image)
            key = _unoriented_key(root)
            pairing = _pair(gram, interior, root)
            if pairing == 0:
                raise ArithmeticError("the double interior lies on a boundary wall")
            if pairing < 0:
                root = tuple(-value for value in root)
            boundary.setdefault(
                key,
                {
                    "root": root,
                    "witness": tuple(_fraction(value) for value in witness),
                    "origin_copy": copy,
                    "origin_seed_facet": index + 1,
                },
            )

    records = [boundary[key] for key in sorted(boundary)]
    roots = tuple(record["root"] for record in records)
    witnesses = tuple(record["witness"] for record in records)
    covectors = tuple(_covector(gram, root) for root in roots)
    witness_squares = tuple(_pair(gram, witness, witness) for witness in witnesses)
    minimum_witness_slack: Fraction | None = None
    for target, witness in enumerate(witnesses):
        if not (
            _pair_from_covector(witness, covectors[target]) == 0
            and witness_squares[target] > 0
        ):
            raise ArithmeticError("a transported facet witness failed its target")
        slacks = tuple(
            _pair_from_covector(witness, covectors[index])
            for index in range(len(roots))
            if index != target
        )
        if not all(value > 0 for value in slacks):
            raise ArithmeticError("a doubled inequality is redundant or misoriented")
        local_minimum = min(slacks)
        minimum_witness_slack = (
            local_minimum
            if minimum_witness_slack is None
            else min(minimum_witness_slack, local_minimum)
        )

    relation_histogram: Counter[str] = Counter()
    ultraparallel_ratios: Counter[str] = Counter()
    for left in range(len(roots)):
        for right in range(left):
            token, ratio = _relation(
                roots[left],
                covectors[left],
                roots[right],
                covectors[right],
            )
            relation_histogram[token] += 1
            if token == "ultraparallel":
                ultraparallel_ratios[str(ratio)] += 1

    sympy_roots = sp.Matrix([[int(value) for value in root] for root in roots])
    interior_covectors = tuple(_pair(gram, interior, root) for root in roots)
    squares = tuple(_pair_from_covector(root, covector) for root, covector in zip(roots, covectors, strict=True))
    divisibilities = tuple(
        math.gcd(*(abs(int(value)) for value in covector if value))
        for covector in covectors
    )
    if not (
        sympy_roots.rank() == 13
        and _pair(gram, interior, interior) > 0
        and all(value > 0 for value in interior_covectors)
        and all(_content(root) == 1 for root in roots)
        and all(_reflection_integral(gram, root) for root in roots)
        and minimum_witness_slack is not None
    ):
        raise ArithmeticError("the doubled boundary failed its global exact gates")

    return {
        "id": f"JK_H12_D2_W{selected + 1}",
        "selected_doubling_wall_one_based": selected + 1,
        "copy_count_and_volume_multiple": 2,
        "hyperbolic_dimension": 12,
        "facet_count": len(roots),
        "full_linear_rank": sympy_roots.rank(),
        "strict_interior": {
            "coordinates": _vector_payload(interior),
            "square": _payload_number(_pair(gram, interior, interior)),
            "minimum_facet_pairing": _payload_number(min(interior_covectors)),
        },
        "strict_facet_witnesses": {
            "count": len(witnesses),
            "minimum_square": _payload_number(min(witness_squares)),
            "minimum_off_target_slack": _payload_number(minimum_witness_slack),
            "construction": "reflect the seed orthogonal-projection witness with its copy",
        },
        "primitive_root_square_histogram": {
            str(key): value for key, value in sorted(Counter(squares).items())
        },
        "primitive_root_divisibility_histogram": {
            str(key): value
            for key, value in sorted(Counter(divisibilities).items())
        },
        "relation_histogram": dict(sorted(relation_histogram.items())),
        "ultraparallel_normalized_square_histogram": dict(
            sorted(ultraparallel_ratios.items(), key=lambda item: Fraction(item[0]))
        ),
        "integral_reflection_count": len(roots),
        "boundary_records": [
            {
                "facet": index + 1,
                "primitive_oriented_coordinates": _vector_payload(record["root"]),
                "origin_copy": record["origin_copy"],
                "origin_seed_facet": record["origin_seed_facet"],
            }
            for index, record in enumerate(records)
        ],
        "finite_volume": True,
        "compact": False,
        "finite_volume_proof": "exact convex union of two finite-volume seed chambers",
        "noncompactness_proof": "index-two subgroup of the nonuniform seed reflection group",
        "all_exact_checks_pass": True,
    }


def build() -> dict[str, Any]:
    source, gram, roots, interior = _load_source()
    relations = _seed_relations(gram, roots)
    doubling = _doubling_walls(relations)
    counts: defaultdict[int, list[int]] = defaultdict(list)
    for selected in doubling:
        counts[_boundary_count(gram, roots, selected)].append(selected + 1)
    expected_counts = {
        296: list(range(31, 34)) + [36, 39, 42, 48, 51, 56, 58, 59, 60],
        317: list(range(61, 77)),
        326: list(range(77, 173)),
    }
    if len(doubling) != 124 or dict(counts) != expected_counts:
        raise ArithmeticError("the JK-H12 doubling-wall census changed")

    redoubling_left, redoubling_right = 30, 35
    left_root = roots[redoubling_left]
    right_root = roots[redoubling_right]
    null_sum = tuple(
        left + right for left, right in zip(left_root, right_root, strict=True)
    )
    primitive_null_sum = _primitive(null_sum)
    translated_left = _reflection_apply(
        gram,
        left_root,
        _reflection_apply(gram, right_root, left_root),
    )
    expected_translated_left = tuple(
        _fraction(left) + 2 * _fraction(delta)
        for left, delta in zip(left_root, null_sum, strict=True)
    )
    if not (
        redoubling_left in doubling
        and redoubling_right in doubling
        and _pair(gram, left_root, left_root) == -4
        and _pair(gram, right_root, right_root) == -4
        and _pair(gram, left_root, right_root) == 4
        and _pair(gram, primitive_null_sum, primitive_null_sum) == 0
        and translated_left == expected_translated_left
    ):
        raise ArithmeticError("the claimed JK-H12 parallel redoubling pair changed")
    representatives = tuple(
        _construct_double(gram, roots, interior, selected, relations)
        for selected in REPRESENTATIVES
    )
    if [row["facet_count"] for row in representatives] != [296, 317, 326]:
        raise ArithmeticError("the representative facet counts changed")

    seed_row_histograms = {}
    for selected in REPRESENTATIVES:
        seed_row_histograms[str(selected + 1)] = dict(
            sorted(Counter(token for index, token in enumerate(relations[selected]) if index != selected).items())
        )
    carrier = source["chamber"]["carrier"]
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-08",
        "all_checks_pass": True,
        "pinned_source": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": SOURCE_SHA256,
            "id": source["chamber"]["id"],
            "finite_volume": source["chamber"]["finite_volume"],
            "seed_facet_count": source["chamber"]["facet_count"],
        },
        "complete_doubling_wall_census": {
            "criterion": "no odd finite Coxeter label incident to the selected facet",
            "doubling_wall_count": len(doubling),
            "non_doubling_wall_count": len(roots) - len(doubling),
            "single_double_facet_count_classes": {
                str(count): walls for count, walls in sorted(counts.items())
            },
            "representative_seed_row_histograms": seed_row_histograms,
        },
        "carrier": {
            "name": carrier["name"],
            "rank": carrier["rank"],
            "signature": carrier["signature"],
            "determinant_absolute": carrier["determinant_absolute"],
            "very_general_period_consequence": carrier[
                "very_general_period_consequence"
            ],
            "same_full_Picard_lattice_for_all_descendants": True,
        },
        "infinite_redoubling_corollary": {
            "parallel_doubling_walls_one_based": [31, 36],
            "generator_squares_and_pairing": [-4, -4, 4],
            "primitive_common_null_direction": _vector_payload(primitive_null_sum),
            "reflection_group": "D_infinity",
            "translation": "T=s_31*s_36",
            "translation_identity": "T^z(a)=a+2*z*(a+b) for every integer z",
            "all_translation_powers_distinct": True,
            "consequence": "for every N>=1 a finite N-chamber path is a finite-volume noncompact H12 Coxeter polytope of volume N times the seed",
            "pairwise_nonisomorphic_reason": "the volumes are distinct positive integer multiples of the seed volume",
            "same_full_Picard_K3_carrier": True,
            "candidate_facet_formula_pending_independent_all_z_audit": "124*N+48",
        },
        "representatives": list(representatives),
        "decision": {
            "new_to_project_exact_polytope_count": 3,
            "infinite_H12_family": True,
            "pairwise_nonisomorphic_reason": "distinct facet counts 296, 317, and 326",
            "classification": "finite_volume_noncompact_H12_Coxeter_polytopes",
            "surface_statement": "weak full-Picard K3 divisor-reflection carriers",
            "worldwide_priority": "open",
        },
        "scope": {
            "all_124_single_doubles_constructible": True,
            "three_representatives_certified_here": True,
            "infinite_parallel_redoubling_family_proved": True,
            "all_parameter_facet_formula_proved_here": False,
            "all_124_labelled_isomorphism_types_classified": False,
            "literal_descendant_Vinberg_replay": False,
            "equivalent_volume_certificate": "exact index-two convex Coxeter gluing from a literal-Vinberg-certified seed",
            "full_nef_or_surface_automorphism_claimed": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if OUTPUT.read_text(encoding="utf-8") != rendered:
            raise ArithmeticError("JK-H12 single-doubling archive changed")
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(rendered, encoding="utf-8")
    print(json.dumps(payload["decision"], sort_keys=True))
    print(OUTPUT)


if __name__ == "__main__":
    main()
