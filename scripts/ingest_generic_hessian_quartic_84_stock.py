#!/usr/bin/env python3
"""Ingest the complete 84-wall chamber of a general Hessian quartic K3.

Dolgachev--Keum, arXiv:math/0103186, embed the rank-16 Picard lattice of
the minimal resolution of a general quartic Hessian in II(1,25).  Their
Lemma 8.3 and Remark 8.4 classify every wall of the induced Conway chamber,
and Lemma 8.2 identifies its exact symmetry group with C2 x S5.

The twenty nodal/line curves form the Desargues (10_3) configuration.  They
generate Pic(X), but have four relations.  We therefore keep formal
20-coordinate divisor expressions and identify two expressions by their
twenty exact intersection numbers.  A nonsingular 16-curve sub-Gram then
gives exact quotient coordinates for every wall.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import combinations, permutations
import json
from math import factorial
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/generic_hessian_quartic_84_wall_stock.json"
SOURCE_TEX = Path("/private/tmp/hessian0103186_src/hess.tex")

SCHEMA = "coxeter.generic_hessian_quartic.complete_84_wall_stock.v1"
SOURCE_URL = "https://export.arxiv.org/e-print/math/0103186"
SOURCE_TEX_SHA256 = "a24cb366240f57dec312d165418efe7709170380c50aea55abf7cfe6b7f3f26d"
SOURCE_TAR_SHA256 = "62cd1e64369ef1c385801e3943fb78f9d8f27fc78fb0f271d72d6f2bdaf13e4e"

Vector = tuple[Fraction, ...]
Permutation5 = tuple[int, ...]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _text(value: Fraction | sp.Rational | int) -> str:
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


N_LABELS = tuple(combinations(range(1, 6), 3))
T_LABELS = tuple(combinations(range(1, 6), 2))
FORMAL_LABELS = tuple(("N", label) for label in N_LABELS) + tuple(
    ("T", label) for label in T_LABELS
)
FORMAL_INDEX = {label: index for index, label in enumerate(FORMAL_LABELS)}


def _curve_name(label: tuple[str, tuple[int, ...]]) -> str:
    return label[0] + "".join(str(value) for value in label[1])


def _formal_gram() -> tuple[tuple[int, ...], ...]:
    """Intersection Gram of the twenty nodal and line curves."""

    gram = [[0 for _ in FORMAL_LABELS] for _ in FORMAL_LABELS]
    for i, (left_kind, left_label) in enumerate(FORMAL_LABELS):
        gram[i][i] = -2
        for j, (right_kind, right_label) in enumerate(FORMAL_LABELS):
            if left_kind == right_kind:
                continue
            node = left_label if left_kind == "N" else right_label
            line = right_label if right_kind == "T" else left_label
            gram[i][j] = int(set(line).issubset(node))
    return tuple(tuple(row) for row in gram)


FORMAL_GRAM = _formal_gram()


def _vector(terms: Iterable[tuple[str, Sequence[int], int | Fraction]]) -> Vector:
    result = [Fraction(0) for _ in FORMAL_LABELS]
    for kind, raw_label, coefficient in terms:
        label = (kind, tuple(sorted(raw_label)))
        result[FORMAL_INDEX[label]] += Fraction(coefficient)
    return tuple(result)


def _add(*vectors: Vector) -> Vector:
    return tuple(sum(entries, Fraction(0)) for entries in zip(*vectors))


def _scale(coefficient: int | Fraction, vector: Vector) -> Vector:
    return tuple(Fraction(coefficient) * entry for entry in vector)


def _unit(label: tuple[str, tuple[int, ...]]) -> Vector:
    return _vector(((label[0], label[1], 1),))


@lru_cache(maxsize=None)
def _class_key(vector: Vector) -> tuple[Fraction, ...]:
    """Canonical key in Pic(X) tensor Q, quotienting the four relations."""

    return tuple(
        sum(
            Fraction(FORMAL_GRAM[i][j]) * vector[j]
            for j in range(len(FORMAL_LABELS))
            if FORMAL_GRAM[i][j] and vector[j]
        )
        for i in range(len(FORMAL_LABELS))
    )


def _pair(left: Vector, right: Vector) -> Fraction:
    right_pairings = _class_key(right)
    return sum(
        coefficient * right_pairings[index]
        for index, coefficient in enumerate(left)
        if coefficient
    )


def _permute(vector: Vector, permutation: Permutation5) -> Vector:
    point_image = {index + 1: permutation[index] for index in range(5)}
    result = [Fraction(0) for _ in FORMAL_LABELS]
    for (kind, label), coefficient in zip(FORMAL_LABELS, vector):
        if not coefficient:
            continue
        target = (kind, tuple(sorted(point_image[value] for value in label)))
        result[FORMAL_INDEX[target]] += coefficient
    return tuple(result)


def _enriques(vector: Vector) -> Vector:
    """The Enriques involution N_A <-> T_(A complement)."""

    universe = set(range(1, 6))
    result = [Fraction(0) for _ in FORMAL_LABELS]
    for (kind, label), coefficient in zip(FORMAL_LABELS, vector):
        if not coefficient:
            continue
        target_kind = "T" if kind == "N" else "N"
        target = (target_kind, tuple(sorted(universe.difference(label))))
        result[FORMAL_INDEX[target]] += coefficient
    return tuple(result)


def _orbit(representative: Vector) -> dict[tuple[Fraction, ...], Vector]:
    result: dict[tuple[Fraction, ...], Vector] = {}
    for permutation in permutations(range(1, 6)):
        image = _permute(representative, permutation)
        result.setdefault(_class_key(image), image)
    return result


def _representatives() -> dict[str, Vector]:
    n_sum = _vector(("N", label, 1) for label in N_LABELS)
    t_sum = _vector(("T", label, 1) for label in T_LABELS)
    eta_h = _add(_scale(Fraction(3, 5), n_sum), _scale(Fraction(2, 5), t_sum))
    eta_s = _add(_scale(Fraction(2, 5), n_sum), _scale(Fraction(3, 5), t_sum))

    # Equation (8.1), multiplied by 3 to the primitive (-6)-root alpha.
    type_1a = _add(
        _scale(2, _add(eta_h, _scale(-1, eta_s))),
        _vector(
            (
                ("T", (2, 3), 2), ("T", (1, 4), 2),
                ("T", (1, 5), 2), ("T", (3, 5), 2),
                ("T", (2, 4), 2),
                ("N", (2, 3, 4), 1), ("N", (1, 2, 4), 1),
                ("N", (1, 3, 5), 1), ("N", (2, 3, 5), 1),
                ("N", (1, 4, 5), 1),
            )
        ),
    )

    # Equation (8.2), multiplied by 2 to the primitive (-4)-root alpha.
    type_2 = _add(
        eta_h,
        _vector(
            (
                ("T", (2, 5), -2),
                ("N", (1, 2, 5), -1), ("N", (2, 4, 5), -1),
                ("N", (2, 3, 5), -1), ("N", (1, 3, 4), -1),
            )
        ),
    )

    # Equation (8.3), multiplied by 6 to the primitive (-24)-vector.
    type_3 = _vector(
        (
            ("N", (3, 4, 5), 3), ("N", (1, 3, 4), 3),
            ("N", (1, 3, 5), 3), ("N", (1, 4, 5), 3),
            ("N", (1, 2, 5), -2), ("N", (2, 3, 4), -2),
            ("N", (2, 3, 5), 2), ("N", (1, 2, 4), 2),
            ("T", (1, 4), 4), ("T", (3, 5), 4),
            ("T", (1, 5), 2), ("T", (3, 4), 2),
        )
    )
    return {"1a": type_1a, "2": type_2, "3": type_3}


def wall_stock() -> tuple[list[Vector], list[str]]:
    representatives = _representatives()
    blocks: list[tuple[str, dict[tuple[Fraction, ...], Vector]]] = [
        ("0N", {_class_key(_unit(("N", label))): _unit(("N", label)) for label in N_LABELS}),
        ("0T", {_class_key(_unit(("T", label))): _unit(("T", label)) for label in T_LABELS}),
        ("1a", _orbit(representatives["1a"])),
        ("1b", _orbit(_enriques(representatives["1a"]))),
        ("2", _orbit(representatives["2"])),
        # S5 combines the paper's 15 cases (3a) and 15 cases (3b).
        ("3a_or_3b", _orbit(representatives["3"])),
    ]
    expected = {"0N": 10, "0T": 10, "1a": 12, "1b": 12, "2": 10, "3a_or_3b": 30}
    if {name: len(block) for name, block in blocks} != expected:
        raise ArithmeticError("the six source families did not reconstruct 84 quotient classes")

    walls: list[Vector] = []
    source_types: list[str] = []
    seen: set[tuple[Fraction, ...]] = set()
    for name, block in blocks:
        for key in sorted(block):
            if key in seen:
                raise ArithmeticError(f"duplicate wall class across source blocks: {name}")
            seen.add(key)
            walls.append(block[key])
            source_types.append(name)
    if len(walls) != 84:
        raise ArithmeticError(f"expected 84 walls, obtained {len(walls)}")
    return walls, source_types


def _wall_action(
    walls: Sequence[Vector], permutation: Permutation5, enriques: int
) -> tuple[int, ...]:
    lookup = {_class_key(wall): index for index, wall in enumerate(walls)}
    result = []
    for wall in walls:
        image = _permute(wall, permutation)
        if enriques:
            image = _enriques(image)
        result.append(lookup[_class_key(image)])
    return tuple(result)


def _cycle_type(permutation: Sequence[int]) -> tuple[int, ...]:
    seen: set[int] = set()
    lengths = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            current = permutation[current] - 1
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def _permutation_cycle_lengths(permutation: Sequence[int]) -> list[int]:
    # Wall actions are zero-based, unlike the five-point permutations above.
    seen: set[int] = set()
    lengths = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            current = permutation[current]
            length += 1
        lengths.append(length)
    return sorted(lengths, reverse=True)


def _class_size(cycle_type: Sequence[int]) -> int:
    multiplicities = Counter(cycle_type)
    denominator = 1
    for length, multiplicity in multiplicities.items():
        denominator *= length**multiplicity * factorial(multiplicity)
    return factorial(5) // denominator


def _basis_data() -> tuple[tuple[int, ...], sp.Matrix, sp.Matrix]:
    formal = sp.Matrix(FORMAL_GRAM)
    pivots = tuple(int(value) for value in formal.rref()[1])
    if len(pivots) != 16:
        raise ArithmeticError(f"formal curve Gram has rank {len(pivots)}, not 16")
    basis_gram = formal.extract(pivots, pivots)
    if basis_gram.det() != -192 or exact_inertia(basis_gram) != (1, 15, 0):
        raise ArithmeticError("the chosen quotient basis changed")
    return pivots, basis_gram, basis_gram.inv()


def _basis_coordinates(
    vector: Vector, pivots: Sequence[int], basis_gram_inverse: sp.Matrix
) -> sp.Matrix:
    pairings = sp.Matrix([sp.Rational(_class_key(vector)[index]) for index in pivots])
    return basis_gram_inverse * pairings


def _divisor_action_matrix(
    permutation: Permutation5,
    enriques: int,
    pivots: Sequence[int],
    basis_gram_inverse: sp.Matrix,
) -> sp.Matrix:
    columns = []
    for pivot in pivots:
        image = _permute(_unit(FORMAL_LABELS[pivot]), permutation)
        if enriques:
            image = _enriques(image)
        columns.append(_basis_coordinates(image, pivots, basis_gram_inverse))
    return sp.Matrix.hstack(*columns)


def _element_class_census(walls: Sequence[Vector]) -> list[dict[str, Any]]:
    pivots, _, basis_gram_inverse = _basis_data()
    representatives: dict[tuple[int, tuple[int, ...]], Permutation5] = {}
    actions: set[tuple[int, ...]] = set()
    for raw in permutations(range(1, 6)):
        permutation = tuple(raw)
        cycle_type = _cycle_type(permutation)
        for enriques in (0, 1):
            representatives.setdefault((enriques, cycle_type), permutation)
            actions.add(_wall_action(walls, permutation, enriques))
    if len(actions) != 240 or len(representatives) != 14:
        raise ArithmeticError("C2 x S5 did not act faithfully with fourteen element classes")

    rows = []
    for (enriques, cycle_type), permutation in sorted(representatives.items()):
        action = _wall_action(walls, permutation, enriques)
        divisor_action = _divisor_action_matrix(
            permutation, enriques, pivots, basis_gram_inverse
        )
        rows.append(
            {
                "central_enriques_factor": enriques,
                "s5_cycle_type": list(cycle_type),
                "class_size": _class_size(cycle_type),
                "order": int(sp.ilcm(*cycle_type, 2 if enriques else 1)),
                "fixed_divisor_rank": 16 - (divisor_action - sp.eye(16)).rank(),
                "wall_cycle_lengths": _permutation_cycle_lengths(action),
                "fixed_walls": sum(value == index for index, value in enumerate(action)),
            }
        )
    return rows


def _coxeter_pair_audit(gram: Sequence[Sequence[Fraction]]) -> dict[str, Any]:
    allowed = {
        Fraction(0): "2",
        Fraction(1, 4): "3",
        Fraction(1, 2): "4",
        Fraction(3, 4): "6",
    }
    counts: Counter[str] = Counter()
    failures = []
    for right in range(len(gram)):
        for left in range(right):
            product = gram[left][right]
            ratio = product * product / (gram[left][left] * gram[right][right])
            if product < 0:
                label = "forbidden_negative_inward_product"
            elif ratio in allowed:
                label = allowed[ratio]
            elif ratio == 1:
                label = "infinity"
            elif ratio > 1:
                label = "ultraparallel"
            else:
                label = "forbidden_non_coxeter_angle"
            counts[label] += 1
            if label.startswith("forbidden"):
                failures.append(
                    {
                        "walls_zero_based": [left, right],
                        "product": _text(product),
                        "normalized_square_ratio": _text(ratio),
                        "reason": label,
                    }
                )
    return {
        "ambient_induced_chamber_is_coxeter": not failures,
        "pair_counts": dict(sorted(counts.items())),
        "failure_count": len(failures),
        "first_failure_witnesses": failures[:20],
    }


def build_payload() -> dict[str, Any]:
    walls, source_types = wall_stock()
    pivots, basis_gram, basis_gram_inverse = _basis_data()
    wall_gram = [[_pair(left, right) for right in walls] for left in walls]
    if sp.Matrix(wall_gram).rank() != 16:
        raise ArithmeticError("84-wall Gram does not have rank 16")

    element_classes = _element_class_census(walls)
    high_classes = [
        row
        for row in element_classes
        if row["order"] > 1 and row["fixed_divisor_rank"] >= 11
    ]
    coordinates = [
        [_text(value) for value in _basis_coordinates(wall, pivots, basis_gram_inverse)]
        for wall in walls
    ]
    payload = {
        "schema": SCHEMA,
        "source": {
            "authors": ["Igor V. Dolgachev", "Jonghae Keum"],
            "title": "Birational automorphisms of quartic Hessian surfaces",
            "arxiv": "math/0103186v1",
            "source_url": SOURCE_URL,
            "source_tar_sha256": SOURCE_TAR_SHA256,
            "source_tex_path": str(SOURCE_TEX),
            "source_tex_sha256": SOURCE_TEX_SHA256,
            "pin_status": "source URL and exact TeX/tar hashes pinned",
        },
        "algebraic_carrier": {
            "object": "minimal nonsingular model of the Hessian quartic of a general cubic surface",
            "type": "complex projective K3 surface",
            "picard_rank": 16,
            "picard_discriminant_absolute": 48,
            "hyperbolic_dimension": 15,
            "curve_configuration": "twenty (-2)-curves in the Desargues (10_3) configuration",
        },
        "completeness": {
            "wall_count": 84,
            "source_partition": {
                "case_0": 20,
                "case_1a": 12,
                "case_1b": 12,
                "case_2": 10,
                "case_3a": 15,
                "case_3b": 15,
            },
            "source_results": [
                "Section 8 defines P(S_H)' as the intersection with a Conway chamber",
                "the paragraph before Lemma 8.2 states that P(S_H)' has finitely many faces and finite volume",
                "Lemma 8.3 classifies every Leech-root wall and Remark 8.4 counts the cases",
                "Lemma 8.2 proves Aut(P(S_H)') is C2 x S5",
            ],
            "bounded_prefix": False,
            "surface_automorphism_gate_imposed": False,
        },
        "exact_ingestion": {
            "formal_curve_labels": [_curve_name(label) for label in FORMAL_LABELS],
            "formal_curve_gram_rank": 16,
            "formal_curve_gram_nullity": 4,
            "quotient_basis_formal_indices_zero_based": list(pivots),
            "quotient_basis_labels": [_curve_name(FORMAL_LABELS[index]) for index in pivots],
            "quotient_basis_gram_determinant": int(basis_gram.det()),
            "quotient_basis_gram_inertia": list(exact_inertia(basis_gram)),
            "wall_gram_rank": 16,
            "wall_gram_inertia": [1, 15, 68],
            "wall_norm_histogram": dict(
                sorted((str(norm), count) for norm, count in Counter(wall_gram[i][i] for i in range(84)).items())
            ),
            "wall_source_types": source_types,
            "wall_coordinates_in_16_curve_quotient_basis": coordinates,
            "wall_gram": [[_text(value) for value in row] for row in wall_gram],
        },
        "symmetry_census": {
            "source_exact_group": "C2 x S5",
            "order": 240,
            "element_conjugacy_class_count": 14,
            "element_classes": element_classes,
            "proper_element_classes_fixed_rank_at_least_11": high_classes,
            "gap_subgroup_conjugacy_class_count": 57,
            "all_subgroup_fold_plan": "enumerate the 57 GAP classes, restrict all 84 inequalities, exact conic-reduce, then apply the Coxeter and literal Vinberg finite-volume gates",
        },
        "ambient_coxeter_gate": _coxeter_pair_audit(wall_gram),
        "lineage_gate": {
            "stock_previously_present_in_project_ledger_before_this_ingestion": False,
            "novel_hyperbolic_polytope_claim": False,
            "next_required": "complete all 57 subgroup folds and compare any cofinite survivors with maximal reflective-lattice and CoxIter/literature lineages",
        },
    }
    payload["all_checks_pass"] = bool(
        payload["exact_ingestion"]["wall_gram_rank"] == 16
        and sum(payload["completeness"]["source_partition"].values()) == 84
        and payload["symmetry_census"]["order"] == 240
        and payload["symmetry_census"]["element_conjugacy_class_count"] == 14
    )
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-source", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()

    if args.check_source:
        if not SOURCE_TEX.exists():
            raise FileNotFoundError(SOURCE_TEX)
        if _sha256(SOURCE_TEX) != SOURCE_TEX_SHA256:
            raise ArithmeticError("pinned Dolgachev--Keum source TeX changed")

    payload = build_payload()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text() != encoded:
            raise ArithmeticError(f"archive mismatch: {OUTPUT}")
        print(f"archive matches: {OUTPUT}")
        return
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(encoded)
    print(json.dumps({
        "output": str(OUTPUT),
        "walls": payload["completeness"]["wall_count"],
        "group_order": payload["symmetry_census"]["order"],
        "high_element_classes": len(payload["symmetry_census"]["proper_element_classes_fixed_rank_at_least_11"]),
        "ambient_coxeter": payload["ambient_coxeter_gate"]["ambient_induced_chamber_is_coxeter"],
        "all_checks_pass": payload["all_checks_pass"],
    }, indent=2))


if __name__ == "__main__":
    main()
