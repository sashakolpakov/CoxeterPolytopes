"""Screen the Frobenius-fixed slice of the characteristic-3 Fermat K3 chamber.

The exact input is Kondo--Shimada's ``compdataAutFQChar3.txt``.  In that
file ``FQWW112``, ``FQWW648`` and ``FQWW5184`` are *covectors*: if ``r`` is
a wall normal and ``G = FQGramS``, the archived row is ``r G``.  Working
with covectors avoids introducing fractions until a final scalar product;
their Gram matrix is therefore ``G^{-1}``.

This first bounded pass treats the arithmetic Frobenius involution
``FQAF``.  It proves the dimension of the fixed Lorentzian space, derives
the wall permutation from the lattice matrix, forms the exact restricted
wall rays, and searches for a decisive local obstruction.  A restricted
wall is certified as a genuine facet only when the orthogonal projection
of the ample class ``FQh0`` to that wall is timelike and pairs strictly
positively with every other distinct restricted wall ray.  Hence two such
facets with a forbidden normalized product are a rigorous Coxeter no-go;
no global facet enumeration or floating-point LP is needed for that
conclusion.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
from functools import lru_cache, reduce
from math import gcd
import json
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "data"
    / "external"
    / "kondo_shimada_fermat_quartic_char3"
    / "compdataAutFQChar3.txt"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_frobenius_fold.json"
)


def _extract_assignment(text: str, name: str) -> str:
    """Extract one bracketed assignment, tolerating the source line escape."""

    marker = text.index(name)
    start = text.index("[", marker)
    depth = 0
    for end in range(start, len(text)):
        character = text[end]
        if character == "[":
            depth += 1
        elif character == "]":
            depth -= 1
            if depth == 0:
                return text[start : end + 1].replace("\\\n", "")
    raise ValueError(f"unterminated assignment {name}")


def _numeric_ast(node: ast.AST) -> Any:
    if isinstance(node, ast.List):
        return [_numeric_ast(element) for element in node.elts]
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return sp.Rational(node.value)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -_numeric_ast(node.operand)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        numerator = _numeric_ast(node.left)
        denominator = _numeric_ast(node.right)
        return sp.Rational(numerator, denominator)
    raise ValueError(f"unsupported source expression: {ast.dump(node)}")


def _parse_array(text: str, name: str) -> Any:
    expression = ast.parse(_extract_assignment(text, name), mode="eval")
    return _numeric_ast(expression.body)


def parse_source(path: Path = SOURCE) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    gram = sp.Matrix(_parse_array(text, "FQGramS"))
    ample = sp.Matrix(1, gram.rows, _parse_array(text, "FQh0"))
    frobenius = sp.Matrix(_parse_array(text, "FQAF"))
    blocks: list[tuple[str, tuple[tuple[int, ...], ...]]] = []
    for name, label in (
        ("FQWW112", "112"),
        ("FQWW648", "648"),
        ("FQWW5184", "5184"),
    ):
        rows = tuple(
            tuple(int(value) for value in row) for row in _parse_array(text, name)
        )
        blocks.append((label, rows))
    wall_types = tuple(
        label for label, rows in blocks for _ in rows
    )
    walls = tuple(row for _, rows in blocks for row in rows)
    return {
        "gram": gram,
        "ample": ample,
        "frobenius": frobenius,
        "walls": walls,
        "wall_types": wall_types,
    }


def _primitive_positive(row: tuple[int, ...], ample: np.ndarray) -> tuple[int, ...]:
    divisor = reduce(gcd, (abs(value) for value in row if value), 0)
    if divisor == 0:
        raise ArithmeticError("a restricted wall covector vanished")
    primitive = tuple(value // divisor for value in row)
    pairing = sum(left * right for left, right in zip(primitive, ample, strict=True))
    if pairing == 0:
        raise ArithmeticError("a wall through the ample interior point was produced")
    if pairing < 0:
        primitive = tuple(-value for value in primitive)
    return primitive


def _permutation_orbits(permutation: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    unseen = set(range(len(permutation)))
    orbits: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        orbit: list[int] = []
        current = start
        while current not in orbit:
            orbit.append(current)
            current = permutation[current]
        if current != start:
            raise ArithmeticError("wall action is not a permutation cycle")
        unseen.difference_update(orbit)
        orbits.append(tuple(orbit))
    return tuple(orbits)


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def _coxeter_label(product: int, left_norm: int, right_norm: int) -> str | None:
    """Classify using numerators in one common metric denominator."""

    if product < 0:
        return None
    if product == 0:
        return "2"
    ratio = Fraction(product * product, left_norm * right_norm)
    finite = {
        Fraction(1, 4): "3",
        Fraction(1, 2): "4",
        Fraction(3, 4): "6",
        Fraction(1): "infinity",
    }
    if ratio in finite:
        return finite[ratio]
    if ratio > 1:
        return "ultraparallel"
    return None


@lru_cache(maxsize=2)
def build_artifact(path: Path = SOURCE) -> dict[str, Any]:
    source = parse_source(path)
    gram: sp.Matrix = source["gram"]
    ample_matrix: sp.Matrix = source["ample"]
    frobenius: sp.Matrix = source["frobenius"]
    walls = source["walls"]
    wall_types = source["wall_types"]

    identity = sp.eye(gram.rows)
    if frobenius * frobenius != identity:
        raise ArithmeticError("FQAF is not an involution")
    if frobenius * gram * frobenius.T != gram:
        raise ArithmeticError("FQAF does not preserve FQGramS")
    if ample_matrix * frobenius != ample_matrix:
        raise ArithmeticError("FQAF does not fix FQh0")

    # Rows act on divisor classes on the right.  Covectors therefore act by
    # F^{-T}; for this involution that matrix is simply F^T.
    covector_action = frobenius.T
    wall_index = {wall: index for index, wall in enumerate(walls)}
    if len(wall_index) != len(walls):
        raise ArithmeticError("published chamber covectors are not distinct")
    permutation: list[int] = []
    for wall in walls:
        image = tuple(int(value) for value in sp.Matrix(1, gram.rows, wall) * covector_action)
        if image not in wall_index:
            raise ArithmeticError("FQAF does not permute the published chamber walls")
        permutation.append(wall_index[image])
    wall_permutation = tuple(permutation)
    if any(wall_permutation[wall_permutation[index]] != index for index in range(len(walls))):
        raise ArithmeticError("derived Frobenius wall permutation is not involutory")
    if any(wall_types[index] != wall_types[image] for index, image in enumerate(wall_permutation)):
        raise ArithmeticError("Frobenius mixes the three published wall types")

    orbits = _permutation_orbits(wall_permutation)
    ample = np.array([int(value) for value in ample_matrix], dtype=np.int64)
    ray_to_source_orbits: dict[tuple[int, ...], list[int]] = {}
    orbit_sums: list[tuple[int, ...]] = []
    for orbit_index, orbit in enumerate(orbits):
        row = tuple(sum(walls[index][column] for index in orbit) for column in range(gram.rows))
        primitive = _primitive_positive(row, ample)
        if primitive not in ray_to_source_orbits:
            orbit_sums.append(primitive)
            ray_to_source_orbits[primitive] = []
        ray_to_source_orbits[primitive].append(orbit_index)

    inverse = gram.inv()
    denominator = int(sp.ilcm(*[sp.Rational(value).q for value in inverse]))
    metric_numerator = np.array(inverse * denominator, dtype=np.int64)
    restricted = np.array(orbit_sums, dtype=np.int64)
    restricted_metric_rows = restricted @ metric_numerator
    norm_numerators = np.sum(restricted_metric_rows * restricted, axis=1)
    ample_pairings = restricted @ ample

    restricted_sympy = sp.Matrix(restricted.tolist())
    fixed_rank = gram.rows - (frobenius - identity).rank()
    restricted_rank = restricted_sympy.rank()
    if restricted_rank != fixed_rank:
        raise ArithmeticError("restricted walls do not span the Frobenius fixed dual space")
    # Never build the roughly 3000 by 3000 singular Gram matrix.  Pivot
    # columns of the transpose select an exact basis of its row span, and a
    # 12 by 12 Gram matrix carries the complete nonzero inertia.
    _, basis_indices = restricted_sympy.T.rref()
    restricted_basis = restricted_sympy[list(basis_indices), :]
    restricted_basis_gram = restricted_basis * inverse * restricted_basis.T
    nonzero_inertia = exact_inertia(restricted_basis_gram)
    inertia = (
        nonzero_inertia[0],
        nonzero_inertia[1],
        len(orbit_sums) - fixed_rank,
    )
    expected_inertia = (1, fixed_rank - 1, len(orbit_sums) - fixed_rank)
    if inertia != expected_inertia:
        raise ArithmeticError(f"unexpected restricted wall inertia {inertia}")

    strict_facets: list[int] = []
    facet_certificates: list[dict[str, Any]] = []
    for index, norm in enumerate(norm_numerators):
        norm = int(norm)
        if norm >= 0:
            continue
        gram_row = restricted_metric_rows[index] @ restricted.T
        # For q_i < 0,
        #   <h - <h,n_i>/q_i n_i, n_j> > 0
        # is equivalent to a_i g_ij - q_i a_j > 0.
        margins = int(ample_pairings[index]) * gram_row - norm * ample_pairings
        if int(margins[index]) != 0:
            raise ArithmeticError("facet-witness self pairing is not zero")
        other_margins = np.delete(margins, index)
        if other_margins.size and int(other_margins.min()) > 0:
            strict_facets.append(index)
            witness_square = Fraction(
                int((ample_matrix * gram * ample_matrix.T)[0]) * norm
                - int(ample_pairings[index]) ** 2 * denominator,
                norm,
            )
            if witness_square <= 0:
                raise ArithmeticError("strict facet witness is not timelike")
            facet_certificates.append(
                {
                    "restricted_ray_index": index + 1,
                    "primitive_normal_covector": list(orbit_sums[index]),
                    "source_wall_orbit_indices": [
                        source_index + 1
                        for source_index in ray_to_source_orbits[orbit_sums[index]]
                    ],
                    "source_wall_indices": sorted(
                        wall_index + 1
                        for source_index in ray_to_source_orbits[orbit_sums[index]]
                        for wall_index in orbits[source_index]
                    ),
                    "normal_square": _fraction_text(Fraction(norm, denominator)),
                    "ample_pairing": int(ample_pairings[index]),
                    "witness_square": _fraction_text(witness_square),
                    "minimum_integer_strict_margin": int(other_margins.min()),
                }
            )

    bad_pairs: list[dict[str, Any]] = []
    label_counts: Counter[str] = Counter()
    for left_position, left in enumerate(strict_facets):
        gram_row = restricted_metric_rows[left] @ restricted.T
        for right in strict_facets[left_position + 1 :]:
            product = int(gram_row[right])
            left_norm = int(norm_numerators[left])
            right_norm = int(norm_numerators[right])
            label = _coxeter_label(product, left_norm, right_norm)
            ratio = Fraction(product * product, left_norm * right_norm)
            if label is None:
                bad_pairs.append(
                    {
                        "restricted_ray_indices": [left + 1, right + 1],
                        "normal_product": _fraction_text(Fraction(product, denominator)),
                        "normalized_squared_product": _fraction_text(ratio),
                        "reason": (
                            "negative normal product (non-acute)"
                            if product < 0
                            else "not a Coxeter cosine square"
                        ),
                    }
                )
            else:
                label_counts[label] += 1

    orbit_type_histogram: Counter[str] = Counter()
    for orbit in orbits:
        orbit_type_histogram[f"{wall_types[orbit[0]]}:size_{len(orbit)}"] += 1

    return {
        "schema": "ss_fermat_quartic_char3_frobenius_fold_v1",
        "source": {
            "paper": "Kondo--Shimada, The automorphism group of a supersingular K3 surface with Artin invariant 1 in characteristic 3",
            "paper_url": "https://arxiv.org/abs/1205.6520",
            "data_url": "https://www.math.sci.hiroshima-u.ac.jp/shimada/preprints/FermatQuartChar3/autFQ/compdataAutFQChar3.txt",
            "assignments": ["FQGramS", "FQh0", "FQAF", "FQWW112", "FQWW648", "FQWW5184"],
        },
        "scope": {
            "element": "arithmetic Frobenius FQAF",
            "element_order": 2,
            "full_chamber_stabilizer_census_completed": False,
            "surface_automorphism_required": False,
            "global_facet_reduction_completed": False,
            "local_strict_facet_gate_is_decisive_when_bad_pair_count_is_positive": True,
        },
        "lattice": {
            "rank": gram.rows,
            "determinant": int(gram.det()),
            "inertia": list(exact_inertia(gram)),
            "ample_square": str((ample_matrix * gram * ample_matrix.T)[0]),
        },
        "frobenius": {
            "fixed_rank": fixed_rank,
            "hyperbolic_dimension": fixed_rank - 1,
            "is_involution": True,
            "is_lattice_isometry": True,
            "fixes_ample_class": True,
            "wall_permutation_is_involution": True,
        },
        "walls": {
            "published_count": len(walls),
            "published_type_counts": dict(Counter(wall_types)),
            "frobenius_orbit_count": len(orbits),
            "frobenius_orbit_type_size_histogram": dict(sorted(orbit_type_histogram.items())),
            "distinct_restricted_ray_count": len(orbit_sums),
            "coincident_restricted_ray_count": len(orbits) - len(orbit_sums),
            "restricted_ray_rank": restricted_rank,
            "restricted_gram_inertia": list(inertia),
            "metric_denominator": denominator,
        },
        "strict_facet_gate": {
            "certified_strict_facet_count": len(strict_facets),
            "certificates": facet_certificates,
            "coxeter_compatible_pair_label_counts": dict(sorted(label_counts.items())),
            "bad_pair_count": len(bad_pairs),
            "bad_pairs": bad_pairs,
            "decision": (
                "not_coxeter"
                if bad_pairs
                else "inconclusive_without_global_facet_reduction"
            ),
        },
        "checks": {
            "source_dimensions_match": len(walls) == 112 + 648 + 5184,
            "frobenius_exact_checks_pass": True,
            "frobenius_preserves_each_wall_type": True,
            "restricted_walls_span_fixed_space": restricted_rank == fixed_rank,
            "restricted_space_is_lorentzian": inertia == expected_inertia,
            "each_listed_facet_has_a_strict_timelike_relative_interior_witness": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    artifact = build_artifact(arguments.source)
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    summary = {
        "hyperbolic_dimension": artifact["frobenius"]["hyperbolic_dimension"],
        "wall_orbits": artifact["walls"]["frobenius_orbit_count"],
        "restricted_rays": artifact["walls"]["distinct_restricted_ray_count"],
        "strict_facets": artifact["strict_facet_gate"]["certified_strict_facet_count"],
        "bad_pairs": artifact["strict_facet_gate"]["bad_pair_count"],
        "decision": artifact["strict_facet_gate"]["decision"],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
