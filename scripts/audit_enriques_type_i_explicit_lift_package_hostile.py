#!/usr/bin/env python3
"""Independent hostile replay of the explicit Type-I Enriques lift package.

The production certificate uses PARI ``qfminim`` and its own reconstruction
helpers.  This audit deliberately does not import that certificate.  It
rebuilds the involution from the raw OSCAR lattice chain, generates the
``E8+A1`` roots by Weyl reflections, and exhausts the L4 ellipsoid with an
exact rational LDL recursion.

Geometric existence is kept separate from the byte-level checks: the OSCAR
file is a lattice/Borcherds context.  The actual Type-I Enriques surface is a
literature input supplied by Brandhorst--Shimada, Section 7.4.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import ast
import hashlib
from math import isqrt
import json
from pathlib import Path
import re
from typing import Any, Iterable

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


ROOT = Path(__file__).resolve().parents[1]
OSCAR_DATA = (
    ROOT
    / "data"
    / "external"
    / "oscar_enriques_type_i"
    / "TauTaubarGenericEnriquesSurfaceNo172.mrdi"
)
OSCAR_README = OSCAR_DATA.with_name("README.md")
SHIMADA_DATA = (
    ROOT
    / "data"
    / "external"
    / "shimada_l10l26note"
    / "L10L26compdata2.txt"
)
PRODUCTION = (
    ROOT
    / "data"
    / "enumerations"
    / "enriques_type_i_explicit_lift_package.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "enriques_type_i_explicit_lift_package_hostile_audit.json"
)

EXPECTED_SHA256 = {
    "oscar_data": "14c9bddd93dbbcbefcda855302ee2f4fccc2f9acfc6e29b1d385e041b36f5fa8",
    "oscar_readme": "80a8a0a5aac90c985f2e48c85a67d4de122e61129ab7e51347a56ee7d6b802a9",
    "shimada": "218c8188e467bde29369cdeda3771d767c3f70105d8fb56def8d698b90e75a6c",
    "production": "d9313e37687e355c9ca8f86d1ee11d74f0810a8c6906a0443f5377d2481fece7",
}
OSCAR_LOADER_URL = (
    "https://raw.githubusercontent.com/oscar-system/Oscar.jl/v1.4.0/"
    "src/AlgebraicGeometry/Surfaces/EnriquesAut.jl"
)
OSCAR_LOADER_SHA256 = (
    "48cebfc2e924b914bfdc4d4d4f9d430a06439f8af35b798c8f7a556395989924"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix(raw: Iterable[Iterable[Any]]) -> sp.Matrix:
    return sp.Matrix([[int(value) for value in row] for row in raw])


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=sp.ZZ)
    return [
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if smith[index, index]
    ]


def _vector_hash(vectors: Iterable[Iterable[int]]) -> str:
    canonical = sorted(tuple(int(value) for value in vector) for vector in vectors)
    return hashlib.sha256(
        json.dumps(canonical, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _find_ambient_gram(document: dict[str, Any]) -> sp.Matrix:
    candidates = []
    for reference in document["_refs"].values():
        value = reference.get("data")
        if (
            isinstance(value, list)
            and len(value) == 26
            and all(isinstance(row, list) and len(row) == 26 for row in value)
        ):
            candidates.append(_matrix(value))
    if len(candidates) != 1:
        raise ArithmeticError("the OSCAR archive has no unique 26-dimensional Gram")
    return candidates[0]


def _first_shimada_record() -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Parse only the three literal integer arrays needed from record 12A."""

    text = SHIMADA_DATA.read_text(encoding="utf-8")
    patterns = {
        "gram": r"L10L26compdata2:=\[ rec\(\s*GramL10 := (\[.*?\]),\s*GramL26",
        "configuration": r"configuration := (\[.*?\]),\s*embS",
        "walls": r'name := "12A",\s*walls := (\[.*?\]\s*\])\s*\),\s*rec\(',
    }
    values = {}
    for name, pattern in patterns.items():
        match = re.search(pattern, text, flags=re.S)
        if match is None:
            raise ArithmeticError(f"could not parse Shimada 12A {name}")
        values[name] = _matrix(ast.literal_eval(match.group(1)))
    return values["gram"], values["configuration"], values["walls"]


def _ceil_div(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


def _ldl_upper(matrix: sp.Matrix) -> tuple[list[sp.Rational], sp.Matrix]:
    """Exact ``matrix = U.T*diag(d)*U`` with unit upper triangular U."""

    if matrix.rows != matrix.cols or matrix != matrix.T:
        raise ArithmeticError("LDL input is not symmetric square")
    size = matrix.rows
    upper = sp.eye(size)
    diagonal: list[sp.Rational] = []
    for column in range(size):
        pivot = sp.Rational(
            matrix[column, column]
            - sum(
                diagonal[index] * upper[index, column] ** 2
                for index in range(column)
            )
        )
        if pivot <= 0:
            raise ArithmeticError("LDL input is not positive definite")
        diagonal.append(pivot)
        for target in range(column + 1, size):
            upper[column, target] = sp.Rational(
                matrix[column, target]
                - sum(
                    diagonal[index]
                    * upper[index, column]
                    * upper[index, target]
                    for index in range(column)
                ),
                pivot,
            )
    return diagonal, upper


def _exact_short_vectors(
    matrix: sp.Matrix, bound: int
) -> tuple[list[tuple[int, ...]], dict[str, Any]]:
    """Complete exact ellipsoid enumeration, independent of PARI qfminim."""

    diagonal, upper = _ldl_upper(matrix)
    vector = [0] * matrix.rows
    vectors: list[tuple[int, ...]] = []
    leaf_count = 0

    def recurse(index: int, remaining: sp.Rational) -> None:
        nonlocal leaf_count
        if index < 0:
            leaf_count += 1
            if any(vector):
                column = sp.Matrix(vector)
                exact = int((column.T * matrix * column)[0])
                if not 0 < exact <= bound:
                    raise ArithmeticError("LDL recursion admitted a vector outside the bound")
                vectors.append(tuple(vector))
            return
        center = sp.Rational(
            sum(
                upper[index, later] * vector[later]
                for later in range(index + 1, matrix.rows)
            )
        )
        allowance = sp.Rational(remaining / diagonal[index])
        if allowance < 0:
            return
        center_numerator, center_denominator = int(center.p), int(center.q)
        radius_squared_numerator = int(allowance.p) * center_denominator**2
        radius_squared_denominator = int(allowance.q)
        radius = isqrt(radius_squared_numerator // radius_squared_denominator)
        low = _ceil_div(-radius - center_numerator, center_denominator)
        high = (radius - center_numerator) // center_denominator
        for value in range(low, high + 1):
            term = diagonal[index] * (sp.Rational(value) + center) ** 2
            if term <= remaining:
                vector[index] = value
                recurse(index - 1, remaining - term)
        vector[index] = 0

    recurse(matrix.rows - 1, sp.Rational(bound))
    vectors = sorted(set(vectors))
    if leaf_count != len(vectors) + 1:
        raise ArithmeticError("LDL recursion produced duplicate leaves")
    return vectors, {
        "engine": "independent exact rational LDL ellipsoid recursion",
        "bound": int(bound),
        "nonzero_vector_count": len(vectors),
        "including_zero_leaf_count": leaf_count,
        "vector_set_sha256": _vector_hash(vectors),
    }


def _weyl_root_closure(cartan: sp.Matrix) -> list[tuple[int, ...]]:
    """Generate all roots from simple roots using exact simple reflections."""

    size = cartan.rows
    roots = set()
    frontier = []
    for index in range(size):
        root = tuple(1 if entry == index else 0 for entry in range(size))
        for signed in (root, tuple(-value for value in root)):
            roots.add(signed)
            frontier.append(signed)
    while frontier:
        root = frontier.pop()
        for simple in range(size):
            pairing = sum(root[index] * int(cartan[index, simple]) for index in range(size))
            reflected = list(root)
            reflected[simple] -= pairing
            reflected_tuple = tuple(reflected)
            if reflected_tuple not in roots:
                column = sp.Matrix(reflected_tuple)
                if int((column.T * cartan * column)[0]) != 2:
                    raise ArithmeticError("simple reflection changed root norm")
                roots.add(reflected_tuple)
                frontier.append(reflected_tuple)
    return sorted(roots)


def analyze() -> dict[str, Any]:
    actual_hashes = {
        "oscar_data": _sha256(OSCAR_DATA),
        "oscar_readme": _sha256(OSCAR_README),
        "shimada": _sha256(SHIMADA_DATA),
        "production": _sha256(PRODUCTION),
    }
    if actual_hashes != EXPECTED_SHA256:
        raise AssertionError(f"a pinned input changed: {actual_hashes}")

    readme = OSCAR_README.read_text(encoding="utf-8")
    if (
        OSCAR_LOADER_SHA256 not in readme
        or "SY, SX, L26, w, u = load(...)" not in readme
        or "not equations for a\nsurface" not in readme
    ):
        raise ArithmeticError("OSCAR source semantics/provenance note changed")

    production = json.loads(PRODUCTION.read_text(encoding="utf-8"))
    raw = json.loads(OSCAR_DATA.read_text(encoding="utf-8"))
    sy = _matrix(raw["data"][0])
    sx = _matrix(raw["data"][1])
    l26_basis = _matrix(raw["data"][2])
    weyl = _matrix(raw["data"][3])
    preprocessing = _matrix(raw["data"][4])
    ambient = _find_ambient_gram(raw)
    if (
        sy.shape != (10, 26)
        or sx.shape != (19, 26)
        or l26_basis != sp.eye(26)
        or weyl.shape != (1, 26)
        or preprocessing.shape != (1, 26)
    ):
        raise ArithmeticError("OSCAR tuple dimensions changed")

    shimada_gram, shimada_configuration, shimada_walls = _first_shimada_record()
    gram_sy = sy * ambient * sy.T
    gram_sx = sx * ambient * sx.T
    if gram_sy != 2 * shimada_gram or sx[:10, :] != sy:
        raise ArithmeticError("the raw OSCAR/Shimada marking does not join exactly")
    if shimada_walls * shimada_gram * shimada_walls.T != shimada_configuration:
        raise ArithmeticError("the raw Shimada wall Gram does not replay")
    archived_walls = _matrix(
        [
            record["downstairs_root_L10_coordinates"]
            for record in production["finite_lift_compiler_calibration"]["lift_records"]
        ]
    )
    if archived_walls != shimada_walls:
        raise ArithmeticError("production lifts are not in Shimada's literal wall order")

    # Row coordinates are forced by z -> z*tau.  The column convention is
    # explicitly tested and rejected so a transposition cannot pass silently.
    d_matrix = sx[10:, :] * ambient * sy.T * shimada_gram.inv()
    if any(value.q != 1 for value in d_matrix):
        raise ArithmeticError("projection/glue matrix is not integral")
    d_matrix = _matrix(_rows(d_matrix))
    tau = sp.Matrix.vstack(
        sp.Matrix.hstack(sp.eye(10), sp.zeros(10, 9)),
        sp.Matrix.hstack(d_matrix, -sp.eye(9)),
    )
    projection = sp.Matrix.vstack(2 * sp.eye(10), d_matrix)
    row_isometry = tau * gram_sx * tau.T == gram_sx
    column_isometry = tau.T * gram_sx * tau == gram_sx
    if not row_isometry or column_isometry or tau * tau != sp.eye(19):
        raise ArithmeticError("row/column involution convention is not uniquely certified")
    if _rows(tau) != production["global_involution"]["tau_matrix_for_row_coordinates"]:
        raise ArithmeticError("independent involution differs from production")
    if _rows(projection) != production["global_involution"]["projection_matrix"]:
        raise ArithmeticError("independent invariant projection differs from production")
    invariant_basis = sp.Matrix.hstack(sp.eye(10), sp.zeros(10, 9))
    anti_basis = sp.Matrix.hstack(-d_matrix, 2 * sp.eye(9))
    if (
        invariant_basis * tau != invariant_basis
        or anti_basis * tau != -anti_basis
        or invariant_basis * gram_sx * anti_basis.T != sp.zeros(10, 9)
    ):
        raise ArithmeticError("integral eigenspaces failed")
    if _smith_diagonal(anti_basis) != [1] * 9:
        raise ArithmeticError("anti eigenspace is not saturated in SX")
    anti_gram = anti_basis * gram_sx * anti_basis.T
    anti_cartan = -anti_gram / 2
    _ldl_upper(anti_cartan)
    if (
        int(anti_gram.det()) != -1024
        or _smith_diagonal(anti_gram) != [2] * 8 + [4]
        or any(anti_gram[index, index] % 4 for index in range(9))
        or any(
            anti_gram[left, right] % 2
            for left in range(9)
            for right in range(left)
        )
    ):
        raise ArithmeticError("anti-lattice determinant, Smith, or norm divisibility failed")
    combined = sp.Matrix.vstack(invariant_basis, anti_basis)
    if abs(int(combined.det())) != 512:
        raise ArithmeticError("invariant/anti glue index changed")

    fiber_roots = _weyl_root_closure(anti_cartan)
    production_fiber = production["finite_lift_compiler_calibration"]["fiber_enumeration"]
    if (
        len(fiber_roots) != 242
        or _vector_hash(fiber_roots) != production_fiber["vector_set_sha256"]
    ):
        raise ArithmeticError("independent Weyl closure differs from qfminim roots")

    lift_pairs: list[list[tuple[int, ...]]] = []
    for wall_index in range(shimada_walls.rows):
        root = tuple(int(value) for value in shimada_walls.row(wall_index))
        candidates = []
        for c_vector in fiber_roots:
            c_row = sp.Matrix([c_vector])
            c_d = c_row * d_matrix
            if all((root[index] - int(c_d[index])) % 2 == 0 for index in range(10)):
                a = tuple(
                    (root[index] - int(c_d[index])) // 2 for index in range(10)
                )
                candidate = a + c_vector
                row = sp.Matrix([candidate])
                if int((row * gram_sx * row.T)[0]) == -2:
                    candidates.append(candidate)
        candidates = sorted(set(candidates))
        archived = production["finite_lift_compiler_calibration"]["lift_records"][wall_index]
        archived_pair = sorted(
            [tuple(archived["lift_0_SX_coordinates"]), tuple(archived["lift_1_SX_coordinates"])]
        )
        if candidates != archived_pair or len(candidates) != 2:
            raise ArithmeticError(f"wall {wall_index + 1} lift fiber changed")
        first, second = (sp.Matrix([candidate]) for candidate in candidates)
        invariant_sum = sp.Matrix([list(root) + [0] * 9])
        if (
            first + second != invariant_sum
            or first * tau != second
            or int((first * gram_sx * second.T)[0]) != 0
            or first * projection != sp.Matrix([root])
        ):
            raise ArithmeticError("lift pair glue/exchange/disjointness failed")
        lift_pairs.append(candidates)

    all_lifts = [lift for pair in lift_pairs for lift in pair]
    histogram: Counter[int] = Counter()
    for left in range(len(all_lifts)):
        for right in range(left + 1, len(all_lifts)):
            value = int(
                (
                    sp.Matrix([all_lifts[left]])
                    * gram_sx
                    * sp.Matrix([all_lifts[right]]).T
                )[0]
            )
            if value < 0:
                raise ArithmeticError("distinct proposed lifts intersect negatively")
            histogram[value] += 1
    if histogram != Counter({0: 246, 1: 28, 2: 2}):
        raise ArithmeticError("lift intersection histogram changed")
    selected = sp.Matrix([list(pair[0]) for pair in lift_pairs])
    selected_c = selected[:, 10:]
    anti_differences = sp.Matrix.vstack(
        *[
            2 * selected[index, :]
            - sp.Matrix([_rows(shimada_walls)[index] + [0] * 9])
            for index in range(12)
        ]
    )
    if (
        anti_differences != selected_c * anti_basis
        or _smith_diagonal(selected_c) != [1] * 9
        or _smith_diagonal(anti_differences) != [1] * 9
    ):
        raise ArithmeticError("selected anti-differences are not the full saturated eigenspace")

    # Independent exact L4 exhaustion.
    h_down = weyl * ambient * sy.T * gram_sy.inv()
    if any(value.q != 1 for value in h_down):
        raise ArithmeticError("Weyl projection is not integral")
    h_down = _matrix(_rows(h_down))
    h_sx = sp.Matrix.hstack(h_down, sp.zeros(1, 9))
    h_square = int((h_sx * gram_sx * h_sx.T)[0])
    h_pairing = gram_sx * h_sx.T
    ellipsoid = -h_square * gram_sx + 2 * h_pairing * h_pairing.T
    short_vectors, short_audit = _exact_short_vectors(ellipsoid, 378)
    production_short = production["ample_orientation_and_L4"]["bounded_enumeration"]
    if (
        len(short_vectors) != 472
        or short_audit["vector_set_sha256"] != production_short["vector_set_sha256"]
    ):
        raise ArithmeticError("independent LDL exhaustion differs from qfminim")
    positive_roots = []
    degrees: Counter[int] = Counter()
    degree_zero_roots = []
    for vector in short_vectors:
        column = sp.Matrix(vector)
        square = int((column.T * gram_sx * column)[0])
        degree = int((column.T * h_pairing)[0])
        if square != -2 or not 1 <= abs(degree) <= 7:
            if square == -2 and degree == 0:
                degree_zero_roots.append(vector)
            raise ArithmeticError("the complete L4 ellipsoid contains an unexpected vector")
        if degree > 0:
            positive_roots.append(vector)
            degrees[degree] += 1
    if degree_zero_roots or degrees != Counter({1: 20, 2: 20, 3: 24, 4: 28, 5: 36, 6: 44, 7: 64}):
        raise ArithmeticError("root-regularity or small-degree root census changed")
    bad_pairs = []
    isotropic_residual_pairs = []
    for wall_index, pair in enumerate(lift_pairs):
        for sheet, lift in enumerate(pair):
            lift_row = sp.Matrix([lift])
            degree = int((lift_row * h_pairing)[0])
            for root in positive_roots:
                root_column = sp.Matrix(root)
                root_degree = int((root_column.T * h_pairing)[0])
                if root_degree >= degree:
                    continue
                pairing = int((lift_row * gram_sx * root_column)[0])
                if pairing < 0:
                    bad_pairs.append((wall_index + 1, sheet, root))
                    residual = lift_row.T - root_column
                    if int((residual.T * gram_sx * residual)[0]) == 0:
                        isotropic_residual_pairs.append((wall_index + 1, sheet, root))
    if bad_pairs:
        raise ArithmeticError("independent L4 search found a negative component")

    # A degree-zero root would have ellipsoid norm 280 and every root with
    # 0 < H.e < 8 has norm at most 378, so the same enumeration is complete
    # for both root-regularity and L4.  If b=e+F with F^2=0, then b.e=-2;
    # hence isotropic residuals are covered without enumerating F separately.
    checks = {
        "raw_OSCAR_and_Shimada_markings_join_exactly": True,
        "row_coordinate_tau_convention_forced_and_column_convention_rejected": True,
        "integral_invariant_projection_and_half_sum_glue_replayed": True,
        "anti_lattice_is_full_saturated_E8_minus2_plus_minus4_and_rootless": True,
        "independent_Weyl_reflection_closure_matches_all_242_qfminim_fiber_roots": True,
        "all_24_lifts_and_their_nonnegative_intersections_replayed": True,
        "independent_exact_LDL_enumeration_matches_all_472_qfminim_vectors": True,
        "L4_and_root_regularity_complete_including_isotropic_residual_logic": True,
        "actual_surface_existence_is_correctly_scoped_as_a_literature_input": True,
        "no_new_polytope_or_surface_is_claimed": True,
    }

    return {
        "schema": "enriques_type_i_explicit_lift_package_hostile_audit_v1",
        "verdict": "PASS/FREEZE",
        "pinned_inputs": {
            "hashes": actual_hashes,
            "oscar_loader_url": OSCAR_LOADER_URL,
            "oscar_loader_sha256_checked_during_hostile_review": OSCAR_LOADER_SHA256,
        },
        "row_column_convention": {
            "coordinate_action": "row z maps to z*tau",
            "row_isometry_identity_tau_G_tauT": row_isometry,
            "column_isometry_identity_tauT_G_tau": column_isometry,
            "column_convention_maximum_error": max(
                abs(int(value)) for value in tau.T * gram_sx * tau - gram_sx
            ),
            "tau_squared_is_identity": True,
            "projection_formula": "z*(I+tau)=(2a+cD,0)",
        },
        "anti_lattice": {
            "rank": anti_gram.rank(),
            "determinant": int(anti_gram.det()),
            "smith_diagonal": _smith_diagonal(anti_gram),
            "basis_primitive_in_SX_smith_diagonal": _smith_diagonal(anti_basis),
            "all_norms_divisible_by_4": True,
            "minus_two_root_count": 0,
            "invariant_plus_anti_index": abs(int(combined.det())),
            "selected_difference_span_smith_in_anti_basis": _smith_diagonal(selected_c),
            "selected_difference_span_primitive_in_SX": _smith_diagonal(anti_differences),
        },
        "independent_fiber_replay": {
            "engine": "simple-reflection closure, not qfminim",
            "root_count": len(fiber_roots),
            "vector_set_sha256": _vector_hash(fiber_roots),
            "exactly_two_lifts_per_wall": True,
            "all_lift_intersection_histogram": {
                str(value): count for value, count in sorted(histogram.items())
            },
        },
        "independent_L4_replay": {
            **short_audit,
            "positive_root_degree_histogram": {
                str(degree): count for degree, count in sorted(degrees.items())
            },
            "degree_zero_root_count": 0,
            "negative_component_witness_count": len(bad_pairs),
            "isotropic_residual_negative_component_witness_count": len(isotropic_residual_pairs),
            "completeness_identity": "Q_H(e)=280+2(H.e)^2 for every (-2)-root e",
            "isotropic_residual_reason": (
                "if b=e+F with e^2=-2, F^2=0, e.F=0, then b.e=-2; "
                "the enumerated negative-component gate detects e"
            ),
        },
        "geometric_provenance_and_scope": {
            "serialized_datum_kind": "lattice/Borcherds context, not surface equations",
            "oscar_1_4_loader_tuple_order": "SY, SX, L26, w, u",
            "oscar_fourth_entry_is_weyl_vector": True,
            "actual_surface_source": (
                "Brandhorst-Shimada, Automorphism Groups of Certain Enriques "
                "Surfaces, Section 7.4, No. 172"
            ),
            "literature_claim_used": (
                "a very general Type-I Enriques surface exists; its K3 cover has "
                "Picard rank 19; it has exactly 12 smooth rational curves; its "
                "nef chamber is the 12A L26/L10(2) chamber"
            ),
            "machine_certificate_alone_constructs_surface_equations": False,
            "actual_surface_identification_supplies_the_bridge_from_SX_in_II_1_25_to_NS_of_a_K3_cover": True,
            "new_abstract_polytope_claimed": False,
            "new_surface_claimed": False,
            "remaining_compiler_problem": (
                "recover the embedding/glue package from downstairs data without "
                "starting from the literature-backed OSCAR SX"
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(analyze(), indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.is_file() or args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived hostile audit differs: {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(hashlib.sha256(rendered.encode("utf-8")).hexdigest())


if __name__ == "__main__":
    main()
