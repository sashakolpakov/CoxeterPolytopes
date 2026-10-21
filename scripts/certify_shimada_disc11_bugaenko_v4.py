"""Certify the Bugaenko fixed slice and its Klein-four K3 stabilizer.

This is a deliberately standalone exact audit.  In the discriminant-11
golden marking, the enhanced involution ``A`` fixes the hyperplane ``x_1=0``.
The induced lattice is Bugaenko's ``<-2 phi> + E8``.  Reversing Bugaenko's
eleven-wall diagram has two integral lifts to the H11 lattice (the two signs
on the anti-invariant line).  We fold both lifts through genuine K3 root
walls using Shimada's SHA-pinned induced-chamber atlas, and check directly
that the two chamber residuals, together with ``A``, form a Klein four group.

The executable certificate proves the arithmetic and chamber statements.
The short geometric argument converting these checks into
``C_Aut(X,tau)(A) = V4`` is recorded in the companion result note.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

from sympy import Matrix, symbols
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.matrix import inertia
from coxeter_polytopes.exact.trace_lattice import (
    INV_SQRT5,
    integral_gram,
    trace_gram,
)
from coxeter_polytopes.geometric.golden_halfspaces import (
    golden_dot,
    golden_quadratic_value,
    inverse_golden_matrix,
)
from coxeter_polytopes.known_examples.bugaenko import (
    bugaenko_basis_gram,
    bugaenko_coxeter_relations,
    bugaenko_h8_certificate,
    bugaenko_simple_roots,
)
from coxeter_polytopes.verify.vinberg import verify_bugaenko_vinberg_completion
from scripts.coupled_discriminant11_candidate import field_gram, multiplication_by_phi
from scripts.explore_shimada_disc11_global_eigen_nef import (
    _Q_GRAM,
    _physical_state_key,
    _restricted_rows,
)
from scripts.fold_h11_seed_via_shimada import (
    GoldenMatrix,
    GoldenVector,
    LocationResult,
    _identity,
    _isometry_verified,
    _matrix_product,
    _reflection_in_functional,
    _row_product,
    _transpose,
    _walk_until_boundary_or_interior,
)
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas
from scripts.shimada_disc11_conjugate import (
    AMPLE_FIELD_COEFFICIENTS,
    DISC11_SOURCE_SHA256,
)


TRANSCRIPT = (
    Path(__file__).parents[1]
    / "data"
    / "known_examples"
    / "bugaenko_vinberg_transcript.json"
)

# Coordinates in the fixed hyperplane are (x_0,x_2,...,x_9).  Bottinelli's
# E8 basis for the archived Bugaenko roots differs from the H11 E8 basis by
# this permutation.  If x_std = x_bug Q, then S_std = Q^t S_bug Q.
BUGAENKO_TO_H11_FIXED = (0, 1, 8, 2, 3, 4, 5, 6, 7)
FIXED_POSITIONS = (0, 2, 3, 4, 5, 6, 7, 8, 9)

# The unique nontrivial automorphism of Bugaenko's eleven-node diagram.
BUGAENKO_REVERSAL = (8, 7, 6, 5, 4, 3, 2, 1, 0, 10, 9)

EXPECTED_FIELD_MATRIX_SHA256 = {
    "A": "5b6e9046cbbaf4cbffaf954c98a0332b169e3d246b851af9adeb43a8468cdefe",
    "B_plus": "842661cdc5ca7c1039e90f449d64f009620bc6ee0a78560f322cb82fa9761abc",
    "B_minus": "f4009609eb740c47a5a5934af8e3ef825703df7ab21d7369ecb4d9e54ecf8bac",
}
EXPECTED_FOLD_WORD_SHA256 = {
    "B_plus": "a7120bb8d2f1466a1156bd6f7d6f827a7472a99ce704c495053cd14dd7cee3cb",
    "B_minus": "68007a8b84a28c6fa2d1fbb66367b6d61be9c6bfb741f15ac80a8487533d8fc1",
}


def _canonical_sha(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _encoded_matrix(matrix: Sequence[Sequence[Golden]]) -> list[list[dict[str, str]]]:
    return [[entry.to_json() for entry in row] for row in matrix]


def _encoded_vector(vector: Sequence[Golden]) -> list[dict[str, str]]:
    return [entry.to_json() for entry in vector]


def _golden_matrix(rows: Sequence[Sequence[Golden | int]]) -> GoldenMatrix:
    return tuple(tuple(Golden.coerce(entry) for entry in row) for row in rows)


def _permutation_matrix(images: Sequence[int]) -> GoldenMatrix:
    size = len(images)
    if sorted(images) != list(range(size)):
        raise ValueError("not a permutation")
    return tuple(
        tuple(Golden(int(images[i] == j)) for j in range(size))
        for i in range(size)
    )


def _diagram_automorphisms(relations: Sequence[Sequence[int | str]]) -> tuple[tuple[int, ...], ...]:
    """Enumerate labeled complete-graph automorphisms by exact backtracking."""

    size = len(relations)
    signatures = tuple(
        tuple(sorted((str(value) for j, value in enumerate(row) if j != i)))
        for i, row in enumerate(relations)
    )
    candidates = tuple(
        tuple(j for j in range(size) if signatures[j] == signatures[i])
        for i in range(size)
    )
    order = tuple(sorted(range(size), key=lambda i: (len(candidates[i]), i)))
    image = [-1] * size
    used: set[int] = set()
    answers: list[tuple[int, ...]] = []

    def visit(depth: int) -> None:
        if depth == size:
            answers.append(tuple(image))
            return
        source = order[depth]
        for target in candidates[source]:
            if target in used:
                continue
            if any(
                image[other] >= 0
                and relations[source][other] != relations[target][image[other]]
                for other in range(size)
            ):
                continue
            image[source] = target
            used.add(target)
            visit(depth + 1)
            used.remove(target)
            image[source] = -1

    visit(0)
    return tuple(sorted(answers))


def _bugaenko_reversal_isometry() -> tuple[GoldenMatrix, GoldenMatrix]:
    roots = _golden_matrix(bugaenko_simple_roots())
    basis = roots[:9]
    target = tuple(roots[BUGAENKO_REVERSAL[i]] for i in range(9))
    inverse = inverse_golden_matrix(basis)
    reversal_bug = _matrix_product(inverse, target)
    permutation = _permutation_matrix(BUGAENKO_TO_H11_FIXED)
    reversal_std = _matrix_product(
        _matrix_product(_transpose(permutation), reversal_bug), permutation
    )
    return reversal_bug, reversal_std


def _field_involution_a() -> GoldenMatrix:
    result = [list(row) for row in _identity(10)]
    result[1][0] = 1 - PHI
    result[1][1] = Golden(-1)
    return tuple(tuple(row) for row in result)


def _lift_reversal(reversal: GoldenMatrix, anti_sign: int) -> GoldenMatrix:
    if anti_sign not in (-1, 1):
        raise ValueError("anti_sign must be +/-1")
    result = [[Golden(0) for _ in range(10)] for _ in range(10)]
    for i, source in enumerate(FIXED_POSITIONS):
        for j, target in enumerate(FIXED_POSITIONS):
            result[source][target] = reversal[i][j]
    result[1][1] = Golden(anti_sign)
    anti_coefficient = PHI - 1
    for j, target in enumerate(FIXED_POSITIONS):
        desired = Golden(anti_sign) * anti_coefficient if j == 0 else Golden(0)
        result[1][target] = (desired - anti_coefficient * reversal[0][j]) / 2
    return tuple(tuple(row) for row in result)


def _fixed_restriction(matrix: GoldenMatrix) -> GoldenMatrix:
    return tuple(
        tuple(matrix[i][j] for j in FIXED_POSITIONS) for i in FIXED_POSITIONS
    )


def _underlying_integer_row_action(matrix: GoldenMatrix) -> Matrix:
    """Expand a right O_F-row action in the interleaved (e,phi e) basis."""

    rank = len(matrix)
    result = [[0 for _ in range(2 * rank)] for _ in range(2 * rank)]
    for i in range(rank):
        for j in range(rank):
            entry = matrix[i][j]
            if not entry.is_integral:
                raise ValueError("field action is not integral over Z[phi]")
            a, b = int(entry.a), int(entry.b)
            result[2 * i][2 * j] = a
            result[2 * i][2 * j + 1] = b
            result[2 * i + 1][2 * j] = b
            result[2 * i + 1][2 * j + 1] = a + b
    return Matrix(result)


def _trace_lattice_gram() -> Matrix:
    alpha = -PHI * INV_SQRT5
    return Matrix(integral_gram(trace_gram(field_gram(), alpha)))


def _discriminant_sign(action: Matrix, gram: Matrix) -> int:
    inverse = gram.inv()
    generator = next(
        (inverse.row(i) for i in range(inverse.rows) if any(not x.is_Integer for x in inverse.row(i))),
        None,
    )
    if generator is None:
        raise ValueError("trace lattice unexpectedly unimodular")
    possible = []
    for sign in (-1, 1):
        difference = generator * action - sign * generator
        if all(entry.is_Integer for entry in difference):
            possible.append(sign)
    if len(possible) != 1:
        raise ValueError(f"discriminant action is not a unique sign: {possible}")
    return possible[0]


def _a_fixed_lattice_audit(gram: Matrix) -> dict[str, Any]:
    """Compute the 2-elementary invariant lattice and Nikulin fixed-locus data."""

    # A fixes exactly the O_F hyperplane x_1=0.  In the interleaved trace
    # basis this removes coordinates 2 and 3.
    fixed_indices = (0, 1, *range(4, 20))
    fixed = gram.extract(fixed_indices, fixed_indices)
    fixed_smith = smith_normal_form(fixed, domain=ZZ)
    fixed_inertia = inertia(
        [[int(fixed[i, j]) for j in range(fixed.cols)] for i in range(fixed.rows)]
    )
    inverse = fixed.inv()
    delta = int(any(not inverse[i, i].is_Integer for i in range(inverse.rows)))

    # The anti-invariant O_F-line is generated by (phi-1,2,0^8).
    anti_basis = Matrix(
        [
            [-1, 1, 2, 0, *([0] * 16)],
            [1, 0, 0, 2, *([0] * 16)],
        ]
    )
    anti = anti_basis * gram * anti_basis.T
    anti_smith = smith_normal_form(anti, domain=ZZ)
    invariants = [abs(int(fixed_smith[i, i])) for i in range(fixed.rows)]
    a_length = sum(value == 2 for value in invariants)
    rank = fixed.rows
    genus = (22 - rank - a_length) // 2
    rational_curves = (rank - a_length) // 2
    return {
        "rank": rank,
        "signature": [fixed_inertia.positive, fixed_inertia.negative],
        "determinant": str(fixed.det()),
        "smith_invariants": invariants,
        "nikulin_invariants_r_a_delta": [rank, a_length, delta],
        "isometry_type": "U + E8(-1) + E7(-1) + A1(-1)",
        "anti_invariant_rank": anti.rows,
        "anti_invariant_determinant": str(anti.det()),
        "anti_invariant_smith_invariants": [
            abs(int(anti_smith[i, i])) for i in range(anti.rows)
        ],
        "fixed_locus": {
            "unique_positive_genus_component_genus": genus,
            "rational_curve_components": rational_curves,
        },
        "verified": (
            rank == 18
            and (fixed_inertia.positive, fixed_inertia.negative) == (1, 17)
            and fixed.det() == -4
            and invariants == [1] * 16 + [2, 2]
            and (rank, a_length, delta) == (18, 2, 1)
            and anti.det() == 44
            and [abs(int(anti_smith[i, i])) for i in range(anti.rows)] == [2, 22]
            and genus == 1
            and rational_curves == 8
        ),
    }


def _fold_with_word(
    atlas: ShimadaDisc11Atlas,
    seed: GoldenMatrix,
    *,
    max_root_reflections: int = 100,
    max_gallery: int = 20_000,
) -> tuple[GoldenMatrix, GoldenVector, LocationResult, tuple[dict[str, Any], ...]]:
    point = _row_product(tuple(AMPLE_FIELD_COEFFICIENTS), seed)
    residual = seed
    word: list[dict[str, Any]] = []
    terminal: LocationResult | None = None
    for index in range(max_root_reflections + 1):
        location = _walk_until_boundary_or_interior(
            atlas, point, max_gallery=max_gallery
        )
        if location.kind == "interior":
            terminal = location
            break
        if location.kind != "root_boundary" or location.boundary_functional is None:
            raise RuntimeError(f"fold stopped at {location.kind} after {index} roots")
        functional = location.boundary_functional
        word.append(
            {
                "index": index,
                "representative": location.state.index,
                "physical_state_sha256": _canonical_sha(_physical_state_key(atlas, location.state)),
                "row_transport_sha256": _canonical_sha(location.state.row_isometry),
                "boundary_label": list(location.boundary_label or ()),
                "boundary_functional": _encoded_vector(functional),
                "auxiliary_gallery": [list(step) for step in location.gallery],
            }
        )
        reflection = _reflection_in_functional(functional)
        point = _row_product(point, reflection)
        residual = _matrix_product(residual, reflection)
    if terminal is None:
        raise RuntimeError("root reflection limit reached")
    return residual, point, terminal, tuple(word)


def _terminal_audit(atlas: ShimadaDisc11Atlas, point: GoldenVector, terminal: LocationResult) -> dict[str, Any]:
    rows = _restricted_rows(atlas, terminal.state)
    signs = tuple(golden_dot(row, point).sign() for row in rows.values())
    return {
        "kind": terminal.kind,
        "representative": terminal.state.index,
        "physical_state_sha256": _canonical_sha(_physical_state_key(atlas, terminal.state)),
        "row_transport_sha256": _canonical_sha(terminal.state.row_isometry),
        "auxiliary_gallery": [list(step) for step in terminal.gallery],
        "all_wall_pairings_strictly_positive": all(sign > 0 for sign in signs),
        "timelike": golden_quadratic_value(_Q_GRAM, point).sign() > 0,
    }


def _matrix_audit(name: str, field_action: GoldenMatrix, gram: Matrix) -> dict[str, Any]:
    integer = _underlying_integer_row_action(field_action)
    variable = symbols("t")
    encoded = _encoded_matrix(field_action)
    phi = Matrix(multiplication_by_phi(10))
    return {
        "name": name,
        "field_row_action_sha256": _canonical_sha(encoded),
        "integral_row_action_sha256": _canonical_sha(
            [[int(entry) for entry in row] for row in integer.tolist()]
        ),
        "integral_over_Zphi": all(entry.is_integral for row in field_action for entry in row),
        "field_lattice_isometry": _isometry_verified(field_action),
        "trace_lattice_isometry": integer * gram * integer.T == gram,
        "commutes_with_tau": integer * phi == phi * integer,
        "determinant": str(integer.det()),
        "order_two": integer * integer == Matrix.eye(20),
        "NS_characteristic_polynomial": str(integer.charpoly(variable).as_expr().factor()),
        "discriminant_action": _discriminant_sign(integer, gram),
        "field_trace": sum((field_action[i][i] for i in range(10)), Golden(0)).to_json(),
        "field_row_action": encoded,
    }


def certify(path: Path, *, include_words: bool = False) -> dict[str, Any]:
    source_digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if source_digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")

    transcript = json.loads(TRANSCRIPT.read_text())
    completion = verify_bugaenko_vinberg_completion(
        transcript, bugaenko_h8_certificate()
    )
    relations = bugaenko_coxeter_relations()
    automorphisms = _diagram_automorphisms(relations)
    reversal_bug, reversal_std = _bugaenko_reversal_isometry()
    roots_bug = _golden_matrix(bugaenko_simple_roots())
    bug_gram = _golden_matrix(bugaenko_basis_gram())
    h11 = _golden_matrix(field_gram())
    fixed_gram = tuple(
        tuple(h11[i][j] for j in FIXED_POSITIONS) for i in FIXED_POSITIONS
    )
    permutation = _permutation_matrix(BUGAENKO_TO_H11_FIXED)
    expected_fixed_gram = _matrix_product(
        _matrix_product(_transpose(permutation), bug_gram), permutation
    )
    reversed_roots = tuple(_row_product(root, reversal_bug) for root in roots_bug)

    a = _field_involution_a()
    lift_plus = _lift_reversal(reversal_std, 1)
    lift_minus = _lift_reversal(reversal_std, -1)
    anti_vector = (PHI - 1, Golden(2), *(Golden(0) for _ in range(8)))
    lift_checks = {
        "fixed_gram_equals_bugaenko_after_permutation": fixed_gram == expected_fixed_gram,
        "diagram_automorphisms": [list(item) for item in automorphisms],
        "diagram_automorphism_group_is_C2": set(automorphisms) == {
            tuple(range(11)), BUGAENKO_REVERSAL
        },
        "reversal_maps_all_eleven_roots": all(
            reversed_roots[i] == roots_bug[BUGAENKO_REVERSAL[i]] for i in range(11)
        ),
        "reversal_isometry": _matrix_product(
            _matrix_product(reversal_bug, bug_gram), _transpose(reversal_bug)
        ) == bug_gram,
        "reversal_order_two": _matrix_product(reversal_bug, reversal_bug) == _identity(9),
        "lifts_integral": all(
            entry.is_integral
            for lift in (lift_plus, lift_minus)
            for row in lift
            for entry in row
        ),
        "lifts_are_H11_isometries": all(
            _isometry_verified(lift) for lift in (lift_plus, lift_minus)
        ),
        "lifts_order_two": all(
            _matrix_product(lift, lift) == _identity(10)
            for lift in (lift_plus, lift_minus)
        ),
        "lifts_commute_with_A": all(
            _matrix_product(lift, a) == _matrix_product(a, lift)
            for lift in (lift_plus, lift_minus)
        ),
        "lifts_restrict_to_reversal": all(
            _fixed_restriction(lift) == reversal_std
            for lift in (lift_plus, lift_minus)
        ),
        "anti_line_signs": (
            _row_product(anti_vector, lift_plus) == anti_vector
            and _row_product(anti_vector, lift_minus)
            == tuple(-entry for entry in anti_vector)
        ),
    }

    atlas = ShimadaDisc11Atlas(path)
    plus, plus_point, plus_terminal, plus_word = _fold_with_word(atlas, lift_plus)
    minus, minus_point, minus_terminal, minus_word = _fold_with_word(atlas, lift_minus)
    plus_terminal_audit = _terminal_audit(atlas, plus_point, plus_terminal)
    minus_terminal_audit = _terminal_audit(atlas, minus_point, minus_terminal)

    gram = _trace_lattice_gram()
    fixed_lattice = _a_fixed_lattice_audit(gram)
    audit_a = _matrix_audit("A", a, gram)
    audit_plus = _matrix_audit("B_plus", plus, gram)
    audit_minus = _matrix_audit("B_minus", minus, gram)
    products = {
        "plus_commutes_with_A": _matrix_product(plus, a) == _matrix_product(a, plus),
        "minus_commutes_with_A": _matrix_product(minus, a) == _matrix_product(a, minus),
        "plus_minus_equals_A": _matrix_product(plus, minus) == a,
        "minus_plus_equals_A": _matrix_product(minus, plus) == a,
        "four_elements_distinct": len({_identity(10), a, plus, minus}) == 4,
        # Folding changes the archived Bugaenko chamber to the one containing
        # the A-fixed slice of the selected H11 chamber.  Thus the residual
        # restrictions need not equal ``reversal_std`` in the archived
        # coordinates.  They are the same nontrivial diagram involution of
        # that conjugate chamber: equality follows from B_plus B_minus=A,
        # while nontriviality follows directly from the matrices.
        "fixed_restrictions_equal": _fixed_restriction(plus)
        == _fixed_restriction(minus),
        "plus_fixed_restriction_nonidentity": _fixed_restriction(plus)
        != _identity(9),
        "minus_fixed_restriction_nonidentity": _fixed_restriction(minus)
        != _identity(9),
    }

    word_summary: dict[str, Any] = {
        "plus_root_count": len(plus_word),
        "minus_root_count": len(minus_word),
        "plus_exact_word_sha256": _canonical_sha(plus_word),
        "minus_exact_word_sha256": _canonical_sha(minus_word),
        "plus_total_auxiliary_steps": sum(len(item["auxiliary_gallery"]) for item in plus_word),
        "minus_total_auxiliary_steps": sum(len(item["auxiliary_gallery"]) for item in minus_word),
    }
    if include_words:
        word_summary["plus_exact_word"] = plus_word
        word_summary["minus_exact_word"] = minus_word

    k3_extensions = {
        item["name"]: (
            "extend by +I on T11" if item["discriminant_action"] == 1
            else "extend by -I on T11"
        )
        for item in (audit_a, audit_plus, audit_minus)
    }
    expected_hashes = {
        "field_matrices": {
            item["name"]: item["field_row_action_sha256"]
            == EXPECTED_FIELD_MATRIX_SHA256[item["name"]]
            for item in (audit_a, audit_plus, audit_minus)
        },
        "fold_words": {
            "B_plus": word_summary["plus_exact_word_sha256"]
            == EXPECTED_FOLD_WORD_SHA256["B_plus"],
            "B_minus": word_summary["minus_exact_word_sha256"]
            == EXPECTED_FOLD_WORD_SHA256["B_minus"],
        },
    }
    verified = (
        completion.chamber.cusps == ()
        and completion.chamber.facet_count == 11
        and all(lift_checks.values())
        and plus_terminal_audit["all_wall_pairings_strictly_positive"]
        and plus_terminal_audit["timelike"]
        and minus_terminal_audit["all_wall_pairings_strictly_positive"]
        and minus_terminal_audit["timelike"]
        and len(plus_word) == len(minus_word) == 88
        and all(
            item[key]
            for item in (audit_a, audit_plus, audit_minus)
            for key in (
                "integral_over_Zphi",
                "field_lattice_isometry",
                "trace_lattice_isometry",
                "commutes_with_tau",
                "order_two",
            )
        )
        and all(products.values())
        and fixed_lattice["verified"]
        and all(expected_hashes["field_matrices"].values())
        and all(expected_hashes["fold_words"].values())
    )
    return {
        "verified": verified,
        "source_sha256": source_digest,
        "bugaenko_fixed_slice": {
            "complete_exact_vinberg_enumeration": len(completion.search.roots) == 11,
            "compact": completion.chamber.cusps == (),
            "facet_count": completion.chamber.facet_count,
            "finite_vertex_count": len(completion.chamber.finite_vertices),
            **lift_checks,
        },
        "folds": {
            "exact_words": word_summary,
            "plus_terminal": plus_terminal_audit,
            "minus_terminal": minus_terminal_audit,
        },
        "residuals": [audit_a, audit_plus, audit_minus],
        "A_fixed_lattice": fixed_lattice,
        "multiplication_table_checks": products,
        "archived_hash_checks": expected_hashes,
        "K3_Hodge_extensions": k3_extensions,
        "theorem": {
            "centralizer_of_A_in_Aut_X_tau": "C2 x C2",
            "order": 4,
            "finiteness_criterion": "Aut(X,tau) finite iff the conjugacy orbit of A is finite",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    parser.add_argument("--include-words", action="store_true")
    args = parser.parse_args()
    print(json.dumps(certify(args.disc11, include_words=args.include_words), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
