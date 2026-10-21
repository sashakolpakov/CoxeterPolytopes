"""Filter Shimada's published discriminant-11 K3 data through the RM marking.

The external ``disc11.txt`` file is not vendored.  It is available in
Shimada's ``threediscs.zip`` computational supplement.  This script parses
only ``GramMatS`` and ``B``: the Neron--Severi Gram matrix and the finite set
of (-2)-vectors found by the automorphism-group computation.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
from collections import Counter
from pathlib import Path

from sympy import Matrix

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.matrix import rank
from coxeter_polytopes.exact.trace_lattice import INV_SQRT5, integral_gram, trace_gram
from coxeter_polytopes.geometric.k3 import (
    K3_TRACE_SCALAR,
    field_inner_product,
    is_integral_reflection_root,
    is_k3_wall_norm,
)
from coxeter_polytopes.lattices import block_diagonal, e8_cartan_gram
from scripts.coupled_discriminant11_candidate import (
    field_gram,
    multiplication_by_phi,
    trace_basis_change,
)


# Two zero-based entries of Shimada's published set B.  They are retained as
# source-coordinate fingerprints so that reordering the external file cannot
# silently change the obstruction.
OPPOSITE_WALL_PLUS = [0, 1, *([0] * 18)]
OPPOSITE_WALL_MINUS = [
    1, 0, 0, 0, -3, -2, -4, -6, -5, -4, -3, -2, *([0] * 8)
]


def matrix_product(left, right):
    return [
        [sum(left[i][k] * right[k][j] for k in range(len(right))) for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def matrix_vector(matrix, vector):
    return [sum(row[j] * vector[j] for j in range(len(vector))) for row in matrix]


def shimada_basis_change() -> list[list[int]]:
    """Map Shimada's NS coordinates to the original interleaved trace basis."""

    # First map Shimada's U_phi + (-T_red) + E8(-1)^2 basis to our
    # U + (-T11) + E8(-1)^2 basis.
    change = [[int(i == j) for j in range(20)] for i in range(20)]
    change[0][0], change[0][1] = 1, -1
    change[1][0], change[1][1] = 0, 1
    change[2][2], change[2][3] = 0, 1
    change[3][2], change[3][3] = 1, -2
    # Shimada numbers the same E8 Dynkin diagram along different arms.
    shimada_to_ours = (7, 0, 1, 2, 3, 4, 5, 6)
    for offset in (4, 12):
        for i in range(8):
            for j in range(8):
                change[offset + i][offset + j] = 0
        for shimada_column, our_row in enumerate(shimada_to_ours):
            change[offset + our_row][offset + shimada_column] = 1
    # Then map our standard block basis to e_i, phi*e_i.
    return matrix_product(trace_basis_change(), change)


def expected_shimada_gram() -> list[list[int]]:
    negative_e8 = [[-entry for entry in row] for row in e8_cartan_gram()]
    permutation = [[0 for _ in range(8)] for _ in range(8)]
    for shimada_column, our_row in enumerate((7, 0, 1, 2, 3, 4, 5, 6)):
        permutation[our_row][shimada_column] = 1
    shimada_e8 = matrix_product(
        transpose(permutation), matrix_product(negative_e8, permutation)
    )
    return block_diagonal(
        [[0, 1], [1, -2]],
        [[-2, -1], [-1, -6]],
        shimada_e8,
        shimada_e8,
    )


def _parse_outer_collection(text: str):
    # Maple inserts a backslash when a very long assignment wraps mid-token.
    value = text.replace("\\", "").strip()
    if value.endswith(";"):
        value = value[:-1]
    if value.startswith("{") and value.endswith("}"):
        value = "[" + value[1:-1] + "]"
    return ast.literal_eval(value)


def read_shimada_data(
    path: Path,
) -> tuple[list[list[int]], list[list[int]], list[list[list[int]]]]:
    """Read the NS form, curve roots, and automorphisms from ``disc11.txt``."""

    wanted = {"GramMatS", "B", "Gamma"}
    parsed = {}
    pending_name = None
    pending_value: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if pending_name is not None:
                pending_value.append(line.strip())
                if line.rstrip().endswith(";"):
                    parsed[pending_name] = _parse_outer_collection("".join(pending_value))
                    pending_name = None
                    pending_value = []
                    if parsed.keys() == wanted:
                        break
                continue
            if ":=" not in line:
                continue
            name, value = line.split(":=", 1)
            if name in wanted:
                if line.rstrip().endswith(";"):
                    parsed[name] = _parse_outer_collection(value)
                    if parsed.keys() == wanted:
                        break
                else:
                    pending_name = name
                    pending_value = [value.strip()]
    missing = wanted - parsed.keys()
    if missing:
        raise ValueError(f"missing assignments in Shimada data: {sorted(missing)}")
    return parsed["GramMatS"], parsed["B"], parsed["Gamma"]


def field_vector_from_shimada(root: list[int]) -> list[Golden]:
    if len(root) != 20:
        raise ValueError("a Shimada NS vector must have 20 coordinates")
    interleaved = matrix_vector(shimada_basis_change(), root)
    return [Golden(interleaved[2 * i], interleaved[2 * i + 1]) for i in range(10)]


def eigen_nef_opposite_wall_obstruction(roots: list[list[int]]) -> dict[str, object]:
    """Certify that the fixed real eigenspace misses the K3 nef cone.

    Shimada's oriented vectors in ``B`` are effective ``(-2)``-curve classes,
    so nef vectors pair nonnegatively with them.  The two selected walls have
    timelike projections in opposite components of the Lorentzian
    identity-eigenspace.  Reverse Cauchy--Schwarz then makes one of the two nef
    inequalities negative on each component of its positive cone.
    """

    root_set = {tuple(root) for root in roots}
    required = {
        "d_plus": tuple(OPPOSITE_WALL_PLUS),
        "d_minus": tuple(OPPOSITE_WALL_MINUS),
    }
    missing = [name for name, root in required.items() if root not in root_set]
    if missing:
        return {"verified": False, "missing_published_roots": missing}

    h = field_gram()
    zero = Golden(0)
    e = [Golden(1), *([zero] * 9)]
    d_plus = field_vector_from_shimada(OPPOSITE_WALL_PLUS)
    d_minus = field_vector_from_shimada(OPPOSITE_WALL_MINUS)

    def h_pair(left, right):
        return field_inner_product(h, left, right)

    def q_pair(left, right):
        return K3_TRACE_SCALAR * h_pair(left, right)

    q_e_e = q_pair(e, e)
    q_plus_plus = q_pair(d_plus, d_plus)
    q_minus_minus = q_pair(d_minus, d_minus)
    q_e_plus = q_pair(e, d_plus)
    q_e_minus = q_pair(e, d_minus)
    verified = (
        q_e_e.sign() > 0
        and q_plus_plus.sign() > 0
        and q_minus_minus.sign() > 0
        and q_e_plus.sign() > 0
        and q_e_minus.sign() < 0
        and q_plus_plus.trace() == -2
        and q_minus_minus.trace() == -2
    )
    return {
        "verified": verified,
        "conclusion": "identity_eigenspace_positive_cone_intersection_with_nef_is_empty",
        "orientation_source": "Shimada2015 Algorithm 6.1 and Proposition 6.3(3)",
        "timecone_argument": "reverse_Lorentz_Cauchy_Schwarz",
        "control_vector": [entry.to_json() for entry in e],
        "q_e_e": q_e_e.to_json(),
        "walls": {
            "d_plus": {
                "shimada_coordinates": OPPOSITE_WALL_PLUS,
                "field_coordinates": [entry.to_json() for entry in d_plus],
                "h_norm": h_pair(d_plus, d_plus).to_json(),
                "q_norm": q_plus_plus.to_json(),
                "q_pairing_with_e": q_e_plus.to_json(),
                "trace_square": str(q_plus_plus.trace()),
            },
            "d_minus": {
                "shimada_coordinates": OPPOSITE_WALL_MINUS,
                "field_coordinates": [entry.to_json() for entry in d_minus],
                "h_norm": h_pair(d_minus, d_minus).to_json(),
                "q_norm": q_minus_minus.to_json(),
                "q_pairing_with_e": q_e_minus.to_json(),
                "trace_square": str(q_minus_minus.trace()),
            },
        },
    }


def analyze_data(
    shimada_gram: list[list[int]],
    roots: list[list[int]],
    automorphisms: list[list[list[int]]] | None = None,
) -> dict[str, object]:
    h = field_gram()
    trace = integral_gram(trace_gram(h, -PHI * INV_SQRT5))
    change = shimada_basis_change()
    transformed = matrix_product(transpose(change), matrix_product(trace, change))
    expected = expected_shimada_gram()

    norms: list[Golden] = []
    relevant: list[tuple[list[int], list[Golden], Golden]] = []
    integral: list[tuple[list[int], list[Golden], Golden]] = []
    trace_root_count = 0
    for root in roots:
        vector = field_vector_from_shimada(root)
        norm = field_inner_product(h, vector, vector)
        norms.append(norm)
        if (-PHI * INV_SQRT5 * norm).trace() == -2:
            trace_root_count += 1
        if is_k3_wall_norm(norm):
            relevant.append((root, vector, norm))
            if is_integral_reflection_root(h, vector):
                integral.append((root, vector, norm))

    distribution = Counter(str(norm) for _, _, norm in relevant)
    integral_distribution = Counter(str(norm) for _, _, norm in integral)
    integral_vectors = [vector for _, vector, _ in integral]
    pair_signs: Counter[str] = Counter()
    for i, left in enumerate(integral_vectors):
        for right in integral_vectors[i + 1 :]:
            sign = field_inner_product(h, left, right).sign()
            pair_signs[{1: "positive", 0: "zero", -1: "negative"}[sign]] += 1

    expected_seed = []
    for scalar in (Golden(1), PHI - 1):
        for coordinate in range(2, 10):
            vector = [Golden(0)] * 10
            vector[coordinate] = scalar
            expected_seed.append(tuple(vector))
    expected_seed.append((Golden(1), Golden(1), *([Golden(0)] * 8)))
    result = {
        "shimada_gram_matches_published_normal_form": shimada_gram == expected,
        "basis_change_recovers_shimada_gram": transformed == expected,
        "published_root_representatives": len(roots),
        "trace_minus_two_verified": trace_root_count,
        "positive_at_both_places": len(relevant),
        "positive_shell_distribution": dict(sorted(distribution.items())),
        "field_integral_reflection_roots": len(integral),
        "field_integral_shell_distribution": dict(sorted(integral_distribution.items())),
        "field_span_rank": rank(integral_vectors),
        "pair_sign_distribution": dict(sorted(pair_signs.items())),
        "seed_matches_two_e8_bases_plus_one_coupled_root": (
            set(map(tuple, integral_vectors)) == set(expected_seed)
        ),
        "field_integral_root_seed": [
            {
                "shimada_coordinates": root,
                "field_coordinates": [entry.to_json() for entry in vector],
                "field_norm": norm.to_json(),
            }
            for root, vector, norm in integral
        ],
        "eigen_nef_opposite_timecone_obstruction": (
            eigen_nef_opposite_wall_obstruction(roots)
        ),
    }

    if automorphisms is not None:
        # Shimada acts on row vectors from the right.  Apply every published
        # generator once to every surviving effective curve.  This is a
        # reproducible seed expansion, not a completeness claim.
        seed_roots = [root for root, _, _ in integral]
        basis_change = Matrix(shimada_basis_change())
        phi_action = basis_change.inv() * Matrix(multiplication_by_phi(10)) * basis_change
        if any(not entry.is_Integer for entry in phi_action):
            raise ValueError("the phi action is not integral in Shimada's basis")
        phi_action = [[int(entry) for entry in row] for row in phi_action.tolist()]
        centralizing_generators = 0
        images = {tuple(root) for root in seed_roots}
        for automorphism in automorphisms:
            if len(automorphism) != 20 or any(len(row) != 20 for row in automorphism):
                raise ValueError("a Shimada automorphism must be a 20 by 20 matrix")
            column_action = transpose(automorphism)
            if matrix_product(column_action, phi_action) == matrix_product(
                phi_action, column_action
            ):
                centralizing_generators += 1
            for root in seed_roots:
                images.add(
                    tuple(
                        sum(root[i] * automorphism[i][j] for i in range(20))
                        for j in range(20)
                    )
                )

        expanded_vectors: list[list[Golden]] = []
        expanded_norms: list[Golden] = []
        expanded_trace_roots = 0
        for image in images:
            vector = field_vector_from_shimada(list(image))
            norm = field_inner_product(h, vector, vector)
            if (-PHI * INV_SQRT5 * norm).trace() == -2:
                expanded_trace_roots += 1
            if is_k3_wall_norm(norm) and is_integral_reflection_root(h, vector):
                expanded_vectors.append(vector)
                expanded_norms.append(norm)

        result.update(
            {
                "published_automorphism_generators": len(automorphisms),
                "published_generators_centralizing_phi": centralizing_generators,
                "one_step_distinct_curve_images": len(images),
                "one_step_trace_minus_two_verified": expanded_trace_roots,
                "one_step_field_integral_reflection_roots": len(expanded_vectors),
                "one_step_field_integral_shell_distribution": dict(
                    sorted(Counter(str(norm) for norm in expanded_norms).items())
                ),
                "one_step_field_span_rank": rank(expanded_vectors),
            }
        )
    return result


def analyze_file(path: Path) -> dict[str, object]:
    gram, roots, automorphisms = read_shimada_data(path)
    result = analyze_data(gram, roots, automorphisms)
    result["source_file"] = str(path)
    result["source_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path, help="path to Shimada's disc11.txt")
    parser.add_argument(
        "--summary-only", action="store_true", help="omit the 17 explicit seed vectors"
    )
    args = parser.parse_args()
    result = analyze_file(args.disc11)
    if args.summary_only:
        result.pop("field_integral_root_seed")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
