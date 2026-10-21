"""Exact reconstruction of Bugaenko's eleven-wall chamber.

The ambient form and primitive roots are translated from the pinned
``VinbergsAlgorithmNF`` Bugaenko regression fixture into the distinguished
embedding where ``phi=(1+sqrt(5))/2``.  The roots are reordered to match the
published horizontal-chain diagram and the final three are multiplied by the
unit ``phi`` so that all eleven roots have norm two.  This gives products
``-1`` on simple edges, ``-phi`` on the two 5-edges, and ``-2*phi`` on the
dotted edge.

This module reconstructs and certifies the *listed polytope*.  It does not yet
replace the separate Vinberg transcript required to prove that the list is the
complete set of simple roots of Bugaenko's ambient reflection group.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Sequence

from coxeter_polytopes.exact.coxeter import classify_root_pair
from coxeter_polytopes.exact.golden import Golden, PHI, ZERO


def bugaenko_basis_gram() -> list[list[Golden]]:
    """The published lattice ``<-2*phi> + E8`` in Bottinelli's E8 basis."""

    e8 = [
        [2, 0, -1, 0, 0, 0, 0, 0],
        [0, 2, 0, -1, 0, 0, 0, 0],
        [-1, 0, 2, -1, 0, 0, 0, 0],
        [0, -1, -1, 2, -1, 0, 0, 0],
        [0, 0, 0, -1, 2, -1, 0, 0],
        [0, 0, 0, 0, -1, 2, -1, 0],
        [0, 0, 0, 0, 0, -1, 2, -1],
        [0, 0, 0, 0, 0, 0, -1, 2],
    ]
    gram = [[ZERO for _ in range(9)] for _ in range(9)]
    gram[0][0] = -2 * PHI
    for i in range(8):
        for j in range(8):
            gram[i + 1][j + 1] = Golden(e8[i][j])
    return gram


def bugaenko_diagram_basis_gram() -> list[list[Golden]]:
    """Gram matrix of the first nine roots in diagram order."""

    gram = [[Golden(2) if i == j else ZERO for j in range(9)] for i in range(9)]
    for i in range(8):
        product = -PHI if i in (0, 7) else Golden(-1)
        gram[i][i + 1] = product
        gram[i + 1][i] = product
    return gram


def bugaenko_primitive_roots() -> list[list[Golden]]:
    """The eleven primitive source roots ``r1,...,r11`` in ambient coordinates."""

    roots: list[list[Golden]] = []
    for i in range(1, 9):
        root = [ZERO] * 9
        root[i] = Golden(-1)
        roots.append(root)
    roots.extend(
        [
            [
                PHI - 1,
                Golden(2),
                Golden(3),
                Golden(4),
                Golden(6),
                Golden(5),
                Golden(4),
                Golden(3),
                Golden(2),
            ],
            [
                Golden(1), Golden(4), Golden(5), Golden(7), Golden(10),
                Golden(8), Golden(6), Golden(4), Golden(2),
            ],
            [
                1 + 2 * PHI,
                4 + 6 * PHI,
                7 + 9 * PHI,
                8 + 12 * PHI,
                12 + 18 * PHI,
                9 + 15 * PHI,
                6 + 12 * PHI,
                4 + 8 * PHI,
                2 + 4 * PHI,
            ],
        ]
    )
    return roots


def _scale(vector: Sequence[Golden], scalar: Golden) -> list[Golden]:
    return [scalar * entry for entry in vector]


def bugaenko_simple_roots() -> list[list[Golden]]:
    """Equal-norm roots in Guglielmetti's eleven-node diagram order."""

    roots = bugaenko_primitive_roots()
    # (g1,...,g11)=(r10,r1,r3,r4,r5,r6,r7,r8,r9,r2,r11).
    # The source roots r9,r10,r11 have norm 2*(2-phi), so multiplication
    # by the unit phi makes their norm equal to two.
    return [
        _scale(roots[9], PHI),
        roots[0],
        roots[2],
        roots[3],
        roots[4],
        roots[5],
        roots[6],
        roots[7],
        _scale(roots[8], PHI),
        roots[1],
        _scale(roots[10], PHI),
    ]


def bugaenko_interior_point() -> list[Golden]:
    """A strict timelike ambient vector; its first nine wall products are -1."""

    coefficients = [
        Golden(Fraction(7, 2)) + Golden(Fraction(9, 2)) * PHI,
        1 + 8 * PHI,
        Golden(Fraction(-3, 2)) + 8 * PHI,
        -3 + 8 * PHI,
        Golden(Fraction(-7, 2)) + 8 * PHI,
        -3 + 8 * PHI,
        Golden(Fraction(-3, 2)) + 8 * PHI,
        1 + 8 * PHI,
        Golden(Fraction(7, 2)) + Golden(Fraction(9, 2)) * PHI,
    ]
    top_roots = bugaenko_simple_roots()[:9]
    return [
        sum(
            (coefficients[j] * top_roots[j][i] for j in range(9)),
            start=ZERO,
        )
        for i in range(9)
    ]


def _mat_vec(matrix: Sequence[Sequence[Golden]], vector: Sequence[Golden]) -> list[Golden]:
    return [sum((row[j] * vector[j] for j in range(len(vector))), start=ZERO) for row in matrix]


def _pair(
    left: Sequence[Golden],
    gram: Sequence[Sequence[Golden]],
    right: Sequence[Golden],
) -> Golden:
    product = _mat_vec(gram, right)
    return sum((left[i] * product[i] for i in range(len(left))), start=ZERO)


def bugaenko_wall_gram() -> list[list[Golden]]:
    gram = bugaenko_basis_gram()
    roots = bugaenko_simple_roots()
    return [[_pair(left, gram, right) for right in roots] for left in roots]


def bugaenko_coxeter_relations() -> list[list[int | str]]:
    wall_gram = bugaenko_wall_gram()
    relations: list[list[int | str]] = []
    for i, row in enumerate(wall_gram):
        relation_row: list[int | str] = []
        for j, product in enumerate(row):
            if i == j:
                relation_row.append(1)
                continue
            relation = classify_root_pair(row[i], product, wall_gram[j][j])
            if relation.kind == "finite":
                assert relation.label is not None
                relation_row.append(relation.label)
            elif relation.kind == "parallel":
                relation_row.append("parallel")
            elif relation.kind == "ultraparallel":
                relation_row.append("ultraparallel")
            else:
                raise ValueError(f"unrecognized Coxeter relation at ({i},{j})")
        relations.append(relation_row)
    return relations


# Active facets of all projective extreme rays.  These are archived certificate
# data, not trusted by the verifier, which independently recomputes them.
BUGAENKO_VERTICES: tuple[tuple[int, ...], ...] = (
    (0, 1, 2, 3, 5, 6, 7, 8),
    (0, 1, 2, 3, 5, 6, 7, 10),
    (0, 1, 2, 3, 5, 6, 8, 10),
    (0, 1, 2, 3, 5, 7, 8, 10),
    (0, 1, 2, 3, 6, 7, 8, 10),
    (0, 1, 2, 4, 5, 6, 7, 9),
    (0, 1, 2, 4, 5, 6, 7, 10),
    (0, 1, 2, 4, 5, 6, 8, 9),
    (0, 1, 2, 4, 5, 6, 8, 10),
    (0, 1, 2, 4, 5, 7, 8, 9),
    (0, 1, 2, 4, 5, 7, 8, 10),
    (0, 1, 2, 4, 6, 7, 8, 9),
    (0, 1, 2, 4, 6, 7, 8, 10),
    (0, 1, 2, 5, 6, 7, 8, 9),
    (0, 1, 3, 4, 5, 6, 7, 9),
    (0, 1, 3, 4, 5, 6, 7, 10),
    (0, 1, 3, 4, 5, 6, 8, 9),
    (0, 1, 3, 4, 5, 6, 8, 10),
    (0, 1, 3, 4, 5, 7, 8, 9),
    (0, 1, 3, 4, 5, 7, 8, 10),
    (0, 1, 3, 4, 6, 7, 8, 9),
    (0, 1, 3, 4, 6, 7, 8, 10),
    (0, 1, 3, 5, 6, 7, 8, 9),
    (0, 2, 3, 4, 5, 6, 7, 9),
    (0, 2, 3, 4, 5, 6, 7, 10),
    (0, 2, 3, 4, 5, 6, 8, 9),
    (0, 2, 3, 4, 5, 6, 8, 10),
    (0, 2, 3, 4, 5, 7, 8, 9),
    (0, 2, 3, 4, 5, 7, 8, 10),
    (0, 2, 3, 4, 6, 7, 8, 9),
    (0, 2, 3, 4, 6, 7, 8, 10),
    (0, 2, 3, 5, 6, 7, 8, 9),
    (1, 2, 3, 4, 5, 6, 7, 9),
    (1, 2, 3, 4, 5, 6, 7, 10),
    (1, 2, 3, 4, 5, 6, 8, 9),
    (1, 2, 3, 4, 5, 6, 8, 10),
    (1, 2, 3, 4, 5, 7, 8, 9),
    (1, 2, 3, 4, 5, 7, 8, 10),
    (1, 2, 3, 4, 6, 7, 8, 9),
    (1, 2, 3, 4, 6, 7, 8, 10),
    (1, 2, 3, 5, 6, 7, 8, 9),
)


def bugaenko_f_vector() -> tuple[int, ...]:
    """Derive the face numbers from the archived vertex-facet incidences."""

    vertex_facets = [set(vertex) for vertex in BUGAENKO_VERTICES]
    faces: dict[frozenset[int], int] = {}
    for mask in range(1 << 11):
        chosen = {i for i in range(11) if mask & (1 << i)}
        incident_vertices = frozenset(
            index for index, facets in enumerate(vertex_facets) if chosen <= facets
        )
        if not incident_vertices:
            continue
        common_facets = set.intersection(
            *(vertex_facets[index] for index in incident_vertices)
        )
        faces[incident_vertices] = 8 - len(common_facets)
    return tuple(sum(dimension == k for dimension in faces.values()) for k in range(9))


def _encoded_vector(vector: Sequence[Golden]) -> list[dict[str, str]]:
    return [entry.to_json() for entry in vector]


def _encoded_matrix(matrix: Sequence[Sequence[Golden]]) -> list[list[dict[str, str]]]:
    return [_encoded_vector(row) for row in matrix]


def bugaenko_h8_certificate() -> dict[str, Any]:
    """Return the exact geometry certificate reconstructed from the diagram."""

    return {
        "schema_version": "0.1.0",
        "certificate_kind": "compact_reflective",
        "field": "Q(sqrt(5))",
        "dimension": 8,
        "lattice_gram": _encoded_matrix(bugaenko_basis_gram()),
        "simple_roots": _encoded_matrix(bugaenko_simple_roots()),
        "wall_gram": _encoded_matrix(bugaenko_wall_gram()),
        "coxeter_relations": bugaenko_coxeter_relations(),
        "completeness": {
            "method": "exact_vinberg_height_enumeration",
            "evidence": {
                "source": "Tumarkin2007",
                "locator": "p. 2, displayed Bugaenko H8 Coxeter diagram",
                "root_source": (
                    "VinbergsAlgorithmNF commit "
                    "ff32ab048ee483c7f3b9df3baad1ac23a22b6bb6, "
                    "test/some_known_lattices/Bugaenko.jl"
                ),
                "lattice": "<-2*phi> + E8 over Z[phi]",
                "coxiter_f_vector": list(bugaenko_f_vector()),
                "transcript": "data/known_examples/bugaenko_vinberg_transcript.json",
                "transcript_verifier": (
                    "coxeter_polytopes.verify.vinberg."
                    "verify_bugaenko_vinberg_completion"
                ),
                "control_vector": "e0 of norm -2*phi",
                "positive_height_solution_counts": [1, 1, 0, 0, 0, 0, 0, 1],
                "status": "verified",
            },
        },
        "geometry": {
            "nonempty": True,
            "interior_point": _encoded_vector(bugaenko_interior_point()),
            "vertices": [list(vertex) for vertex in BUGAENKO_VERTICES],
            "cusps": [],
        },
    }
