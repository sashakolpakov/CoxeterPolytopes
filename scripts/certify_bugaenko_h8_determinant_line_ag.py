#!/usr/bin/env python3
"""Compile Bugaenko's compact H8 chamber into a determinant-line AG carrier.

The Coxeter lattice, roots, chamber vertices, and Vinberg completeness flag are
read directly from ``data/known_examples/bugaenko_compact_h8.json``.  No
hard-coded Bugaenko lattice module is imported.  The compact norm-two mirror
enumeration is also replayed from those data and compared with the independent
archived K3 mirror checksum.

This is an arrangement certificate.  It does not assert that the resulting
Heegner divisor is principal or is the divisor of a Borcherds product.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any, Iterable, Sequence

from coxeter_polytopes.exact.golden import Golden, PHI, ZERO
from coxeter_polytopes.exact.matrix import determinant, inertia, rectangular_pairing


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "known_examples" / "bugaenko_compact_h8.json"
MIRROR_ARCHIVE = (
    ROOT / "data" / "trace_lattices" / "bugaenko_k3_eigen_nef.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "bugaenko_h8_determinant_line_ag.json"
)

Vector = tuple[Golden, ...]
Matrix = list[list[Golden]]
Residue = tuple[int, int]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _golden_matrix(payload: Sequence[Sequence[dict[str, str]]]) -> Matrix:
    return [[Golden.from_json(entry) for entry in row] for row in payload]


def _golden_vectors(payload: Sequence[Sequence[dict[str, str]]]) -> list[Vector]:
    return [tuple(Golden.from_json(entry) for entry in row) for row in payload]


def _inertia_record(matrix: Sequence[Sequence[Golden]], embedding: int) -> dict[str, int]:
    value = inertia(matrix, embedding=embedding)
    return {
        "positive": value.positive,
        "negative": value.negative,
        "zero": value.zero,
    }


def _pair(
    left: Sequence[Golden],
    gram: Sequence[Sequence[Golden]],
    right: Sequence[Golden],
) -> Golden:
    return sum(
        (
            left[i] * gram[i][j] * right[j]
            for i in range(len(left))
            for j in range(len(right))
        ),
        start=ZERO,
    )


def _identity(size: int) -> Matrix:
    return [
        [Golden(int(left == right)) for right in range(size)]
        for left in range(size)
    ]


def _transpose(matrix: Sequence[Sequence[Golden]]) -> Matrix:
    return [list(column) for column in zip(*matrix, strict=True)]


def _matrix_product(
    left: Sequence[Sequence[Golden]], right: Sequence[Sequence[Golden]]
) -> Matrix:
    return [
        [
            sum(
                (left[i][k] * right[k][j] for k in range(len(right))),
                start=ZERO,
            )
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def _block_diagonal(
    left: Sequence[Sequence[Golden]], right_entry: Golden
) -> Matrix:
    size = len(left)
    return [
        [
            left[i][j]
            if i < size and j < size
            else right_entry
            if i == size and j == size
            else ZERO
            for j in range(size + 1)
        ]
        for i in range(size + 1)
    ]


def _is_unit(value: Golden) -> bool:
    return value.is_integral and abs(value.norm()) == 1


def _reflection_matrix(root: Sequence[Golden], gram: Matrix) -> Matrix:
    """Matrix of x -> x + (x,r)r for a square-minus-two root."""

    size = len(root)
    covector = [
        sum((root[k] * gram[k][j] for k in range(size)), start=ZERO)
        for j in range(size)
    ]
    return [
        [Golden(int(i == j)) + root[i] * covector[j] for j in range(size)]
        for i in range(size)
    ]


def _canonical_mirror(root: Sequence[Golden]) -> Vector:
    vector = tuple(root)
    first = next((entry for entry in vector if entry), None)
    if first is None:
        raise ArithmeticError("zero vector does not define a mirror")
    return vector if first.sign(0) > 0 else tuple(-entry for entry in vector)


def _vector_key(vector: Sequence[Golden]) -> tuple[tuple[int, int, int, int], ...]:
    return tuple(
        (entry.a.numerator, entry.a.denominator, entry.b.numerator, entry.b.denominator)
        for entry in vector
    )


def _encoded_vector(vector: Sequence[Golden]) -> list[dict[str, str]]:
    return [entry.to_json() for entry in vector]


def _mirror_digest(mirrors: Iterable[Vector]) -> str:
    payload = [
        _encoded_vector(vector) for vector in sorted(mirrors, key=_vector_key)
    ]
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _spherical_mirrors_at_vertex(
    active_facets: tuple[int, ...],
    roots: Sequence[Vector],
    positive_gram: Matrix,
) -> tuple[Vector, ...]:
    generators = tuple(roots[index] for index in active_facets)
    local_gram = [
        [_pair(left, positive_gram, right) for right in generators]
        for left in generators
    ]
    local_inertia = inertia(local_gram)
    if (local_inertia.positive, local_inertia.negative, local_inertia.zero) != (
        len(generators),
        0,
        0,
    ):
        raise ArithmeticError(f"vertex {active_facets} is not spherical")

    rank = len(generators)
    coefficient_roots = {
        tuple(Golden(int(i == j)) for i in range(rank)) for j in range(rank)
    }
    frontier = list(coefficient_roots)
    while frontier:
        coefficients = frontier.pop()
        for wall_index in range(rank):
            product = sum(
                (
                    coefficients[i] * local_gram[i][wall_index]
                    for i in range(rank)
                ),
                start=ZERO,
            )
            reflected = list(coefficients)
            reflected[wall_index] -= product
            candidate = tuple(reflected)
            if candidate not in coefficient_roots:
                coefficient_roots.add(candidate)
                frontier.append(candidate)

    if any(_pair(root, local_gram, root) != 2 for root in coefficient_roots):
        raise ArithmeticError("finite parabolic closure produced a non-root")
    ambient_roots = {
        tuple(
            sum(
                (coefficients[i] * generators[i][j] for i in range(rank)),
                start=ZERO,
            )
            for j in range(len(generators[0]))
        )
        for coefficients in coefficient_roots
    }
    if len(ambient_roots) != len(coefficient_roots):
        raise ArithmeticError("vertex roots are linearly dependent")
    mirrors = {_canonical_mirror(root) for root in ambient_roots}
    if len(ambient_roots) != 2 * len(mirrors):
        raise ArithmeticError("finite root system is not closed under sign")
    return tuple(sorted(mirrors, key=_vector_key))


def _closed_chamber_mirrors(
    vertices: Sequence[Sequence[int]],
    roots: Sequence[Vector],
    positive_gram: Matrix,
    dimension: int,
) -> tuple[tuple[Vector, ...], Counter[int], list[int]]:
    local_sets: list[set[Vector]] = []
    local_counts: list[int] = []
    mirrors: set[Vector] = set()
    for vertex in vertices:
        local = set(
            _spherical_mirrors_at_vertex(tuple(vertex), roots, positive_gram)
        )
        local_sets.append(local)
        local_counts.append(len(local))
        mirrors.update(local)

    ordered = tuple(sorted(mirrors, key=_vector_key))
    vertex_facets = [set(vertex) for vertex in vertices]
    contact_dimensions: Counter[int] = Counter()
    for mirror in ordered:
        incident = [index for index, local in enumerate(local_sets) if mirror in local]
        if not incident:
            raise ArithmeticError("enumerated mirror has no incident vertex")
        common_facets = set.intersection(
            *(vertex_facets[index] for index in incident)
        )
        contact_dimensions[dimension - len(common_facets)] += 1
    return ordered, contact_dimensions, local_counts


def _residue_mul(left: Residue, right: Residue) -> Residue:
    """Multiply pairs a+b*phi using phi^2=phi+1, before reduction."""

    a, b = left
    c, d = right
    return a * c + b * d, a * d + b * c + b * d


def _discriminant_enumeration() -> dict[str, Any]:
    """Enumerate q=(2*a^2+b^2)/(8*phi) on O/(2) + O/(4)."""

    records: list[dict[str, Any]] = []
    isotropic: list[dict[str, list[int]]] = []
    histogram: Counter[tuple[int, int]] = Counter()
    for a in itertools.product(range(2), repeat=2):
        for b in itertools.product(range(4), repeat=2):
            a_square = _residue_mul(a, a)
            b_square = _residue_mul(b, b)
            numerator = (
                2 * a_square[0] + b_square[0],
                2 * a_square[1] + b_square[1],
            )
            residue = numerator[0] % 8, numerator[1] % 8
            record = {
                "a_mod_2": list(a),
                "b_mod_4": list(b),
                "q_numerator_mod_8": list(residue),
            }
            records.append(record)
            histogram[residue] += 1
            if residue == (0, 0):
                isotropic.append(
                    {"a_mod_2": list(a), "b_mod_4": list(b)}
                )

    encoded = json.dumps(records, sort_keys=True, separators=(",", ":")).encode()
    if len(records) != 64 or isotropic != [
        {"a_mod_2": [0, 0], "b_mod_4": [0, 0]}
    ]:
        raise ArithmeticError("unexpected isotropic discriminant elements")
    return {
        "element_count": len(records),
        "underlying_abelian_group": "(Z/2)^2 + (Z/4)^2",
        "quadratic_module": {
            "module": "O/(2) * u + O/(4) * v",
            "q_u": "1/(4*phi) mod O",
            "q_v": "1/(8*phi) mod O",
            "formula": "q(a*u+b*v)=(2*a^2+b^2)/(8*phi) mod O",
        },
        "enumeration_sha256": hashlib.sha256(encoded).hexdigest(),
        "q_numerator_mod_8_histogram": {
            f"{left},{right}": count
            for (left, right), count in sorted(histogram.items())
        },
        "isotropic_elements": isotropic,
        "isotropic_element_count": len(isotropic),
        "totally_isotropic_O_submodules": [
            {"generators": [], "order": 1, "overlattice_index": 1}
        ],
        "totally_isotropic_O_submodule_count": 1,
        "reflection_stable_even_overlattices": ["Lhat_0"],
        "reflection_stable_even_overlattice_count": 1,
    }


def analyze() -> dict[str, Any]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    mirror_archive = json.loads(MIRROR_ARCHIVE.read_text(encoding="utf-8"))
    dimension = int(source["dimension"])
    if dimension != 8:
        raise ArithmeticError("this calibration expects Bugaenko H8")
    if source["completeness"]["evidence"]["status"] != "verified":
        raise ArithmeticError("the source chamber is not marked Vinberg-complete")

    source_gram = _golden_matrix(source["lattice_gram"])
    lattice_gram = [[-entry for entry in row] for row in source_gram]
    roots = _golden_vectors(source["simple_roots"])
    source_wall_gram = _golden_matrix(source["wall_gram"])
    normalized_wall_gram = [[-entry for entry in row] for row in source_wall_gram]
    root_columns = [
        [roots[column][row] for column in range(len(roots))]
        for row in range(len(lattice_gram))
    ]
    computed_wall_gram = rectangular_pairing(
        root_columns, lattice_gram, root_columns
    )
    if computed_wall_gram != normalized_wall_gram:
        raise ArithmeticError("normalized root Gram does not match source wall Gram")

    lattice_determinant = determinant(lattice_gram)
    if lattice_determinant != 2 * PHI:
        raise ArithmeticError(f"unexpected normalized determinant {lattice_determinant}")
    e8_block = [row[1:] for row in lattice_gram[1:]]
    if determinant(e8_block) != 1:
        raise ArithmeticError("the E8(-1) block is not unimodular")
    determinant_line_square = 2 * ((-1) ** dimension) * lattice_determinant
    if determinant_line_square != 4 * PHI:
        raise ArithmeticError("determinant-line square is not 4*phi")
    extended_gram = _block_diagonal(lattice_gram, determinant_line_square)
    extended_determinant = determinant(extended_gram)
    if extended_determinant != 8 * PHI * PHI:
        raise ArithmeticError("unexpected determinant-line extension determinant")

    lattice_signatures = {
        "identity": _inertia_record(lattice_gram, 0),
        "conjugate": _inertia_record(lattice_gram, 1),
    }
    extended_signatures = {
        "identity": _inertia_record(extended_gram, 0),
        "conjugate": _inertia_record(extended_gram, 1),
    }
    if lattice_signatures != {
        "identity": {"positive": 1, "negative": 8, "zero": 0},
        "conjugate": {"positive": 0, "negative": 9, "zero": 0},
    }:
        raise ArithmeticError("normalized lattice has the wrong signatures")
    if extended_signatures != {
        "identity": {"positive": 2, "negative": 8, "zero": 0},
        "conjugate": {"positive": 0, "negative": 10, "zero": 0},
    }:
        raise ArithmeticError("determinant-line extension has the wrong signatures")
    if not all(
        entry.is_integral for row in extended_gram for entry in row
    ) or not all((extended_gram[index][index] / 2).is_integral for index in range(10)):
        raise ArithmeticError("determinant-line extension is not even integral")

    root_records: list[dict[str, Any]] = []
    reflection_digests: list[str] = []
    for index, root in enumerate(roots):
        square = _pair(root, lattice_gram, root)
        if square != -2:
            raise ArithmeticError(f"root {index} has square {square}, not -2")
        coordinate_witness = next(
            ((position, entry) for position, entry in enumerate(root) if _is_unit(entry)),
            None,
        )
        pairings = [
            sum(
                (lattice_gram[row][column] * root[column] for column in range(9)),
                start=ZERO,
            )
            for row in range(9)
        ]
        divisibility_witness = next(
            ((position, entry) for position, entry in enumerate(pairings) if _is_unit(entry)),
            None,
        )
        if coordinate_witness is None or divisibility_witness is None:
            raise ArithmeticError(f"root {index} lacks a unit primitivity/divisibility witness")
        reflection = _reflection_matrix(root, lattice_gram)
        if any(not entry.is_integral for row in reflection for entry in row):
            raise ArithmeticError(f"reflection {index} is not integral")
        if _matrix_product(reflection, reflection) != _identity(9):
            raise ArithmeticError(f"reflection {index} is not an involution")
        if _matrix_product(
            _matrix_product(_transpose(reflection), lattice_gram), reflection
        ) != lattice_gram:
            raise ArithmeticError(f"reflection {index} is not an isometry")
        if determinant(reflection) != -1:
            raise ArithmeticError(f"reflection {index} has wrong determinant")
        encoded_reflection = json.dumps(
            [[entry.to_json() for entry in row] for row in reflection],
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
        reflection_digests.append(hashlib.sha256(encoded_reflection).hexdigest())
        root_records.append(
            {
                "facet": index,
                "square_in_Lhat": "-2",
                "primitive": True,
                "primitive_unit_coordinate": {
                    "basis_index": coordinate_witness[0],
                    "value": str(coordinate_witness[1]),
                },
                "divisibility_ideal": "O",
                "divisibility_unit_pairing": {
                    "basis_index": divisibility_witness[0],
                    "value": str(divisibility_witness[1]),
                },
                "heegner_type_signature_2_8": {
                    "norm": "-2",
                    "discriminant_coset": "0",
                },
                "bruinier_sign_type": {"m": 1, "mu": "0"},
            }
        )

    discriminant = _discriminant_enumeration()

    relation_matrix = source["coxeter_relations"]
    vertices = source["geometry"]["vertices"]
    vertex_pairs = {
        pair
        for vertex in vertices
        for pair in itertools.combinations(sorted(vertex), 2)
    }
    finite_pairs: set[tuple[int, int]] = set()
    pair_records: list[dict[str, Any]] = []
    label_counts: Counter[int] = Counter()
    expected_pairings = {2: ZERO, 3: Golden(1), 5: PHI}
    ultraparallel_pairs: list[list[int]] = []
    for left in range(len(roots)):
        for right in range(left + 1, len(roots)):
            label = relation_matrix[left][right]
            pairing = computed_wall_gram[left][right]
            if label == "ultraparallel":
                if pairing != 2 * PHI:
                    raise ArithmeticError("wrong ultraparallel Gram entry")
                ultraparallel_pairs.append([left, right])
                pair_records.append(
                    {
                        "facets": [left, right],
                        "relation": "ultraparallel",
                        "pairing": str(pairing),
                        "is_ridge": False,
                        "real_slice_codimension_two_meeting": False,
                    }
                )
                continue
            numerical_label = int(label)
            if pairing != expected_pairings[numerical_label]:
                raise ArithmeticError(f"wrong Gram entry for I2({label})")
            finite_pairs.add((left, right))
            label_counts[numerical_label] += 1
            pair_records.append(
                {
                    "facets": [left, right],
                    "relation": f"I2({numerical_label})",
                    "pairing": str(pairing),
                    "is_ridge": (left, right) in vertex_pairs,
                    "distinguished_stack_inertia_order": 2 * numerical_label,
                }
            )
    ridge_count = source["completeness"]["evidence"]["coxiter_f_vector"][
        dimension - 2
    ]
    if finite_pairs != vertex_pairs or len(finite_pairs) != ridge_count:
        raise ArithmeticError("finite wall pairs do not equal the archived ridges")
    if label_counts != Counter({2: 44, 3: 8, 5: 2}):
        raise ArithmeticError("unexpected Bugaenko ridge-label profile")
    if ultraparallel_pairs != [[9, 10]]:
        raise ArithmeticError("unexpected ultraparallel pair")

    mirrors, contact_dimensions, local_counts = _closed_chamber_mirrors(
        vertices, roots, source_gram, dimension
    )
    facet_mirrors = {_canonical_mirror(root) for root in roots}
    mirror_set = set(mirrors)
    if len(facet_mirrors) != 11 or not facet_mirrors <= mirror_set:
        raise ArithmeticError("facet mirrors are missing from the compact mirror set")
    if any(_pair(root, source_gram, root) != 2 for root in mirrors):
        raise ArithmeticError("compact mirror enumeration contains a non-root")
    if any(
        not all(entry.is_integral for entry in root)
        for root in mirrors
    ):
        raise ArithmeticError("compact mirror enumeration contains a nonintegral root")

    mirror_digest = _mirror_digest(mirrors)
    archived_closed = mirror_archive["closed_chamber"]
    computed_closed = {
        "distinct_k3_mirrors_meeting": len(mirrors),
        "mirror_representative_sha256": mirror_digest,
        "facet_mirrors": len(facet_mirrors),
        "nonfacet_mirrors_touching_lower_faces": len(mirrors) - len(facet_mirrors),
        "mirrors_cutting_interior": 0,
        "contact_dimension_distribution": {
            str(key): value for key, value in sorted(contact_dimensions.items())
        },
    }
    if computed_closed != archived_closed:
        raise ArithmeticError("source-only compact mirror replay disagrees with archive")
    if mirror_archive["compact_vertex_count"] != len(vertices):
        raise ArithmeticError("mirror archive has a different compact vertex set")
    if not mirror_archive["vinberg_completion_verified"]:
        raise ArithmeticError("mirror archive lacks Vinberg completion")

    nonfacet_contact_dimensions = Counter(contact_dimensions)
    nonfacet_contact_dimensions[dimension - 1] -= len(facet_mirrors)
    if nonfacet_contact_dimensions[dimension - 1] != 0:
        raise ArithmeticError("a nonfacet mirror contains a codimension-one face")
    del nonfacet_contact_dimensions[dimension - 1]

    payload: dict[str, Any] = {
        "schema_version": "0.1.0",
        "certificate_kind": "bugaenko_h8_canonical_determinant_line_AG_arrangement",
        "source": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": _sha256(SOURCE),
            "only_polytope_lattice_input": True,
            "dimension": dimension,
            "facet_count": len(roots),
            "compact_vertex_count": len(vertices),
            "cusp_count": len(source["geometry"]["cusps"]),
            "vinberg_completeness_status": source["completeness"]["evidence"]["status"],
            "f_vector": source["completeness"]["evidence"]["coxiter_f_vector"],
        },
        "normalization": {
            "operation": "multiply the archived bilinear form by -1",
            "lattice": "L=<2*phi> + E8(-1) over O=Z[phi]",
            "rank": len(lattice_gram),
            "even_integral": True,
            "determinant": str(lattice_determinant),
            "signatures": lattice_signatures,
        },
        "determinant_line": {
            "coordinate_free_module": "K_L=det_O(L)",
            "bilinear_form": "b_K(wedge v_i,wedge w_i)=2*(-1)^8*det(b_L(v_i,w_j))",
            "computational_generator": "e=b_0 wedge ... wedge b_8",
            "generator_square": str(determinant_line_square),
            "basis_independence": "under A in GL(L), e scales by det(A) and b_K by det(A)^2",
            "sign_identity_embedding": determinant_line_square.sign(0),
            "sign_conjugate_embedding": determinant_line_square.sign(1),
        },
        "canonical_extension": {
            "lattice": "Lhat_0=L orthogonal_sum K_L=<2*phi>+E8(-1)+<4*phi>",
            "rank": len(extended_gram),
            "even_integral": True,
            "signatures": extended_signatures,
            "determinant": str(extended_determinant),
            "determinant_ideal": "(8*phi^2)=(8)",
            "absolute_discriminant_module_order": 64,
            "anisotropic_over_Q_sqrt_5": True,
            "anisotropy_proof": "an F-isotropic vector would remain isotropic at the negative-definite conjugate real place",
        },
        "discriminant_gluing": discriminant,
        "root_extensions": {
            "extension_rule": "r_i maps to (r_i,0); every wall reflection fixes K_L pointwise",
            "wall_gram_matches_normalized_source": True,
            "root_count": len(root_records),
            "all_primitive_square_minus_two_divisibility_one": True,
            "all_reflections_integral_involutive_isometries": True,
            "reflection_matrix_sha256": reflection_digests,
            "discriminant_action": "identity",
            "discriminant_action_proof": "for x in L^vee, s_r(x)-x=(x,r)r lies in L; K_L is fixed",
            "roots": root_records,
        },
        "canonical_real_slice": {
            "positive_line_generator": "e in K_L with e^2=4*phi",
            "map": "j_e([x])=span^+(e,x)",
            "source": "H(L_identity)=H^8",
            "target": "D^+(Lhat_0), type IV of complex dimension 8",
            "facet_mirror_pullback": "j_e^{-1}(D_(r_i,0))=r_i^perp",
            "arithmetic_stack": "[D^+(Lhat_0)/Gamma], for an arithmetic Gamma containing the extended reflections",
            "proper_projective_coarse_quotient": True,
            "surface_extraction": "on a neat projective cover, intersect six sufficiently general ample divisors; restricted Heegner divisors become curves",
            "surface_is_canonical_without_auxiliary_choices": False,
        },
        "facet_heegner_types": {
            "facet_count": len(roots),
            "all_same_numerical_type": True,
            "signature_2_8_type": {"norm": "-2", "divisibility": 1, "coset": "0"},
            "bruinier_sign_type": {"m": 1, "mu": "0"},
            "global_orbit_or_irreducible_component_equality_claimed": False,
        },
        "stacky_coxeter_labels": {
            "ridge_count": len(finite_pairs),
            "finite_pair_equals_ridge_set": True,
            "distinguished_dihedral_inertia_counts": {
                f"I2({label})_order_{2 * label}": count
                for label, count in sorted(label_counts.items())
            },
            "generic_inertia_contains_not_necessarily_equals_distinguished_group": True,
            "ultraparallel_pairs": ultraparallel_pairs,
            "pairs": pair_records,
        },
        "hidden_slice_gate": {
            "shell": {"norm_in_L": "-2", "bruinier_type": {"m": 1, "mu": "0"}},
            "method": "exact finite-parabolic reflection closure at all 41 compact vertices",
            "completeness_proof": "a root mirror cannot cut the interior of a complete reflection chamber; if it meets the compact closed chamber, its contact face has a vertex and it lies in that vertex's finite standard parabolic root system",
            "classification": {
                "facet_mirrors": len(facet_mirrors),
                "nonfacet_lower_boundary_contacts": len(mirrors) - len(facet_mirrors),
                "interior_cutting_mirrors": 0,
            },
            "nonfacet_contact_dimension_distribution": {
                str(key): value
                for key, value in sorted(nonfacet_contact_dimensions.items())
            },
            "all_contact_dimension_distribution": {
                str(key): value for key, value in sorted(contact_dimensions.items())
            },
            "distinct_mirror_count": len(mirrors),
            "mirror_representative_sha256": mirror_digest,
            "local_mirror_count_distribution": {
                str(key): value
                for key, value in sorted(Counter(local_counts).items())
            },
            "sum_local_mirror_counts": sum(local_counts),
            "closed_slice_is_facet_only": False,
            "open_chamber_interior_is_uncut_by_extra_shell_mirrors": True,
            "interpretation": "the 309 extra mirrors touch only lower-dimensional boundary faces; they obstruct a facet-only closed-slice divisor but do not refine the chamber interior",
        },
        "mirror_archive_cross_check": {
            "path": str(MIRROR_ARCHIVE.relative_to(ROOT)),
            "sha256": _sha256(MIRROR_ARCHIVE),
            "source_only_replay_matches_archive_exactly": True,
        },
        "AG_arrangement_conclusion": {
            "new_realization": "canonical determinant-line type-IV arrangement carrier for the known compact Bugaenko H8 chamber",
            "associated_projective_AG_object_exists": True,
            "associated_surface_exists_after_neat_cover_and_complete_intersection": True,
            "Coxeter_angles_recovered_stack_theoretically": True,
            "full_numerical_Heegner_type_is_facet_only_on_closed_P": False,
            "full_numerical_Heegner_type_cuts_open_P": False,
            "automorphic_principality_or_Borcherds_product_claimed": False,
            "new_abstract_Coxeter_polytope_claimed": False,
        },
        "all_checks_pass": True,
    }
    return payload


def main() -> None:
    payload = analyze()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
