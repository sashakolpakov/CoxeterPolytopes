#!/usr/bin/env python3
"""Certify an explicit good K3-cover lift package for Enriques Type I.

This reconstruction joins three independently pinned objects:

* Shimada's Type-I ``12A`` roots in ``L10``;
* OSCAR 1.4.0's serialized Picard lattice for
  ``generic_enriques_surface(172)``; and
* the finite ``L1``--``L4`` gates of the variable Enriques lift theorem.

The result is an explicit positive calibration, not a new abstract Coxeter
polytope or a new Enriques surface.  In particular, the script recovers all
twenty-four K3 curve lifts, the global Enriques involution on the rank-19
Picard lattice, its saturated rootless anti-invariant lattice, a small ample
orientation, and the complete bounded-root irreducibility certificate.

PARI/GP's exact ``qfminim`` is used twice as a complete positive-definite
short-vector enumerator.  Every returned vector is rechecked with integer
matrix arithmetic before it contributes to the certificate.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from dataclasses import asdict
import hashlib
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any, Iterable

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

try:
    from scripts.certify_enriques_variable_lift_package import parse_gap_records
except ModuleNotFoundError:  # direct ``python scripts/...py`` execution
    from certify_enriques_variable_lift_package import parse_gap_records
from coxeter_polytopes.exact.matrix import inertia


ROOT = Path(__file__).resolve().parents[1]
OSCAR_DATA = (
    ROOT
    / "data"
    / "external"
    / "oscar_enriques_type_i"
    / "TauTaubarGenericEnriquesSurfaceNo172.mrdi"
)
SHIMADA_DATA = (
    ROOT
    / "data"
    / "external"
    / "shimada_l10l26note"
    / "L10L26compdata2.txt"
)
FULL_NEF_ARCHIVE = (
    ROOT
    / "data"
    / "known_examples"
    / "enriques_finite_automorphism_nef_chambers.json"
)
VARIABLE_THEOREM_ARCHIVE = (
    ROOT / "data" / "enumerations" / "enriques_variable_lift_package.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "enriques_type_i_explicit_lift_package.json"
)

OSCAR_UPSTREAM_SHA256 = (
    "73e8e8612a9223156fb6024c62d51b32dc8a739ec9cd1cb09129a65acb73caa1"
)
OSCAR_ARCHIVE_SHA256 = (
    "14c9bddd93dbbcbefcda855302ee2f4fccc2f9acfc6e29b1d385e041b36f5fa8"
)
SHIMADA_SHA256 = (
    "218c8188e467bde29369cdeda3771d767c3f70105d8fb56def8d698b90e75a6c"
)
FULL_NEF_ARCHIVE_SHA256 = (
    "245715b7352219a7ed1bfc34a43b5b37921e289d39829df51015c347b97072bc"
)
VARIABLE_THEOREM_ARCHIVE_SHA256 = (
    "27115e714ef4acd5a8bb7801a12795668f9c9b1f0309816f336454c591583a6a"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _matrix(raw: Iterable[Iterable[Any]]) -> sp.Matrix:
    return sp.Matrix([[int(value) for value in row] for row in raw])


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=sp.ZZ)
    return [
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if smith[index, index] != 0
    ]


def _vector_hash(vectors: Iterable[Iterable[int]]) -> str:
    canonical = sorted(tuple(int(value) for value in vector) for vector in vectors)
    encoded = json.dumps(canonical, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _gp_matrix(matrix: sp.Matrix) -> str:
    if any(value.q != 1 for value in matrix):
        raise ArithmeticError("PARI form is not integral")
    return "[" + ";".join(
        ",".join(str(int(matrix[row, column])) for column in range(matrix.cols))
        for row in range(matrix.rows)
    ) + "]"


def _pari_short_vectors(
    form: sp.Matrix,
    bound: int,
) -> tuple[list[tuple[int, ...]], dict[str, Any]]:
    """Enumerate every nonzero integral vector of ``form`` up to ``bound``.

    PARI returns one vector from each sign pair and separately reports the
    signed total.  We expand both signs and verify every norm exactly.
    """

    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP is required for the exact root enumeration")
    if form.rows != form.cols:
        raise ArithmeticError("short-vector form is not square")
    form_inertia = inertia(_rows(form))
    if asdict(form_inertia) != {
        "positive": form.rows,
        "negative": 0,
        "zero": 0,
    }:
        raise ArithmeticError("short-vector form is not positive definite")

    program = (
        f"A={_gp_matrix(form)};"
        f"R=qfminim(A,{int(bound)});"
        "print(R[1]);print(R[2]);print(matsize(R[3]));"
        "V=R[3];for(j=1,matsize(V)[2],print(Vec(V[,j])));\n"
    )
    completed = subprocess.run(
        [gp, "-fq"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if len(lines) < 3:
        raise ArithmeticError(f"truncated qfminim output: {completed.stdout!r}")
    signed_total = int(lines[0])
    maximum_norm = int(lines[1])
    dimensions = ast.literal_eval(lines[2])
    representatives = [tuple(int(value) for value in ast.literal_eval(line)) for line in lines[3:]]
    if dimensions != [form.rows, len(representatives)]:
        raise ArithmeticError("qfminim matrix dimensions changed")

    vectors = representatives + [tuple(-value for value in vector) for vector in representatives]
    if signed_total != len(vectors) or len(set(vectors)) != len(vectors):
        raise ArithmeticError("qfminim sign expansion is inconsistent")
    exact_norms = []
    for vector in vectors:
        column = sp.Matrix(vector)
        norm = int((column.T * form * column)[0])
        if not 0 < norm <= bound:
            raise ArithmeticError("qfminim returned a vector outside its bound")
        exact_norms.append(norm)
    if maximum_norm != max(exact_norms, default=0):
        raise ArithmeticError("qfminim maximum norm disagrees with exact replay")

    vectors = sorted(vectors)
    return vectors, {
        "engine": "PARI/GP qfminim exact positive-definite enumeration",
        "bound": int(bound),
        "signed_vector_count": signed_total,
        "sign_representative_count": len(representatives),
        "maximum_norm": maximum_norm,
        "vector_set_sha256": _vector_hash(vectors),
        "all_vectors_rechecked_by_exact_integer_arithmetic": True,
    }


def _find_ambient_gram(document: dict[str, Any]) -> sp.Matrix:
    candidates = []
    for reference in document["_refs"].values():
        data = reference.get("data")
        if (
            isinstance(data, list)
            and len(data) == 26
            and all(isinstance(row, list) and len(row) == 26 for row in data)
        ):
            candidates.append(_matrix(data))
    if len(candidates) != 1:
        raise ArithmeticError(f"expected one ambient Gram, found {len(candidates)}")
    return candidates[0]


def _is_signed_e8_plus_a1_cartan(cartan: sp.Matrix) -> bool:
    """Recognize a sign-switch of the ``E8+A1`` Cartan matrix.

    The non-isolated component must be a tree whose unique trivalent vertex
    has arms of lengths 1, 2, and 4.  On a tree every choice of edge signs is
    removed by diagonal sign switches, so this recognizes the integral
    isometry type rather than only the unsigned graph.
    """

    if cartan.shape != (9, 9) or cartan != cartan.T:
        return False
    if any(cartan[index, index] != 2 for index in range(9)):
        return False
    if any(
        cartan[left, right] not in (-1, 0, 1)
        for left in range(9)
        for right in range(9)
        if left != right
    ):
        return False
    adjacency = {
        vertex: [
            other
            for other in range(9)
            if other != vertex and cartan[vertex, other] != 0
        ]
        for vertex in range(9)
    }
    isolated = [vertex for vertex, neighbors in adjacency.items() if not neighbors]
    branches = [vertex for vertex, neighbors in adjacency.items() if len(neighbors) == 3]
    if len(isolated) != 1 or len(branches) != 1:
        return False
    if sum(len(neighbors) for neighbors in adjacency.values()) // 2 != 7:
        return False
    branch = branches[0]
    arm_lengths = []
    visited = {branch, isolated[0]}
    for neighbor in adjacency[branch]:
        previous = branch
        current = neighbor
        length = 1
        visited.add(current)
        while len(adjacency[current]) == 2:
            following = next(vertex for vertex in adjacency[current] if vertex != previous)
            previous, current = current, following
            length += 1
            visited.add(current)
        if len(adjacency[current]) != 1:
            return False
        arm_lengths.append(length)
    return sorted(arm_lengths) == [1, 2, 4] and len(visited) == 9


def _type_i_archive_record() -> dict[str, Any]:
    document = json.loads(FULL_NEF_ARCHIVE.read_text(encoding="utf-8"))
    records = [record for record in document["types"] if record["nikulin_kondo_type"] == "I"]
    if len(records) != 1:
        raise ArithmeticError("full-nef archive has no unique Type-I record")
    record = records[0]
    expected = {
        "facet_count": 12,
        "classification": "finite_volume_noncompact_H9_Coxeter_polytope",
    }
    if any(record[key] != value for key, value in expected.items()):
        raise ArithmeticError("Type-I downstream archive summary changed")
    if record["exact_cone"]["hyperideal_vertex_count"] != 0:
        raise ArithmeticError("Type-I cone has a hyperideal extreme ray")
    if not record["literal_vinberg"]["finite_volume"] or record["literal_vinberg"]["compact"]:
        raise ArithmeticError("Type-I literal Vinberg decision changed")
    return record


def analyze() -> dict[str, Any]:
    raw_oscar = OSCAR_DATA.read_bytes()
    if _sha256(OSCAR_DATA) != OSCAR_ARCHIVE_SHA256:
        raise AssertionError("archived OSCAR datum changed")
    if not raw_oscar.endswith(b"\n") or raw_oscar.endswith(b"\n\n"):
        raise AssertionError("unexpected OSCAR archive newline normalization")
    if hashlib.sha256(raw_oscar[:-1]).hexdigest() != OSCAR_UPSTREAM_SHA256:
        raise AssertionError("OSCAR datum no longer matches the upstream byte stream")
    if _sha256(SHIMADA_DATA) != SHIMADA_SHA256:
        raise AssertionError("Shimada datum changed")
    if _sha256(FULL_NEF_ARCHIVE) != FULL_NEF_ARCHIVE_SHA256:
        raise AssertionError("full-nef transfer archive changed")
    if _sha256(VARIABLE_THEOREM_ARCHIVE) != VARIABLE_THEOREM_ARCHIVE_SHA256:
        raise AssertionError("variable lift theorem archive changed")

    oscar = json.loads(raw_oscar)
    if oscar.get("_ns", {}).get("Oscar") != [
        "https://github.com/oscar-system/Oscar.jl",
        "1.4.0",
    ]:
        raise ArithmeticError("unexpected OSCAR serialization namespace")
    if len(oscar["data"]) != 5:
        raise ArithmeticError("OSCAR Type-I tuple no longer has five entries")

    sy = _matrix(oscar["data"][0])
    sx = _matrix(oscar["data"][1])
    l26_basis = _matrix(oscar["data"][2])
    weyl = _matrix(oscar["data"][3])
    preprocessing_vector = _matrix(oscar["data"][4])
    gram_l26 = _find_ambient_gram(oscar)
    if sy.shape != (10, 26) or sx.shape != (19, 26):
        raise ArithmeticError("unexpected OSCAR lattice ranks")
    if l26_basis != sp.eye(26):
        raise ArithmeticError("serialized L26 basis is not the identity")
    if weyl.shape != (1, 26) or preprocessing_vector.shape != (1, 26):
        raise ArithmeticError("serialized Borcherds vectors have wrong size")
    if gram_l26 != gram_l26.T or abs(int(gram_l26.det())) != 1:
        raise ArithmeticError("ambient L26 form is not symmetric unimodular")
    if asdict(inertia(_rows(gram_l26))) != {
        "positive": 1,
        "negative": 25,
        "zero": 0,
    }:
        raise ArithmeticError("ambient lattice is not II(1,25)")

    records = parse_gap_records(SHIMADA_DATA)
    record = records["12A"]
    gram_l10 = sp.Matrix(record["GramL10"])
    walls = sp.Matrix(record["walls"])
    configuration = sp.Matrix(record["configuration"])
    if walls * gram_l10 * walls.T != configuration:
        raise ArithmeticError("Shimada Type-I wall Gram reconstruction failed")

    gram_sy = sy * gram_l26 * sy.T
    gram_sx = sx * gram_l26 * sx.T
    if gram_sy != 2 * gram_l10 or sx[:10, :] != sy:
        raise ArithmeticError("OSCAR and Shimada do not use the same L10 basis")
    if int(gram_sx.det()) != 4 or asdict(inertia(_rows(gram_sx))) != {
        "positive": 1,
        "negative": 18,
        "zero": 0,
    }:
        raise ArithmeticError("OSCAR K3-cover Picard lattice changed")
    if _smith_diagonal(sx) != [1] * 19:
        raise ArithmeticError("SX is not primitive in the Borcherds ambient lattice")

    # If z=(a,c) is an SX coordinate row, its invariant projection is
    # r=2a+cD.  This determines the Enriques involution on all of SX.
    d_matrix = sx[10:, :] * gram_l26 * sy.T * gram_l10.inv()
    if any(value.q != 1 for value in d_matrix):
        raise ArithmeticError("projection/glue matrix D is not integral")
    d_matrix = sp.Matrix(_rows(d_matrix))
    tau = sp.Matrix.vstack(
        sp.Matrix.hstack(sp.eye(10), sp.zeros(10, 9)),
        sp.Matrix.hstack(d_matrix, -sp.eye(9)),
    )
    projection = sp.Matrix.vstack(2 * sp.eye(10), d_matrix)
    if tau * tau != sp.eye(19) or tau * gram_sx * tau.T != gram_sx:
        raise ArithmeticError("reconstructed Enriques involution is not an isometry")
    if (sp.eye(19) + tau)[:, :10] != projection or (sp.eye(19) + tau)[:, 10:] != sp.zeros(19, 9):
        raise ArithmeticError("invariant projection formula changed")

    invariant_basis = sp.Matrix.hstack(sp.eye(10), sp.zeros(10, 9))
    anti_basis = sp.Matrix.hstack(-d_matrix, 2 * sp.eye(9))
    if invariant_basis * tau != invariant_basis or anti_basis * tau != -anti_basis:
        raise ArithmeticError("wrong invariant/anti-invariant eigenspaces")
    if invariant_basis * gram_sx * anti_basis.T != sp.zeros(10, 9):
        raise ArithmeticError("invariant and anti-invariant lattices are not orthogonal")
    if _smith_diagonal(anti_basis) != [1] * 9:
        raise ArithmeticError("anti-invariant lattice is not saturated in SX")
    anti_gram = anti_basis * gram_sx * anti_basis.T
    if asdict(inertia(_rows(anti_gram))) != {
        "positive": 0,
        "negative": 9,
        "zero": 0,
    }:
        raise ArithmeticError("anti-invariant algebraic lattice is not negative definite")
    if any(anti_gram[index, index] % 4 for index in range(9)) or any(
        anti_gram[left, right] % 2
        for left in range(9)
        for right in range(left)
    ):
        raise ArithmeticError("anti-lattice parity check failed")
    anti_cartan = -anti_gram / 2
    if any(value.q != 1 for value in anti_cartan):
        raise ArithmeticError("anti-lattice is not twice an integral root lattice")
    anti_cartan = sp.Matrix(_rows(anti_cartan))
    if not _is_signed_e8_plus_a1_cartan(anti_cartan):
        raise ArithmeticError("anti-lattice is not a signed E8+A1 Cartan form")
    if int(anti_gram.det()) != -1024 or _smith_diagonal(anti_gram) != [2] * 8 + [4]:
        raise ArithmeticError("anti-lattice discriminant data changed")
    combined = sp.Matrix.vstack(invariant_basis, anti_basis)
    invariant_anti_index = abs(int(combined.det()))
    if invariant_anti_index != 512:
        raise ArithmeticError("invariant/anti-invariant glue index changed")
    if abs(int(gram_sy.det() * anti_gram.det())) // invariant_anti_index**2 != abs(int(gram_sx.det())):
        raise ArithmeticError("invariant/anti-invariant determinant formula failed")

    # The lift equation reduces to the roots of E8+A1 in the c-coordinates.
    fiber_roots, fiber_enumeration = _pari_short_vectors(anti_cartan, 2)
    if len(fiber_roots) != 242 or any(
        int((sp.Matrix(root).T * anti_cartan * sp.Matrix(root))[0]) != 2
        for root in fiber_roots
    ):
        raise ArithmeticError("the E8+A1 root census changed")

    lift_pairs: list[list[tuple[int, ...]]] = []
    lift_records = []
    for wall_index in range(walls.rows):
        root = tuple(int(value) for value in walls.row(wall_index))
        candidates = []
        for c_vector in fiber_roots:
            c_row = sp.Matrix([c_vector])
            c_d = c_row * d_matrix
            if all((root[column] - int(c_d[column])) % 2 == 0 for column in range(10)):
                a = tuple(
                    (root[column] - int(c_d[column])) // 2
                    for column in range(10)
                )
                candidate = a + c_vector
                candidate_row = sp.Matrix([candidate])
                if int((candidate_row * gram_sx * candidate_row.T)[0]) != -2:
                    raise ArithmeticError("lift candidate is not a minus-two root")
                candidates.append(candidate)
        candidates = sorted(set(candidates))
        if len(candidates) != 2:
            raise ArithmeticError(f"wall {wall_index + 1} has {len(candidates)} lifts")
        invariant_sum = root + (0,) * 9
        if tuple(candidates[0][index] + candidates[1][index] for index in range(19)) != invariant_sum:
            raise ArithmeticError("the two lifts do not sum to the invariant pullback")
        first = sp.Matrix([candidates[0]])
        second = sp.Matrix([candidates[1]])
        if first * tau != second or second * tau != first:
            raise ArithmeticError("the lift pair is not exchanged by tau")
        if int((first * gram_sx * second.T)[0]) != 0:
            raise ArithmeticError("the exchanged lift pair is not disjoint")
        if first * projection != sp.Matrix([root]):
            raise ArithmeticError("lift has the wrong downstairs projection")
        lift_pairs.append(candidates)
        lift_records.append(
            {
                "wall_one_based": wall_index + 1,
                "downstairs_root_L10_coordinates": list(root),
                "lift_0_SX_coordinates": list(candidates[0]),
                "lift_1_SX_coordinates": list(candidates[1]),
                "fiber_root_c_coordinates": [
                    list(candidates[0][10:]),
                    list(candidates[1][10:]),
                ],
                "both_squares": [-2, -2],
                "mutual_intersection": 0,
            }
        )

    all_lifts = [lift for pair in lift_pairs for lift in pair]
    if len(set(all_lifts)) != 24:
        raise ArithmeticError("the twenty-four lifted roots are not distinct")
    intersection_histogram: Counter[int] = Counter()
    for left in range(len(all_lifts)):
        left_row = sp.Matrix([all_lifts[left]])
        for right in range(left + 1, len(all_lifts)):
            value = int((left_row * gram_sx * sp.Matrix([all_lifts[right]]).T)[0])
            if value < 0:
                raise ArithmeticError("two distinct proposed curve lifts intersect negatively")
            intersection_histogram[value] += 1
    if intersection_histogram != Counter({0: 246, 1: 28, 2: 2}):
        raise ArithmeticError("twenty-four-curve intersection histogram changed")

    for left in range(walls.rows):
        downstairs = int(configuration[left, left])
        if downstairs != -2:
            raise ArithmeticError("downstairs wall is not a root")
        for right in range(left + 1, walls.rows):
            m = int(configuration[left, right])
            for left_sheet in range(2):
                b_left = sp.Matrix([lift_pairs[left][left_sheet]])
                same = int((b_left * gram_sx * sp.Matrix([lift_pairs[right][0]]).T)[0])
                cross = int((b_left * gram_sx * sp.Matrix([lift_pairs[right][1]]).T)[0])
                if same < 0 or cross < 0 or same + cross != m:
                    raise ArithmeticError("signed two-sheet intersection gate failed")

    # Any choice of one sheet gives the same saturated anti span.  We freeze
    # the lexicographically first lift for a deterministic package.
    selected_lifts = sp.Matrix([list(pair[0]) for pair in lift_pairs])
    selected_c = selected_lifts[:, 10:]
    anti_differences = []
    for index in range(walls.rows):
        invariant_sum = sp.Matrix([_rows(walls)[index] + [0] * 9])
        anti_differences.append(2 * selected_lifts[index, :] - invariant_sum)
    anti_differences_matrix = sp.Matrix.vstack(*anti_differences)
    if anti_differences_matrix != selected_c * anti_basis:
        raise ArithmeticError("anti-difference factorization failed")
    if selected_c.rank() != 9 or _smith_diagonal(selected_c) != [1] * 9:
        raise ArithmeticError("selected differences do not generate the full anti lattice")
    if _smith_diagonal(anti_differences_matrix) != [1] * 9:
        raise ArithmeticError("selected anti span is not saturated in SX")

    # The projection of OSCAR's Borcherds Weyl vector is a much smaller ample
    # orientation than the generic sum-of-rays witness used downstairs.
    h_down = weyl * gram_l26 * sy.T * gram_sy.inv()
    if any(value.q != 1 for value in h_down):
        raise ArithmeticError("Weyl-vector projection is not integral")
    h_down = sp.Matrix(_rows(h_down))
    if _rows(h_down) != [[38, 23, 46, 70, 57, 45, 34, 24, 15, 7]]:
        raise ArithmeticError("small Type-I ample orientation changed")
    h_sx = sp.Matrix.hstack(h_down, sp.zeros(1, 9))
    h_square = int((h_sx * gram_sx * h_sx.T)[0])
    downstairs_h_square = int((h_down * gram_l10 * h_down.T)[0])
    if h_square != 140 or downstairs_h_square != 70 or h_sx * tau != h_sx:
        raise ArithmeticError("ample orientation square or invariance failed")
    wall_degrees = [int(value) for value in walls * gram_l10 * h_down.T]
    if wall_degrees != [1] * 10 + [8, 8]:
        raise ArithmeticError("small ample wall degrees changed")

    h_pairing = gram_sx * h_sx.T
    irreducibility_bound = 2 * h_square + 2 * 7**2
    short_form = -h_square * gram_sx + 2 * h_pairing * h_pairing.T
    short_vectors, short_enumeration = _pari_short_vectors(
        short_form,
        irreducibility_bound,
    )
    if len(short_vectors) != 472:
        raise ArithmeticError("bounded SX short-vector census changed")
    positive_small_roots = []
    positive_degree_histogram: Counter[int] = Counter()
    for vector in short_vectors:
        column = sp.Matrix(vector)
        square = int((column.T * gram_sx * column)[0])
        degree = int((column.T * h_pairing)[0])
        if square != -2 or not 1 <= abs(degree) <= 7:
            raise ArithmeticError("bounded enumeration contains an unexpected vector")
        if degree > 0:
            positive_small_roots.append(vector)
            positive_degree_histogram[degree] += 1
    if positive_degree_histogram != Counter({1: 20, 2: 20, 3: 24, 4: 28, 5: 36, 6: 44, 7: 64}):
        raise ArithmeticError("bounded root degree histogram changed")

    l4_records = []
    for wall_index, pair in enumerate(lift_pairs):
        for sheet, lift in enumerate(pair):
            lift_row = sp.Matrix([lift])
            degree = int((lift_row * h_pairing)[0])
            if degree != wall_degrees[wall_index]:
                raise ArithmeticError("lift and downstairs wall degrees disagree")
            tested = [
                vector
                for vector in positive_small_roots
                if int((sp.Matrix(vector).T * h_pairing)[0]) < degree
            ]
            bad = [
                vector
                for vector in tested
                if int((lift_row * gram_sx * sp.Matrix(vector))[0]) < 0
            ]
            if bad:
                raise ArithmeticError("L4 found a smaller-degree decomposing root")
            l4_records.append(
                {
                    "wall_one_based": wall_index + 1,
                    "sheet": sheet,
                    "H_degree": degree,
                    "smaller_positive_degree_roots_tested": len(tested),
                    "negative_pairing_witness_count": 0,
                    "L4_pass": True,
                }
            )

    type_i_archive = _type_i_archive_record()
    anti_minimum_divisibility = 4
    checks = {
        "L1_integral_half_sum_glue_and_nonnegative_signed_intersections": True,
        "L2_saturated_negative_definite_rank_9_anti_span": True,
        "L2_anti_span_rootless_by_norm_divisibility": anti_minimum_divisibility == 4,
        "L3_exact_rank_19_Picard_lattice_and_invariant_positive_orientation": True,
        "L3_H_root_regular_by_complete_degree_zero_root_search": True,
        "L4_all_24_lifts_pass_bounded_root_irreducibility": True,
        "downstairs_whole_cone_and_literal_Vinberg_finite_volume_gates_imported_from_pinned_archive": True,
    }
    if not all(checks.values()):
        raise ArithmeticError("a good-lift-package gate failed")

    return {
        "schema": "enriques_type_i_explicit_lift_package_v1",
        "scope_and_lineage": {
            "result": (
                "explicit good L1-L4 K3-cover lift package for the known Type-I "
                "finite-volume noncompact H9/12 full-nef Enriques chamber"
            ),
            "new_abstract_coxeter_polytope_claimed": False,
            "new_algebraic_surface_claimed": False,
            "project_gain": (
                "first exact worked positive package for the variable Enriques lift theorem; "
                "it exposes the finite data the downstairs compiler must recover"
            ),
            "remaining_projection_problem": (
                "derive and exhaust such anti-invariant/glue packages from the downstairs "
                "Coxeter datum without starting from OSCAR's K3-cover Picard lattice"
            ),
        },
        "pinned_inputs": {
            "oscar_version": "1.4.0",
            "oscar_function": "generic_enriques_surface(172)",
            "oscar_archived_path": str(OSCAR_DATA.relative_to(ROOT)),
            "oscar_archived_sha256": OSCAR_ARCHIVE_SHA256,
            "oscar_upstream_sha256_before_terminal_newline": OSCAR_UPSTREAM_SHA256,
            "oscar_upstream_url": (
                "https://raw.githubusercontent.com/oscar-system/Oscar.jl/v1.4.0/"
                "data/TauTaubarGenericEnriquesSurfaces/"
                "TauTaubarGenericEnriquesSurfaceNo172.mrdi"
            ),
            "shimada_path": str(SHIMADA_DATA.relative_to(ROOT)),
            "shimada_sha256": SHIMADA_SHA256,
            "full_nef_archive_path": str(FULL_NEF_ARCHIVE.relative_to(ROOT)),
            "full_nef_archive_sha256": FULL_NEF_ARCHIVE_SHA256,
            "variable_theorem_archive_path": str(VARIABLE_THEOREM_ARCHIVE.relative_to(ROOT)),
            "variable_theorem_archive_sha256": VARIABLE_THEOREM_ARCHIVE_SHA256,
        },
        "surface_and_lattices": {
            "surface": "very general Enriques Type I / OSCAR number 172",
            "downstairs_lattice": "L10=U+E8(-1)",
            "downstairs_gram": _rows(gram_l10),
            "K3_cover_Picard_rank": 19,
            "K3_cover_Picard_determinant": int(gram_sx.det()),
            "K3_cover_Picard_inertia": asdict(inertia(_rows(gram_sx))),
            "K3_cover_Picard_gram": _rows(gram_sx),
            "Borcherds_ambient": "II(1,25); this is not the K3 cohomology lattice",
            "invariant_lattice": {
                "isometry_type": "U(2)+E8(-2)",
                "rank": 10,
                "determinant_absolute": abs(int(gram_sy.det())),
                "gram": _rows(gram_sy),
            },
            "anti_invariant_algebraic_lattice": {
                "isometry_type": "E8(-2)+<-4>",
                "basis_in_SX_coordinates": _rows(anti_basis),
                "gram": _rows(anti_gram),
                "rank": 9,
                "determinant_absolute": abs(int(anti_gram.det())),
                "smith_diagonal": _smith_diagonal(anti_gram),
                "primitive_in_SX_smith_diagonal": _smith_diagonal(anti_basis),
                "negative_definite": True,
                "all_norms_divisible_by": anti_minimum_divisibility,
                "contains_minus_two_root": False,
            },
            "invariant_plus_anti_index_in_SX": invariant_anti_index,
            "determinant_glue_identity": "2^10 * 2^10 / 512^2 = 4",
            "period_domain_complex_dimension": 1,
        },
        "global_involution": {
            "SX_row_coordinate_formula": "tau(a,c)=(a+cD,-c)",
            "D": _rows(d_matrix),
            "tau_matrix_for_row_coordinates": _rows(tau),
            "tau_squared_is_identity": True,
            "tau_is_Picard_isometry": True,
            "projection_formula": "r=2a+cD",
            "projection_matrix": _rows(projection),
        },
        "finite_lift_compiler_calibration": {
            "fiber_positive_form": _rows(anti_cartan),
            "fiber_root_system": "E8+A1",
            "fiber_root_count": len(fiber_roots),
            "fiber_enumeration": fiber_enumeration,
            "each_downstairs_wall_has_exactly_two_lifts": True,
            "lift_records": lift_records,
            "all_24_lifts_distinct": True,
            "all_24_pairwise_intersections_nonnegative": True,
            "all_24_unordered_pair_intersection_histogram": {
                str(value): count for value, count in sorted(intersection_histogram.items())
            },
            "all_4096_sheet_selections_pass_pairwise_nonnegativity": True,
            "frozen_selection": "lexicographically first lift of every wall",
            "frozen_selected_c_matrix": _rows(selected_c),
            "frozen_selected_c_smith_diagonal": _smith_diagonal(selected_c),
            "frozen_anti_difference_gram": _rows(
                anti_differences_matrix * gram_sx * anti_differences_matrix.T
            ),
            "frozen_differences_generate_full_saturated_anti_lattice": True,
        },
        "ample_orientation_and_L4": {
            "source": "orthogonal projection of OSCAR's serialized Borcherds Weyl vector",
            "downstairs_H_coordinates": _rows(h_down)[0],
            "downstairs_H_square": downstairs_h_square,
            "K3_pullback_H_square": h_square,
            "wall_degrees": wall_degrees,
            "strict_on_all_12_downstairs_walls": True,
            "root_regular": True,
            "positive_definite_enumeration_form": "-H^2*Gram(SX)+2*(H-pairing column)^2",
            "irreducibility_enumeration_bound": irreducibility_bound,
            "bounded_enumeration": short_enumeration,
            "all_bounded_vectors_are_minus_two_roots_of_absolute_H_degree_1_to_7": True,
            "positive_root_degree_histogram": {
                str(degree): count
                for degree, count in sorted(positive_degree_histogram.items())
            },
            "degree_zero_root_count": 0,
            "L4_records": l4_records,
            "all_24_lifts_irreducible_by_L4": True,
        },
        "downstairs_polytope_gate": {
            "shimada_record": "12A",
            "hyperbolic_dimension": 9,
            "facet_count": 12,
            "finite_vertex_count": type_i_archive["exact_cone"]["finite_vertex_count"],
            "ideal_vertex_count": type_i_archive["exact_cone"]["ideal_vertex_count"],
            "hyperideal_vertex_count": type_i_archive["exact_cone"]["hyperideal_vertex_count"],
            "literal_vinberg_elliptic_edge_count": type_i_archive["literal_vinberg"]["elliptic_rank_8_edge_count"],
            "literal_vinberg_affine_audit_count": type_i_archive["literal_vinberg"]["connected_affine_component_audit_count"],
            "finite_volume": True,
            "compact": False,
            "complete_full_nef_wall_set": True,
        },
        "good_lift_package_checks": checks,
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
            raise SystemExit(f"archived certificate differs: {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(hashlib.sha256(rendered.encode("utf-8")).hexdigest())


if __name__ == "__main__":
    main()
