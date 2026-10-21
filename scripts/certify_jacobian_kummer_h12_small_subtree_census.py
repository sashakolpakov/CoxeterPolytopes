#!/usr/bin/env python3
"""Certify the JK-H12 one-edge and three-vertex finite-subtree census.

The input chamber is the exact 172-facet Jacobian--Kummer H12 chamber.  We
work with the 112 pairwise Allcock-disjoint doubling walls W61,...,W172.
The script independently reconstructs the 384-element action on the folded
chamber from Dolgachev--Shimada's published wall data, exhausts its orbits on
single walls and unordered wall pairs, and rebuilds every orbit-representative
boundary from integral reflection matrices.

Finite volume is supplied by Allcock's finite-subtree theorem.  The script
also checks the facts that theorem does not encode in the archive: exact
boundary-hyperplane deduplication, acute Coxeter pair labels, primitive
integral roots, a strict timelike interior, and an index-one Smith core.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter, deque
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = (
    ROOT / "data/external/dolgachev_shimada_15_nodal_2019/15nodalcompdata.txt"
)
FOLD = ROOT / "data/enumerations/jacobian_kummer_h12_interim_facet_gram.json"
CARRIER = (
    ROOT / "data/enumerations/jacobian_kummer_h12_full_picard_k3_carrier.json"
)
CORPUS_AUDIT = (
    ROOT / "data/enumerations/jacobian_kummer_h12_path_cube_hostile_audit.json"
)
OUTPUT = (
    ROOT / "data/enumerations/jacobian_kummer_h12_small_subtree_census.json"
)

SCHEMA = "coxeter.jacobian_kummer.h12_small_subtree_census.v1"
STATUS_DATE = "2026-08-08"
PINNED_SHA256 = {
    "published_Dolgachev_Shimada_data": (
        "413d7363faf4574595e9c28329737eb4ad8010e0aebae53afb76762ce9f8966a"
    ),
    "H12_fold_orbits": (
        "53c39642b81bbd7f5be172408bfecb857bd4ed57467c91bdbaed63f9e734fb92"
    ),
    "full_Picard_K3_carrier": (
        "f2c57b9fb0423560c9e238541b0ac1a81ca141917c864929501ebdfa7f294ff2"
    ),
    "bounded_CoxIter_corpus_audit": (
        "5c81e4c311e96a64a86ec377d10aa40b54dcd94f685432fee11d9bc9f8ea5663"
    ),
}

Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]
Permutation = tuple[int, ...]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(payload: Any) -> str:
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(rendered.encode()).hexdigest()


def _extract_assignment(text: str, name: str) -> str:
    marker = f"{name}:="
    start = text.index(marker) + len(marker)
    while text[start].isspace():
        start += 1
    if text[start] != "[":
        raise ArithmeticError(f"published assignment {name} is not a list")
    depth = 0
    quote: str | None = None
    for end in range(start, len(text)):
        character = text[end]
        if quote is not None:
            if character == quote and text[end - 1] != "\\":
                quote = None
            continue
        if character in {'"', "'"}:
            quote = character
        elif character == "[":
            depth += 1
        elif character == "]":
            depth -= 1
            if depth == 0:
                return text[start : end + 1]
    raise ArithmeticError(f"unterminated published assignment {name}")


def _content(vector: Iterable[int]) -> int:
    values = tuple(abs(value) for value in vector if value)
    return math.gcd(*values) if values else 0


def _primitive_key(vector: Sequence[int]) -> Vector:
    content = _content(vector)
    if content == 0:
        raise ArithmeticError("zero vector has no hyperplane key")
    primitive = tuple(int(value) // content for value in vector)
    first = next(value for value in primitive if value)
    return primitive if first > 0 else tuple(-value for value in primitive)


def _identity(rank: int) -> Matrix:
    return tuple(
        tuple(1 if row == column else 0 for column in range(rank))
        for row in range(rank)
    )


def _matvec(matrix: Matrix, vector: Sequence[int]) -> Vector:
    return tuple(
        sum(entry * value for entry, value in zip(row, vector, strict=True))
        for row in matrix
    )


def _matmul(left: Matrix, right: Matrix) -> Matrix:
    rank = len(left)
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column] for middle in range(rank))
            for column in range(rank)
        )
        for row in range(rank)
    )


def _pair(gram: Matrix, left: Sequence[int], right: Sequence[int]) -> int:
    return sum(
        left[row] * gram[row][column] * right[column]
        for row in range(len(left))
        for column in range(len(right))
    )


def _covector(gram: Matrix, vector: Sequence[int]) -> Vector:
    return tuple(
        sum(gram[row][column] * vector[column] for column in range(len(vector)))
        for row in range(len(vector))
    )


def _reflection_matrix(gram: Matrix, root: Vector) -> Matrix:
    square = _pair(gram, root, root)
    covector = _covector(gram, root)
    rank = len(root)
    rows: list[tuple[int, ...]] = []
    for row in range(rank):
        values = []
        for column in range(rank):
            coefficient = Fraction(2 * covector[column], square)
            value = Fraction(int(row == column)) - coefficient * root[row]
            if value.denominator != 1:
                raise ArithmeticError("a selected wall reflection is nonintegral")
            values.append(int(value))
        rows.append(tuple(values))
    matrix = tuple(rows)
    if _matmul(matrix, matrix) != _identity(rank):
        raise ArithmeticError("a selected wall reflection is not involutive")
    return matrix


def _compose_permutations(left: Permutation, right: Permutation) -> Permutation:
    """Permutation induced by the row-vector matrix product M_left*M_right."""

    return tuple(right[left[index]] for index in range(len(left)))


def _wall_permutation(matrix: sp.Matrix, walls: tuple[sp.Matrix, ...]) -> Permutation:
    index = {tuple(wall): position for position, wall in enumerate(walls)}
    images = tuple(tuple(wall * matrix) for wall in walls)
    if any(image not in index for image in images):
        raise ArithmeticError("a published stabilizer generator loses a wall")
    permutation = tuple(index[image] for image in images)
    if len(set(permutation)) != len(walls):
        raise ArithmeticError("a published wall action is not a permutation")
    return permutation


def _load_inputs() -> tuple[
    dict[str, Any], Matrix, tuple[Vector, ...], Vector, dict[str, Any]
]:
    paths = {
        "published_Dolgachev_Shimada_data": PUBLISHED,
        "H12_fold_orbits": FOLD,
        "full_Picard_K3_carrier": CARRIER,
        "bounded_CoxIter_corpus_audit": CORPUS_AUDIT,
    }
    hashes = {name: _sha256(path) for name, path in paths.items()}
    if hashes != PINNED_SHA256:
        raise ArithmeticError(f"a pinned small-subtree input changed: {hashes}")

    carrier = json.loads(CARRIER.read_text(encoding="utf-8"))
    chamber = carrier["chamber"]
    gram = tuple(
        tuple(int(value) for value in row) for row in chamber["carrier"]["gram"]
    )
    roots = tuple(
        tuple(int(value) for value in record["oriented_coordinates"])
        for record in chamber["oriented_wall_system"]["wall_records"]
    )
    interior = tuple(
        int(value)
        for value in chamber["oriented_wall_system"]["strict_interior_coordinates"]
    )
    finite_volume = chamber["finite_volume"]
    if not (
        carrier["all_checks_pass"]
        and chamber["dimension"] == 12
        and chamber["facet_count"] == 172
        and finite_volume["literal_Vinberg"]["classification"]
        == "finite_volume_noncompact"
        and finite_volume["independent_exact_GMP_cone"]["classification"]
        == "finite_volume_noncompact"
        and len(gram) == 13
        and len(roots) == 172
        and _pair(gram, interior, interior) > 0
        and all(_content(root) == 1 for root in roots)
        and all(_pair(gram, root, root) < 0 for root in roots)
        and all(_pair(gram, interior, root) > 0 for root in roots)
    ):
        raise ArithmeticError("the pinned full-Picard seed failed exact gates")

    corpus = json.loads(CORPUS_AUDIT.read_text(encoding="utf-8"))
    bounded = corpus["bounded_coxiter_corpus_replay"]
    if not (
        corpus["all_checks_pass"]
        and bounded["all_exact_checks_pass"]
        and bounded["ordinary_H12_facet_histogram"]
        == {"13": 6, "14": 3, "16": 1, "18": 1}
        and bounded["index_two_H12_facet_histogram"]
        == {"18": 1, "20": 1, "25": 1, "35": 1}
    ):
        raise ArithmeticError("the bounded CoxIter denominator changed")
    return carrier, gram, roots, interior, bounded


def _folded_chamber_actions() -> tuple[Permutation, ...]:
    text = PUBLISHED.read_text(encoding="utf-8").replace("\\\n", "")
    generators = tuple(
        sp.Matrix(matrix)
        for matrix in ast.literal_eval(
            _extract_assignment(text, "GeneratorsOS16D16")
        )
    )
    wall_text = _extract_assignment(text, "D16WallRecs")
    vector_texts = re.findall(r"\bvector\s*:=\s*\[([^\]]+)\]", wall_text)
    walls = tuple(
        sp.Matrix(
            1,
            17,
            [sp.Rational(entry.strip()) for entry in vector.split(",")],
        )
        for vector in vector_texts
    )
    fold = json.loads(FOLD.read_text(encoding="utf-8"))
    involution = sp.Matrix(fold["group_class"]["representative_matrix"])
    if not (len(generators) == 10 and len(walls) == 316):
        raise ArithmeticError("the published D16 stabilizer parse changed")

    generator_actions = tuple(_wall_permutation(matrix, walls) for matrix in generators)
    involution_action = _wall_permutation(involution, walls)
    identity = tuple(range(316))
    ambient_group: set[Permutation] = {identity}
    queue: deque[Permutation] = deque((identity,))
    while queue:
        element = queue.popleft()
        for generator in generator_actions:
            image = _compose_permutations(element, generator)
            if image not in ambient_group:
                ambient_group.add(image)
                queue.append(image)
    centralizer = tuple(
        element
        for element in ambient_group
        if _compose_permutations(element, involution_action)
        == _compose_permutations(involution_action, element)
    )

    all_orbits = tuple(
        frozenset(int(index) - 1 for index in orbit)
        for orbit in fold["wall_orbits"]
    )
    orbit_index = {orbit: index for index, orbit in enumerate(all_orbits)}
    facet_orbit_indices = tuple(
        int(index) - 1 for index in fold["facet_orbit_indices"]
    )
    facet_position = {
        orbit: position for position, orbit in enumerate(facet_orbit_indices)
    }
    facet_actions: set[Permutation] = set()
    for element in centralizer:
        orbit_action = tuple(
            orbit_index[frozenset(element[wall] for wall in orbit)]
            for orbit in all_orbits
        )
        if any(orbit_action[index] not in facet_position for index in facet_orbit_indices):
            raise ArithmeticError("a centralizer element loses a folded facet")
        facet_actions.add(
            tuple(
                facet_position[orbit_action[index]] for index in facet_orbit_indices
            )
        )
    if not (
        len(ambient_group) == 23_040
        and len(centralizer) == 768
        and len(facet_actions) == 384
        and all(len(set(action)) == 172 for action in facet_actions)
    ):
        raise ArithmeticError("the exact folded chamber action changed")

    return tuple(sorted(facet_actions))


def _orbits(
    elements: Iterable[Any], action_image: Any
) -> tuple[tuple[Any, ...], ...]:
    remaining = set(elements)
    rows: list[tuple[Any, ...]] = []
    while remaining:
        representative = min(remaining)
        orbit = tuple(sorted(action_image(representative)))
        if not set(orbit) <= remaining | (set(elements) - remaining):
            raise ArithmeticError("an orbit escaped its prescribed universe")
        rows.append(orbit)
        remaining.difference_update(orbit)
    return tuple(rows)


def _relation_label(pairing: int, left_square: int, right_square: int) -> str:
    if pairing < 0 or left_square >= 0 or right_square >= 0:
        raise ArithmeticError("a descendant pair is not acute-spacelike")
    numerator = pairing * pairing
    denominator = left_square * right_square
    if pairing == 0:
        return "m2"
    if 4 * numerator == denominator:
        return "m3"
    if 2 * numerator == denominator:
        return "m4"
    if 4 * numerator == 3 * denominator:
        return "m6"
    if numerator == denominator:
        return "parallel"
    if numerator > denominator:
        ratio = Fraction(numerator, denominator)
        return f"ultraparallel:{ratio.numerator}/{ratio.denominator}"
    raise ArithmeticError(
        f"non-Coxeter normalized square {numerator}/{denominator}"
    )


def _boundary(
    gram: Matrix,
    seed_roots: tuple[Vector, ...],
    seed_interior: Vector,
    selected: tuple[int, ...],
    core_smith_diagonal: list[int],
    *,
    compute_diagram: bool = True,
) -> dict[str, Any]:
    reflections = {
        index: _reflection_matrix(gram, seed_roots[index]) for index in selected
    }
    if len(selected) == 1:
        transformations = (_identity(13), reflections[selected[0]])
    elif len(selected) == 2:
        left, right = selected
        transformations = (
            _identity(13),
            reflections[left],
            _matmul(reflections[left], reflections[right]),
        )
    else:
        raise ArithmeticError("this census only supports one-edge and length-two paths")

    transformation_set = set(transformations)
    internal: set[Vector] = set()
    for transformation in transformations:
        for index in selected:
            if _matmul(transformation, reflections[index]) in transformation_set:
                internal.add(
                    _primitive_key(_matvec(transformation, seed_roots[index]))
                )
    if len(internal) != len(selected):
        raise ArithmeticError("the finite path has the wrong internal-mirror count")

    interior = tuple(
        sum(_matvec(transformation, seed_interior)[coordinate]
            for transformation in transformations)
        for coordinate in range(13)
    )
    boundary: dict[Vector, Vector] = {}
    for transformation in transformations:
        for seed_root in seed_roots:
            image = _matvec(transformation, seed_root)
            key = _primitive_key(image)
            if key in internal:
                continue
            pairing = _pair(gram, interior, image)
            if pairing == 0:
                raise ArithmeticError("the summed path interior lies on a boundary")
            oriented = image if pairing > 0 else tuple(-value for value in image)
            boundary.setdefault(key, oriented)

    ordered_keys = tuple(sorted(boundary))
    roots = tuple(boundary[key] for key in ordered_keys)
    core_keys = {_primitive_key(root) for root in seed_roots[:60]}
    if not core_keys <= set(ordered_keys):
        raise ArithmeticError("the fixed W1--W60 Picard core left the boundary")

    if core_smith_diagonal != [1] * 13:
        raise ArithmeticError("the W1--W60 core no longer spans Picard with index one")

    covectors = tuple(_covector(gram, root) for root in roots)
    squares = tuple(
        sum(root[index] * covector[index] for index in range(13))
        for root, covector in zip(roots, covectors, strict=True)
    )
    integral_reflections = 0
    for root, square, covector in zip(roots, squares, covectors, strict=True):
        if not (
            _content(root) == 1
            and square < 0
            and all(Fraction(2 * value, square).denominator == 1 for value in covector)
        ):
            raise ArithmeticError("a boundary root lost primitive integral reflectivity")
        integral_reflections += 1

    interior_pairings = tuple(_pair(gram, interior, root) for root in roots)
    if not (
        _pair(gram, interior, interior) > 0
        and all(value > 0 for value in interior_pairings)
    ):
        raise ArithmeticError("a descendant failed its timelike/rank gates")

    record: dict[str, Any] = {
        "copy_count_and_volume_multiple": len(transformations),
        "internal_mirror_count": len(internal),
        "facet_count": len(roots),
        "full_linear_rank": 13,
        "strict_interior_square": _pair(gram, interior, interior),
        "minimum_strict_interior_facet_pairing": min(interior_pairings),
        "primitive_integral_reflection_root_count": integral_reflections,
        "fixed_boundary_core_one_based": "W1--W60",
        "fixed_boundary_core_Smith_diagonal": core_smith_diagonal,
        "marked_boundary_root_sha256": _canonical_sha256(
            [list(root) for root in roots]
        ),
    }
    if not compute_diagram:
        return record

    vertex_profiles = [Counter() for _root in roots]
    relation_profiles = [Counter() for _root in roots]
    coxeter_matrix_profiles = [Counter() for _root in roots]
    edge_labels: Counter[str] = Counter()
    for left in range(len(roots)):
        for right in range(left):
            pairing = sum(
                roots[left][index] * covectors[right][index] for index in range(13)
            )
            label = _relation_label(pairing, squares[left], squares[right])
            relation = label.split(":", 1)[0]
            coxeter_matrix_label = (
                "infinity"
                if relation in {"parallel", "ultraparallel"}
                else relation
            )
            edge_labels[label] += 1
            vertex_profiles[left][label] += 1
            vertex_profiles[right][label] += 1
            relation_profiles[left][relation] += 1
            relation_profiles[right][relation] += 1
            coxeter_matrix_profiles[left][coxeter_matrix_label] += 1
            coxeter_matrix_profiles[right][coxeter_matrix_label] += 1

    exact_profiles = tuple(
        sorted(tuple(sorted(profile.items())) for profile in vertex_profiles)
    )
    relation_profile_multiset = tuple(
        sorted(tuple(sorted(profile.items())) for profile in relation_profiles)
    )
    coxeter_matrix_profile_multiset = tuple(
        sorted(tuple(sorted(profile.items())) for profile in coxeter_matrix_profiles)
    )
    exact_edge_histogram = tuple(sorted(edge_labels.items()))
    coarse = Counter()
    for label, count in edge_labels.items():
        coarse[label.split(":", 1)[0]] += count
    record.update(
        {
            "coarse_relation_histogram": dict(sorted(coarse.items())),
            "exact_edge_label_histogram_sha256": _canonical_sha256(
                exact_edge_histogram
            ),
            "exact_vertex_edge_profile_multiset_sha256": _canonical_sha256(
                exact_profiles
            ),
            "distance_forgotten_vertex_relation_profile_multiset_sha256": _canonical_sha256(
                relation_profile_multiset
            ),
            "ordinary_Coxeter_matrix_vertex_profile_multiset_sha256": _canonical_sha256(
                coxeter_matrix_profile_multiset
            ),
            "exact_vertex_profile_class_count": len(set(exact_profiles)),
            "maximum_exact_vertex_profile_multiplicity": max(
                Counter(exact_profiles).values()
            ),
            "_exact_vertex_profiles": exact_profiles,
            "_distance_forgotten_vertex_profiles": relation_profile_multiset,
            "_ordinary_Coxeter_matrix_vertex_profiles": coxeter_matrix_profile_multiset,
        }
    )
    return record


def _public(record: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in record.items() if not key.startswith("_")}


def build() -> dict[str, Any]:
    carrier, gram, roots, interior, bounded_corpus = _load_inputs()
    actions = _folded_chamber_actions()

    # Check the exact carrier Gram action on all 172 walls.  Primitive roots
    # may be flipped in sign, so compare pairings after strict-interior orientation.
    root_gram = tuple(tuple(_pair(gram, left, right) for right in roots) for left in roots)
    if any(
        root_gram[left][right] != root_gram[action[left]][action[right]]
        for action in actions
        for left in range(172)
        for right in range(172)
    ):
        raise ArithmeticError("a reconstructed chamber action fails the exact root Gram")

    stock = tuple(range(60, 172))
    stock_set = set(stock)
    if any({action[index] for index in stock} != stock_set for action in actions):
        raise ArithmeticError("the 384-element action does not preserve W61--W172")

    # Recheck Allcock's two hypotheses directly from the seed Gram.
    seed_squares = tuple(_pair(gram, root, root) for root in roots)
    for left in stock:
        for right in range(172):
            label = _relation_label(
                root_gram[left][right], seed_squares[left], seed_squares[right]
            ) if left != right else "self"
            if label in {"m3"}:
                raise ArithmeticError("a selected wall is not doubling")
        for right in stock:
            if left < right:
                label = _relation_label(
                    root_gram[left][right], seed_squares[left], seed_squares[right]
                )
                if label not in {"parallel"} and not label.startswith("ultraparallel:"):
                    raise ArithmeticError("two selected doubling walls meet in H12")

    wall_orbits = _orbits(
        stock,
        lambda index: {action[index] for action in actions},
    )
    if [len(orbit) for orbit in wall_orbits] != [16, 96]:
        raise ArithmeticError("the selected single-wall orbit census changed")

    core_matrix = sp.Matrix(roots[:60])
    core_smith = smith_normal_form(core_matrix, domain=ZZ)
    core_smith_diagonal = [
        abs(int(core_smith[index, index]))
        for index in range(min(core_smith.rows, core_smith.cols))
        if core_smith[index, index]
    ]
    if core_smith_diagonal != [1] * 13:
        raise ArithmeticError("the fixed W1--W60 core is not index one")

    unordered_pairs = {
        (left, right)
        for left in stock
        for right in stock
        if left < right
    }
    pair_orbits = _orbits(
        unordered_pairs,
        lambda pair: {
            tuple(sorted((action[pair[0]], action[pair[1]])))
            for action in actions
        },
    )
    if not (
        len(unordered_pairs) == 6_216
        and len(pair_orbits) == 35
        and sum(len(orbit) for orbit in pair_orbits) == 6_216
    ):
        raise ArithmeticError("the complete unordered-pair orbit census changed")

    # Exhaust every one-edge double directly, not merely its two orbit reps.
    single_rows_by_wall: dict[int, dict[str, Any]] = {}
    single_profiles: dict[int, tuple[Any, ...]] = {}
    single_representatives = {orbit[0] for orbit in wall_orbits}
    for wall in stock:
        record = _boundary(
            gram,
            roots,
            interior,
            (wall,),
            core_smith_diagonal,
            compute_diagram=wall in single_representatives,
        )
        if wall in single_representatives:
            single_profiles[wall] = record["_exact_vertex_profiles"]
        single_rows_by_wall[wall] = {
            "selected_wall_one_based": wall + 1,
            **_public(record),
        }
    single_types = []
    for orbit in wall_orbits:
        representative = orbit[0]
        row = single_rows_by_wall[representative]
        inherited_fields = {
            key: row[key]
            for key in (
                "coarse_relation_histogram",
                "exact_edge_label_histogram_sha256",
                "exact_vertex_edge_profile_multiset_sha256",
                "distance_forgotten_vertex_relation_profile_multiset_sha256",
                "ordinary_Coxeter_matrix_vertex_profile_multiset_sha256",
                "exact_vertex_profile_class_count",
                "maximum_exact_vertex_profile_multiplicity",
            )
        }
        for wall in orbit:
            single_rows_by_wall[wall].update(inherited_fields)
            single_rows_by_wall[wall]["diagram_type_representative_wall_one_based"] = (
                representative + 1
            )
        single_types.append(
            {
                "representative_wall_one_based": representative + 1,
                "seed_wall_orbit_size": len(orbit),
                "facet_count": row["facet_count"],
                "coarse_relation_histogram": row["coarse_relation_histogram"],
                "exact_vertex_edge_profile_multiset_sha256": row[
                    "exact_vertex_edge_profile_multiset_sha256"
                ],
            }
        )
    if not (
        [row["facet_count"] for row in single_types] == [317, 326]
        and len(set(single_profiles.values())) == 2
    ):
        raise ArithmeticError("the complete single-double type theorem changed")
    single_rows = [single_rows_by_wall[wall] for wall in stock]

    path_rows = []
    path_profiles: list[tuple[Any, ...]] = []
    distance_forgotten_path_profiles: list[tuple[Any, ...]] = []
    coxeter_matrix_path_profiles: list[tuple[Any, ...]] = []
    for orbit_number, orbit in enumerate(pair_orbits, 1):
        left, right = orbit[0]
        generator_label = _relation_label(
            root_gram[left][right], seed_squares[left], seed_squares[right]
        )
        record = _boundary(
            gram,
            roots,
            interior,
            (left, right),
            core_smith_diagonal,
        )
        path_profiles.append(record["_exact_vertex_profiles"])
        distance_forgotten_path_profiles.append(
            record["_distance_forgotten_vertex_profiles"]
        )
        coxeter_matrix_path_profiles.append(
            record["_ordinary_Coxeter_matrix_vertex_profiles"]
        )
        path_rows.append(
            {
                "orbit_number": orbit_number,
                "representative_walls_one_based": [left + 1, right + 1],
                "raw_unordered_pair_orbit_size": len(orbit),
                "generator_pair_relation": generator_label,
                **_public(record),
            }
        )

    facet_histogram = Counter(row["facet_count"] for row in path_rows)
    raw_pair_facet_histogram = Counter()
    for row in path_rows:
        raw_pair_facet_histogram[row["facet_count"]] += row[
            "raw_unordered_pair_orbit_size"
        ]
    generator_histogram = Counter(
        "parallel"
        if row["generator_pair_relation"] == "parallel"
        else "ultraparallel"
        for row in path_rows
    )
    if not (
        facet_histogram == Counter({462: 4, 471: 7, 480: 24})
        and raw_pair_facet_histogram == Counter({462: 120, 471: 1_536, 480: 4_560})
        and generator_histogram == Counter({"parallel": 6, "ultraparallel": 29})
        and len(set(path_profiles)) == 35
        and len(set(distance_forgotten_path_profiles)) == 35
        and len(set(coxeter_matrix_path_profiles)) == 35
        and len(
            {
                row["exact_vertex_edge_profile_multiset_sha256"]
                for row in path_rows
            }
        )
        == 35
    ):
        raise ArithmeticError("the 35-type path theorem changed")

    new_counts = {462, 471, 480}
    corpus_counts = {
        int(value) for value in bounded_corpus["ordinary_H12_facet_histogram"]
    } | {
        int(value) for value in bounded_corpus["index_two_H12_facet_histogram"]
    }
    if new_counts & corpus_counts:
        raise ArithmeticError("a new path header appeared in the bounded corpus")

    checks = {
        "all_pinned_inputs_match": True,
        "published_D16_group_order_is_23040": True,
        "trace_nine_centralizer_order_is_768": True,
        "folded_chamber_action_order_is_384": len(actions) == 384,
        "all_folded_actions_preserve_the_exact_172_root_Gram": True,
        "W61_through_W172_are_112_pairwise_disjoint_doubling_walls": True,
        "all_112_single_doubles_rebuilt_directly": len(single_rows) == 112,
        "single_double_census_has_exactly_two_diagram_types": len(single_types) == 2,
        "all_6216_unordered_generator_pairs_partitioned": sum(
            len(orbit) for orbit in pair_orbits
        )
        == 6_216,
        "all_35_path_orbit_representatives_rebuilt_directly": len(path_rows) == 35,
        "all_35_exact_Coxeter_vertex_profile_invariants_are_distinct": len(
            set(path_profiles)
        )
        == 35,
        "all_35_remain_distinct_after_forgetting_ultraparallel_distances": len(
            set(distance_forgotten_path_profiles)
        )
        == 35,
        "all_35_remain_distinct_as_ordinary_Coxeter_matrices": len(
            set(coxeter_matrix_path_profiles)
        )
        == 35,
        "every_descendant_has_primitive_integral_roots_and_Smith_1_to_13_core": all(
            row["primitive_integral_reflection_root_count"] == row["facet_count"]
            and row["fixed_boundary_core_Smith_diagonal"] == [1] * 13
            for row in single_rows + path_rows
        ),
        "every_descendant_passes_exact_acute_Coxeter_pair_gate_directly_or_by_seed_action": True,
        "new_462_471_480_headers_absent_from_pinned_CoxIter_lanes": True,
        "worldwide_priority_explicitly_open": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"small-subtree census failed: {checks}")

    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": True,
        "pinned_inputs": {
            name: {"sha256": digest}
            for name, digest in PINNED_SHA256.items()
        },
        "seed": {
            "id": "JK-H12/172",
            "hyperbolic_dimension": 12,
            "facet_count": 172,
            "finite_volume": True,
            "compact": False,
            "full_Picard_lattice": "U(2)+D5(-1)^2+<-4>",
            "K3_moduli_dimension": 7,
            "source_carrier_theorem": carrier["theorem"],
        },
        "selected_stock": {
            "walls_one_based": "W61--W172",
            "wall_count": 112,
            "all_doubling": True,
            "all_pairwise_disjoint_in_Allcock_sense": True,
            "folded_chamber_action_order": len(actions),
            "wall_orbit_sizes": [len(orbit) for orbit in wall_orbits],
        },
        "complete_one_edge_census": {
            "raw_selected_wall_count": len(single_rows),
            "exact_diagram_type_count": len(single_types),
            "type_representatives": single_types,
            "all_wall_records": single_rows,
            "finite_volume_noncompact": True,
            "reflection_index_and_volume_multiple": 2,
            "project_status": (
                "classification refinement of the already recorded H12/317 and "
                "H12/326 examples; no new facet count"
            ),
        },
        "complete_three_vertex_coloured_path_census": {
            "raw_unordered_generator_pair_count": len(unordered_pairs),
            "exact_seed_action_orbit_count": len(pair_orbits),
            "pair_orbit_size_histogram": {
                str(key): value
                for key, value in sorted(Counter(map(len, pair_orbits)).items())
            },
            "generator_relation_orbit_histogram": dict(sorted(generator_histogram.items())),
            "diagram_type_count": len(path_rows),
            "all_35_pairwise_nonisomorphic": True,
            "nonisomorphism_certificate": (
                "the sorted multiset of edge-label histograms at vertices is "
                "different for every representative even after all parallel and "
                "ultraparallel pairs are collapsed to the ordinary Coxeter-matrix "
                "label infinity"
            ),
            "diagram_type_facet_histogram": {
                str(key): value for key, value in sorted(facet_histogram.items())
            },
            "raw_pair_facet_histogram": {
                str(key): value
                for key, value in sorted(raw_pair_facet_histogram.items())
            },
            "orbit_representatives": path_rows,
            "finite_volume_noncompact": True,
            "reflection_index_and_volume_multiple": 3,
            "new_to_project_facet_counts": [462, 471, 480],
        },
        "allcock_and_direct_certification": {
            "allcock_hypotheses_checked_directly": [
                "the seed is a finite-volume Coxeter polytope",
                "every selected wall is doubling",
                "the selected walls are pairwise disjoint (parallel or ultraparallel)",
            ],
            "finite_subtree_consequence": (
                "each two- or three-vertex union is a Coxeter polytope, with "
                "index and volume multiplier equal to its vertex count"
            ),
            "noncompactness": (
                "a finite-index subgroup of the nonuniform seed lattice is nonuniform"
            ),
            "direct_extra_gates": [
                "exact reflected boundary-root reconstruction and hyperplane deduplication",
                "strict timelike summed interior",
                "every representative boundary pair has an exact Coxeter/parallel/ultraparallel label; seed-action equivariance transports the two single types",
                "every primitive boundary reflection is integral",
                "the exterior W1--W60 core has Smith diagonal 1^13",
            ],
        },
        "bounded_lineage": {
            "CoxIter_commit": bounded_corpus["commit"],
            "ordinary_H12_file_count": bounded_corpus["ordinary_file_count"],
            "index_two_file_count": bounded_corpus["index_two_file_count"],
            "headers_462_471_480_absent": True,
            "evidence_grade": "B1: bounded same-header absence only",
            "worldwide_priority": "open",
            "warning": (
                "Allcock already proved exponential H12 abundance; these exact "
                "surface-derived diagrams have not received a publication-complete audit"
            ),
        },
        "carrier_scope": {
            "same_full_Picard_K3_family": True,
            "weak_divisor_reflection_carrier_only": True,
            "full_nef_or_effectivity_claim": False,
            "surface_automorphism_claim": False,
        },
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    payload = build()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.exists():
            raise SystemExit(f"archive is missing: {arguments.output}")
        if arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archive differs from exact small-subtree rebuild")
        print("archive_check=PASS")
    else:
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("one_edge_types=2")
    print("three_vertex_path_types=35")
    print("new_facet_counts=462,471,480")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
