"""Close all subgroup-fixed sections of the rank-18 and rank-19 K3 stocks.

The two source configurations are complete published nef-wall stocks.  Their
exact weighted wall-Gram groups are C2 and S3.  This certifier enumerates all
element and subgroup conjugacy classes, lifts every diagram symmetry to the
integral Neron--Severi lattice, and retains every proper subgroup with fixed
rank at least eleven.

There is only one retained class: the transposition class in the rank-19 S3.
It is audited directly (strict facets, exact Coxeter relations, exact cone
rays, and literal Vinberg--CoxIter finite-covolume criterion), then matched
entrywise to the already frozen quotient-nef fold.  All other new subgroup
types fail the rank gate.  A bounded exact diagram comparison covers every
positive in the local lineage ledger and every H10/13 record in the pinned
ordinary and index-two-doubling CoxIter snapshots.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import tarfile
from typing import Any, Iterable, Sequence

import networkx as nx
import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.certify_jacobian_kummer_h10_maximal_reflection_lineage import (
    _parse_coxiter,
    _permutation as _relation_permutation,
)
from scripts.certify_k3_high_rank_all_involution_folds import (
    _cycles,
    _permutation_from_cycles,
    exact_weighted_gram_automorphisms,
)
from scripts.certify_k3_high_rank_all_subgroup_completion import (
    _group_type,
    _json_digest,
    _matrix_payload,
    _one_based_orbits,
    _orbit_gram,
    _sha256,
    _subgroup_conjugates,
)
from scripts.certify_ss_fermat_char3_survivor_fixed_lattices import (
    _saturated_fixed_basis,
)
from scripts.inventory_alexeev_20_2_1_15_9_1_wall_symmetries import (
    _gap_subgroup_classes,
    _greedy_generators,
)
from scripts.search_apery_fermi_chamber_folds import (
    Permutation,
    Subgroup,
    _fixed_rank,
    _generated_subgroup,
    _orbits,
    _permutation_order,
)
from scripts.search_k3_high_rank_involution_folds import (
    _conjugacy_classes,
    _group_checks,
    _identity,
)
from scripts.verify_k3_high_rank_nef_chambers import build_gram, exact_inertia
from scripts.verify_k3_quotient_nef_slices import _coxeter_label, analyze_case


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/known_examples/k3_high_rank_nef_chambers.json"
QUOTIENT_SOURCE = ROOT / "data/known_examples/k3_quotient_nef_slices.json"
INVOLUTION_ARCHIVE = ROOT / "data/candidates/k3_high_rank_all_involution_folds.json"
F2_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h10_maximal_reflection_lineage.json"
LEDGER = ROOT / "results/higher_dimensional_polytope_lineage_ledger.md"
OUTPUT = ROOT / "data/enumerations/k3_rank18_19_all_subgroup_completion.json"

COXITER_ROOT = Path("/private/tmp/coxiter-current/CoxIter-master")
COXITER_ZIP = Path("/private/tmp/coxiter-current.zip")
DOUBLING_ARCHIVE = COXITER_ROOT / "graphs/testing_index2.tar.gz"

EXPECTED_HASHES = {
    "data/known_examples/k3_high_rank_nef_chambers.json": (
        "0dcc7ec9c10f3437cd39713cd0d6036dc4fa8ad0c978ffdad4c2b34748269120"
    ),
    "data/known_examples/k3_quotient_nef_slices.json": (
        "20dbe5535bc9832947ca3c6a4b7dde5afe60432a06913c383f0e63718ea3caea"
    ),
    "data/candidates/k3_high_rank_all_involution_folds.json": (
        "460754e9f576ce1b0323f3b36d43301510ebc8b065e74f1701e83c2a954c5000"
    ),
    "data/enumerations/jacobian_kummer_h10_maximal_reflection_lineage.json": (
        "3cb8b0e9120d13479f9b9e3a32cfb63922e03f41f2f0fc1dc501b0a9f6654c29"
    ),
}
COXITER_ZIP_SHA256 = "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
DOUBLING_SHA256 = "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"

EXPECTED = {
    18: {
        "wall_count": 19,
        "group_name": "C2",
        "group_order": 2,
        "element_histogram": {1: 1, 2: 1},
        "element_classes": 2,
        "subgroup_classes": 2,
    },
    19: {
        "wall_count": 24,
        "group_name": "S3",
        "group_order": 6,
        "element_histogram": {1: 1, 2: 3, 3: 2},
        "element_classes": 3,
        "subgroup_classes": 4,
    },
}
MIN_FIXED_RANK = 11

# These are precisely the positive IDs in the first table of the frozen local
# lineage ledger, including the rank-19 quotient-nef calibration added by this
# completion.  Dimension/facet count is already an exact obstruction to a
# square facet-Gram permutation outside the single H10/13 row.
FROZEN_LEDGER_CHAMBERS = (
    ("F3-C2-H10-196", 10, 196),
    ("F3-V4-H10-123", 10, 123),
    ("A15-C2-H10-34", 10, 34),
    ("SS2-2D-H10-47", 10, 47),
    ("K3-S3-C2-H10-13", 10, 13),
    ("F3-C2-H11-206", 11, 206),
    ("A20-C2b-H11-20", 11, 20),
    ("F3-V4-H11-328", 11, 328),
    ("SS2-3B-H11-22", 11, 22),
    ("D16-D5-H10-30", 10, 30),
    ("D16-D4-C7-H10-27", 10, 27),
    ("JK-H10-0007", 10, 84),
    ("K3-D4-H12", 12, 17),
    ("JK-H12", 12, 172),
    ("A20-C2a-H12-19", 12, 19),
    ("D16-D5-H14", 14, 60),
    ("K3-D6-H14", 14, 17),
    ("X3-H19", 19, 36),
    ("BOR-D4-H21-210", 21, 210),
)


def _ambient_reconstruction(
    wall_gram: sp.Matrix, basis_indices_one_based: Sequence[int]
) -> tuple[sp.Matrix, sp.Matrix, tuple[int, ...]]:
    basis_indices = tuple(int(label) - 1 for label in basis_indices_one_based)
    ambient_gram = wall_gram.extract(basis_indices, basis_indices)
    wall_rows = wall_gram[:, list(basis_indices)] * ambient_gram.inv()
    if not (
        ambient_gram.det() != 0
        and wall_rows[list(basis_indices), :] == sp.eye(len(basis_indices))
        and wall_rows * ambient_gram * wall_rows.T == wall_gram
        and all(value.is_Integer for value in wall_rows)
    ):
        raise ArithmeticError("published curve basis did not reconstruct the wall stock")
    return ambient_gram, wall_rows, basis_indices


def _lift_actions(
    group: tuple[Permutation, ...],
    wall_rows: sp.Matrix,
    ambient_gram: sp.Matrix,
    basis_indices: tuple[int, ...],
) -> dict[Permutation, sp.Matrix]:
    actions: dict[Permutation, sp.Matrix] = {}
    for permutation in group:
        action = wall_rows[[permutation[index] for index in basis_indices], :]
        expected = wall_rows[[permutation[index] for index in range(wall_rows.rows)], :]
        if not (
            wall_rows * action == expected
            and action * ambient_gram * action.T == ambient_gram
            and all(value.is_Integer for value in action)
            and abs(int(action.det())) == 1
        ):
            raise ArithmeticError("a diagram action failed its integral NS lift")
        actions[permutation] = action
    return actions


def _relation_matrix(facet_gram: sp.Matrix) -> list[list[str]]:
    matrix: list[list[str]] = []
    for left in range(facet_gram.rows):
        row = []
        for right in range(facet_gram.cols):
            if left == right:
                row.append("1")
            else:
                label = _coxeter_label(
                    sp.Rational(facet_gram[left, left]),
                    sp.Rational(facet_gram[right, right]),
                    sp.Rational(facet_gram[left, right]),
                )
                row.append("0" if label == "infinity" else label)
        matrix.append(row)
    return matrix


def _weighted_gram_permutation(
    source: sp.Matrix, target: sp.Matrix
) -> list[int] | None:
    if source.shape != target.shape:
        return None

    def graph(matrix: sp.Matrix) -> nx.Graph:
        answer = nx.Graph()
        for index in range(matrix.rows):
            answer.add_node(index, norm=str(sp.Rational(matrix[index, index])))
        for left in range(matrix.rows):
            for right in range(left + 1, matrix.cols):
                answer.add_edge(
                    left,
                    right,
                    pairing=str(sp.Rational(matrix[left, right])),
                )
        return answer

    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph(source),
        graph(target),
        node_match=nx.algorithms.isomorphism.categorical_node_match("norm", None),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match(
            "pairing", None
        ),
    )
    mapping = next(matcher.isomorphisms_iter(), None)
    if mapping is None:
        return None
    permutation = [int(mapping[index]) for index in range(source.rows)]
    if any(
        source[left, right] != target[permutation[left], permutation[right]]
        for left in range(source.rows)
        for right in range(source.cols)
    ):
        raise ArithmeticError("weighted-Gram matcher returned a false permutation")
    return permutation


def _parse_coxiter_text(text: str) -> tuple[int, int, list[list[str]]]:
    lines = [
        content
        for line in text.splitlines()
        if (content := line.split("#", 1)[0].strip())
    ]
    facets, dimension = map(int, lines[0].split()[:2])
    matrix = [
        ["1" if left == right else "2" for right in range(facets)]
        for left in range(facets)
    ]
    vertex_line = next(
        (line for line in lines[1:] if line.startswith("vertices labels:")), None
    )
    labels = (
        [str(index + 1) for index in range(facets)]
        if vertex_line is None
        else vertex_line.split(":", 1)[1].split()
    )
    if len(labels) != facets or len(set(labels)) != facets:
        raise ArithmeticError("invalid CoxIter labels")
    index = {label: position for position, label in enumerate(labels)}
    for line in lines[1:]:
        if line.startswith("vertices labels:"):
            continue
        tokens = line.split()
        if len(tokens) != 3:
            raise ArithmeticError("unsupported CoxIter record line")
        left, right = index[tokens[0]], index[tokens[1]]
        matrix[left][right] = matrix[right][left] = tokens[2]
    return facets, dimension, matrix


def _coxiter_corpus_audit(
    candidate_relations: list[list[str]],
    coxiter_root: Path,
    doubling_archive: Path,
) -> dict[str, Any]:
    ordinary_paths = sorted(coxiter_root.rglob("*.coxiter"))
    ordinary_candidates = []
    for path in ordinary_paths:
        header = next(
            line.split("#", 1)[0].strip()
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
            if line.split("#", 1)[0].strip()
        )
        facets, dimension = map(int, header.split()[:2])
        if (dimension, facets) != (10, 13):
            continue
        parsed_facets, parsed_dimension, relations = _parse_coxiter(path)
        if (parsed_dimension, parsed_facets) != (10, 13):
            raise ArithmeticError("ordinary CoxIter header changed during parsing")
        ordinary_candidates.append(
            {
                "path": path.relative_to(coxiter_root).as_posix(),
                "sha256": _sha256(path),
                "exact_relation_graph_permutation": _relation_permutation(
                    candidate_relations, relations
                ),
            }
        )

    doubling_candidates = []
    doubling_count = 0
    with tarfile.open(doubling_archive, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            doubling_count += 1
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError("could not read a CoxIter doubling member")
            raw = handle.read()
            facets, dimension, relations = _parse_coxiter_text(
                raw.decode("utf-8", errors="replace")
            )
            if (dimension, facets) != (10, 13):
                continue
            doubling_candidates.append(
                {
                    "path": member.name,
                    "sha256": hashlib.sha256(raw).hexdigest(),
                    "exact_relation_graph_permutation": _relation_permutation(
                        candidate_relations, relations
                    ),
                }
            )
    return {
        "commit": COXITER_COMMIT,
        "ordinary_file_count": len(ordinary_paths),
        "ordinary_H10_13_records": ordinary_candidates,
        "ordinary_exact_matches": [
            row for row in ordinary_candidates if row["exact_relation_graph_permutation"]
        ],
        "doubling_file_count": doubling_count,
        "doubling_H10_13_records": doubling_candidates,
        "doubling_exact_matches": [
            row for row in doubling_candidates if row["exact_relation_graph_permutation"]
        ],
    }


def _element_class_records(
    group: tuple[Permutation, ...],
    actions: dict[Permutation, sp.Matrix],
    ambient_rank: int,
    wall_count: int,
) -> list[dict[str, Any]]:
    records = []
    for class_index, conjugacy_class in enumerate(_conjugacy_classes(group, group)):
        representative = conjugacy_class[0]
        cyclic = _generated_subgroup((representative,))
        fixed_rank = _fixed_rank((actions[representative],), ambient_rank)
        records.append(
            {
                "class_index": class_index,
                "element_order": _permutation_order(representative),
                "conjugacy_class_size": len(conjugacy_class),
                "representative_images_zero_based": list(representative),
                "fixed_rank": fixed_rank,
                "hyperbolic_dimension": fixed_rank - 1,
                "cyclic_wall_orbits_one_based": _one_based_orbits(
                    _orbits(cyclic, wall_count)
                ),
            }
        )
    return records


def _subgroup_class_records(
    group: tuple[Permutation, ...],
    actions: dict[Permutation, sp.Matrix],
    ambient_rank: int,
    wall_gram: sp.Matrix,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows = _gap_subgroup_classes(group, _greedy_generators(group))
    records = []
    internal_high = []
    for class_index, row in enumerate(rows):
        subgroup: Subgroup = row["representative_subgroup"]
        generators: tuple[Permutation, ...] = row["representative_generators"]
        fixed_rank = (
            _fixed_rank(tuple(actions[element] for element in generators), ambient_rank)
            if generators
            else ambient_rank
        )
        orbits = _orbits(subgroup, wall_gram.rows)
        orbit_gram = _orbit_gram(wall_gram, orbits)
        conjugates = _subgroup_conjugates(group, subgroup)
        if not (
            orbit_gram.rank() == fixed_rank
            and len(conjugates) == int(row["conjugacy_class_size"])
        ):
            raise ArithmeticError("subgroup fixed-rank/conjugacy audit failed")
        proper_high = len(subgroup) > 1 and fixed_rank >= MIN_FIXED_RANK
        public = {
            "class_index": class_index,
            "group_type": _group_type(subgroup),
            "order": len(subgroup),
            "conjugacy_class_size": int(row["conjugacy_class_size"]),
            "cyclic": any(
                _permutation_order(element) == len(subgroup) for element in subgroup
            ),
            "fixed_rank": fixed_rank,
            "hyperbolic_dimension": fixed_rank - 1,
            "proper_fixed_rank_at_least_11": proper_high,
            "element_order_histogram": {
                str(order): count
                for order, count in sorted(
                    Counter(_permutation_order(element) for element in subgroup).items()
                )
            },
            "representative_generators_zero_based": [
                list(element) for element in generators
            ],
            "representative_elements_zero_based": [
                list(element) for element in sorted(subgroup)
            ],
            "representative_wall_orbits_one_based": _one_based_orbits(orbits),
            "orbit_gram_rank": orbit_gram.rank(),
            "orbit_gram_inertia": list(exact_inertia(orbit_gram)),
            "orbit_normal_square_histogram": {
                str(square): count
                for square, count in sorted(
                    Counter(
                        int(orbit_gram[index, index])
                        for index in range(orbit_gram.rows)
                    ).items()
                )
            },
        }
        records.append(public)
        if proper_high:
            internal_high.append(
                {
                    **public,
                    "_subgroup": subgroup,
                    "_conjugates": conjugates,
                    "_orbit_gram": orbit_gram,
                }
            )
    return records, internal_high


def _strict_facet_audits(audit: dict[str, Any]) -> list[dict[str, Any]]:
    orbit_gram = sp.Matrix(audit["orbit_gram"])
    orbit_cycles = [tuple(orbit) for orbit in audit["orbit_cycles"]]
    basis_indices = [
        orbit_cycles.index(tuple(orbit)) for orbit in audit["vertex_basis_orbits"]
    ]
    pairing = orbit_gram.extract(range(orbit_gram.rows), basis_indices)
    vertices = audit["vertices"]
    records = []
    for facet in [index for index in range(orbit_gram.rows) if orbit_gram[index, index] < 0]:
        equality_vectors = [
            sp.Matrix(vertex["basis_coordinates"])
            for vertex in vertices
            if facet + 1 in vertex["active_facets"]
        ]
        equality_rank = (
            sp.Matrix.hstack(*equality_vectors).rank() if equality_vectors else 0
        )
        witness = next(
            (
                vertex
                for vertex in vertices
                if (pairing.extract([facet], range(pairing.cols))
                    * sp.Matrix(vertex["basis_coordinates"]))[0]
                > 0
            ),
            None,
        )
        if equality_rank != audit["hyperbolic_dimension"] or witness is None:
            raise ArithmeticError("a retained wall lacks a strict facet witness")
        witness_pairing = (
            pairing.extract([facet], range(pairing.cols))
            * sp.Matrix(witness["basis_coordinates"])
        )[0]
        records.append(
            {
                "facet_orbit_index": facet + 1,
                "source_wall_orbit": list(orbit_cycles[facet]),
                "equality_ray_span_rank": equality_rank,
                "strict_witness": {
                    "type": witness["type"],
                    "basis_coordinates": witness["basis_coordinates"],
                    "active_facets": witness["active_facets"],
                    "strict_pairing": str(sp.Rational(witness_pairing)),
                },
                "irredundant_with_strict_relative_interior": True,
            }
        )
    return records


def _fixed_lattice_and_mirrors(
    action: sp.Matrix,
    ambient_gram: sp.Matrix,
    wall_rows: sp.Matrix,
    cycles: Sequence[Sequence[int]],
) -> dict[str, Any]:
    basis, kernel = _saturated_fixed_basis((action,), ambient_gram.rows)
    fixed_gram = basis * ambient_gram * basis.T
    fixed_smith = smith_normal_form(fixed_gram, domain=ZZ)
    mirror_rows = []
    facet_coordinates = []
    for orbit_index, orbit in enumerate(cycles):
        normal = sum(
            (wall_rows[int(label) - 1, :] for label in orbit),
            sp.zeros(1, ambient_gram.rows),
        )
        solution = sp.linsolve((basis.T, normal.T))
        if solution is sp.EmptySet:
            raise ArithmeticError("fixed orbit normal is absent from the fixed lattice")
        coordinates = sp.Matrix(tuple(next(iter(solution)))).T
        if not (
            coordinates * basis == normal
            and all(value.is_Integer for value in coordinates)
        ):
            raise ArithmeticError("orbit normal is not integral in the fixed lattice")
        square = int((coordinates * fixed_gram * coordinates.T)[0])
        reflection = (
            sp.eye(fixed_gram.rows)
            - sp.Rational(2, square)
            * fixed_gram
            * coordinates.T
            * coordinates
        )
        record = {
            "orbit_index": orbit_index + 1,
            "source_walls": [int(label) for label in orbit],
            "fixed_basis_coordinates": [int(value) for value in coordinates],
            "square": square,
            "reflection_integral": all(value.is_Integer for value in reflection),
            "role": "facet" if square < 0 else "redundant_timelike",
        }
        mirror_rows.append(record)
        if square < 0:
            facet_coordinates.append(coordinates)
    facet_matrix = sp.Matrix.vstack(*facet_coordinates)
    facet_smith = smith_normal_form(facet_matrix, domain=ZZ)
    facet_span_diagonal = [
        abs(int(facet_smith[index, index]))
        for index in range(fixed_gram.rows)
    ]
    if not (
        fixed_gram.rows == 11
        and exact_inertia(fixed_gram) == (1, 10, 0)
        and int(fixed_gram.det()) == 512
        and all(row["reflection_integral"] for row in mirror_rows if row["role"] == "facet")
        and facet_matrix.rank() == 11
        and facet_span_diagonal == [1] * 11
    ):
        raise ArithmeticError("rank-19 fixed lattice or mirror audit failed")
    return {
        "basis_rows_in_published_NS_basis": _matrix_payload(basis),
        "gram": _matrix_payload(fixed_gram),
        "rank": fixed_gram.rows,
        "inertia": list(exact_inertia(fixed_gram)),
        "determinant": int(fixed_gram.det()),
        "smith_diagonal_absolute": [
            abs(int(fixed_smith[index, index]))
            for index in range(fixed_gram.rows)
        ],
        "even": all(int(fixed_gram[index, index]) % 2 == 0 for index in range(fixed_gram.rows)),
        "smith_kernel_audit": kernel,
        "mirror_records": mirror_rows,
        "facet_mirrors_generate_saturated_fixed_lattice": facet_span_diagonal
        == [1] * 11,
        "facet_mirror_coordinate_smith_diagonal": facet_span_diagonal,
    }


def _ledger_comparison(ledger_text: str) -> list[dict[str, Any]]:
    rows = []
    for chamber_id, dimension, facets in FROZEN_LEDGER_CHAMBERS:
        if f"`{chamber_id}`" not in ledger_text:
            raise ArithmeticError(f"frozen ledger ID is missing: {chamber_id}")
        if dimension != 10:
            decision = "not_isomorphic_dimension_differs"
        elif facets != 13:
            decision = "not_isomorphic_facet_count_differs"
        else:
            decision = "exact_match_to_the_same_frozen_rank19_quotient_nef_chamber"
        rows.append(
            {
                "id": chamber_id,
                "hyperbolic_dimension": dimension,
                "facet_count": facets,
                "decision": decision,
            }
        )
    return rows


def build(
    coxiter_root: Path = COXITER_ROOT,
    coxiter_zip: Path = COXITER_ZIP,
    doubling_archive: Path = DOUBLING_ARCHIVE,
) -> dict[str, Any]:
    local_paths = (SOURCE, QUOTIENT_SOURCE, INVOLUTION_ARCHIVE, F2_ARCHIVE)
    for path in local_paths:
        relative = str(path.relative_to(ROOT))
        if _sha256(path) != EXPECTED_HASHES[relative]:
            raise ArithmeticError(f"unexpected source SHA-256: {relative}")
    if _sha256(coxiter_zip) != COXITER_ZIP_SHA256:
        raise ArithmeticError("pinned CoxIter zip hash changed")
    if _sha256(doubling_archive) != DOUBLING_SHA256:
        raise ArithmeticError("pinned CoxIter doubling archive hash changed")

    source = json.loads(SOURCE.read_text())
    quotient = json.loads(QUOTIENT_SOURCE.read_text())
    involution_archive = json.loads(INVOLUTION_ARCHIVE.read_text())
    f2_archive = json.loads(F2_ARCHIVE.read_text())
    source_cases = {int(case["rho"]): case for case in source["cases"]}
    quotient_cases = {
        int(case["ambient_rho"]): case for case in quotient["cases"]
    }

    cases = []
    all_high = []
    cached_case_data: dict[int, dict[str, Any]] = {}
    for rho in (18, 19):
        case = source_cases[rho]
        expected = EXPECTED[rho]
        wall_gram = build_gram(case)
        ambient_gram, wall_rows, basis_indices = _ambient_reconstruction(
            wall_gram, case["basis_indices"]
        )
        group = exact_weighted_gram_automorphisms(wall_gram)
        actions = _lift_actions(group, wall_rows, ambient_gram, basis_indices)
        order_histogram = Counter(_permutation_order(element) for element in group)
        if not (
            wall_gram.rows == expected["wall_count"]
            and wall_gram.rank() == rho
            and exact_inertia(wall_gram) == (1, rho - 1, wall_gram.rows - rho)
            and int(ambient_gram.det()) == int(case["expected_basis_determinant"])
            and exact_inertia(ambient_gram) == (1, rho - 1, 0)
            and len(group) == expected["group_order"]
            and dict(order_histogram) == expected["element_histogram"]
            and all(_group_checks(group).values())
        ):
            raise ArithmeticError(f"rank-{rho} wall/group inventory changed")
        element_classes = _element_class_records(
            group, actions, rho, wall_gram.rows
        )
        subgroup_classes, high = _subgroup_class_records(
            group, actions, rho, wall_gram
        )
        if not (
            len(element_classes) == expected["element_classes"]
            and len(subgroup_classes) == expected["subgroup_classes"]
        ):
            raise ArithmeticError(f"rank-{rho} exact class count changed")
        for record in high:
            record["rho"] = rho
            all_high.append(record)
        cases.append(
            {
                "rho": rho,
                "ambient_lattice": case["lattice"],
                "wall_count": wall_gram.rows,
                "ambient_basis_indices_one_based": case["basis_indices"],
                "ambient_basis_gram": _matrix_payload(ambient_gram),
                "ambient_basis_determinant": int(ambient_gram.det()),
                "weighted_diagram_group": {
                    "name": expected["group_name"],
                    "order": len(group),
                    "element_order_histogram": {
                        str(order): count
                        for order, count in sorted(order_histogram.items())
                    },
                    "checks": _group_checks(group),
                    "all_actions_integral_unimodular_Lorentz": True,
                    "generators_zero_based": [
                        list(element) for element in _greedy_generators(group)
                    ],
                    "elements_zero_based": [list(element) for element in group],
                },
                "element_conjugacy_classes": element_classes,
                "subgroup_conjugacy_classes": subgroup_classes,
                "summary": {
                    "element_conjugacy_class_count": len(element_classes),
                    "subgroup_conjugacy_class_count": len(subgroup_classes),
                    "proper_fixed_rank_at_least_11_count": len(high),
                },
            }
        )
        cached_case_data[rho] = {
            "wall_gram": wall_gram,
            "ambient_gram": ambient_gram,
            "wall_rows": wall_rows,
            "basis_indices": basis_indices,
            "group": group,
            "actions": actions,
        }

    if not (
        len(all_high) == 1
        and all_high[0]["rho"] == 19
        and all_high[0]["group_type"] == "C2"
        and all_high[0]["fixed_rank"] == 11
        and all_high[0]["conjugacy_class_size"] == 3
    ):
        raise ArithmeticError("the rank-18/19 high subgroup queue changed")

    high = all_high[0]
    case19 = cached_case_data[19]
    published_cycles = quotient_cases[19]["involution_cycles"]
    published_permutation = _permutation_from_cycles(24, published_cycles)
    if not any(
        published_permutation in conjugate for conjugate in high["_conjugates"]
    ):
        raise ArithmeticError("published transposition left the unique C2 class")
    direct_audit = analyze_case(quotient_cases[19], source_cases[19])
    if not direct_audit["all_checks_pass"]:
        raise ArithmeticError("direct rank-19 geometry replay failed")
    strict_facets = _strict_facet_audits(direct_audit)
    orbit_gram = sp.Matrix(direct_audit["orbit_gram"])
    facet_indices = [
        index for index in range(orbit_gram.rows) if orbit_gram[index, index] < 0
    ]
    facet_gram = orbit_gram.extract(facet_indices, facet_indices)
    relation_matrix = _relation_matrix(facet_gram)

    fixed_lattice = _fixed_lattice_and_mirrors(
        case19["actions"][published_permutation],
        case19["ambient_gram"],
        case19["wall_rows"],
        published_cycles,
    )

    archived_case19 = next(
        case
        for case in involution_archive["cases"]
        if int(case["ambient_rho"]) == 19
    )
    archived_published = next(
        fold for fold in archived_case19["folds"] if fold["equals_archived_published_action"]
    )
    archived_orbit_gram = sp.Matrix(archived_published["orbit_gram"])
    archived_facets = [
        index
        for index in range(archived_orbit_gram.rows)
        if archived_orbit_gram[index, index] < 0
    ]
    archived_facet_gram = archived_orbit_gram.extract(archived_facets, archived_facets)
    archive_permutation = _weighted_gram_permutation(facet_gram, archived_facet_gram)
    if archive_permutation is None:
        raise ArithmeticError("the survivor no longer matches its frozen facet Gram")

    conjugate_fold_matches = []
    for fold in archived_case19["folds"]:
        conjugate_orbit_gram = sp.Matrix(fold["orbit_gram"])
        conjugate_facets = [
            index
            for index in range(conjugate_orbit_gram.rows)
            if conjugate_orbit_gram[index, index] < 0
        ]
        conjugate_facet_gram = conjugate_orbit_gram.extract(
            conjugate_facets, conjugate_facets
        )
        permutation = _weighted_gram_permutation(facet_gram, conjugate_facet_gram)
        if permutation is None:
            raise ArithmeticError("an S3-conjugate transposition has a different facet Gram")
        conjugate_fold_matches.append(
            {
                "fold_id": fold["fold_id"],
                "source_facet_to_conjugate_facet_zero_based": permutation,
                "entrywise_exact_facet_gram_match": True,
                "surface_conjugator_realized": True,
            }
        )

    f2_relations = f2_archive["vinberg_maximal_integral_reflection_chamber"][
        "normalized_coxeter_relation_matrix"
    ]
    f2_permutation = _relation_permutation(relation_matrix, f2_relations)
    if f2_permutation is not None:
        raise ArithmeticError("rank-19 H10/13 unexpectedly became the f2 chamber")

    corpus = _coxiter_corpus_audit(
        relation_matrix, coxiter_root, doubling_archive
    )
    if not (
        corpus["ordinary_file_count"] == 1214
        and len(corpus["ordinary_H10_13_records"]) == 4
        and corpus["ordinary_exact_matches"] == []
        and corpus["doubling_file_count"] == 2959
        and len(corpus["doubling_H10_13_records"]) == 4
        and corpus["doubling_exact_matches"] == []
    ):
        raise ArithmeticError("bounded CoxIter H10/13 comparison changed")

    ledger_rows = _ledger_comparison(LEDGER.read_text())
    if sum(
        row["decision"].startswith("exact_match") for row in ledger_rows
    ) != 1:
        raise ArithmeticError("frozen ledger deduplication is not unique")

    vertices = direct_audit["vertices"]
    finite = [vertex for vertex in vertices if vertex["type"] == "finite"]
    ideal = [vertex for vertex in vertices if vertex["type"] == "ideal"]
    hyperideal = [vertex for vertex in vertices if vertex["type"] == "hyperideal"]
    literal = direct_audit["literal_coxiter_prop_1_10"]
    geometry_checks = {
        "unique_proper_fixed_rank_at_least_11_class": len(all_high) == 1,
        "all_13_facets_have_strict_witnesses": (
            len(strict_facets) == 13
            and all(
                row["irredundant_with_strict_relative_interior"]
                for row in strict_facets
            )
        ),
        "exact_Coxeter_angles": direct_audit["coxeter_edge_counts"]
        == {"3": 9, "4": 2, "infinity": 2},
        "all_folded_reflections_integral": all(
            direct_audit["fold_reflection_checks"].values()
        )
        and all(
            row["reflection_integral"]
            for row in fixed_lattice["mirror_records"]
            if row["role"] == "facet"
        ),
        "exact_cone_has_34_finite_4_ideal_0_hyperideal_rays": (
            len(finite), len(ideal), len(hyperideal)
        )
        == (34, 4, 0),
        "literal_208_edges_each_have_two_endpoints": (
            literal["spherical_rank_n_minus_one_count"] == 208
            and literal["extension_histogram"] == {"2": 208}
            and literal["bad_edges"] == []
        ),
        "literal_Vinberg_finite_volume_noncompact": (
            literal["finite_covolume"] is True
            and direct_audit["finite_volume"] is True
            and direct_audit["compact"] is False
        ),
        "direct_geometry_all_checks_pass": direct_audit["all_checks_pass"] is True,
        "exact_match_to_frozen_published_facet_Gram": archive_permutation
        is not None,
        "all_three_transposition_folds_are_exact_Gram_isomorphic": len(
            conjugate_fold_matches
        )
        == 3,
        "distinct_from_f2_H10_13_relation_graph": f2_permutation is None,
        "no_match_in_eight_pinned_CoxIter_H10_13_records": (
            corpus["ordinary_exact_matches"] == []
            and corpus["doubling_exact_matches"] == []
        ),
    }
    if not all(geometry_checks.values()):
        raise ArithmeticError("rank-19 survivor geometry/lineage check failed")

    public_high = {
        key: value for key, value in high.items() if not key.startswith("_")
    }
    payload: dict[str, Any] = {
        "schema_version": "1.0.0",
        "status_date": "2026-08-06",
        "artifact_id": "k3_rank18_19_all_subgroup_completion",
        "certificate_kind": (
            "exhaustive_exact_subgroup_fixed_section_geometry_and_bounded_lineage"
        ),
        "source_hashes": EXPECTED_HASHES,
        "external_snapshot": {
            "coxiter_commit": COXITER_COMMIT,
            "coxiter_zip_path": str(coxiter_zip),
            "coxiter_zip_sha256": COXITER_ZIP_SHA256,
            "doubling_archive_path": str(doubling_archive),
            "doubling_archive_sha256": DOUBLING_SHA256,
        },
        "scope": {
            "ambient_picard_ranks": [18, 19],
            "complete_wall_counts": {"18": 19, "19": 24},
            "group_denominator": "all subgroups modulo conjugacy in the full exact weighted wall-Gram group",
            "minimum_proper_fixed_rank_for_geometry": MIN_FIXED_RANK,
            "surface_realization_is_a_separate_gate": True,
            "bounded_lineage": (
                "every positive in the frozen local ledger plus every H10/13 item in the pinned ordinary and index-two-doubling CoxIter snapshots"
            ),
        },
        "cases": cases,
        "proper_fixed_rank_at_least_11_queue": [public_high],
        "survivor_geometry": {
            "id": "K3-S3-C2-H10-13",
            "source": "rank-19 K3 nef chamber, unique C2 subgroup class in S3",
            "surface_realization": {
                "status": "published_surface_automorphism_and_quotient_nef_cone",
                "reason": (
                    "Kondo realizes the full rank-19 S3 action; all three transpositions are surface-conjugate."
                ),
                "new_realization": False,
            },
            "fixed_lattice": fixed_lattice,
            "wall_orbits_one_based": published_cycles,
            "redundant_timelike_orbits": direct_audit["redundant_timelike_orbits"],
            "facet_gram": _matrix_payload(facet_gram),
            "facet_gram_sha256": _json_digest(_matrix_payload(facet_gram)),
            "normalized_coxeter_relation_matrix": relation_matrix,
            "coxeter_edge_counts": direct_audit["coxeter_edge_counts"],
            "strict_facet_audits": strict_facets,
            "exact_cone": {
                "vertices": vertices,
                "finite_vertex_count": len(finite),
                "ideal_vertex_count": len(ideal),
                "hyperideal_vertex_count": len(hyperideal),
            },
            "optimized_vinberg": direct_audit["vinberg_volume_criterion"],
            "literal_coxiter_prop_1_10": literal,
            "finite_volume": direct_audit["finite_volume"],
            "compact": direct_audit["compact"],
            "checks": geometry_checks,
        },
        "exact_deduplication": {
            "frozen_published_rank19_fold": {
                "artifact": "data/candidates/k3_high_rank_all_involution_folds.json",
                "source_facet_to_archive_facet_zero_based": archive_permutation,
                "entrywise_exact_facet_gram_match": True,
            },
            "three_realized_S3_conjugates": conjugate_fold_matches,
            "frozen_positive_ledger": ledger_rows,
            "classical_f2_H10_13": {
                "artifact": "data/enumerations/jacobian_kummer_h10_maximal_reflection_lineage.json",
                "exact_relation_graph_permutation": f2_permutation,
                "isomorphic": False,
                "independent_face_obstruction": (
                    "candidate has 34 finite vertices and 4 cusps; f2 has 29 finite vertices and 3 cusps"
                ),
            },
            "pinned_coxiter_H10_13": corpus,
            "bounded_conclusion": (
                "The only exact frozen-ledger match is the already archived rank-19 quotient-nef chamber itself. The diagram is distinct from f2 and absent from all eight H10/13 records in the two pinned CoxIter snapshots; no worldwide priority inference is made."
            ),
        },
        "summary": {
            "weighted_group_orders": {"18": 2, "19": 6},
            "element_conjugacy_class_counts": {"18": 2, "19": 3},
            "subgroup_conjugacy_class_counts": {"18": 2, "19": 4},
            "total_subgroup_conjugacy_classes": 6,
            "proper_fixed_rank_at_least_11_classes": 1,
            "new_noninvolution_high_classes": 0,
            "certified_finite_volume_noncompact_classes": 1,
            "new_hyperbolic_Coxeter_polytopes": 0,
            "new_algebraic_realizations": 0,
            "closure_statement": (
                "The complete rank-18 and rank-19 K3 weighted-wall subgroup denominator is closed. C3 and S3 in rank 19 lie below H10, and the sole high C2 class is the already frozen realized quotient-nef H10/13 chamber."
            ),
        },
    }
    payload["all_checks_pass"] = all(geometry_checks.values())
    payload["certificate_payload_sha256"] = _json_digest(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coxiter-root", type=Path, default=COXITER_ROOT)
    parser.add_argument("--coxiter-zip", type=Path, default=COXITER_ZIP)
    parser.add_argument("--doubling-archive", type=Path, default=DOUBLING_ARCHIVE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build(args.coxiter_root, args.coxiter_zip, args.doubling_archive)
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text() != encoded:
            raise SystemExit(f"archive mismatch: {args.output}")
        print(f"archive verified: {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
        print(args.output)
    print(json.dumps(payload["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
