#!/usr/bin/env python3
"""Independent hostile audit of the 15-edge no-parallel Coble lane.

The audit deliberately does not import the production lane or its cone,
overlattice, Vinberg, or isomorphism helpers.  It reruns nauty directly,
uses exact rational congruence for inertia, enumerates the determinant-four
half-cosets explicitly, calls cddlib directly, and constructs marked
unimodular isometries to the archived Type-II full wall system.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any, Sequence

from flint import fmpz_mat
import networkx as nx
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "data/enumerations/coble_nodot_15edge_lane.json"
SOURCE_FULL = ROOT / "data/enumerations/coble_one_boundary_full_ns_h10.json"
SOURCE_PROMOTION = ROOT / "data/enumerations/coble_one_boundary_raw_k3_promotion.json"
REPORT = ROOT / "results/coble_nodot_15edge_lane.md"
OUTPUT = ROOT / "data/enumerations/coble_nodot_15edge_lane_hostile_audit.json"

PINNED = {
    TARGET: "c9e2da6bcc08736fbc9848d64ca3fd273c23d535394ef727b671f369e7b20dae",
    SOURCE_FULL: "1e90fcd5f9b1ba8538937b91d2f26ea401f6952e8484a597624939edaf2a873e",
    SOURCE_PROMOTION: "1c8b8699da8f8747fdad22d520743bb4640b41b8e988c16844ecfe5e27c1958c",
}
SCHEMA = "coxeter.coble_nodot_15edge_lane_hostile_audit.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _integer_rows(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.q != 1 for value in matrix):
        raise ArithmeticError("expected an integral matrix")
    return [[int(value) for value in matrix.row(index)] for index in range(matrix.rows)]


def _rational_rows(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [int(value) if value.q == 1 else str(value) for value in matrix.row(index)]
        for index in range(matrix.rows)
    ]


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _exact_inertia_by_congruence(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Compute inertia over Q by exact symmetric Schur complements."""

    active = sp.Matrix(matrix)
    positive = negative = 0
    while active.rows:
        diagonal = next(
            (index for index in range(active.rows) if active[index, index] != 0),
            None,
        )
        if diagonal is not None:
            if diagonal:
                active.row_swap(0, diagonal)
                active.col_swap(0, diagonal)
            pivot = active[0, 0]
            if pivot > 0:
                positive += 1
            else:
                negative += 1
            if active.rows == 1:
                active = sp.zeros(0, 0)
            else:
                column = active[1:, 0]
                active = active[1:, 1:] - column * column.T / pivot
            continue

        pair = next(
            (
                (left, right)
                for left in range(active.rows)
                for right in range(left)
                if active[left, right] != 0
            ),
            None,
        )
        if pair is None:
            return positive, negative, active.rows
        left, right = pair
        order = [left, right, *[index for index in range(active.rows) if index not in pair]]
        active = active.extract(order, order)
        block = active[:2, :2]
        positive += 1
        negative += 1
        if active.rows == 2:
            active = sp.zeros(0, 0)
        else:
            cross = active[2:, :2]
            active = active[2:, 2:] - cross * block.inv() * cross.T
    return positive, negative, 0


def _support_matrix(graph: nx.Graph) -> sp.Matrix:
    return sp.Matrix(
        [
            [
                2 if left == right else -1 if graph.has_edge(left, right) else 0
                for right in range(12)
            ]
            for left in range(12)
        ]
    )


def _enumerate_supports(geng: str) -> tuple[Counter[int], list[dict[str, Any]]]:
    process = subprocess.Popen(
        [geng, "-q", "-c", "-D3", "12", "15:15"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert process.stdout is not None
    rank_histogram: Counter[int] = Counter()
    survivors: list[dict[str, Any]] = []
    graph_count = 0
    for raw in process.stdout:
        graph6 = raw.strip().decode()
        graph = nx.from_graph6_bytes(raw.strip())
        matrix = _support_matrix(graph)
        exact_rank = fmpz_mat(_integer_rows(matrix)).rank()
        rank_histogram[exact_rank] += 1
        graph_count += 1
        if exact_rank != 10:
            continue
        inertia = _exact_inertia_by_congruence(matrix)
        if inertia == (9, 1, 2):
            survivors.append(
                {"graph6": graph6, "positive_cartan": matrix, "inertia": inertia}
            )
    stderr = process.stderr.read().decode() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code:
        raise RuntimeError(stderr or f"geng exited with status {return_code}")
    if graph_count != sum(rank_histogram.values()):
        raise ArithmeticError("the direct geng stream count is inconsistent")
    return rank_histogram, survivors


def _quotient(matrix: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix]:
    flint = fmpz_mat(_integer_rows(matrix))
    _, transform = flint.transpose().hnf(transform=True)
    basis = sp.Matrix(transform.transpose().tolist())
    if abs(int(basis.det())) != 1:
        raise ArithmeticError("the direct quotient basis is not unimodular")
    transformed = basis.T * matrix * basis
    if any(
        transformed[left, right]
        for left in range(12)
        for right in range(12)
        if left >= 10 or right >= 10
    ):
        raise ArithmeticError("the direct quotient failed to isolate the radical")
    quotient = transformed[:10, :10]
    roots = basis.inv()[:10, :]
    if roots.T * quotient * roots != matrix:
        raise ArithmeticError("the direct quotient failed Gram replay")
    return quotient, roots


def _extension_cosets(quotient: sp.Matrix) -> list[tuple[int, ...] | None]:
    determinant = abs(int(quotient.det()))
    if determinant == 1:
        return [None]
    if determinant != 4:
        raise ArithmeticError("the hostile audit encountered an unexpected determinant")
    result = []
    for bits in itertools.product((0, 1), repeat=10):
        if not any(bits):
            continue
        vector = sp.Matrix(bits)
        if any(value % 2 != 0 for value in quotient * vector):
            continue
        if (vector.T * quotient * vector)[0] % 8 != 0:
            continue
        result.append(bits)
    return result


def _primitive(values: Sequence[Fraction]) -> tuple[int, ...]:
    denominator = math.lcm(*(value.denominator for value in values))
    integers = [value.numerator * (denominator // value.denominator) for value in values]
    divisor = math.gcd(*(abs(value) for value in integers))
    return tuple(value // divisor for value in integers)


def _cdd_rays(inequality: sp.Matrix, cdd: str) -> list[tuple[int, ...]]:
    lines = [
        "H-representation",
        "begin",
        f"{inequality.rows} {inequality.cols + 1} rational",
        *[
            "0 " + " ".join(str(inequality[row, column]) for column in range(inequality.cols))
            for row in range(inequality.rows)
        ],
        "end",
    ]
    completed = subprocess.run(
        [cdd, "--rep"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        check=False,
        timeout=600,
    )
    if completed.returncode:
        raise RuntimeError(completed.stderr[-4000:] or completed.stdout[-4000:])
    output = completed.stdout.splitlines()
    begin = output.index("begin")
    count, width, kind = output[begin + 1].split()
    if kind != "rational" or int(width) != inequality.cols + 1:
        raise ArithmeticError("cdd returned an unexpected direct V-representation")
    rays = []
    for line in output[begin + 2 : begin + 2 + int(count)]:
        tokens = [Fraction(token) for token in line.split()]
        if tokens[0] != 0:
            raise ArithmeticError("the homogeneous cone acquired an affine vertex")
        rays.append(_primitive(tokens[1:]))
    if len(rays) != len(set(rays)):
        raise ArithmeticError("cdd returned duplicate oriented rays")
    return rays


def _cone_audit(
    gram: sp.Matrix, normals: sp.Matrix, cdd: str
) -> tuple[dict[str, Any], list[sp.Matrix]]:
    inequality = normals.T * gram
    rays = [sp.Matrix(ray) for ray in _cdd_rays(inequality, cdd)]
    types: Counter[str] = Counter()
    facet_vectors: list[list[sp.Matrix]] = [[] for _ in range(normals.cols)]
    hyperideal: list[tuple[int, ...]] = []
    ideal: list[sp.Matrix] = []
    for vector in rays:
        slacks = inequality * vector
        if any(value < 0 for value in slacks):
            raise ArithmeticError("cdd returned an infeasible direct ray")
        square = int((vector.T * gram * vector)[0])
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        types[kind] += 1
        if kind == "ideal":
            ideal.append(vector)
        if kind == "hyperideal":
            hyperideal.append(tuple(int(value) for value in vector))
        for index, slack in enumerate(slacks):
            if slack == 0:
                facet_vectors[index].append(vector)
    face_ranks = [
        sp.Matrix.hstack(*vectors).rank() if vectors else 0 for vectors in facet_vectors
    ]
    primitive = sorted(tuple(int(value) for value in vector) for vector in rays)
    return {
        "ray_count": len(rays),
        "finite": types["finite"],
        "ideal": types["ideal"],
        "hyperideal": types["hyperideal"],
        "all_facets_irredundant": all(rank == gram.rows - 1 for rank in face_ranks),
        "ray_sha256": _digest(primitive),
        "hyperideal_ray_sha256": _digest(sorted(hyperideal)),
    }, ideal


def _full_graph(matrix: sp.Matrix) -> nx.Graph:
    graph = nx.Graph()
    for index in range(matrix.rows):
        graph.add_node(index, norm=int(matrix[index, index]))
    for left in range(matrix.rows):
        for right in range(left):
            graph.add_edge(left, right, pairing=int(matrix[left, right]))
    return graph


def _canonical_graph6(graph6: str, labelg: str) -> str:
    completed = subprocess.run(
        [labelg, "-q"], input=graph6 + "\n", text=True, capture_output=True, check=True
    )
    return completed.stdout.strip()


def _full_type_ii_match(
    gram: sp.Matrix,
    walls: sp.Matrix,
    known_gram: sp.Matrix,
    known_walls: sp.Matrix,
    known_wall_gram: sp.Matrix,
) -> dict[str, Any]:
    wall_gram = walls.T * gram * walls
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        _full_graph(wall_gram),
        _full_graph(known_wall_gram),
        node_match=lambda left, right: left["norm"] == right["norm"],
        edge_match=lambda left, right: left["pairing"] == right["pairing"],
    )
    mapping = next(matcher.isomorphisms_iter(), None)
    if mapping is None:
        raise ArithmeticError("a finite branch has no full Type-II Gram isomorphism")
    permutation = [mapping[index] for index in range(walls.cols)]
    matched = known_walls[:, permutation]
    pivots = list(walls.rref()[1])
    isometry = matched[:, pivots] * walls[:, pivots].inv()
    smith = smith_normal_form(walls, domain=ZZ)
    span_diagonal = [abs(int(smith[index, index])) for index in range(gram.rows)]
    characteristic = [
        all(
            (gram[column, column] - (gram * walls[:, index])[column]) % 2 == 0
            for column in range(gram.rows)
        )
        for index in range(walls.cols)
    ]
    checks = {
        "permutation_replays_all_900_Gram_entries": wall_gram
        == known_wall_gram.extract(permutation, permutation),
        "marked_isometry_integral": all(value.q == 1 for value in isometry),
        "marked_isometry_unimodular": abs(int(isometry.det())) == 1,
        "marked_isometry_maps_all_walls": isometry * walls == matched,
        "marked_isometry_preserves_lattice_form": isometry.T * known_gram * isometry
        == gram,
        "walls_span_whole_I_1_10": span_diagonal == [1] * gram.rows,
        "unique_characteristic_wall_zero_based": [
            index for index, value in enumerate(characteristic) if value
        ]
        == [0],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"the marked Type-II return failed: {checks}")
    return {
        "facet_permutation_to_archived_Type_II_zero_based": permutation,
        "facet_permutation_sha256": _digest(permutation),
        "marked_isometry_determinant": int(isometry.det()),
        "marked_isometry_sha256": _digest(_rational_rows(isometry)),
        "wall_span_smith_diagonal": span_diagonal,
        "checks": checks,
    }


def build(geng: str | None = None, labelg: str | None = None, cdd: str | None = None) -> dict[str, Any]:
    for path, expected in PINNED.items():
        if _sha256(path) != expected:
            raise ArithmeticError(f"a hostile-audit input changed: {path}")
    geng = geng or shutil.which("geng")
    labelg = labelg or shutil.which("labelg")
    cdd = cdd or shutil.which("cddexec_gmp")
    if not geng or not labelg or not cdd:
        raise RuntimeError("geng, labelg, and cddexec_gmp are required")

    target = json.loads(TARGET.read_text(encoding="utf-8"))
    source = json.loads(SOURCE_FULL.read_text(encoding="utf-8"))["cases"]["II"]
    promotion = json.loads(SOURCE_PROMOTION.read_text(encoding="utf-8"))
    raw_30 = next(case for case in promotion["cases"] if case["raw_id"] == "raw_30")
    known_gram = sp.Matrix(source["lattice_gram"])
    known_walls = sp.Matrix(source["wall_vectors"]).T
    known_wall_gram = sp.Matrix(source["negative_wall_gram"])

    rank_histogram, survivors = _enumerate_supports(geng)
    known_root_indices = [
        index for index in range(known_wall_gram.rows) if known_wall_gram[index, index] == -2
    ]
    known_root_gram = known_wall_gram.extract(known_root_indices, known_root_indices)
    known_root_graph = nx.Graph()
    known_root_graph.add_nodes_from(range(12))
    for left in range(12):
        for right in range(left):
            if known_root_gram[left, right] != 0:
                known_root_graph.add_edge(left, right)
    known_root_canonical = _canonical_graph6(
        nx.to_graph6_bytes(known_root_graph, header=False).decode().strip(), labelg
    )

    branch_records = []
    branch_classes: Counter[str] = Counter()
    for class_index, survivor in enumerate(survivors):
        quotient, root_coordinates = _quotient(survivor["positive_cartan"])
        smith = smith_normal_form(survivor["positive_cartan"], domain=ZZ)
        smith_diagonal = [abs(int(smith[index, index])) for index in range(12)]
        canonical = _canonical_graph6(survivor["graph6"], labelg)
        extensions = _extension_cosets(quotient)
        extension_records = []
        for extension_index, residue in enumerate(extensions):
            transition = sp.eye(10)
            if residue is not None:
                vector = sp.Matrix(residue)
                pivot = next(index for index, value in enumerate(vector) if value)
                transition[:, pivot] = vector / 2
            gram = -(transition.T * quotient * transition)
            roots = transition.inv() * root_coordinates
            if not (
                all(value.q == 1 for value in gram)
                and all(value.q == 1 for value in roots)
                and gram.det() == -1
                and all(gram[index, index] % 2 == 0 for index in range(10))
                and _exact_inertia_by_congruence(gram) == (1, 9, 0)
                and roots.T * gram * roots == -survivor["positive_cartan"]
            ):
                raise ArithmeticError("an independently reconstructed E10 branch failed")
            cone, ideal = _cone_audit(gram, roots, cdd)
            classification = (
                "finite_volume_noncompact" if cone["hyperideal"] == 0 and cone["ideal"] else
                "compact" if cone["hyperideal"] == 0 else
                "infinite_volume_hyperideal"
            )
            branch_classes[classification] += 1
            full = None
            if classification == "finite_volume_noncompact":
                ambient = sp.diag(gram, sp.Matrix([[-1]]))
                boundary = sp.zeros(11, 1)
                boundary[10] = 1
                full_walls = sp.Matrix.hstack(
                    boundary,
                    *[roots[:, index].col_join(sp.Matrix([0])) for index in range(12)],
                    *[cusp.col_join(sp.Matrix([-1])) for cusp in ideal],
                )
                full_cone, _ = _cone_audit(ambient, full_walls, cdd)
                match = _full_type_ii_match(
                    ambient, full_walls, known_gram, known_walls, known_wall_gram
                )
                full = {"wall_count": full_walls.cols, "cone": full_cone, **match}
            extension_records.append(
                {
                    "extension_index_zero_based": extension_index,
                    "isotropic_half_coset_binary": list(residue) if residue is not None else None,
                    "root_cone": cone,
                    "classification": classification,
                    "full_Type_II_return": full,
                }
            )
        branch_records.append(
            {
                "class_index_zero_based": class_index,
                "graph6": survivor["graph6"],
                "canonical_graph6": canonical,
                "exact_inertia": list(survivor["inertia"]),
                "singular_gram_smith_diagonal": smith_diagonal,
                "radical_quotient_determinant": int(quotient.det()),
                "is_Type_II_root_graph": canonical == known_root_canonical,
                "even_unimodular_extension_count": len(extensions),
                "extension_branches": extension_records,
            }
        )

    report_text = REPORT.read_text(encoding="utf-8")
    scope_checks = {
        "target_explicitly_denies_strict_L4": target["claims"][
            "strict_L4_for_full_D_1_split_nodot"
        ]
        is False,
        "target_limits_denominator_to_lane": target["scope"][
            "not_the_full_D_1_split_nodot_denominator"
        ]
        is True,
        "return_inferred_only_after_full_Gram_match": target["claims"][
            "full_split_compiler_only_inferred_after_literal_known_gram_match"
        ]
        is True,
        "report_denies_full_denominator": "not the whole nodot denominator" in report_text,
        "report_denies_new_abstract_polytope": "not a new abstract polytope" in report_text,
        "report_calls_result_scoped": "scoped exact inverse" in report_text,
        "pinned_raw_30_promotion_passes": raw_30["all_checks_pass"] is True,
    }
    finite = [
        branch
        for record in branch_records
        for branch in record["extension_branches"]
        if branch["classification"] == "finite_volume_noncompact"
    ]
    rejected = [
        branch
        for record in branch_records
        for branch in record["extension_branches"]
        if branch["classification"] == "infinite_volume_hyperideal"
    ]
    checks = {
        "direct_geng_5741": sum(rank_histogram.values()) == 5741,
        "exact_rank_histogram": rank_histogram == {12: 5588, 11: 151, 10: 2},
        "exact_two_Lorentzian_corank_two_Grams": len(survivors) == 2,
        "expected_graph6_pair": [record["graph6"] for record in branch_records]
        == ["K??E@_KK_YWC", "K?AA@_ghECHG"],
        "determinants_and_extension_counts": [
            (record["radical_quotient_determinant"], record["even_unimodular_extension_count"])
            for record in branch_records
        ]
        == [(-4, 2), (-1, 1)],
        "two_finite_Type_II_returns": len(finite) == 2
        and all(branch["full_Type_II_return"] is not None for branch in finite),
        "finite_root_cones_38_17_0": all(
            (branch["root_cone"]["finite"], branch["root_cone"]["ideal"], branch["root_cone"]["hyperideal"])
            == (38, 17, 0)
            for branch in finite
        ),
        "full_cones_1333_207_0": all(
            (
                branch["full_Type_II_return"]["cone"]["finite"],
                branch["full_Type_II_return"]["cone"]["ideal"],
                branch["full_Type_II_return"]["cone"]["hyperideal"],
            )
            == (1333, 207, 0)
            for branch in finite
        ),
        "determinant_one_has_24_hyperideal_rays": len(rejected) == 1
        and rejected[0]["root_cone"]["hyperideal"] == 24,
        "volume_ledger_matches_target": dict(sorted(branch_classes.items()))
        == target["volume_ledger"],
        "scope_and_novelty_conservative": all(scope_checks.values()),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"hostile 15-edge audit failed: {checks}")
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "method_independence": {
            "imports_production_lane": False,
            "imports_production_helpers": False,
            "support_enumerator": "direct nauty geng",
            "inertia": "exact rational congruence, no numerical spectral prefilter",
            "overclosures": "direct binary half-coset enumeration",
            "cones": "direct cddexec_gmp H-to-V conversion and independent parser",
            "full_return": "explicit facet permutation and unimodular marked isometry to frozen Type II",
        },
        "input_sha256": {str(path.relative_to(ROOT)): value for path, value in PINNED.items()},
        "direct_support_rank_histogram": {
            str(rank): count for rank, count in sorted(rank_histogram.items())
        },
        "known_Type_II_root_canonical_graph6": known_root_canonical,
        "class_records": branch_records,
        "volume_ledger": dict(sorted(branch_classes.items())),
        "scope_audit": scope_checks,
        "checks": checks,
        "conclusion": (
            "PASS: the 15-edge maximum-degree-three no-parallel lane has exactly "
            "two Lorentzian corank-two Grams and three E10 branches; the two "
            "determinant-four branches return the same archived Type-II H10/30 "
            "chamber, while the determinant-one branch has 24 hyperideal rays. "
            "The result is a scoped inverse calibration, not strict L4 or a new polytope."
        ),
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geng", default=None)
    parser.add_argument("--labelg", default=None)
    parser.add_argument("--cdd", default=None)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build(args.geng, args.labelg, args.cdd)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit("hostile audit archive mismatch")
        print(hashlib.sha256(rendered.encode()).hexdigest())
    elif args.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
