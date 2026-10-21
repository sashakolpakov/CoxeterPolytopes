"""Audit full weighted-wall symmetries of every published disc16 chamber.

Shimada archives the lattice/chamber stabilizer ``Aut(D)``.  The fishing
factory deliberately permits a broader input: any exact finite symmetry of
the complete oriented wall configuration, whether or not it is induced by a
surface automorphism.  This script therefore enumerates the automorphism
group of each exact weighted wall Gram graph in all 4,539 published chambers,
reconstructs every graph symmetry as a rational Lorentzian isometry, compares
it with ``Aut(D)``, and inventories every subgroup-fixed section above H9.

This remains only a rank/symmetry inventory.  Coxeter angles, irredundancy,
and volume are separate gates.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp
import networkx as nx

from scripts.certify_shimada_disc16_d5_full_subgroup_fishing import (
    _subgroup_conjugacy_classes,
)
from scripts.search_apery_fermi_chamber_folds import (
    Permutation,
    Subgroup,
    _all_subgroups,
    _compose,
    _fixed_rank,
    _identity,
    _permutation_order,
    _rational_text,
)
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas, _row_product


ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path("/private/tmp/threediscs/disc16.txt")
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_weighted_wall_symmetry_inventory.json"
)
ARCHIVED_INVENTORY = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_atlas_fixed_section_inventory.json"
)
EXPECTED_SOURCE_SHA256 = (
    "a19b797d6c314bc85e1ac4991fe3b190fe29a18a4746103bc6fa9b05e26a7db6"
)
AMBIENT_RANK = 20
MIN_FIXED_RANK = 11


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[str]]:
    return [
        [_rational_text(sp.Rational(entry)) for entry in row]
        for row in matrix.tolist()
    ]


def _weighted_wall_graph_automorphisms_fast(
    integral_scaled_dual_gram: sp.Matrix,
    walls: tuple[tuple[int, ...], ...],
) -> tuple[Permutation, ...]:
    """Enumerate the same exact group with complete-graph colour refinement.

    The shared generic helper performs a fresh 1-by-20 matrix product for
    every pair and gives VF2 only the norm colours.  Across 4,539 chambers
    that is needlessly expensive.  Here the full wall Gram is multiplied
    once, and exact one-dimensional Weisfeiler--Leman refinement includes
    every zero and nonzero pairing before VF2 starts.  Singleton refinement
    cells prove triviality immediately; otherwise VF2 still enumerates the
    complete colour-preserving graph automorphism group.
    """

    wall_matrix = sp.Matrix(walls)
    gram = wall_matrix * integral_scaled_dual_gram * wall_matrix.T
    size = gram.rows
    labels = tuple(
        tuple(int(gram[left, right]) for right in range(size))
        for left in range(size)
    )
    norm_values = sorted(set(labels[index][index] for index in range(size)))
    norm_colour = {value: index for index, value in enumerate(norm_values)}
    colours = tuple(norm_colour[labels[index][index]] for index in range(size))

    while True:
        signatures = []
        for left in range(size):
            neighbourhood = Counter(
                (labels[left][right], colours[right])
                for right in range(size)
                if right != left
            )
            signatures.append(
                (
                    colours[left],
                    tuple(
                        sorted(
                            (pairing, colour, count)
                            for (pairing, colour), count in neighbourhood.items()
                        )
                    ),
                )
            )
        unique = {signature: index for index, signature in enumerate(sorted(set(signatures)))}
        refined = tuple(unique[signature] for signature in signatures)
        if len(set(refined)) == len(set(colours)):
            colours = refined
            break
        colours = refined

    if len(set(colours)) == size:
        return (_identity(size),)

    graph = nx.Graph()
    for index in range(size):
        graph.add_node(index, refined_colour=colours[index])
    for left in range(size):
        for right in range(left):
            product = labels[left][right]
            if product != 0:
                graph.add_edge(left, right, product=product)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph,
        graph,
        node_match=nx.algorithms.isomorphism.categorical_node_match(
            "refined_colour", None
        ),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match(
            "product", None
        ),
    )
    return tuple(
        sorted(
            {
                tuple(int(mapping[index]) for index in range(size))
                for mapping in matcher.isomorphisms_iter()
            }
        )
    )


def _wall_action(
    walls: tuple[tuple[int, ...], ...], matrix: sp.Matrix
) -> Permutation:
    index = {wall: position for position, wall in enumerate(walls)}
    images = tuple(tuple(sp.Rational(entry) for entry in sp.Matrix(1, 20, wall) * matrix) for wall in walls)
    if any(image not in index for image in images):
        raise ArithmeticError("an archived symmetry does not preserve the wall stock")
    permutation = tuple(index[image] for image in images)
    if len(set(permutation)) != len(walls):
        raise ArithmeticError("a wall action is not a permutation")
    return permutation


def _archived_stabilizer_fast(
    atlas: ShimadaDisc11Atlas,
    chamber_index: int,
    expected_order: int,
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Generate Aut(D) without parsing unused adjacency records."""

    generators = tuple(
        atlas.assignment(atlas._chamber_name(chamber_index, "generatorsaut"))
    )
    identity = tuple(
        tuple(int(row == column) for column in range(AMBIENT_RANK))
        for row in range(AMBIENT_RANK)
    )
    seen = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            product = tuple(
                tuple(
                    sum(
                        current[row][middle] * generator[middle][column]
                        for middle in range(AMBIENT_RANK)
                    )
                    for column in range(AMBIENT_RANK)
                )
                for row in range(AMBIENT_RANK)
            )
            if product not in seen:
                seen.add(product)
                frontier.append(product)
                if len(seen) > expected_order:
                    raise ArithmeticError("Aut(D) generators exceed the archived order")
    if len(seen) != expected_order:
        raise ArithmeticError("Aut(D) generators do not attain the archived order")
    return tuple(sorted(seen))


def _linear_actions_from_graph_group(
    walls: tuple[tuple[int, ...], ...],
    graph_group: tuple[Permutation, ...],
    dual_gram: sp.Matrix,
    ambient_gram: sp.Matrix,
) -> dict[Permutation, sp.Matrix]:
    wall_matrix = sp.Matrix(walls)
    if wall_matrix.rank() != AMBIENT_RANK:
        raise ArithmeticError("a complete wall stock does not span the dual space")
    independent_rows = tuple(wall_matrix.T.rref()[1])
    basis = wall_matrix[list(independent_rows), :]
    if basis.rows != basis.cols or basis.det() == 0:
        raise ArithmeticError("failed to select an independent wall basis")
    basis_inverse = basis.inv()
    actions: dict[Permutation, sp.Matrix] = {}
    for permutation in graph_group:
        image_basis = wall_matrix[
            [permutation[index] for index in independent_rows], :
        ]
        action = basis_inverse * image_basis
        expected = wall_matrix[[permutation[index] for index in range(len(walls))], :]
        if wall_matrix * action != expected:
            raise ArithmeticError("a weighted-graph symmetry does not lift linearly")
        if action * dual_gram * action.T != dual_gram:
            raise ArithmeticError("a lifted graph symmetry is not Lorentzian")
        primal = ambient_gram * action * dual_gram
        if primal * ambient_gram * primal.T != ambient_gram:
            raise ArithmeticError("the conjugate primal action is not Lorentzian")
        actions[permutation] = action
    if len(actions) != len(graph_group):
        raise ArithmeticError("two graph symmetries induced the same action")
    return actions


def _subgroup_row(
    subgroup: Subgroup,
    actions: dict[Permutation, sp.Matrix],
) -> dict[str, Any]:
    matrices = tuple(actions[element] for element in subgroup)
    fixed_rank = _fixed_rank(matrices, AMBIENT_RANK)
    return {
        "order": len(subgroup),
        "fixed_rank": fixed_rank,
        "hyperbolic_dimension": fixed_rank - 1,
        "above_h9": fixed_rank >= MIN_FIXED_RANK,
        "element_order_histogram": {
            str(order): count
            for order, count in sorted(
                Counter(_permutation_order(element) for element in subgroup).items()
            )
        },
        "elements": [list(element) for element in sorted(subgroup)],
    }


def _regular_action(
    wall_group: tuple[Permutation, ...],
) -> tuple[tuple[Permutation, ...], dict[Permutation, Permutation]]:
    """Replace long wall permutations by the group's short regular action."""

    index = {element: position for position, element in enumerate(wall_group)}
    regular_group: list[Permutation] = []
    lookup: dict[Permutation, Permutation] = {}
    for left in wall_group:
        regular = tuple(
            index[_compose(left, right)] for right in wall_group
        )
        regular_group.append(regular)
        lookup[regular] = left
    if len(lookup) != len(wall_group):
        raise ArithmeticError("the weighted group regular action is not faithful")
    return tuple(regular_group), lookup


def build(source: Path = SOURCE, progress_every: int = 100) -> dict[str, Any]:
    source_hash = _sha256(source)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise ValueError(f"unexpected disc16.txt SHA-256 {source_hash}")
    atlas = ShimadaDisc11Atlas(source)
    ambient_gram = sp.Matrix(atlas.gram)
    dual_gram = ambient_gram.inv()
    dual_scale = sp.ilcm(
        *[sp.Rational(entry).q for entry in dual_gram]
    )
    integral_scaled_dual_gram = dual_scale * dual_gram
    if any(not entry.is_Integer for entry in integral_scaled_dual_gram):
        raise ArithmeticError("failed to clear the dual Gram denominators")
    archived_inventory_hash = _sha256(ARCHIVED_INVENTORY)
    archived_inventory = json.loads(ARCHIVED_INVENTORY.read_text(encoding="utf-8"))
    if (
        archived_inventory["source"]["sha256"] != source_hash
        or archived_inventory["source"]["chamber_count"] != atlas.chamber_count
    ):
        raise ArithmeticError("the archived Aut(D) inventory has incompatible scope")
    archived_rows = {
        row["chamber_index"]: row for row in archived_inventory["chambers"]
    }

    compact_rows: list[dict[str, Any]] = []
    nontrivial_rows: list[dict[str, Any]] = []
    graph_order_histogram: Counter[int] = Counter()
    equality_count = 0
    strict_extension_count = 0
    total_subgroups = 0
    total_classes = 0
    proper_high_classes = 0
    proper_high_dimension_histogram: Counter[int] = Counter()

    for chamber_index in range(atlas.chamber_count):
        orbit_count = int(
            atlas.assignment(atlas._chamber_name(chamber_index, "numborbitswalls"))
        )
        wall_orbits = tuple(
            atlas.assignment(atlas._chamber_name(chamber_index, "wallorbit", orbit))
            for orbit in range(orbit_count)
        )
        walls = tuple(wall for orbit in wall_orbits for wall in orbit)
        if len(set(walls)) != len(walls):
            raise ArithmeticError(f"chamber {chamber_index} has duplicate walls")
        graph_group = _weighted_wall_graph_automorphisms_fast(
            integral_scaled_dual_gram, walls
        )
        graph_order_histogram[len(graph_group)] += 1

        archived_order = int(
            atlas.assignment(atlas._chamber_name(chamber_index, "orderaut"))
        )
        archived_matrices = tuple(
            sp.Matrix(atlas.dual_action(matrix))
            for matrix in _archived_stabilizer_fast(
                atlas, chamber_index, archived_order
            )
        )
        archived_group = tuple(
            sorted(_wall_action(walls, matrix) for matrix in archived_matrices)
        )
        graph_set = set(graph_group)
        archived_set = set(archived_group)
        if not archived_set.issubset(graph_set):
            raise ArithmeticError("Aut(D) is not contained in the weighted graph group")
        equal = archived_set == graph_set
        equality_count += int(equal)
        strict_extension_count += int(not equal)
        compact_rows.append(
            {
                "chamber_index": chamber_index,
                "wall_count": len(walls),
                "archived_stabilizer_order": len(archived_group),
                "weighted_wall_group_order": len(graph_group),
                "groups_equal": equal,
            }
        )

        if len(graph_group) > 1 and equal:
            try:
                archived_row = archived_rows[chamber_index]
            except KeyError as error:
                raise ArithmeticError(
                    "a nontrivial equal group is absent from the Aut(D) inventory"
                ) from error
            total_subgroups += archived_row["subgroup_count"]
            total_classes += archived_row["subgroup_conjugacy_class_count"]
            for row in archived_row["classes"]:
                if row["order"] > 1 and row["above_h9"]:
                    proper_high_classes += 1
                    proper_high_dimension_histogram[row["hyperbolic_dimension"]] += 1
            nontrivial_rows.append(
                {
                    "chamber_index": chamber_index,
                    "wall_count": len(walls),
                    "archived_stabilizer_order": len(archived_group),
                    "weighted_wall_group_order": len(graph_group),
                    "groups_equal": True,
                    "class_source": "reused_exact_archived_Aut(D)_inventory",
                    "archived_inventory_chamber_index": chamber_index,
                    "subgroup_count": archived_row["subgroup_count"],
                    "subgroup_conjugacy_class_count": archived_row[
                        "subgroup_conjugacy_class_count"
                    ],
                    "classes": archived_row["classes"],
                }
            )
        elif len(graph_group) > 1:
            actions = _linear_actions_from_graph_group(
                walls, graph_group, dual_gram, ambient_gram
            )
            regular_group, regular_lookup = _regular_action(graph_group)
            subgroups = _all_subgroups(regular_group)
            classes = _subgroup_conjugacy_classes(regular_group, subgroups)
            class_rows: list[dict[str, Any]] = []
            for class_index, conjugacy_class in enumerate(classes):
                representative_index = conjugacy_class[0]
                wall_subgroup = frozenset(
                    regular_lookup[element]
                    for element in subgroups[representative_index]
                )
                row = _subgroup_row(wall_subgroup, actions)
                if row["order"] > 1 and row["above_h9"]:
                    proper_high_classes += 1
                    proper_high_dimension_histogram[row["hyperbolic_dimension"]] += 1
                class_rows.append(
                    {
                        "class_index": class_index,
                        "representative_subgroup_index": representative_index,
                        "subgroup_indices": list(conjugacy_class),
                        "conjugacy_class_size": len(conjugacy_class),
                        **row,
                    }
                )
            total_subgroups += len(subgroups)
            total_classes += len(classes)
            group_rows = []
            for group_index, permutation in enumerate(graph_group):
                action = actions[permutation]
                primal = ambient_gram * action * dual_gram
                group_rows.append(
                    {
                        "group_index": group_index,
                        "permutation": list(permutation),
                        "order": _permutation_order(permutation),
                        "dual_action": _matrix_payload(action),
                        "primal_action_integral": all(
                            entry.is_Integer for entry in primal
                        ),
                    }
                )
            nontrivial_rows.append(
                {
                    "chamber_index": chamber_index,
                    "wall_count": len(walls),
                    "archived_stabilizer_order": len(archived_group),
                    "weighted_wall_group_order": len(graph_group),
                    "groups_equal": False,
                    "class_source": "fresh_full_weighted_wall_group_enumeration",
                    "archived_group_elements": [list(element) for element in archived_group],
                    "weighted_group_elements": group_rows,
                    "subgroup_count": len(subgroups),
                    "subgroup_conjugacy_class_count": len(classes),
                    "classes": class_rows,
                }
            )
        if progress_every and (chamber_index + 1) % progress_every == 0:
            print(
                "weighted-wall progress",
                chamber_index + 1,
                "/",
                atlas.chamber_count,
                "nontrivial=",
                len(nontrivial_rows),
                "extensions=",
                strict_extension_count,
                flush=True,
            )

    return {
        "schema": "coxeter.shimada_disc16_weighted_wall_symmetry_inventory.v1",
        "source": {
            "path": str(source),
            "sha256": source_hash,
            "chamber_count": atlas.chamber_count,
            "ambient_rank": AMBIENT_RANK,
            "archived_AutD_inventory_path": str(ARCHIVED_INVENTORY),
            "archived_AutD_inventory_sha256": archived_inventory_hash,
        },
        "scope": {
            "all_published_chambers": True,
            "exact_weighted_wall_graph_automorphisms": True,
            "all_graph_symmetries_lifted_to_rational_lorentzian_isometries": True,
            "all_subgroups_of_nontrivial_graph_groups": True,
            "minimum_fixed_rank": MIN_FIXED_RANK,
            "wall_angle_gate_performed": False,
            "volume_gate_performed": False,
        },
        "summary": {
            "chamber_count": atlas.chamber_count,
            "weighted_group_order_histogram": {
                str(order): count
                for order, count in sorted(graph_order_histogram.items())
            },
            "nontrivial_weighted_group_chamber_count": len(nontrivial_rows),
            "archived_and_weighted_groups_equal_count": equality_count,
            "weighted_group_strict_extension_count": strict_extension_count,
            "total_subgroup_count_in_nontrivial_weighted_groups": total_subgroups,
            "total_subgroup_conjugacy_class_count": total_classes,
            "proper_high_dimensional_subgroup_class_count": proper_high_classes,
            "proper_high_dimension_histogram": {
                str(dimension): count
                for dimension, count in sorted(
                    proper_high_dimension_histogram.items()
                )
            },
        },
        "all_chambers": compact_rows,
        "nontrivial_weighted_group_chambers": nontrivial_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument(
        "--check-archive",
        action="store_true",
        help="rebuild the complete payload and require byte-identical JSON",
    )
    arguments = parser.parse_args()
    payload = build(arguments.source, arguments.progress_every)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file():
            raise SystemExit(f"archive is absent: {arguments.output}")
        if arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("weighted-wall symmetry archive replay mismatch")
        print("weighted-wall symmetry archive replay passed")
        print(json.dumps(payload["summary"], indent=2, sort_keys=True))
        return
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(rendered, encoding="utf-8")
    print(json.dumps(payload["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
