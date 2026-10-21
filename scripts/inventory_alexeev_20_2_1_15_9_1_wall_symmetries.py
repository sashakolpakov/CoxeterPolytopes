"""Exact weighted-diagram subgroup inventory for two Alexeev chambers.

The pinned arXiv ancillary file supplies complete primitive-root Gram
matrices for the reflective 2-elementary lattices ``(20,2,1)`` and
``(15,9,1)``.  This script reconstructs an integral ambient lattice basis
from the leading nonsingular root block, enumerates the full automorphism
group of each *weighted* Gram matrix, lifts every wall permutation to an
integral Lorentz isometry, and asks GAP for all subgroup conjugacy classes.

This is an inventory/rank gate only.  A diagram subgroup is never silently
identified with a geometric automorphism group, and no fixed section is
called Coxeter or finite-volume until the separate exact geometry screen.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any

import networkx as nx
import sympy as sp

from scripts.certify_alexeev_high_dim_reflection_examples import _extract_matrix
from scripts.search_apery_fermi_chamber_folds import (
    Permutation,
    Subgroup,
    _compose,
    _fixed_rank,
    _generated_subgroup,
    _identity,
    _matrix_payload,
    _permutation_order,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path("/private/tmp/alexeev_2209_source/all_cases_in_2el_paper.txt")
SOURCE_SHA256 = "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457"
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_20_2_1_15_9_1_wall_symmetry_inventory.json"
)
SCHEMA = "coxeter.alexeev_20_2_1_15_9_1_wall_symmetry_inventory.v1"
MIN_FIXED_RANK = 11
CONFIGS: dict[tuple[int, int, int], dict[str, Any]] = {
    (20, 2, 1): {
        "wall_count": 50,
        "rank": 20,
        "determinant": -4,
        "diagram_group_name": "S5",
        "diagram_group_order": 120,
        "geometric_context": "Picard lattice of the unique K3 surface with these invariants; this inventory nevertheless uses the full diagram group only",
    },
    (15, 9, 1): {
        "wall_count": 98,
        "rank": 15,
        "determinant": 512,
        "diagram_group_name": "PGL(2,7)",
        "diagram_group_order": 336,
        "geometric_context": "Picard lattice of a Hilbert square Hilb^2(X); its PGL(2,7) diagram group is larger than the S4 diagram group of the associated K3 lattice",
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _key_text(key: tuple[int, int, int]) -> str:
    return "(" + ",".join(map(str, key)) + ")"


def _refined_weighted_gram_automorphisms(gram: sp.Matrix) -> tuple[Permutation, ...]:
    """Enumerate the exact complete edge-coloured graph automorphism group."""

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
        unique = {
            signature: index for index, signature in enumerate(sorted(set(signatures)))
        }
        refined = tuple(unique[signature] for signature in signatures)
        if refined == colours:
            break
        colours = refined

    graph = nx.Graph()
    for index in range(size):
        graph.add_node(index, colour=colours[index])
    for left in range(size):
        for right in range(left):
            product = labels[left][right]
            if product:
                graph.add_edge(left, right, product=product)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph,
        graph,
        node_match=nx.algorithms.isomorphism.categorical_node_match("colour", None),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match("product", None),
    )
    group = tuple(
        sorted(
            {
                tuple(int(mapping[index]) for index in range(size))
                for mapping in matcher.isomorphisms_iter()
            }
        )
    )
    if not group or group[0] != _identity(size):
        raise ArithmeticError("weighted-graph enumeration omitted the identity")
    return group


def _ambient_reconstruction(
    wall_gram: sp.Matrix, rank: int
) -> tuple[sp.Matrix, sp.Matrix]:
    ambient_gram = wall_gram[:rank, :rank]
    wall_rows = wall_gram[:, :rank] * ambient_gram.inv()
    if not (
        ambient_gram.det() != 0
        and wall_rows[:rank, :] == sp.eye(rank)
        and wall_rows * ambient_gram * wall_rows.T == wall_gram
        and all(value.is_Integer for value in wall_rows)
        and all(
            sp.gcd_list([abs(int(value)) for value in wall_rows.row(index)]) == 1
            for index in range(wall_rows.rows)
        )
    ):
        raise ArithmeticError("the leading root block did not reconstruct an integral wall basis")
    return ambient_gram, wall_rows


def _lift_actions(
    group: tuple[Permutation, ...], wall_rows: sp.Matrix, ambient_gram: sp.Matrix
) -> dict[Permutation, sp.Matrix]:
    rank = ambient_gram.rows
    actions = {}
    for permutation in group:
        action = wall_rows[[permutation[index] for index in range(rank)], :]
        expected = wall_rows[[permutation[index] for index in range(wall_rows.rows)], :]
        if not (
            wall_rows * action == expected
            and action * ambient_gram * action.T == ambient_gram
            and all(value.is_Integer for value in action)
            and abs(int(action.det())) == 1
        ):
            raise ArithmeticError("a weighted diagram automorphism failed its integral lift")
        actions[permutation] = action
    return actions


def _greedy_generators(group: tuple[Permutation, ...]) -> tuple[Permutation, ...]:
    generated: Subgroup = frozenset({_identity(len(group[0]))})
    generators = []
    for element in group:
        if element in generated:
            continue
        generators.append(element)
        generated = _generated_subgroup(generators)
        if len(generated) == len(group):
            break
    if generated != frozenset(group):
        raise ArithmeticError("greedy permutations do not generate the full diagram group")
    return tuple(generators)


def _gap_subgroup_classes(
    group: tuple[Permutation, ...], generators: tuple[Permutation, ...]
) -> tuple[dict[str, Any], ...]:
    degree = len(group[0])
    # GAP's ``Group([])`` has no permutation degree to infer, whereas the
    # exact weighted-graph inventory legitimately returns a trivial group
    # for several of Alexeev's small chambers.  Its subgroup census is
    # canonical and needs no external enumeration.
    if len(group) == 1:
        if group != (_identity(degree),) or generators:
            raise ArithmeticError("malformed trivial weighted diagram group")
        return (
            {
                "order": 1,
                "conjugacy_class_size": 1,
                "representative_generators": tuple(),
                "representative_subgroup": frozenset({_identity(degree)}),
            },
        )
    gap_generators = ",".join(
        "PermList(" + str([value + 1 for value in generator]) + ")"
        for generator in generators
    )
    program = f"""
SizeScreen([1000000,1000000]);;
G:=Group([{gap_generators}]);;
cc:=ConjugacyClassesSubgroups(G);;
Print("GROUP|",Size(G),"|",Length(cc),"\\n");;
for c in cc do
  h:=Representative(c);;
  Print("CLASS|",Size(h),"|",Size(c),"|",List(GeneratorsOfGroup(h),x->ListPerm(x,{degree})),"\\n");;
od;
QUIT;
"""
    completed = subprocess.run(
        ["gap", "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    header = None
    rows = []
    for line in completed.stdout.splitlines():
        if line.startswith("GROUP|"):
            _, order, count = line.split("|")
            header = (int(order), int(count))
        elif line.startswith("CLASS|"):
            _, order, class_size, generator_text = line.split("|", 3)
            gap_rows = ast.literal_eval(generator_text)
            representative_generators = tuple(
                tuple(int(value) - 1 for value in row) for row in gap_rows
            )
            subgroup = (
                _generated_subgroup(representative_generators)
                if representative_generators
                else frozenset({_identity(degree)})
            )
            if len(subgroup) != int(order) or not subgroup.issubset(group):
                raise ArithmeticError("GAP returned an invalid subgroup representative")
            rows.append(
                {
                    "order": int(order),
                    "conjugacy_class_size": int(class_size),
                    "representative_generators": representative_generators,
                    "representative_subgroup": subgroup,
                }
            )
    if header != (len(group), len(rows)):
        raise ArithmeticError(
            f"GAP subgroup header mismatch: {header}, group {len(group)}, rows {len(rows)}"
        )
    return tuple(rows)


def _inventory_one(source_text: str, key: tuple[int, int, int]) -> dict[str, Any]:
    config = CONFIGS[key]
    wall_gram = _extract_matrix(source_text, key)
    rank = int(config["rank"])
    if not (
        wall_gram.shape == (config["wall_count"], config["wall_count"])
        and wall_gram == wall_gram.T
        and wall_gram.rank() == rank
        and all(
            wall_gram[left, right] >= 0
            for left in range(wall_gram.rows)
            for right in range(wall_gram.cols)
            if left != right
        )
    ):
        raise ArithmeticError(f"Alexeev matrix {_key_text(key)} failed basic checks")
    ambient_gram, wall_rows = _ambient_reconstruction(wall_gram, rank)
    if not (
        int(ambient_gram.det()) == config["determinant"]
        and exact_inertia(ambient_gram) == (1, rank - 1, 0)
    ):
        raise ArithmeticError("the reconstructed ambient lattice has wrong invariants")

    group = _refined_weighted_gram_automorphisms(wall_gram)
    print(json.dumps({"key": list(key), "weighted_group_order": len(group)}), flush=True)
    if len(group) != config["diagram_group_order"]:
        raise ArithmeticError(
            f"expected diagram group order {config['diagram_group_order']}, got {len(group)}"
        )
    actions = _lift_actions(group, wall_rows, ambient_gram)
    generators = _greedy_generators(group)
    subgroup_classes = _gap_subgroup_classes(group, generators)
    print(
        json.dumps({"key": list(key), "subgroup_class_count": len(subgroup_classes)}),
        flush=True,
    )

    classes = []
    dimension_histogram: Counter[int] = Counter()
    cyclic_dimension_histogram: Counter[int] = Counter()
    high_count = 0
    for class_index, row in enumerate(subgroup_classes):
        subgroup: Subgroup = row.pop("representative_subgroup")
        representative_generators = row.pop("representative_generators")
        matrices = tuple(actions[element] for element in representative_generators)
        fixed_rank = _fixed_rank(matrices, rank) if matrices else rank
        cyclic = any(_permutation_order(element) == len(subgroup) for element in subgroup)
        proper_high = len(subgroup) > 1 and fixed_rank >= MIN_FIXED_RANK
        if proper_high:
            high_count += 1
            dimension_histogram[fixed_rank - 1] += 1
            if cyclic:
                cyclic_dimension_histogram[fixed_rank - 1] += 1
        classes.append(
            {
                "class_index": class_index,
                "order": len(subgroup),
                "conjugacy_class_size": int(row["conjugacy_class_size"]),
                "cyclic": cyclic,
                "fixed_rank": fixed_rank,
                "hyperbolic_dimension": fixed_rank - 1,
                "proper_fixed_section_above_h9": proper_high,
                "element_order_histogram": {
                    str(order): count
                    for order, count in sorted(
                        Counter(_permutation_order(element) for element in subgroup).items()
                    )
                },
                "representative_generators": [list(element) for element in representative_generators],
                "representative_elements": [list(element) for element in sorted(subgroup)]
                if proper_high
                else None,
            }
        )

    orbit_sizes = Counter()
    for index in range(wall_gram.rows):
        orbit_sizes[len({element[index] for element in group})] += 1
    return {
        "key": list(key),
        "ambient": {
            "rank": rank,
            "hyperbolic_dimension": rank - 1,
            "gram": _matrix_payload(ambient_gram),
            "determinant": int(ambient_gram.det()),
            "inertia": list(exact_inertia(ambient_gram)),
            "leading_walls_are_integral_lattice_basis": True,
        },
        "walls": {
            "count": wall_gram.rows,
            "diagonal_histogram": {
                str(value): count
                for value, count in sorted(
                    Counter(int(wall_gram[index, index]) for index in range(wall_gram.rows)).items()
                )
            },
            "coordinates_in_leading_root_basis": _matrix_payload(wall_rows),
            "gram_reconstructed_exactly": True,
        },
        "weighted_diagram_group": {
            "name_from_Alexeev": config["diagram_group_name"],
            "order": len(group),
            "greedy_generator_count": len(generators),
            "generators": [list(element) for element in generators],
            "elements": [list(element) for element in group],
            "wall_orbit_size_histogram_by_vertex_count": {
                str(size): count for size, count in sorted(orbit_sizes.items())
            },
            "all_actions_integral_unimodular_Lorentz": True,
        },
        "subgroup_classes": classes,
        "summary": {
            "subgroup_conjugacy_class_count": len(classes),
            "proper_fixed_rank_at_least_11_count": high_count,
            "proper_high_dimension_histogram": {
                str(key): value for key, value in sorted(dimension_histogram.items())
            },
            "proper_cyclic_high_dimension_histogram": {
                str(key): value for key, value in sorted(cyclic_dimension_histogram.items())
            },
        },
        "geometric_distinction": {
            "diagram_group_used_as_search_denominator": True,
            "surface_or_hyperkahler_automorphism_membership_used_as_gate": False,
            "context_from_Alexeev": config["geometric_context"],
            "individual_subgroup_geometric_realization_status": "not certified by this inventory",
        },
    }


def build(source: Path = SOURCE) -> dict[str, Any]:
    digest = _sha256(source)
    if digest != SOURCE_SHA256:
        raise ArithmeticError(f"unexpected Alexeev ancillary SHA-256 {digest}")
    source_text = source.read_text(encoding="utf-8")
    examples = [_inventory_one(source_text, key) for key in CONFIGS]
    return {
        "schema": SCHEMA,
        "source": {
            "path": str(source),
            "sha256": digest,
            "arxiv": "2209.09110v4",
            "attachment": "all_cases_in_2el_paper.txt",
        },
        "scope": {
            "matrix_keys": [list(key) for key in CONFIGS],
            "excluded_parallel_assignment": [22, 2, 0],
            "group": "full exact automorphism group of the oriented weighted primitive-root Gram matrix",
            "minimum_fixed_rank": MIN_FIXED_RANK,
            "finite_volume_claimed": False,
        },
        "examples": examples,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    payload = build(args.source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"summaries": [row["summary"] for row in payload["examples"]]}, sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
