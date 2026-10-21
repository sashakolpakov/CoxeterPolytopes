"""Exact prime-order fold census for the general 15-nodal quartic K3.

The input is Ichiro Shimada's GAP-format computation file accompanying
Dolgachev--Shimada, *15-nodal quartic surfaces. Part II*.  It supplies the
rank-16 lattice ``S15``, all 314 walls of the induced chamber ``D15``, and
all 720 elements of ``O(S15,D15)`` together with its natural identification
with ``S6``.

This program uses the *full lattice chamber stabilizer*.  A surface
automorphism realization of an element is neither assumed nor required.
It verifies the published action on every wall, then screens one
representative of every prime-order cyclic conjugacy class in ``S6``:

* the three involution classes ``2 1^4``, ``2^2 1^2``, and ``2^3``;
* the two order-three classes ``3 1^3`` and ``3^2``; and
* the order-five class ``5 1``.

For each class it forms orbit-sum normals in the fixed Lorentzian subspace,
removes redundant inequalities with exact rational certificates, proves
that every retained inequality is a genuine hyperbolic facet, and checks
all normalized products against the Coxeter-angle list.  Every Coxeter
survivor is passed to the exact Vinberg edge/end-point and full parabolic
neighbourhood criterion.  Therefore ``finite_volume`` is never inferred
from acuteness or Coxeter labels alone.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from functools import lru_cache
from itertools import permutations
import json
from math import lcm
from pathlib import Path
import re
from typing import Any

import numpy as np
import sympy as sp

import scripts.search_apery_fermi_chamber_folds as apery_fold
from scripts.search_apery_fermi_chamber_folds import (
    _compose,
    _extract_assignment,
    _identity,
    _matrix_payload,
    _permutation_order,
    _rational_text,
    _screen_subgroup,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "data"
    / "external"
    / "dolgachev_shimada_15_nodal_2019"
    / "15nodalcompdata.txt"
)
STOCK_OUTPUT = (
    ROOT / "data" / "known_examples" / "general_15_nodal_quartic_chamber_stock.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "general_15_nodal_quartic_prime_order_folds.json"
)

Permutation = tuple[int, ...]

TARGET_CYCLE_TYPES = (
    (2, 1, 1, 1, 1),
    (2, 2, 1, 1),
    (2, 2, 2),
    (3, 1, 1, 1),
    (3, 3),
    (5, 1),
)


def _preprocess_gap(text: str) -> str:
    """Remove GAP's escaped physical newlines without changing tokens."""

    return text.replace("\\\n", "")


def _parse_rational_vector(payload: str) -> tuple[sp.Rational, ...]:
    payload = payload.strip()
    if payload.startswith("[") and payload.endswith("]"):
        payload = payload[1:-1]
    return tuple(sp.Rational(entry.strip()) for entry in payload.split(","))


def _cycle_type(permutation: Permutation) -> tuple[int, ...]:
    unseen = set(range(len(permutation)))
    lengths: list[int] = []
    while unseen:
        start = min(unseen)
        current = start
        length = 0
        while current in unseen:
            unseen.remove(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def _parse_source(path: Path = SOURCE) -> dict[str, Any]:
    text = _preprocess_gap(path.read_text(encoding="utf-8"))
    gram = sp.Matrix(ast.literal_eval(_extract_assignment(text, "GramS15")))
    ample = sp.Matrix(
        1, 16, _parse_rational_vector(_extract_assignment(text, "ample15"))
    )
    matrices = tuple(
        tuple(tuple(int(entry) for entry in row) for row in matrix)
        for matrix in ast.literal_eval(_extract_assignment(text, "OS15D15"))
    )
    s6_permutations = tuple(
        tuple(int(image) - 1 for image in permutation)
        for permutation in ast.literal_eval(
            _extract_assignment(text, "OS15D15permutation")
        )
    )

    records = _extract_assignment(text, "D15WallRecs")
    numbers = tuple(
        int(value) for value in re.findall(r"\bno\s*:=\s*(\d+)", records)
    )
    orbit_numbers = tuple(
        int(value) for value in re.findall(r"\borbit\s*:=\s*(\d+)", records)
    )
    wall_types = tuple(
        re.findall(r'\binnout\s*:=\s*"([^"]+)"', records)
    )
    wall_vectors = tuple(
        _parse_rational_vector(payload)
        for payload in re.findall(r"\bvector\s*:=\s*\[([^\]]+)\]", records)
    )
    published_norms = tuple(
        sp.Rational(value.strip())
        for value in re.findall(r"\bn\s*:=\s*([^,\)]+)", records)
    )
    published_ample_pairings = tuple(
        sp.Rational(value.strip())
        for value in re.findall(r"\ba\s*:=\s*([^,\)]+)", records)
    )

    record_lengths = {
        len(numbers),
        len(orbit_numbers),
        len(wall_types),
        len(wall_vectors),
        len(published_norms),
        len(published_ample_pairings),
    }
    if record_lengths != {314} or numbers != tuple(range(1, 315)):
        raise ValueError("D15WallRecs is not the published ordered 314-wall list")
    if gram.shape != (16, 16) or any(len(vector) != 16 for vector in wall_vectors):
        raise ValueError("unexpected S15 lattice or wall-vector rank")
    if len(matrices) != 720 or any(
        len(matrix) != 16 or any(len(row) != 16 for row in matrix)
        for matrix in matrices
    ):
        raise ValueError("OS15D15 is not a list of 720 rank-16 matrices")
    if len(s6_permutations) != 720:
        raise ValueError("OS15D15permutation is not a list of 720 permutations")

    walls = tuple(sp.Matrix(1, 16, vector) for vector in wall_vectors)
    calculated_norms = tuple(
        sp.Rational((wall * gram * wall.T)[0]) for wall in walls
    )
    calculated_ample_pairings = tuple(
        sp.Rational((ample * gram * wall.T)[0]) for wall in walls
    )
    if calculated_norms != published_norms:
        raise ArithmeticError("a published D15 wall norm failed exact replay")
    if calculated_ample_pairings != published_ample_pairings:
        raise ArithmeticError("a published D15 ample pairing failed exact replay")

    return {
        "gram": gram,
        "ample": ample,
        "matrices": matrices,
        "s6_permutations": s6_permutations,
        "wall_vectors": wall_vectors,
        "wall_types": wall_types,
        "wall_orbit_numbers": orbit_numbers,
        "wall_norms": calculated_norms,
        "wall_ample_pairings": calculated_ample_pairings,
    }


def _scaled_integer_rows(
    rows: tuple[tuple[sp.Rational, ...], ...]
) -> tuple[int, np.ndarray]:
    scale = lcm(*(int(value.q) for row in rows for value in row))
    array = np.array(
        [[int(value * scale) for value in row] for row in rows], dtype=np.int64
    )
    return scale, array


def _derive_full_wall_action(source: dict[str, Any]) -> tuple[Permutation, ...]:
    """Verify all 720 published isometries and return their wall actions.

    Integer NumPy multiplication is used only after clearing the wall and
    ample denominators.  A conservative a-priori bound is checked before the
    products, so every int64 operation here is exact (there is no floating
    arithmetic and no possibility of overflow).
    """

    gram = np.array(source["gram"].tolist(), dtype=np.int64)
    wall_scale, walls = _scaled_integer_rows(source["wall_vectors"])
    ample_scale, ample_rows = _scaled_integer_rows(
        (tuple(sp.Rational(value) for value in source["ample"]),)
    )
    ample = ample_rows[0]
    wall_index = {tuple(int(value) for value in row): index for index, row in enumerate(walls)}
    if len(wall_index) != 314:
        raise ArithmeticError("the scaled D15 wall vectors are not distinct")

    maximum_matrix = max(
        abs(entry) for matrix in source["matrices"] for row in matrix for entry in row
    )
    maximum_wall = int(np.max(np.abs(walls)))
    maximum_gram = int(np.max(np.abs(gram)))
    # Bounds for W M and M G M^T, respectively.
    if max(
        16 * maximum_wall * maximum_matrix,
        16 * 16 * maximum_matrix * maximum_gram * maximum_matrix,
    ) >= 2**62:
        raise OverflowError("published integer action exceeds the exact int64 guard")

    actions: list[Permutation] = []
    for raw_matrix in source["matrices"]:
        matrix = np.array(raw_matrix, dtype=np.int64)
        if not np.array_equal(matrix @ gram @ matrix.T, gram):
            raise ArithmeticError("a published OS15D15 matrix is not an S15 isometry")
        if not np.array_equal(ample @ matrix, ample):
            raise ArithmeticError("a published OS15D15 matrix does not fix ample15")
        images = walls @ matrix
        try:
            action = tuple(
                wall_index[tuple(int(value) for value in image)] for image in images
            )
        except KeyError as error:
            raise ArithmeticError(
                "a published OS15D15 matrix does not permute all D15 walls"
            ) from error
        if len(set(action)) != 314:
            raise ArithmeticError("a published wall action is not a permutation")
        actions.append(action)

    expected_s6 = set(permutations(range(6)))
    if set(source["s6_permutations"]) != expected_s6:
        raise ArithmeticError("published S6 labels are not all 720 permutations")
    if len(set(actions)) != 720:
        raise ArithmeticError("the 720-element action on D15 walls is not faithful")
    identity_index = source["s6_permutations"].index(tuple(range(6)))
    if actions[identity_index] != _identity(314):
        raise ArithmeticError("the S6 identity does not act identically on D15 walls")
    if wall_scale != 4 or ample_scale != 2:
        raise ArithmeticError("unexpected denominator scale in the published D15 data")

    # Verify the published S6 labelling against the complete wall action.
    # It is enough to multiply every element on the right by the five
    # adjacent transpositions, which generate S6.
    action_by_s6 = dict(zip(source["s6_permutations"], actions, strict=True))
    adjacent_transpositions: list[Permutation] = []
    for left in range(5):
        permutation = list(range(6))
        permutation[left], permutation[left + 1] = (
            permutation[left + 1],
            permutation[left],
        )
        adjacent_transpositions.append(tuple(permutation))
    for permutation in source["s6_permutations"]:
        for generator in adjacent_transpositions:
            product = _compose(generator, permutation)
            if action_by_s6[product] != _compose(
                action_by_s6[generator], action_by_s6[permutation]
            ):
                raise ArithmeticError("published S6 labels fail the wall-action group law")
    return tuple(actions)


def _power_sequence(permutation: Permutation) -> tuple[Permutation, ...]:
    order = _permutation_order(permutation)
    identity = _identity(len(permutation))
    powers = [identity]
    for _ in range(1, order):
        powers.append(_compose(permutation, powers[-1]))
    if _compose(permutation, powers[-1]) != identity:
        raise ArithmeticError("cyclic wall action failed its exact order check")
    return tuple(powers)


def _screen_cycle_type(
    source: dict[str, Any],
    wall_actions: tuple[Permutation, ...],
    cycle_type: tuple[int, ...],
    *,
    run_vinberg: bool,
) -> dict[str, Any]:
    indices = [
        index
        for index, permutation in enumerate(source["s6_permutations"])
        if _cycle_type(permutation) == cycle_type
    ]
    if not indices:
        raise ArithmeticError(f"missing S6 conjugacy class {cycle_type}")
    representative_index = min(
        indices, key=lambda index: source["s6_permutations"][index]
    )
    wall_generator = wall_actions[representative_index]
    order = _permutation_order(wall_generator)
    expected_order = max(cycle_type)
    if order != expected_order:
        raise ArithmeticError("wall action and S6 label have different orders")

    matrix_generator = sp.Matrix(source["matrices"][representative_index])
    matrix_powers = [sp.eye(16)]
    for _ in range(1, order):
        matrix_powers.append(matrix_powers[-1] * matrix_generator)
    if matrix_powers[-1] * matrix_generator != sp.eye(16):
        raise ArithmeticError("published matrix representative has the wrong order")
    wall_powers = _power_sequence(wall_generator)

    local_source = {
        "gram": source["gram"],
        "ample": source["ample"],
        "wall_vectors": source["wall_vectors"],
        "stabilizer": tuple(matrix_powers),
    }
    if run_vinberg:
        screen = _screen_subgroup(
            local_source,
            wall_powers,
            frozenset(wall_powers),
            0,
        )
    else:
        # The common screen calls Vinberg only after every exact Coxeter gate
        # has passed.  Replace that final call locally so that the high-rank
        # theorem census cannot be trapped by an irrelevant lower-dimensional
        # combination explosion.  The H3 control below is run literally.
        original_vinberg = apery_fold._vinberg_payload
        apery_fold._vinberg_payload = lambda _gram, _dimension: {"deferred": True}
        try:
            screen = _screen_subgroup(
                local_source,
                wall_powers,
                frozenset(wall_powers),
                0,
            )
        finally:
            apery_fold._vinberg_payload = original_vinberg
    screen["id"] = "cycle_" + "_".join(str(length) for length in cycle_type)
    screen["s6_cycle_type"] = list(cycle_type)
    screen["s6_conjugacy_class_size"] = len(indices)
    screen["published_representative_index"] = representative_index + 1
    screen["s6_representative"] = [
        image + 1 for image in source["s6_permutations"][representative_index]
    ]
    screen["surface_automorphism_realization_required"] = False
    screen["non_coxeter_ratio_histogram"] = dict(
        sorted(
            Counter(
                pair["normalized_squared_product"]
                for pair in screen["non_coxeter_pairs"]
            ).items()
        )
    )
    if screen["all_facet_pairs_are_coxeter"]:
        if run_vinberg:
            screen["vinberg"]["status"] = "completed_exactly"
        else:
            screen["vinberg"] = {
                "status": "deferred_lower_dimensional",
                "decision": None,
                "finite_volume": None,
                "compact": None,
                "reason": (
                    "Coxeter gate passed, but this dimension is outside the "
                    "complete >H9 theorem target; no volume conclusion is made"
                ),
            }
    return screen


@lru_cache(maxsize=4)
def build_artifacts(path: Path = SOURCE) -> tuple[dict[str, Any], dict[str, Any]]:
    source = _parse_source(path)
    wall_actions = _derive_full_wall_action(source)

    cycle_histogram = Counter(
        _cycle_type(permutation) for permutation in source["s6_permutations"]
    )
    order_histogram = Counter(
        _permutation_order(permutation) for permutation in source["s6_permutations"]
    )
    cycle_type_fixed_ranks: dict[tuple[int, ...], int] = {}
    for cycle_type in cycle_histogram:
        representative_index = next(
            index
            for index, permutation in enumerate(source["s6_permutations"])
            if _cycle_type(permutation) == cycle_type
        )
        matrix = sp.Matrix(source["matrices"][representative_index])
        cycle_type_fixed_ranks[cycle_type] = 16 - (matrix.T - sp.eye(16)).rank()
    walls = tuple(sp.Matrix(1, 16, vector) for vector in source["wall_vectors"])
    stock = {
        "schema": "general_15_nodal_quartic_exact_chamber_stock_v1",
        "source": {
            "paper_title": "15-nodal quartic surfaces. Part II: The automorphism group",
            "authors": ["Igor Dolgachev", "Ichiro Shimada"],
            "paper_doi": "10.1007/s12215-019-00464-7",
            "paper_arxiv": "https://arxiv.org/abs/1908.05390",
            "data_url": "https://www.math.sci.hiroshima-u.ac.jp/shimada/preprints/15NodalQuart/15nodal/compdata/15nodalcompdata.txt",
            "explanation_url": "https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/15NodalQuart/15nodal/compdatatex/15nodalCompData.pdf",
            "source_file": "15nodalcompdata.txt",
            "source_sha256": "413d7363faf4574595e9c28329737eb4ad8010e0aebae53afb76762ce9f8966a",
            "explanation_file": "15nodalCompData.pdf",
            "explanation_sha256": "bdbcad07b4854d1d7992455523a2924aaf74d0cc2d80887a3f6a297f1e69daad",
        },
        "lattice": {
            "name": "S15",
            "description": "Neron-Severi lattice of a general 15-nodal quartic K3 resolution",
            "rank": 16,
            "determinant": int(source["gram"].det()),
            "inertia": list(exact_inertia(source["gram"])),
            "gram": _matrix_payload(source["gram"]),
            "ample15": [_rational_text(value) for value in source["ample"]],
            "ample_square": _rational_text(
                sp.Rational((source["ample"] * source["gram"] * source["ample"].T)[0])
            ),
        },
        "walls": [
            {
                "index": index + 1,
                "vector": [_rational_text(value) for value in source["wall_vectors"][index]],
                "square": _rational_text(source["wall_norms"][index]),
                "ample_pairing": _rational_text(source["wall_ample_pairings"][index]),
                "kind": source["wall_types"][index],
                "published_orbit_number": source["wall_orbit_numbers"][index],
            }
            for index, _wall in enumerate(walls)
        ],
        "wall_summary": {
            "count": 314,
            "kind_histogram": dict(sorted(Counter(source["wall_types"]).items())),
            "published_orbit_histogram": {
                str(number): count
                for number, count in sorted(Counter(source["wall_orbit_numbers"]).items())
            },
            "square_histogram": {
                _rational_text(square): count
                for square, count in sorted(Counter(source["wall_norms"]).items())
            },
        },
        "chamber_isometry_group": {
            "name": "O(S15,D15)",
            "order": 720,
            "published_identification": "S6",
            "element_order_histogram": {
                str(order): count for order, count in sorted(order_histogram.items())
            },
            "cycle_type_histogram": {
                "_".join(str(length) for length in cycle_type): count
                for cycle_type, count in sorted(cycle_histogram.items(), reverse=True)
            },
            "cycle_type_fixed_space": {
                "_".join(str(length) for length in cycle_type): {
                    "fixed_rank": cycle_type_fixed_ranks[cycle_type],
                    "hyperbolic_dimension": cycle_type_fixed_ranks[cycle_type] - 1,
                }
                for cycle_type in sorted(cycle_histogram, reverse=True)
            },
            "elements": [
                {
                    "published_index": index + 1,
                    "s6_permutation": [image + 1 for image in source["s6_permutations"][index]],
                    "cycle_type": list(_cycle_type(source["s6_permutations"][index])),
                    "wall_permutation": [image + 1 for image in wall_actions[index]],
                }
                for index in range(720)
            ],
            "checks": {
                "all_720_s6_labels_are_distinct": True,
                "all_720_matrices_preserve_gram_exactly": True,
                "all_720_matrices_fix_ample15_exactly": True,
                "all_720_matrices_permute_all_314_walls_exactly": True,
                "wall_action_is_faithful": True,
                "published_s6_labels_obey_group_law_on_wall_action": True,
            },
            "surface_automorphism_realization_used_as_search_gate": False,
        },
    }

    screens = [
        _screen_cycle_type(
            source,
            wall_actions,
            cycle_type,
            # One small exact control demonstrates the complete Vinberg
            # pipeline.  H9 and H5 survivors are recorded without a volume
            # claim because they do not bear on the >H9 theorem.
            run_vinberg=cycle_type == (5, 1),
        )
        for cycle_type in TARGET_CYCLE_TYPES
    ]
    high = [screen for screen in screens if screen["hyperbolic_dimension"] > 9]
    coxeter = [screen for screen in screens if screen["all_facet_pairs_are_coxeter"]]
    vinberg_runs = [
        screen
        for screen in coxeter
        if screen["vinberg"].get("status") == "completed_exactly"
    ]
    finite_volume = [
        screen
        for screen in vinberg_runs
        if screen["vinberg"]["finite_volume"]
    ]
    enumeration = {
        "schema": "general_15_nodal_quartic_prime_order_fold_census_v1",
        "method": {
            "symmetry_scope": "all prime-order cyclic conjugacy classes of O(S15,D15) = S6",
            "surface_automorphism_required": False,
            "normal_restriction": "sum all D15 wall normals in each cyclic subgroup orbit",
            "wall_reduction": "exact positive rational dependences or strict timelike facet witnesses",
            "coxeter_gate": "exact normalized squared products",
            "vinberg_gate": "exact edge endpoints and full affine-parabolic neighbourhoods on every Coxeter survivor",
            "complete_simple_walls_for_vinberg": True,
        },
        "searched_cycle_types": [list(cycle_type) for cycle_type in TARGET_CYCLE_TYPES],
        "dimension_above_nine_exhaustion_theorem": {
            "statement": (
                "Every nontrivial subgroup of O(S15,D15)=S6 whose fixed "
                "hyperbolic space has dimension >9 is generated by one transposition"
            ),
            "fixed_space_inclusion": "Fix(H) is contained in Fix(g) for every g in H",
            "only_nonidentity_cycle_type_with_fixed_rank_above_10": [2, 1, 1, 1, 1],
            "transposition_fixed_rank": 12,
            "distinct_transposition_product": (
                "two distinct transpositions multiply to a 3-cycle or a double "
                "transposition, whose fixed ranks are at most 10"
            ),
            "conclusion": (
                "the transposition C2 screen is exhaustive for all nontrivial "
                "subgroup folds in dimensions >9"
            ),
        },
        "subgroups": screens,
        "summary": {
            "searched_conjugacy_class_count": len(screens),
            "dimension_histogram": {
                str(dimension): count
                for dimension, count in sorted(
                    Counter(screen["hyperbolic_dimension"] for screen in screens).items(),
                    reverse=True,
                )
            },
            "dimension_above_nine_count": len(high),
            "dimension_above_nine_cycle_types": [
                screen["s6_cycle_type"] for screen in high
            ],
            "dimension_above_nine_coxeter_count": sum(
                screen["all_facet_pairs_are_coxeter"] for screen in high
            ),
            "coxeter_survivor_count": len(coxeter),
            "vinberg_run_count": len(vinberg_runs),
            "finite_volume_count": len(finite_volume),
            "finite_volume_cycle_types": [
                screen["s6_cycle_type"] for screen in finite_volume
            ],
            "all_exact_reduction_checks_pass": all(
                all(screen["checks"].values()) for screen in screens
            ),
        },
    }
    return stock, enumeration


def _write_or_check(path: Path, payload: dict[str, Any], check: bool) -> None:
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if check:
        if not path.is_file() or path.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(rendered, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--stock-output", type=Path, default=STOCK_OUTPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    stock, enumeration = build_artifacts(arguments.source)
    _write_or_check(arguments.stock_output, stock, arguments.check_archive)
    _write_or_check(arguments.output, enumeration, arguments.check_archive)
    print(json.dumps(enumeration["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
