#!/usr/bin/env python3
"""Certify a JK-H12 D_infinity x C2^3 path-cube family."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Sequence

import scripts.certify_jacobian_kummer_h12_single_doublings as base
import scripts.scout_jacobian_kummer_h12_path_cubes as scout


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/jacobian_kummer_h12_path_cube_family.json"
SCHEMA = "coxeter.jacobian_kummer.h12_path_cube_family.v1"

SCOUT_SHA256 = "147d5cdd90ee2a5e209ba2b167453b63c9efa34dc91885b5e264f4ff34183dea"
SINGLE_BUILDER_SHA256 = "0bd5672b54e091da680192be8873fa7397a4e9c162c592dff6cf1ebd6e15386b"
SINGLE_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_single_doublings.json"
SINGLE_ARCHIVE_SHA256 = "01eb89baef2cedcf78d3a20dd1c55931db9a9a2ca83bc5f4b7c1a1329cbdca98"
CUBE_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_doubling_cube.json"
CUBE_ARCHIVE_SHA256 = "0267e3d56b5b0acaa7015e4cadac1d357d387532831e89485c1846eae28be3f1"
HOSTILE_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_path_cube_hostile_audit.json"
HOSTILE_ARCHIVE_SHA256 = "5c81e4c311e96a64a86ec377d10aa40b54dcd94f685432fee11d9bc9f8ea5663"

# W68/W63 and W31/W58/W59.  The path starts across W68, so N=1 and N=2
# recover the frozen 941- and 1827-facet cube members.
CANONICAL_PAIR = (67, 62)
CANONICAL_CUBE = (30, 57, 58)
DIRECT_RANGE = range(1, 9)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _vector_sum(vectors: Sequence[scout.Vector]) -> scout.Vector:
    return tuple(sum(vector[index] for vector in vectors) for index in range(len(vectors[0])))


def _root_digest(roots: Sequence[Sequence[int]]) -> str:
    return sha256(
        json.dumps([list(root) for root in roots], separators=(",", ":")).encode()
    ).hexdigest()


def _pair_int(gram: base.Matrix, left: scout.Vector, right: scout.Vector) -> int:
    return sum(
        left[row] * gram[row][column] * right[column]
        for row in range(len(left))
        for column in range(len(right))
    )


def _full_rank_mod_prime(rows: Sequence[scout.Vector], prime: int = 1_000_000_007) -> bool:
    """A rank-13 reduction modulo one prime proves rational rank 13."""

    pivots: dict[int, list[int]] = {}
    for source in rows:
        row = [value % prime for value in source]
        for column in sorted(pivots):
            if row[column]:
                coefficient = row[column]
                pivot = pivots[column]
                row = [
                    (value - coefficient * pivot_value) % prime
                    for value, pivot_value in zip(row, pivot, strict=True)
                ]
        pivot_column = next((index for index, value in enumerate(row) if value), None)
        if pivot_column is None:
            continue
        inverse = pow(row[pivot_column], -1, prime)
        row = [(value * inverse) % prime for value in row]
        pivots[pivot_column] = row
        if len(pivots) == len(source):
            return True
    return False


def _boundary(
    gram: base.Matrix,
    roots: tuple[scout.Vector, ...],
    interior: scout.Vector,
    reflections: dict[int, scout.Matrix],
    transforms: tuple[scout.Matrix, ...],
    generators: tuple[int, ...],
) -> dict[str, Any]:
    transform_set = set(transforms)
    internal: set[tuple[int, ...]] = set()
    for transform in transforms:
        for wall in generators:
            if scout._matmul(transform, reflections[wall]) in transform_set:
                internal.add(scout._primitive_key(scout._matvec(transform, roots[wall])))

    union_interior = _vector_sum(
        tuple(scout._matvec(transform, interior) for transform in transforms)
    )
    boundary: dict[tuple[int, ...], tuple[scout.Vector, int]] = {}
    for transform in transforms:
        for seed_index, seed_root in enumerate(roots):
            image = scout._matvec(transform, seed_root)
            key = scout._primitive_key(image)
            if key in internal:
                continue
            pairing = _pair_int(gram, union_interior, image)
            if pairing == 0:
                raise ArithmeticError("the product interior lies on a boundary")
            oriented = image if pairing > 0 else tuple(-value for value in image)
            boundary.setdefault(key, (oriented, seed_index))

    records = tuple(boundary[key] for key in sorted(boundary))
    ordered = tuple(record[0] for record in records)
    if not (
        _pair_int(gram, union_interior, union_interior) > 0
        and all(_pair_int(gram, union_interior, root) > 0 for root in ordered)
        and _full_rank_mod_prime(ordered)
    ):
        raise ArithmeticError("a path-cube boundary failed exact global gates")

    seed_covectors = tuple(base._covector(gram, root) for root in roots)
    seed_squares = tuple(
        base._pair_from_covector(root, covector)
        for root, covector in zip(roots, seed_covectors, strict=True)
    )
    seed_divisibilities = tuple(
        math.gcd(*(abs(int(value)) for value in covector if value))
        for covector in seed_covectors
    )
    squares = tuple(seed_squares[record[1]] for record in records)
    divisibilities = tuple(seed_divisibilities[record[1]] for record in records)
    return {
        "copy_count_and_volume_multiple": len(transforms),
        "internal_mirror_count": len(internal),
        "facet_count": len(ordered),
        "boundary_root_digest": _root_digest(ordered),
        "full_linear_rank": 13,
        "full_rank_certificate": "rank 13 modulo 1000000007",
        "all_boundary_roots_primitive": True,
        "primitivity_reason": "images under integral involutions of primitive seed roots",
        "all_boundary_reflections_integral": True,
        "reflection_integrality_reason": "conjugates of integral seed reflections by integral lattice isometries",
        "strict_interior_square": base._payload_number(
            _pair_int(gram, union_interior, union_interior)
        ),
        "primitive_root_square_histogram": {
            str(key): value for key, value in sorted(Counter(squares).items())
        },
        "primitive_root_divisibility_histogram": {
            str(key): value for key, value in sorted(Counter(divisibilities).items())
        },
    }


def _direct_census(
    gram: base.Matrix,
    roots: tuple[scout.Vector, ...],
    reflections: dict[int, scout.Matrix],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    expected = [886 * copies + 55 for copies in DIRECT_RANGE]
    for pair, cube in scout.CANDIDATES:
        cube_actions = scout._cube_transforms(reflections, cube)
        phases: list[list[int]] = []
        for start in (0, 1):
            values = []
            for copies in DIRECT_RANGE:
                transforms = scout._product(
                    scout._path_transforms(reflections, pair, copies, start),
                    cube_actions,
                )
                values.append(
                    scout._boundary_count(
                        roots, reflections, transforms, (*pair, *cube)
                    )
                )
            if values != expected:
                raise ArithmeticError("a direct path-cube count changed")
            phases.append(values)
        rows.append(
            {
                "parallel_pair_one_based": [index + 1 for index in pair],
                "orthogonal_cube_one_based": [index + 1 for index in cube],
                "phase_0_counts_N_1_through_8": phases[0],
                "phase_1_counts_N_1_through_8": phases[1],
                "both_phases_match_886N_plus_55_through_N_8": True,
            }
        )
    return rows


def build() -> dict[str, Any]:
    pinned = {
        "scout_script": _sha(Path(scout.__file__)),
        "single_builder": _sha(Path(base.__file__)),
        "single_archive": _sha(SINGLE_ARCHIVE),
        "cube_archive": _sha(CUBE_ARCHIVE),
        "all_parameter_hostile_archive": _sha(HOSTILE_ARCHIVE),
    }
    if pinned != {
        "scout_script": SCOUT_SHA256,
        "single_builder": SINGLE_BUILDER_SHA256,
        "single_archive": SINGLE_ARCHIVE_SHA256,
        "cube_archive": CUBE_ARCHIVE_SHA256,
        "all_parameter_hostile_archive": HOSTILE_ARCHIVE_SHA256,
    }:
        raise ArithmeticError("a pinned path-cube input changed")

    single = json.loads(SINGLE_ARCHIVE.read_text(encoding="utf-8"))
    cube_payload = json.loads(CUBE_ARCHIVE.read_text(encoding="utf-8"))
    hostile = json.loads(HOSTILE_ARCHIVE.read_text(encoding="utf-8"))
    if not (
        single["all_checks_pass"]
        and cube_payload["all_checks_pass"]
        and hostile["all_checks_pass"]
        and hostile["universal_family"]["facet_formula"]
        == "f(N)=886*N+55 for every integer N>=1"
        and hostile["universal_family"][
            "one_isometry_class_among_the_eight_listed_configurations_for_each_N"
        ]
    ):
        raise ArithmeticError("a prerequisite certificate is not green")

    source, gram, root_rows, interior = base._load_source()
    roots = tuple(tuple(root) for root in root_rows)
    relations = base._seed_relations(gram, roots)
    doubling = set(base._doubling_walls(relations))
    selected_stock = {
        index for pair, cube in scout.CANDIDATES for index in (*pair, *cube)
    }
    if not selected_stock <= doubling:
        raise ArithmeticError("a proposed path-cube generator is not doubling")
    reflections = {
        index: scout._reflection_matrix(gram, roots[index])
        for index in selected_stock
    }

    generator_rows = []
    for pair, cube in scout.CANDIDATES:
        if not (
            base._pair(gram, roots[pair[0]], roots[pair[0]]) == -4
            and base._pair(gram, roots[pair[1]], roots[pair[1]]) == -4
            and base._pair(gram, roots[pair[0]], roots[pair[1]]) == 4
            and all(
                base._pair(gram, roots[left], roots[right]) == 0
                for left in cube
                for right in cube
                if left != right
            )
            and all(
                base._pair(gram, roots[path], roots[finite]) == 0
                for path in pair
                for finite in cube
            )
        ):
            raise ArithmeticError("a D_infinity x C2^3 Gram factor changed")
        generator_rows.append(
            {
                "parallel_pair_one_based": [index + 1 for index in pair],
                "orthogonal_cube_one_based": [index + 1 for index in cube],
                "path_gram": [[-4, 4], [4, -4]],
                "cube_gram": [[-4 if i == j else 0 for j in range(3)] for i in range(3)],
                "cross_gram_zero": True,
                "group_factor": "D_infinity x C2^3",
            }
        )

    census = _direct_census(gram, roots, reflections)
    canonical_cube = scout._cube_transforms(reflections, CANONICAL_CUBE)
    canonical_rows = []
    for copies in DIRECT_RANGE:
        transforms = scout._product(
            scout._path_transforms(
                reflections, CANONICAL_PAIR, copies, start=0
            ),
            canonical_cube,
        )
        row = _boundary(
            gram,
            roots,
            interior,
            reflections,
            transforms,
            (*CANONICAL_PAIR, *CANONICAL_CUBE),
        )
        row["N"] = copies
        row["expected_facet_count_886N_plus_55"] = 886 * copies + 55
        if row["facet_count"] != row["expected_facet_count_886N_plus_55"]:
            raise ArithmeticError("the canonical boundary formula changed")
        canonical_rows.append(row)

    frozen_cube_rows = {
        tuple(row["selected_walls_one_based"]): row
        for row in cube_payload["constructions"]
    }
    if not (
        canonical_rows[0]["boundary_root_digest"]
        == frozen_cube_rows[(31, 58, 59)]["boundary_root_digest"]
        and canonical_rows[1]["boundary_root_digest"]
        == frozen_cube_rows[(31, 58, 59, 68)]["boundary_root_digest"]
    ):
        raise ArithmeticError("the N=1/N=2 cube compatibility changed")

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-08",
        "pinned_inputs": pinned,
        "seed": {
            "id": "JK-H12/172",
            "dimension": 12,
            "facet_count": 172,
            "finite_volume": True,
            "compact": False,
            "full_picard_lattice": "U(2)+D5(-1)^2+<-4>",
            "k3_moduli_dimension": 7,
            "source_archive_sha256": base.SOURCE_SHA256,
        },
        "generator_configurations": generator_rows,
        "complete_configuration_count": len(generator_rows),
        "direct_N_1_through_8_census": census,
        "canonical_configuration": {
            "parallel_pair_one_based": [68, 63],
            "orthogonal_cube_one_based": [31, 58, 59],
            "phase": 0,
            "boundaries": canonical_rows,
            "N_1_matches_frozen_H12_941_cube": True,
            "N_2_matches_frozen_H12_1827_cube": True,
        },
        "family": {
            "group_factor": "D_infinity x C2^3",
            "exists_for_every_integer_N_at_least_1": True,
            "copy_index_and_volume_multiplier": "8*N",
            "finite_volume": True,
            "compact": False,
            "facet_formula": "886*N+55",
            "facet_formula_directly_replayed_range": "1<=N<=8, both phases, all eight configurations",
            "facet_formula_all_parameter_status": "frozen by separate symbolic all-integer hostile audit",
            "all_eight_configurations_and_both_phases_one_isometry_class_for_fixed_N": True,
            "pairwise_nonisometric_infinite_family": True,
        },
        "independent_all_parameter_audit": {
            "archive_sha256": HOSTILE_ARCHIVE_SHA256,
            "all_integer_collision_equation_count": 3_526_884,
            "facet_formula_proved_for_every_N_at_least_1": True,
            "configuration_and_phase_conjugacy_proved": True,
            "tests": "5 passed in 93.27s",
        },
        "carrier_scope": {
            "same_full_picard_k3_family": True,
            "all_boundary_roots_integral_divisor_classes": True,
            "weak_divisor_reflection_carrier_only": True,
            "full_nef_effectivity_or_surface_action_claim": False,
        },
        "novelty_scope": {
            "new_to_project_ledger": True,
            "worldwide_priority": "open",
            "doubling_mechanism": "classical",
        },
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build()
    encoded = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_bytes(encoded)
    if args.check_archive and OUTPUT.read_bytes() != encoded:
        raise ArithmeticError("the frozen path-cube archive changed")
    if not args.write and not args.check_archive:
        print(json.dumps(payload["family"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
