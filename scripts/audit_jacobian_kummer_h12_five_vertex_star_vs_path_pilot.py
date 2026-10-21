#!/usr/bin/env python3
"""Compare the complete JK-H12 K1,4 census with the five P5 pilots."""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scripts.certify_jacobian_kummer_h12_small_subtree_census as base
import scripts.certify_jacobian_kummer_h12_five_chamber_pilot as pilot
import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star


STAR_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_five_vertex_star_census.py"
STAR_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_star_census.json"
PILOT_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_five_chamber_pilot.py"
PILOT_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_chamber_pilot.json"
OUTPUT = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_star_vs_path_pilot_audit.json"
SCHEMA = "coxeter.jacobian_kummer.h12_five_vertex_star_vs_path_pilot_audit.v1"
STATUS_DATE = "2026-08-08"
PINNED_SHA256 = {
    "star_builder": "bbb99d59726e0b6bcf0d6b3748823ca9bc8f1efe5ecbe5c29a16f4a5db6a53fb",
    "star_archive": "21fbd0946b910360cdd4096392e4e629907f971be85359d5b87e312a78d4f309",
    "path_pilot_builder": "d1e678dc2eea7c67998fff1d5fcef34a8ad45281fac87424f74f3876ea4c9f6d",
    "path_pilot_archive": "809302f4aa10fcf723c68061c34a4d857c9d5c64a4c03e04227c4c54c917ebec",
}

Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _path_boundary_roots(
    gram: Matrix,
    roots: tuple[Vector, ...],
    interior: Vector,
    colors: tuple[int, int, int, int],
) -> np.ndarray:
    transformations, reflections = pilot._path_transformations(gram, roots, colors)
    transform_index = {
        transformation: index for index, transformation in enumerate(transformations)
    }
    internal: set[Vector] = set()
    adjacency = set()
    for chamber, transformation in enumerate(transformations):
        for color in colors:
            neighbor = base._matmul(transformation, reflections[color])
            if neighbor not in transform_index:
                continue
            other = transform_index[neighbor]
            if chamber < other:
                adjacency.add((chamber, other, color))
            internal.add(base._primitive_key(base._matvec(transformation, roots[color])))
    expected = {(index, index + 1, color) for index, color in enumerate(colors)}
    if adjacency != expected or len(internal) != 4:
        raise ArithmeticError("the exact pilot path adjacency changed")
    union_interior = tuple(
        sum(
            base._matvec(transformation, interior)[coordinate]
            for transformation in transformations
        )
        for coordinate in range(star.RANK)
    )
    boundary: dict[Vector, Vector] = {}
    for transformation in transformations:
        for seed_root in roots:
            image = base._matvec(transformation, seed_root)
            key = base._primitive_key(image)
            if key in internal:
                continue
            pairing = base._pair(gram, union_interior, image)
            if pairing == 0:
                raise ArithmeticError("a pilot interior lies on a boundary")
            oriented = image if pairing > 0 else tuple(-value for value in image)
            boundary.setdefault(key, oriented)
    ordered = np.array([boundary[key] for key in sorted(boundary)], dtype=np.int64)
    if not {base._primitive_key(root) for root in roots[:60]} <= set(boundary):
        raise ArithmeticError("the pilot lost the fixed Picard core")
    return ordered


def build() -> dict[str, Any]:
    paths = {
        "star_builder": STAR_SCRIPT,
        "star_archive": STAR_ARCHIVE,
        "path_pilot_builder": PILOT_SCRIPT,
        "path_pilot_archive": PILOT_ARCHIVE,
    }
    actual = {name: _sha256(path) for name, path in paths.items()}
    if actual != PINNED_SHA256:
        raise ArithmeticError(f"a pinned cross-shape prerequisite changed: {actual}")
    star_archive = json.loads(STAR_ARCHIVE.read_text(encoding="utf-8"))
    pilot_archive = json.loads(PILOT_ARCHIVE.read_text(encoding="utf-8"))
    if not (star_archive["all_checks_pass"] and pilot_archive["all_checks_pass"]):
        raise ArithmeticError("a cross-shape prerequisite is not green")

    star_degree_hashes: defaultdict[int, set[str]] = defaultdict(set)
    star_type_hashes = set()
    for row in star_archive["orbit_representatives"]:
        star_degree_hashes[row["facet_count"]].add(
            row["ordinary_degree_profile_multiset_sha256"]
        )
        star_type_hashes.add(row["ordinary_Coxeter_type_invariant_sha256"])
    if len(star_type_hashes) != 16_872:
        raise ArithmeticError("the frozen star nonisomorphism count changed")

    _carrier, gram, roots, interior, _bounded = base._load_inputs()
    engine = star.ExactStarBoundaryEngine(gram, roots, interior)
    rows = []
    degree_disjoint = True
    for (b_count, colors), archived in zip(
        pilot.PILOTS, pilot_archive["pilots"], strict=True
    ):
        boundary = _path_boundary_roots(gram, roots, interior, colors)
        labels, degree_key, ordinary = engine.ordinary_data(boundary)
        del labels, degree_key
        boundary_hash = star._json_sha256(boundary.tolist())
        if not (
            len(boundary) == archived["facet_count"] == 752 + 9 * b_count
            and boundary_hash == archived["marked_boundary_root_sha256"]
            and ordinary["ordinary_Coxeter_relation_histogram"]
            == archived["ordinary_Coxeter_relation_histogram"]
        ):
            raise ArithmeticError("an exact P5 pilot reconstruction changed")
        degree_hash = ordinary["ordinary_degree_profile_multiset_sha256"]
        collision = degree_hash in star_degree_hashes[len(boundary)]
        degree_disjoint &= not collision
        rows.append(
            {
                "B_edge_occurrence_count": b_count,
                "facet_count": len(boundary),
                "P5_edge_colors_one_based": [color + 1 for color in colors],
                "marked_boundary_root_sha256": boundary_hash,
                "ordinary_degree_profile_multiset_sha256": degree_hash,
                "equal_degree_profile_multiset_found_in_complete_K14_census": collision,
                "K14_type_count_with_same_facet_count": star_archive["theorem"][
                    "type_facet_histogram"
                ][str(len(boundary))],
            }
        )
    if not degree_disjoint:
        raise ArithmeticError("a P5 pilot degree profile collided with the K1,4 census")
    checks = {
        "complete_16872_type_K14_archive_pinned_and_green": True,
        "all_five_P5_pilots_rebuilt_from_exact_rank13_seed": True,
        "all_five_rebuilt_boundary_hashes_match_the_frozen_pilot": True,
        "all_five_ordinary_relation_histograms_match_the_frozen_pilot": True,
        "SHA256_digests_of_exact_degree_profile_multisets_separate_P5_from_K14": degree_disjoint,
        "certificate_grade_16877_five_chamber_types_are_pairwise_nonisomorphic": degree_disjoint,
        "claim_remains_project_local_with_worldwide_priority_open": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"cross-shape audit failed: {checks}")
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": True,
        "pinned_inputs": PINNED_SHA256,
        "result": {
            "complete_K14_type_count": 16_872,
            "explicit_P5_pilot_type_count": 5,
            "pairwise_nonisomorphic_combined_five_chamber_type_count": 16_877,
            "cross_shape_separator": (
                "SHA-256 digest of the exact ordinary degree-profile multiset"
            ),
            "cross_shape_evidence_grade": "certificate-grade hash separation",
            "finite_volume": True,
            "compact": False,
            "reflection_index_and_volume_multiplier": 5,
            "worldwide_priority": "open",
        },
        "P5_rows": rows,
        "claim_boundary": {
            "complete_P5_census_claimed": False,
            "complete_degree_3_fork_census_claimed": False,
            "worldwide_novelty_claimed": False,
            "project_internal_nonisomorphism_claimed": True,
        },
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = build()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archive differs from the exact star-vs-path rebuild")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("combined_five_chamber_types=16877")
    print(f"payload_sha256={_sha256(arguments.output)}")


if __name__ == "__main__":
    main()
