"""Certify the D35 counterexample to one-witness restricted-facet traversal.

The fixture stores one decorated physical D35 state and two exact points in
the relative interior of the same restricted facet.  Both witnesses satisfy
the complete primitive wall system, but exact transverse crossings terminate
in different induced Shimada chambers.  This is a small regression for the
facet-star/carrier repair; it does not rerun the obsolete bounded traversal.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden
from scripts.analyze_shimada_disc11_enhanced_quotient import _canonical_sha
from scripts.explore_shimada_disc11_global_eigen_nef import (
    _cross_internal_facet,
    _face_is_exact,
    _physical_state_key,
    _restricted_rows,
)
from scripts.shimada_disc11_atlas import ShimadaChamberState, ShimadaDisc11Atlas
from scripts.shimada_disc11_conjugate import DISC11_SOURCE_SHA256


DEFAULT_FIXTURE = Path(
    "data/candidates/shimada_disc11_d35_facet_subdivision.json"
)


def certify_d35_facet_subdivision(
    disc11: Path,
    fixture_path: Path = DEFAULT_FIXTURE,
) -> dict[str, object]:
    fixture = json.loads(fixture_path.read_text())
    source_digest = hashlib.sha256(disc11.read_bytes()).hexdigest()
    if source_digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")
    if fixture["source_sha256"] != source_digest:
        raise ValueError("the D35 fixture belongs to a different atlas")

    state = ShimadaChamberState(
        int(fixture["state_index"]),
        tuple(tuple(int(entry) for entry in row) for row in fixture["row_isometry"]),
    )
    interior = tuple(Golden.from_json(entry) for entry in fixture["interior"])
    wall_class = tuple(
        (int(label[0]), int(label[1])) for label in fixture["wall_class"]
    )
    atlas = ShimadaDisc11Atlas(disc11)
    state_key = _physical_state_key(atlas, state)
    rows = _restricted_rows(atlas, state)

    records = []
    for expected in fixture["witnesses"]:
        witness = tuple(Golden.from_json(entry) for entry in expected["witness"])
        witness_digest = _canonical_sha([entry.to_json() for entry in witness])
        face_exact = _face_is_exact(rows, wall_class, witness, interior)
        neighbor, _, gallery, denominator = _cross_internal_facet(
            atlas, state, interior, witness, wall_class
        )
        target_key = _physical_state_key(atlas, neighbor)
        representatives = [item[0] for item in gallery] + [neighbor.index]
        record = {
            "witness_sha256": witness_digest,
            "full_exact_relative_interior_witness": face_exact,
            "crossing_denominator": denominator,
            "gallery_representatives": representatives,
            "target_representative": neighbor.index,
            "target_physical_state_sha256": _canonical_sha(target_key),
        }
        if witness_digest != expected["witness_sha256"]:
            raise AssertionError("a D35 witness digest changed")
        if not face_exact:
            raise AssertionError("a D35 fixture point is not a strict exact witness")
        for key in (
            "gallery_representatives",
            "target_representative",
            "target_physical_state_sha256",
        ):
            if record[key] != expected[key]:
                raise AssertionError(f"D35 crossing field {key!r} changed")
        records.append(record)

    targets = {record["target_physical_state_sha256"] for record in records}
    verified = (
        _canonical_sha(state_key) == fixture["physical_state_sha256"]
        and len(records) == 2
        and len(targets) == 2
        and all(record["full_exact_relative_interior_witness"] for record in records)
    )
    return {
        "analysis_id": fixture["analysis_id"],
        "verified": verified,
        "source_sha256": source_digest,
        "physical_state_sha256": _canonical_sha(state_key),
        "state_representative": state.index,
        "wall_class": [list(label) for label in wall_class],
        "witness_crossings": records,
        "distinct_opposite_physical_states": len(targets),
        "conclusion": (
            "one_exact_restricted_facet_has_at_least_two_open_opposite_"
            "induced_cells;_one_witness_per_facet_is_not_graph_complete"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    parser.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE)
    args = parser.parse_args()
    print(
        json.dumps(
            certify_d35_facet_subdivision(args.disc11, args.fixture),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
