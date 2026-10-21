#!/usr/bin/env python3
"""Certify the silver Vinberg shell ``x=4+3sqrt(2)`` exactly."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    _silver_pair,
    run_silver_k3_vinberg_four_plus_three_sqrt2_extension,
    silver_k3_alexeev_field_gram,
)
from scripts.independent_silver_k3_vinberg_4_plus_3sqrt2 import (
    EXPECTED_REPLAY_SHA256,
    audit_block_silver_four_plus_three_shell,
)


PRIOR_CERTIFICATE_SHA256 = (
    "e2283e6330c93632035febe1b84ef1d499a7b9507931008940b956f1470e10d5"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "391b00b000fb57e8073c873a8db0685d42138dd5cf795f3150bfe6ab430b62c1"
)
ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT / "data/candidates/silver_k3_norm2_vinberg_extension_4_plus_3sqrt2.json"
)


def _phase(message: str) -> None:
    print(f"[x=4+3sqrt(2) certificate] {message}", file=sys.stderr, flush=True)


def _jsonable(value: Any) -> Any:
    if isinstance(value, Silver):
        return value.to_json()
    if isinstance(value, Fraction):
        return (
            str(value.numerator)
            if value.denominator == 1
            else f"{value.numerator}/{value.denominator}"
        )
    if isinstance(value, Inertia):
        return {
            "positive": value.positive,
            "negative": value.negative,
            "zero": value.zero,
        }
    if isinstance(value, tuple):
        return [_jsonable(entry) for entry in value]
    if isinstance(value, list):
        return [_jsonable(entry) for entry in value]
    if isinstance(value, Mapping):
        return {str(key): _jsonable(entry) for key, entry in value.items()}
    return value


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _silver_histogram(values: Sequence[Silver]) -> tuple[dict[str, object], ...]:
    return tuple(
        {"value": value, "multiplicity": multiplicity}
        for value, multiplicity in sorted(Counter(values).items())
    )


def _inertia_histogram(values: Sequence[Inertia]) -> tuple[dict[str, object], ...]:
    return tuple(
        {"inertia": inertia_value, "multiplicity": multiplicity}
        for inertia_value, multiplicity in sorted(
            Counter(values).items(),
            key=lambda item: (
                item[0].positive,
                item[0].negative,
                item[0].zero,
            ),
        )
    )


def build_manifest() -> dict[str, object]:
    _phase("building the exact production extension")
    extension = run_silver_k3_vinberg_four_plus_three_sqrt2_extension()
    _phase(
        f"production extension complete; walls={len(extension.roots)}, "
        f"new roots={len(extension.new_level.roots)}"
    )
    level = extension.new_level
    _phase("running the independent block-factorized replay")
    replay = audit_block_silver_four_plus_three_shell()
    _phase("independent replay complete; assembling the literal certificate")
    first_new = extension.prior_wall_count
    new_indices = tuple(range(first_new, len(extension.roots)))
    prior_pairings = tuple(
        tuple(extension.wall_gram[index][old] for old in range(first_new))
        for index in new_indices
    )
    same_height_pairings = tuple(
        extension.wall_gram[index][other]
        for index in new_indices
        for other in new_indices
        if other < index
    )
    lattice_gram = silver_k3_alexeev_field_gram()
    rejection_witnesses = tuple(
        (
            root,
            tuple(
                (
                    index,
                    _silver_pair(root, lattice_gram, extension.roots[index]),
                )
                for index in blockers
            ),
        )
        for root, blockers in level.previous_wall_rejection_witnesses
    )
    blocker_size_histogram = tuple(
        sorted(Counter(len(blockers) for _root, blockers in rejection_witnesses).items())
    )
    blocker_wall_histogram = tuple(
        sorted(
            Counter(
                index
                for _root, blockers in rejection_witnesses
                for index, _pairing in blockers
            ).items()
        )
    )
    integral_count = level.exact_weight_solutions - level.nonintegral_reconstructions
    accepted_count = len(level.roots)

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "silver_k3_intrinsic_norm_two_vinberg_four_plus_three_shell",
        "claim": (
            f"At x=4+3sqrt(2), exact two-embedding enumeration gives "
            f"{level.exact_weight_solutions} weight solutions. Of these, "
            f"{level.nonintegral_reconstructions} reconstruct nonintegrally, "
            f"{level.previous_wall_rejections} integral roots are blocked by "
            f"lower-height walls, and {accepted_count} roots enter the "
            "current Vinberg prefix."
        ),
        "terminology": (
            f"The {accepted_count} accepted roots are simple relative to the "
            "complete accepted lower-height prefix and their simultaneous "
            "same-height cohort; the resulting wall list is not claimed to "
            "be the complete chamber."
        ),
        "prior_empty_five_plus_two_shell_certificate_sha256": (
            PRIOR_CERTIFICATE_SHA256
        ),
        "ordered_prior_prefix": {
            "wall_count": first_new,
            "certificate_sha256": PRIOR_CERTIFICATE_SHA256,
            "reproducer": (
                "run_silver_k3_vinberg_five_plus_two_sqrt2_extension().roots"
            ),
            "ordered_roots_sha256": _canonical_sha256(
                _jsonable(extension.roots[:first_new])
            ),
        },
        "ordered_prefix_root_count": len(extension.roots),
        "ordered_prefix_root_order": (
            "333-root certified lower-height prefix followed by the "
            "393-root production weight-recursion order"
        ),
        "ordered_prefix_roots": _jsonable(extension.roots),
        "ordered_prefix_roots_sha256": _canonical_sha256(
            _jsonable(extension.roots)
        ),
        "completed_through_control_coordinate": _jsonable(
            extension.completed_through_control_coordinate
        ),
        "next_admissible_control_coordinate": _jsonable(
            extension.next_admissible_control_coordinate
        ),
        "candidate_coordinates": _jsonable(extension.candidate_coordinates),
        "new_level": {
            "control_coordinate": _jsonable(level.control_coordinate),
            "height": _jsonable(level.height),
            "target": _jsonable(level.target),
            "alphabet": _jsonable(asdict(level.alphabet)),
            "visited_by_depth": list(level.visited_by_depth),
            "schur_pruned_by_depth": list(level.schur_pruned_by_depth),
            "exact_weight_solutions": level.exact_weight_solutions,
            "nonintegral_reconstructions": level.nonintegral_reconstructions,
            "integral_reconstructions": integral_count,
            "previous_wall_rejections": level.previous_wall_rejections,
            "blocker_set_size_histogram": _jsonable(blocker_size_histogram),
            "blocking_wall_occurrence_histogram": _jsonable(
                blocker_wall_histogram
            ),
            "production_rejection_witnesses_sha256": _canonical_sha256(
                _jsonable(rejection_witnesses)
            ),
            "accepted_root_indices": list(new_indices),
            "accepted_root_count": accepted_count,
            "accepted_root_order": "production_weight_recursion_order",
            "accepted_roots": _jsonable(level.roots),
            "accepted_roots_sha256": _canonical_sha256(_jsonable(level.roots)),
            "first_accepted_root": (
                _jsonable(level.roots[0]) if level.roots else None
            ),
            "last_accepted_root": (
                _jsonable(level.roots[-1]) if level.roots else None
            ),
            "pairings_to_prior_walls_sha256": _canonical_sha256(
                _jsonable(prior_pairings)
            ),
            "all_pairings_to_prior_walls_nonpositive": all(
                entry.sign(0) <= 0
                for row in prior_pairings
                for entry in row
            ),
            "same_height_pairings_sha256": _canonical_sha256(
                _jsonable(same_height_pairings)
            ),
            "same_height_pair_count": len(same_height_pairings),
            "same_height_negative_pair_count": (
                extension.same_height_negative_pair_count
            ),
            "same_height_zero_pair_count": extension.same_height_zero_pair_count,
            "all_distinct_same_height_pairings_nonpositive": all(
                pairing.sign(0) <= 0 for pairing in same_height_pairings
            ),
        },
        "independent_block_factorized_replay": {
            "block_partition": _jsonable(replay.block_partition),
            "component_box_counts": list(replay.component_box_counts),
            "component_bounded_counts": list(replay.component_bounded_counts),
            "norm_join_probe_count": replay.norm_join_probe_count,
            "exact_weight_solutions": replay.exact_weight_solutions,
            "nonintegral_reconstructions": replay.nonintegral_reconstructions,
            "rejected_root_count": replay.rejected_root_count,
            "accepted_root_count": len(replay.accepted_roots),
            "same_height_negative_pair_count": (
                replay.same_height_negative_pair_count
            ),
            "same_height_zero_pair_count": replay.same_height_zero_pair_count,
            "all_prefix_pairings_nonpositive": (
                replay.all_prefix_pairings_nonpositive
            ),
            "prefix_wall_count": replay.prefix_wall_count,
            "prefix_identity_inertia": _jsonable(
                replay.prefix_identity_inertia
            ),
            "prefix_conjugate_inertia": _jsonable(
                replay.prefix_conjugate_inertia
            ),
            "edge_indices": list(replay.edge_indices),
            "finite_endpoint_root_indices": list(
                replay.finite_endpoint_root_indices
            ),
            "new_finite_endpoint_root_indices": list(
                replay.new_finite_endpoint_root_indices
            ),
            "endpoint_index_semantics": "certified_production_prefix_order",
            "native_replay_finite_endpoint_root_indices": list(
                replay.native_finite_endpoint_root_indices
            ),
            "native_replay_new_finite_endpoint_root_indices": list(
                replay.native_new_finite_endpoint_root_indices
            ),
            "next_admissible_control_coordinate": _jsonable(
                replay.next_admissible_control_coordinate
            ),
            "all_reconstructed_roots_sha256": (
                replay.all_reconstructed_roots_sha256
            ),
            "rejection_witnesses_sha256": replay.rejection_witnesses_sha256,
            "production_accepted_root_set_matches": (
                replay.production_accepted_root_set_matches
            ),
            "production_rejection_witness_set_matches": (
                replay.production_rejection_witness_set_matches
            ),
            "replay_sha256": replay.replay_sha256,
            "expected_replay_sha256": EXPECTED_REPLAY_SHA256,
        },
        "prior_wall_count": extension.prior_wall_count,
        "wall_count": len(extension.roots),
        "wall_gram_sha256": _canonical_sha256(_jsonable(extension.wall_gram)),
        "wall_identity_inertia": _jsonable(extension.wall_identity_inertia),
        "wall_conjugate_inertia": _jsonable(extension.wall_conjugate_inertia),
        "wall_inertia_rank_witness": {
            "spanning_wall_indices": list(range(10)),
            "identity_inertia": {"positive": 9, "negative": 1, "zero": 0},
            "conjugate_inertia": {"positive": 10, "negative": 0, "zero": 0},
            "ambient_rank_upper_bound": 10,
        },
        "target_edge_audit": {
            "edge_indices": list(extension.edge_indices),
            "known_finite_endpoint_root_indices": list(
                extension.finite_endpoint_root_indices
            ),
            "new_finite_endpoint_root_indices": list(
                extension.new_finite_endpoint_root_indices
            ),
            "new_root_edge_pairings_sha256": _canonical_sha256(
                _jsonable(extension.new_root_edge_pairings)
            ),
            "new_root_edge_extension_determinants_sha256": _canonical_sha256(
                _jsonable(extension.new_root_edge_extension_determinants)
            ),
            "new_root_edge_extension_determinant_histogram": _jsonable(
                _silver_histogram(extension.new_root_edge_extension_determinants)
            ),
            "new_root_edge_extension_identity_inertia_histogram": _jsonable(
                _inertia_histogram(
                    extension.new_root_edge_extension_identity_inertias
                )
            ),
            "new_root_edge_extension_conjugate_inertia_histogram": _jsonable(
                _inertia_histogram(
                    extension.new_root_edge_extension_conjugate_inertias
                )
            ),
            "ideal_endpoint_audit_status": extension.ideal_endpoint_audit_status,
            "second_endpoint_status": extension.second_endpoint_status,
        },
        "endpoint_completion_criterion": {
            "finite": (
                "A distinct accepted wall r for which the principal Gram "
                "matrix on E union {r} is positive definite of rank nine."
            ),
            "ideal": (
                "A full componentwise-affine zero-neighbourhood of total "
                "parabolic rank eight incident to E, checked against the "
                "complete simple-root list."
            ),
            "terminal": (
                "After simple-root completeness, every elliptic rank-eight "
                "edge must have exactly two finite-or-ideal endpoints and "
                "every connected affine subdiagram must pass Vinberg's full "
                "neighbourhood criterion."
            ),
        },
        "global_diagram_audit": (
            f"not_run_for_the_incomplete_{len(extension.roots)}_wall_prefix"
        ),
        "stopping_status": extension.stopping_status,
        "scope": {
            "proved": [
                "Every norm-two Vinberg candidate at x=4+3sqrt(2) was exhausted exactly.",
                f"Of {level.exact_weight_solutions} exact weights, {integral_count} reconstruct integrally; {level.previous_wall_rejections} are blocked and {accepted_count} survive.",
                f"All {len(same_height_pairings)} distinct same-height pairings are nonpositive, so all {accepted_count} accepted roots enter simultaneously.",
                f"All {len(extension.roots)} roots in the current prefix are archived literally in order; the final {accepted_count} are the accepted shell roots and their set is reproduced by the block-factorized replay.",
                f"The targeted A2+F4+B2 finite endpoint set was recomputed through the {len(extension.roots)}-wall prefix.",
                "The next admissible control coordinate is x=5+3sqrt(2).",
            ],
            "not_proved": [
                "Any ideal endpoint for the targeted edge before the simple-root list is complete.",
                "A full parabolic-neighbourhood audit for the incomplete wall prefix.",
                "Simple-root completeness or any global Vinberg stopping criterion.",
                "Compactness, finite covolume, or infinite covolume of the eventual chamber.",
                "Exhaustion of primitive reflective vectors of norms other than two.",
            ],
        },
    }
    manifest["certificate_sha256"] = _canonical_sha256(manifest)
    return manifest


def verify_manifest(manifest: Mapping[str, object]) -> dict[str, object]:
    expected = build_manifest()
    if manifest != expected:
        raise ValueError(
            "the supplied x=4+3sqrt(2) silver K3 manifest does not reproduce"
        )
    digest = manifest.get("certificate_sha256")
    if not isinstance(digest, str):
        raise ValueError("the x=4+3sqrt(2) silver K3 digest is missing")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and digest != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("the x=4+3sqrt(2) silver K3 certificate is not frozen")
    new_level = manifest["new_level"]
    target = manifest["target_edge_audit"]
    replay = manifest["independent_block_factorized_replay"]
    if not all(isinstance(item, Mapping) for item in (new_level, target, replay)):
        raise ValueError("the x=4+3sqrt(2) silver K3 summary is malformed")
    return _manifest_summary(manifest)


def _manifest_summary(manifest: Mapping[str, object]) -> dict[str, object]:
    digest = manifest.get("certificate_sha256")
    if not isinstance(digest, str):
        raise ValueError("the x=4+3sqrt(2) silver K3 digest is missing")
    new_level = manifest["new_level"]
    target = manifest["target_edge_audit"]
    replay = manifest["independent_block_factorized_replay"]
    if not all(isinstance(item, Mapping) for item in (new_level, target, replay)):
        raise ValueError("the x=4+3sqrt(2) silver K3 summary is malformed")
    return {
        "certificate_sha256": digest,
        "exact_weight_solutions": new_level["exact_weight_solutions"],
        "accepted_root_count": new_level["accepted_root_count"],
        "wall_count": manifest["wall_count"],
        "known_finite_endpoint_root_indices": target[
            "known_finite_endpoint_root_indices"
        ],
        "independent_replay_sha256": replay["replay_sha256"],
        "next_admissible_control_coordinate": manifest[
            "next_admissible_control_coordinate"
        ],
        "stopping_status": manifest["stopping_status"],
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--recompute", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    recomputed = args.recompute or args.write
    manifest = (
        build_manifest()
        if recomputed
        else json.loads(ARCHIVE.read_text(encoding="utf-8"))
    )
    if args.write:
        ARCHIVE.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if recomputed:
        digest = manifest.get("certificate_sha256")
        if not isinstance(digest, str):
            raise ValueError("the recomputed x=4+3sqrt(2) digest is missing")
        if (
            EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
            and digest != EXPECTED_CERTIFICATE_SHA256
        ):
            raise ValueError("the recomputed x=4+3sqrt(2) certificate changed")
        result = _manifest_summary(manifest)
    else:
        result = verify_manifest(manifest)
    if args.full:
        print(json.dumps(manifest, indent=2, sort_keys=True))
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
