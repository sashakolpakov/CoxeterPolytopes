#!/usr/bin/env python3
"""Certify the 56-wall silver Vinberg shell ``x=2+2sqrt(2)`` exactly."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    _silver_pair,
    run_silver_k3_vinberg_two_plus_two_sqrt2_extension,
    silver_k3_alexeev_field_gram,
)
from scripts.independent_silver_k3_vinberg_2_plus_2sqrt2 import (
    EXPECTED_REPLAY_SHA256,
    audit_reverse_silver_two_plus_two_shell,
)


PRIOR_CERTIFICATE_SHA256 = (
    "172301cae3723d15b80425c7303ab12638574a84e69a5cc1ba28dd19179cf4a1"
)
EXPECTED_CERTIFICATE_SHA256 = "4cbc5d390cccf160cac882918c229900b3b5540331d12db82d1686ab2beb6e49"
ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data/candidates/silver_k3_norm2_vinberg_extension_2_plus_2sqrt2.json"
)


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
    extension = run_silver_k3_vinberg_two_plus_two_sqrt2_extension()
    level = extension.new_level
    reverse = audit_reverse_silver_two_plus_two_shell()
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

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "silver_k3_intrinsic_norm_two_vinberg_fifty_six_wall_shell",
        "claim": (
            "At x=2+2sqrt(2), exact two-embedding enumeration gives 277 "
            "weight solutions. Of the 148 integral roots, 92 are blocked by "
            "lower-height walls and 56 mutually acute roots enter the current "
            "Vinberg prefix. None is a second finite endpoint of the targeted "
            "A2+F4+B2 edge."
        ),
        "terminology": (
            "The 56 roots are simple relative to the complete accepted "
            "lower-height prefix and their simultaneous same-height cohort; "
            "the 93-wall list is not claimed to be the complete chamber."
        ),
        "prior_three_wall_shell_certificate_sha256": PRIOR_CERTIFICATE_SHA256,
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
            "integral_reconstructions": (
                level.exact_weight_solutions - level.nonintegral_reconstructions
            ),
            "previous_wall_rejections": level.previous_wall_rejections,
            "blocker_set_size_histogram": _jsonable(blocker_size_histogram),
            "blocking_wall_occurrence_histogram": _jsonable(
                blocker_wall_histogram
            ),
            "production_rejection_witnesses_sha256": _canonical_sha256(
                _jsonable(rejection_witnesses)
            ),
            "accepted_root_indices": list(new_indices),
            "accepted_root_count": len(level.roots),
            "accepted_roots_sha256": _canonical_sha256(_jsonable(level.roots)),
            "first_accepted_root": _jsonable(level.roots[0]),
            "last_accepted_root": _jsonable(level.roots[-1]),
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
        "independent_reverse_coordinate_replay": {
            "coordinate_permutation": list(reverse.coordinate_permutation),
            "visited_by_depth": list(reverse.visited_by_depth),
            "schur_pruned_by_depth": list(reverse.schur_pruned_by_depth),
            "exact_weight_solutions": reverse.exact_weight_solutions,
            "nonintegral_reconstructions": reverse.nonintegral_reconstructions,
            "rejected_root_count": reverse.rejected_root_count,
            "accepted_root_count": len(reverse.accepted_roots),
            "same_height_negative_pair_count": (
                reverse.same_height_negative_pair_count
            ),
            "same_height_zero_pair_count": reverse.same_height_zero_pair_count,
            "all_reconstructed_roots_sha256": (
                reverse.all_reconstructed_roots_sha256
            ),
            "rejection_witnesses_sha256": reverse.rejection_witnesses_sha256,
            "production_accepted_root_set_matches": (
                reverse.production_accepted_root_set_matches
            ),
            "production_rejection_witness_set_matches": (
                reverse.production_rejection_witness_set_matches
            ),
            "replay_sha256": reverse.replay_sha256,
            "expected_replay_sha256": EXPECTED_REPLAY_SHA256,
        },
        "prior_wall_count": extension.prior_wall_count,
        "wall_count": len(extension.roots),
        "wall_gram_sha256": _canonical_sha256(_jsonable(extension.wall_gram)),
        "wall_identity_inertia": _jsonable(extension.wall_identity_inertia),
        "wall_conjugate_inertia": _jsonable(extension.wall_conjugate_inertia),
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
            "new_root_edge_extension_determinants": _jsonable(
                extension.new_root_edge_extension_determinants
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
        "global_diagram_audit": "not_run_for_the_incomplete_93_wall_prefix",
        "stopping_status": extension.stopping_status,
        "scope": {
            "proved": [
                "Every norm-two Vinberg candidate at x=2+2sqrt(2) was exhausted exactly.",
                "Of 277 exact weight solutions, 129 reconstruct nonintegrally, 92 integral roots are blocked by lower-height walls, and 56 roots survive.",
                "All 1,540 distinct same-height pairings are nonpositive, so all 56 roots enter this bounded prefix simultaneously and independently of candidate order.",
                "None of the 56 roots gives a finite endpoint of the targeted A2+F4+B2 edge.",
                "The next admissible control coordinate is x=3+2sqrt(2).",
            ],
            "not_proved": [
                "A second finite or ideal endpoint for the targeted edge.",
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
            "the supplied x=2+2sqrt(2) silver K3 manifest does not reproduce"
        )
    digest = manifest.get("certificate_sha256")
    if not isinstance(digest, str):
        raise ValueError("the x=2+2sqrt(2) silver K3 digest is missing")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and digest != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("the x=2+2sqrt(2) silver K3 certificate is not frozen")
    new_level = manifest["new_level"]
    target = manifest["target_edge_audit"]
    reverse = manifest["independent_reverse_coordinate_replay"]
    if not all(isinstance(item, Mapping) for item in (new_level, target, reverse)):
        raise ValueError("the x=2+2sqrt(2) silver K3 summary is malformed")
    return {
        "certificate_sha256": digest,
        "exact_weight_solutions": new_level["exact_weight_solutions"],
        "accepted_root_count": new_level["accepted_root_count"],
        "wall_count": manifest["wall_count"],
        "known_finite_endpoint_root_indices": target[
            "known_finite_endpoint_root_indices"
        ],
        "independent_replay_sha256": reverse["replay_sha256"],
        "next_admissible_control_coordinate": manifest[
            "next_admissible_control_coordinate"
        ],
        "stopping_status": manifest["stopping_status"],
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args(argv)
    manifest = (
        build_manifest()
        if args.recompute
        else json.loads(ARCHIVE.read_text(encoding="utf-8"))
    )
    result = verify_manifest(manifest)
    if args.full:
        print(json.dumps(manifest, indent=2, sort_keys=True))
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
