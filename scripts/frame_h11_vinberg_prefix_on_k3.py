"""Transport an exact H11 Vinberg prefix to the discriminant-11 K3.

The bound-20 transcript contains 95 certified simple norm-two walls and may
be followed by a resumable exact extension transcript.  This script
constructs a strict point in that chamber, reverses to
the time component containing the K3 ample class, and folds it through exact
field reflections until Shimada's complete nef atlas locates it in the
selected eigen-nef section.  The accumulated integral field isometry then
transports every supplied simple wall to paired rational-curve walls of the
enhanced K3 surface.

This proves a facet *lower bound* for the K3 eigen-nef chamber.  Neither a
bounded base transcript nor an extension proves that its wall list is
complete.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.geometric.golden_halfspaces import golden_quadratic_value
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse, _mat_vec
from coxeter_polytopes.search.h11_vinberg import D, _standard_e8_gram, h11_basis_gram
from coxeter_polytopes.verify.chamber import verify_hyperideal_extreme_ray
from scripts.diagnose_h11_vinberg_prefix import HYPERIDEAL_WITNESS
from scripts.explore_shimada_disc11_global_eigen_nef import (
    _ray_key,
    explore_global_cell,
)
from scripts.fold_h11_seed_via_shimada import (
    _field_pair,
    _isometry_verified,
    _row_product,
    fold_point_to_eigen_nef,
)
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas
from scripts.shimada_disc11_conjugate import (
    AMPLE_FIELD_COEFFICIENTS,
    DISC11_SOURCE_SHA256,
)


DEFAULT_TRANSCRIPT = Path("data/candidates/h11_vinberg_prefix_bound20.json")


def _canonical_sha(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _golden_vector(value) -> tuple[Golden, ...]:
    return tuple(Golden.from_json(entry) for entry in value)


def _transcript_roots(transcript: dict[str, object]) -> tuple[tuple[Golden, ...], ...]:
    height_zero = transcript["height_zero"]
    levels = transcript["levels"]
    assert isinstance(height_zero, dict) and isinstance(levels, list)
    roots = [
        _golden_vector(value) for value in height_zero["simple_roots"]
    ]
    for level in levels:
        assert isinstance(level, dict)
        roots.extend(_golden_vector(value) for value in level["roots"])
    return tuple(roots)


def _extension_roots(extension: dict[str, object]) -> tuple[tuple[Golden, ...], ...]:
    scan = extension["extension_scan"]
    assert isinstance(scan, dict)
    levels = scan["levels"]
    assert isinstance(levels, list)
    roots = []
    for level in levels:
        assert isinstance(level, dict)
        roots.extend(_golden_vector(value) for value in level["roots"])
    return tuple(roots)


def _vinberg_strict_point() -> tuple[Golden, ...]:
    inverse_e8 = _fraction_inverse(_standard_e8_gram())
    rho = _mat_vec(inverse_e8, [Golden(1)] * 8)
    control = (Golden(2), Golden(-1), *(Golden(0) for _ in range(8)))
    inward = (Golden(0), Golden(1), *rho)
    # The archived Vinberg roots use the chamber around +control.  The K3
    # ample point lies in the opposite time component, so reverse the strict
    # perturbation and reverse every simple-root orientation below.
    return tuple(-(control[i] + inward[i] / 10) for i in range(10))


def _positive_functional(root: tuple[Golden, ...]) -> tuple[Golden, ...]:
    gram = h11_basis_gram()
    return tuple(
        -sum((gram[i][j] * root[j] for j in range(10)), Golden(0))
        for i in range(10)
    )


def frame_vinberg_prefix(
    disc11: Path,
    transcript_path: Path = DEFAULT_TRANSCRIPT,
    extension_path: Path | None = None,
) -> dict[str, object]:
    disc_digest = hashlib.sha256(disc11.read_bytes()).hexdigest()
    if disc_digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")
    transcript_bytes = transcript_path.read_bytes()
    transcript_digest = hashlib.sha256(transcript_bytes).hexdigest()
    transcript = json.loads(transcript_bytes)
    base_roots = _transcript_roots(transcript)
    if len(base_roots) != 95:
        raise ValueError("the bound-20 H11 transcript no longer has 95 walls")
    extension_digest = None
    extension = None
    added_roots: tuple[tuple[Golden, ...], ...] = ()
    if extension_path is not None:
        extension_bytes = extension_path.read_bytes()
        extension_digest = hashlib.sha256(extension_bytes).hexdigest()
        extension = json.loads(extension_bytes)
        if extension.get("transcript_kind") != "exact_norm_two_vinberg_height_extension":
            raise ValueError("unexpected H11 extension transcript kind")
        base = extension.get("base_prefix")
        if not isinstance(base, dict) or base.get("sha256") != transcript_digest:
            raise ValueError("H11 extension is not pinned to the supplied base prefix")
        added_roots = _extension_roots(extension)
        if extension.get("output", {}).get("new_wall_count") != len(added_roots):
            raise ValueError("H11 extension wall count is inconsistent")
    roots = (*base_roots, *added_roots)

    gram = h11_basis_gram()
    initial_point = _vinberg_strict_point()
    control_same_component = (
        Golden(-2),
        Golden(1),
        *(Golden(0) for _ in range(8)),
    )
    control_square = _field_pair(control_same_component, control_same_component)
    point_square = _field_pair(initial_point, initial_point)
    control_point_pair = _field_pair(control_same_component, initial_point)
    control_to_point_cosh_squared = (
        control_point_pair * control_point_pair / (control_square * point_square)
    )
    first_unprocessed_wall_cosh_squared_lower_bound = 1 + 100 * D
    point_before_unprocessed_walls = (
        first_unprocessed_wall_cosh_squared_lower_bound
        - control_to_point_cosh_squared
    ).sign() > 0
    initial_roots = tuple(tuple(-entry for entry in root) for root in roots)
    initial_strict = all(_field_pair(initial_point, root).sign() < 0 for root in initial_roots)
    initial_norms = all(_field_pair(root, root) == 2 for root in initial_roots)
    initial_acute = all(
        _field_pair(left, right).sign() <= 0
        for i, left in enumerate(initial_roots)
        for right in initial_roots[i + 1 :]
    )

    atlas = ShimadaDisc11Atlas(disc11)
    fold = fold_point_to_eigen_nef(
        atlas,
        initial_point,
        max_root_reflections=200,
        max_gallery=2000,
    )
    if fold.terminal.kind != "interior":
        raise RuntimeError("the Vinberg point did not fold into the K3 eigen-nef chamber")
    transformed_roots = tuple(
        _row_product(root, fold.row_action) for root in initial_roots
    )
    encoded_weyl_action = [
        [entry.to_json() for entry in row] for row in fold.row_action
    ]
    encoded_transformed_roots = [
        [entry.to_json() for entry in root] for root in transformed_roots
    ]
    transformed_strict = all(
        _field_pair(fold.point, root).sign() < 0 for root in transformed_roots
    )
    ample_strict = all(
        _field_pair(tuple(AMPLE_FIELD_COEFFICIENTS), root).sign() < 0
        for root in transformed_roots
    )
    transformed_integral = all(
        entry.is_integral for root in transformed_roots for entry in root
    )
    transformed_norms = all(
        _field_pair(root, root) == 2 for root in transformed_roots
    )
    distinct_rays = {
        _ray_key(_positive_functional(root)) for root in transformed_roots
    }

    # Cross-check the eleven walls discovered independently in the first ten
    # physical induced chambers.  They must occur among the transported
    # Vinberg prefix, but the prefix contains many walls whose Shimada cells
    # lie much farther away.
    local_traversal = explore_global_cell(disc11, max_states=10)
    local_rays = {
        tuple(
            Golden.from_json(entry)
            for entry in record["normalized_restricted_functional"]
        )
        for record in local_traversal["visible_root_rows"]
    }
    local_matches = local_rays & distinct_rays

    hyperideal = verify_hyperideal_extreme_ray(
        transformed_roots,
        gram,
        HYPERIDEAL_WITNESS,
        dimension=9,
    )
    hyperideal_pairings_consistent = (
        all(value.sign() <= 0 for value in hyperideal.wall_pairings)
        or all(value.sign() >= 0 for value in hyperideal.wall_pairings)
    )
    verified = (
        initial_strict
        and initial_norms
        and initial_acute
        and point_before_unprocessed_walls
        and _isometry_verified(fold.row_action)
        and all(entry.is_integral for row in fold.row_action for entry in row)
        and transformed_strict
        and ample_strict
        and transformed_integral
        and transformed_norms
        and len(distinct_rays) == len(roots)
        and local_matches == local_rays
        and hyperideal.norm.sign() > 0
        and hyperideal_pairings_consistent
    )
    return {
        "verified": verified,
        "disc11_sha256": disc_digest,
        "vinberg_transcript_sha256": transcript_digest,
        "vinberg_extension_sha256": extension_digest,
        "vinberg_identity_bound": transcript["control_coordinate_scan"]["identity_bound"],
        "extension_last_processed_control_coordinate": (
            extension["status"]["last_processed_control_coordinate"]
            if extension is not None
            else None
        ),
        "archived_simple_walls": len(roots),
        "distinct_transported_curve_wall_rays": len(distinct_rays),
        "initial_vinberg_point_strict": initial_strict,
        "initial_simple_roots_norm_two": initial_norms,
        "initial_simple_system_acute": initial_acute,
        "control_to_strict_point_cosh_squared": (
            control_to_point_cosh_squared.to_json()
        ),
        "first_unprocessed_wall_cosh_squared_lower_bound": (
            first_unprocessed_wall_cosh_squared_lower_bound.to_json()
        ),
        "strict_point_lies_before_every_unprocessed_wall": (
            point_before_unprocessed_walls
        ),
        "weyl_field_reflections": sum(
            record["result"] == "root_boundary" for record in fold.records
        ),
        "auxiliary_shimada_crossings": sum(
            record["auxiliary_gallery_length"] for record in fold.records
        ),
        "terminal_physical_representative": fold.terminal.state.index,
        "terminal_strict_interior": transformed_strict,
        "fixed_ample_strict_on_all_transported_walls": ample_strict,
        "weyl_action_integral_over_Zphi": all(
            entry.is_integral for row in fold.row_action for entry in row
        ),
        "weyl_action_isometry": _isometry_verified(fold.row_action),
        "weyl_action_sha256": _canonical_sha(encoded_weyl_action),
        "weyl_fold_records_sha256": _canonical_sha(list(fold.records)),
        "transported_roots_integral_norm_two": (
            transformed_integral and transformed_norms
        ),
        "transported_roots_sha256": _canonical_sha(encoded_transformed_roots),
        "independent_local_wall_rays": len(local_rays),
        "independent_local_rays_matched_in_prefix": len(local_matches),
        "prefix_hyperideal_witness": {
            "active_facets": list(hyperideal.active_facets),
            "norm": hyperideal.norm.to_json(),
            "wall_pairings_have_one_consistent_orientation": (
                hyperideal_pairings_consistent
            ),
        },
        "conclusion": (
            "the_discriminant_11_K3_eigen_nef_chamber_has_at_least_"
            f"{len(roots)}_paired_rational_curve_facets_but_prefix_"
            "completeness_is_open"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    parser.add_argument("--transcript", type=Path, default=DEFAULT_TRANSCRIPT)
    parser.add_argument("--extension", type=Path)
    args = parser.parse_args()
    print(
        json.dumps(
            frame_vinberg_prefix(args.disc11, args.transcript, args.extension),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
