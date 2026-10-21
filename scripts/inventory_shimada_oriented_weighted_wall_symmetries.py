"""Resumable exact oriented-wall symmetry inventory for Shimada atlases.

This is the determinant-parameterized successor to the original disc16
weighted-wall inventory.  It keeps three finite groups/actions separate:

* Shimada's archived integral chamber stabilizer ``Aut(D)``;
* the full automorphism group of the *oriented* exact wall Gram matrix; and
* the rational Lorentzian action lifted from that wall permutation group.

The wall rows in Shimada's files describe mirrors and are not globally
half-space oriented.  Each wall orbit is therefore oriented toward the
published Weyl side.  If the Weyl vector lies on the orbit, the opposite
side of the explicitly archived adjacent chamber fixes the sign.  Archived
``Aut(D)`` must then preserve the resulting oriented rows exactly.

Every nontrivial group is exhaustively subgroup-enumerated modulo inner
conjugacy.  Per-chamber certificates make the long census restartable.  This
script is only the symmetry/rank denominator: facet, Coxeter-angle, and
finite-volume gates are performed separately.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.certify_shimada_disc16_d5_full_subgroup_fishing import (
    _subgroup_conjugacy_classes,
)
from scripts.inventory_shimada_disc16_weighted_wall_symmetries import (
    _archived_stabilizer_fast,
    _linear_actions_from_graph_group,
    _regular_action,
    _sha256,
    _wall_action,
    _weighted_wall_graph_automorphisms_fast,
)
from scripts.search_apery_fermi_chamber_folds import (
    Permutation,
    Subgroup,
    _all_subgroups,
    _fixed_rank,
    _matrix_payload,
    _permutation_order,
)
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas, _row_product


ROOT = Path(__file__).resolve().parents[1]
AMBIENT_RANK = 20
MIN_FIXED_RANK = 11
CASE_SCHEMA = "coxeter.shimada_oriented_weighted_wall_chamber.v1"
INVENTORY_SCHEMA = "coxeter.shimada_oriented_weighted_wall_inventory.v1"
CONFIG_VERSION = 1

CONFIGS: dict[int, dict[str, Any]] = {
    11: {
        "source": Path("/private/tmp/threediscs/disc11.txt"),
        "source_sha256": (
            "fa39e5b59f4a4b2b1b7fe7632e5737d1734ee3eb7e7023ae26bf2b83387afb68"
        ),
        "chamber_count": 1_098,
        "gram_determinant": -11,
    },
    15: {
        "source": Path("/private/tmp/threediscs/disc15.txt"),
        "source_sha256": (
            "f2fdcf80bd980c363fbdb68957a85878284bdb5e7d75c4dac92604aefefc4081"
        ),
        "chamber_count": 2_051,
        "gram_determinant": -15,
    },
    16: {
        "source": Path("/private/tmp/threediscs/disc16.txt"),
        "source_sha256": (
            "a19b797d6c314bc85e1ac4991fe3b190fe29a18a4746103bc6fa9b05e26a7db6"
        ),
        "chamber_count": 4_539,
        "gram_determinant": -16,
    },
}


def _default_case_directory(discriminant: int) -> Path:
    return (
        ROOT
        / "data"
        / "enumerations"
        / f"shimada_disc{discriminant}_oriented_weighted_wall_chambers"
    )


def _default_output(discriminant: int) -> Path:
    return (
        ROOT
        / "data"
        / "enumerations"
        / f"shimada_disc{discriminant}_oriented_weighted_wall_symmetry_inventory.json"
    )


def _case_id(chamber_index: int) -> str:
    return f"D{chamber_index:04d}"


def _case_path(directory: Path, chamber_index: int) -> Path:
    return directory / f"{_case_id(chamber_index)}.json"


def _denominator_lcm(matrix: sp.Matrix) -> int:
    return int(sp.ilcm(*(sp.Rational(value).q for value in matrix)))


def _zero_orbit_orientation(
    atlas: ShimadaDisc11Atlas,
    chamber_index: int,
    orbit_index: int,
    representative_wall: tuple[int, ...],
    dual_gram: sp.Matrix,
) -> tuple[int, str]:
    """Orient a wall orbit whose published Weyl pairing is zero."""

    adjacency = atlas.chamber(chamber_index).adjacencies[orbit_index]
    if adjacency.kind == "_minustwo":
        return 1, "effective_minustwo_orientation_convention"
    if adjacency.kind == "_isom":
        adjacent_weyl = atlas.assignment(
            atlas._chamber_name(chamber_index, "adjacentweyl", orbit_index)
        )
    elif adjacency.kind in {"_backto", "_newcham"} and adjacency.target is not None:
        adjacent_weyl = atlas.chamber(adjacency.target).weyl
    else:
        raise ArithmeticError("a zero-Weyl wall orbit has no archived opposite side")
    adjacent_dual = sp.Matrix(
        1,
        AMBIENT_RANK,
        _row_product(adjacent_weyl, atlas.projection_to_s_dual),
    )
    adjacent_pairing = sp.Rational(
        (adjacent_dual * dual_gram * sp.Matrix(representative_wall))[0]
    )
    if adjacent_pairing == 0:
        raise ArithmeticError("both archived Weyl vectors lie on the wall orbit")
    return (
        -1 if adjacent_pairing > 0 else 1,
        f"opposite_adjacent_{adjacency.kind}_Weyl_side",
    )


def _oriented_walls(
    atlas: ShimadaDisc11Atlas,
    chamber_index: int,
    dual_gram: sp.Matrix,
) -> tuple[tuple[tuple[int, ...], ...], dict[str, Any]]:
    """Return exact inward-oriented rows, preserving each archived orbit."""

    chamber = atlas.chamber(chamber_index)
    published_weyl_dual = sp.Matrix(
        1,
        AMBIENT_RANK,
        _row_product(chamber.weyl, atlas.projection_to_s_dual),
    )
    oriented: list[tuple[int, ...]] = []
    orbit_records = []
    sign_histogram: Counter[int] = Counter()
    method_histogram: Counter[str] = Counter()
    zero_orbit_count = 0
    for orbit_index, orbit in enumerate(chamber.wall_orbits):
        pairings = tuple(
            sp.Rational((published_weyl_dual * dual_gram * sp.Matrix(wall))[0])
            for wall in orbit
        )
        signs = {0 if value == 0 else 1 if value > 0 else -1 for value in pairings}
        if len(signs) != 1:
            raise ArithmeticError("a wall orbit has mixed published Weyl-side signs")
        raw_sign = next(iter(signs))
        if raw_sign > 0:
            sign, method = 1, "positive_published_Weyl_pairing"
        elif raw_sign < 0:
            sign, method = -1, "negative_published_Weyl_pairing_reversed"
        else:
            zero_orbit_count += 1
            sign, method = _zero_orbit_orientation(
                atlas, chamber_index, orbit_index, orbit[0], dual_gram
            )
        signed_orbit = tuple(
            tuple(sign * int(value) for value in wall) for wall in orbit
        )
        oriented.extend(signed_orbit)
        sign_histogram[sign] += len(orbit)
        method_histogram[method] += len(orbit)
        orbit_records.append(
            {
                "orbit_index": orbit_index,
                "size": len(orbit),
                "orientation_sign": sign,
                "orientation_method": method,
                "published_Weyl_pairing": str(pairings[0]),
                "adjacency_kind": chamber.adjacencies[orbit_index].kind,
            }
        )
    if len(set(oriented)) != len(oriented):
        raise ArithmeticError("oriented wall stock contains duplicate rows")
    return tuple(oriented), {
        "wall_count": len(oriented),
        "wall_orbit_count": len(chamber.wall_orbits),
        "zero_Weyl_orbit_count": zero_orbit_count,
        "orientation_sign_histogram": {
            str(key): value for key, value in sorted(sign_histogram.items())
        },
        "orientation_method_histogram": dict(sorted(method_histogram.items())),
        "orbit_records": orbit_records,
    }


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


def _subgroup_classes(
    wall_group: tuple[Permutation, ...],
    actions: dict[Permutation, sp.Matrix],
    *,
    containing_group: set[Permutation] | None = None,
) -> tuple[int, list[dict[str, Any]]]:
    regular_group, regular_lookup = _regular_action(wall_group)
    subgroups = _all_subgroups(regular_group)
    classes = _subgroup_conjugacy_classes(regular_group, subgroups)
    rows = []
    for class_index, conjugacy_class in enumerate(classes):
        representative_index = conjugacy_class[0]
        wall_subgroup: Subgroup = frozenset(
            regular_lookup[element] for element in subgroups[representative_index]
        )
        row = _subgroup_row(wall_subgroup, actions)
        if containing_group is not None:
            row["representative_contained_in_archived_AutD"] = set(
                wall_subgroup
            ).issubset(containing_group)
        rows.append(
            {
                "class_index": class_index,
                "representative_subgroup_index": representative_index,
                "subgroup_indices": list(conjugacy_class),
                "conjugacy_class_size": len(conjugacy_class),
                **row,
            }
        )
    return len(subgroups), rows


def _action_records(
    graph_group: tuple[Permutation, ...],
    actions: dict[Permutation, sp.Matrix],
    ambient_gram: sp.Matrix,
    dual_gram: sp.Matrix,
) -> list[dict[str, Any]]:
    records = []
    for group_index, permutation in enumerate(graph_group):
        dual = actions[permutation]
        primal = ambient_gram * dual * dual_gram
        records.append(
            {
                "group_index": group_index,
                "permutation": list(permutation),
                "order": _permutation_order(permutation),
                "dual_row_action": _matrix_payload(dual),
                "dual_action_denominator_lcm": _denominator_lcm(dual),
                "primal_row_action": _matrix_payload(primal),
                "primal_action_denominator_lcm": _denominator_lcm(primal),
                "primal_action_integral": all(value.is_Integer for value in primal),
            }
        )
    return records


def inventory_chamber(
    atlas: ShimadaDisc11Atlas,
    discriminant: int,
    source_hash: str,
    chamber_index: int,
    ambient_gram: sp.Matrix,
    dual_gram: sp.Matrix,
    integral_scaled_dual_gram: sp.Matrix,
) -> dict[str, Any]:
    chamber = atlas.chamber(chamber_index)
    walls, orientation = _oriented_walls(atlas, chamber_index, dual_gram)
    graph_group = _weighted_wall_graph_automorphisms_fast(
        integral_scaled_dual_gram, walls
    )
    actions = _linear_actions_from_graph_group(
        walls, graph_group, dual_gram, ambient_gram
    )

    archived_primal = _archived_stabilizer_fast(
        atlas, chamber_index, chamber.order_aut
    )
    archived_dual = tuple(sp.Matrix(atlas.dual_action(matrix)) for matrix in archived_primal)
    archived_group = tuple(sorted(_wall_action(walls, matrix) for matrix in archived_dual))
    if not set(archived_group).issubset(graph_group):
        raise ArithmeticError("oriented Aut(D) is not contained in the wall group")
    if any(actions[permutation] not in archived_dual for permutation in archived_group):
        raise ArithmeticError("an archived wall permutation lifted to another action")

    if len(archived_group) > 1:
        archived_subgroup_count, archived_classes = _subgroup_classes(
            archived_group, actions
        )
    else:
        archived_subgroup_count, archived_classes = 0, []
    groups_equal = set(graph_group) == set(archived_group)
    if len(graph_group) > 1 and groups_equal:
        # The regular actions are identical because both tuples are sorted.
        # Re-enumerating all subgroups is especially expensive for the
        # order-48 disc16 chamber and produces byte-for-byte the same class
        # data apart from the explicit containment flag.
        weighted_subgroup_count = archived_subgroup_count
        weighted_classes = [
            record | {"representative_contained_in_archived_AutD": True}
            for record in archived_classes
        ]
    elif len(graph_group) > 1:
        weighted_subgroup_count, weighted_classes = _subgroup_classes(
            graph_group, actions, containing_group=set(archived_group)
        )
    else:
        weighted_subgroup_count, weighted_classes = 0, []

    return {
        "schema": CASE_SCHEMA,
        "configuration_version": CONFIG_VERSION,
        "source": {
            "discriminant": discriminant,
            "sha256": source_hash,
            "ambient_rank": AMBIENT_RANK,
        },
        "chamber_index": chamber_index,
        "wall_count": len(walls),
        "orientation": orientation,
        "archived_AutD": {
            "order": len(archived_group),
            "wall_permutations": [list(element) for element in archived_group],
            "subgroup_count": archived_subgroup_count,
            "subgroup_conjugacy_class_count": len(archived_classes),
            "classes": archived_classes,
        },
        "oriented_weighted_wall_group": {
            "order": len(graph_group),
            "strictly_extends_archived_AutD": not groups_equal,
            "all_actions": _action_records(
                graph_group, actions, ambient_gram, dual_gram
            ),
            "subgroup_count": weighted_subgroup_count,
            "subgroup_conjugacy_class_count": len(weighted_classes),
            "classes": weighted_classes,
        },
        "checks": {
            "all_wall_orbits_uniformly_oriented": True,
            "archived_AutD_preserves_oriented_walls": True,
            "archived_AutD_contained_in_full_oriented_wall_group": True,
            "all_weighted_wall_permutations_lifted_exactly": True,
            "rational_and_integral_actions_recorded_separately": True,
        },
    }


def _valid_case(
    path: Path,
    discriminant: int,
    source_hash: str,
    chamber_index: int,
) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not (
        payload.get("schema") == CASE_SCHEMA
        and payload.get("configuration_version") == CONFIG_VERSION
        and payload.get("source", {}).get("discriminant") == discriminant
        and payload.get("source", {}).get("sha256") == source_hash
        and payload.get("chamber_index") == chamber_index
        and all(payload.get("checks", {}).values())
    ):
        return None
    return payload


def _aggregate(
    discriminant: int,
    source: Path,
    source_hash: str,
    atlas: ShimadaDisc11Atlas,
    case_directory: Path,
) -> dict[str, Any]:
    cases = []
    pending = []
    for chamber_index in range(atlas.chamber_count):
        payload = _valid_case(
            _case_path(case_directory, chamber_index),
            discriminant,
            source_hash,
            chamber_index,
        )
        if payload is None:
            pending.append(chamber_index)
        else:
            cases.append(payload)

    archived_order_histogram: Counter[int] = Counter()
    weighted_order_histogram: Counter[int] = Counter()
    proper_archived_high_histogram: Counter[int] = Counter()
    proper_weighted_high_histogram: Counter[int] = Counter()
    archived_total_subgroups = 0
    archived_total_classes = 0
    weighted_total_subgroups = 0
    weighted_total_classes = 0
    strict_extensions = 0
    zero_orbits = 0
    rational_action_denominators: Counter[int] = Counter()
    for row in cases:
        archived = row["archived_AutD"]
        weighted = row["oriented_weighted_wall_group"]
        archived_order_histogram[archived["order"]] += 1
        weighted_order_histogram[weighted["order"]] += 1
        strict_extensions += int(weighted["strictly_extends_archived_AutD"])
        zero_orbits += row["orientation"]["zero_Weyl_orbit_count"]
        archived_total_subgroups += archived["subgroup_count"]
        archived_total_classes += archived["subgroup_conjugacy_class_count"]
        weighted_total_subgroups += weighted["subgroup_count"]
        weighted_total_classes += weighted["subgroup_conjugacy_class_count"]
        for record in weighted["all_actions"]:
            rational_action_denominators[record["primal_action_denominator_lcm"]] += 1
        for record in archived["classes"]:
            if record["order"] > 1 and record["fixed_rank"] >= MIN_FIXED_RANK:
                proper_archived_high_histogram[record["hyperbolic_dimension"]] += 1
        for record in weighted["classes"]:
            if record["order"] > 1 and record["fixed_rank"] >= MIN_FIXED_RANK:
                proper_weighted_high_histogram[record["hyperbolic_dimension"]] += 1

    complete = not pending
    return {
        "schema": INVENTORY_SCHEMA,
        "configuration_version": CONFIG_VERSION,
        "source": {
            "path": str(source),
            "sha256": source_hash,
            "discriminant": discriminant,
            "chamber_count": atlas.chamber_count,
            "ambient_rank": AMBIENT_RANK,
            "gram_determinant": int(sp.Matrix(atlas.gram).det()),
        },
        "scope": {
            "all_published_chambers": complete,
            "wall_halfspaces_exactly_oriented_before_graph_automorphisms": True,
            "archived_AutD_and_full_oriented_wall_groups_both_inventoried": True,
            "all_subgroups_modulo_inner_conjugacy": True,
            "rational_Lorentz_lifts_allowed": True,
            "minimum_fixed_rank_for_high_section_ledger": MIN_FIXED_RANK,
            "facet_angle_gate_performed": False,
            "finite_volume_gate_performed": False,
            "resumable_per_chamber_certificates": True,
        },
        "summary": {
            "completed_chamber_count": len(cases),
            "pending_chamber_count": len(pending),
            "archived_AutD_order_histogram": {
                str(key): value for key, value in sorted(archived_order_histogram.items())
            },
            "oriented_weighted_group_order_histogram": {
                str(key): value for key, value in sorted(weighted_order_histogram.items())
            },
            "weighted_group_strict_extension_count": strict_extensions,
            "zero_published_Weyl_wall_orbit_count": zero_orbits,
            "archived_total_subgroup_count_nontrivial_groups": archived_total_subgroups,
            "archived_total_subgroup_conjugacy_class_count": archived_total_classes,
            "weighted_total_subgroup_count_nontrivial_groups": weighted_total_subgroups,
            "weighted_total_subgroup_conjugacy_class_count": weighted_total_classes,
            "proper_archived_fixed_rank_at_least_11_dimension_histogram": {
                str(key): value
                for key, value in sorted(proper_archived_high_histogram.items())
            },
            "proper_weighted_fixed_rank_at_least_11_dimension_histogram": {
                str(key): value
                for key, value in sorted(proper_weighted_high_histogram.items())
            },
            "weighted_primal_action_denominator_histogram": {
                str(key): value for key, value in sorted(rational_action_denominators.items())
            },
        },
        "pending_chamber_indices": pending,
        "chambers": cases,
        "complete": complete,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--discriminant", type=int, choices=sorted(CONFIGS), default=15)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--case-directory", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-chambers", type=int)
    parser.add_argument("--recompute", action="store_true")
    parser.add_argument("--progress-every", type=int, default=25)
    arguments = parser.parse_args()

    config = CONFIGS[arguments.discriminant]
    source = arguments.source or config["source"]
    case_directory = arguments.case_directory or _default_case_directory(
        arguments.discriminant
    )
    output = arguments.output or _default_output(arguments.discriminant)
    source_hash = _sha256(source)
    if source_hash != config["source_sha256"]:
        raise SystemExit(
            f"disc{arguments.discriminant} source SHA-256 changed: {source_hash}"
        )
    atlas = ShimadaDisc11Atlas(source)
    ambient_gram = sp.Matrix(atlas.gram)
    if not (
        atlas.chamber_count == config["chamber_count"]
        and int(ambient_gram.det()) == config["gram_determinant"]
    ):
        raise SystemExit("the pinned atlas chamber count or Gram determinant changed")
    dual_gram = ambient_gram.inv()
    dual_scale = _denominator_lcm(dual_gram)
    integral_scaled_dual_gram = dual_scale * dual_gram
    if any(not value.is_Integer for value in integral_scaled_dual_gram):
        raise ArithmeticError("failed to clear dual Gram denominators")

    case_directory.mkdir(parents=True, exist_ok=True)
    processed = 0
    for chamber_index in range(atlas.chamber_count):
        path = _case_path(case_directory, chamber_index)
        if not arguments.recompute and _valid_case(
            path, arguments.discriminant, source_hash, chamber_index
        ) is not None:
            continue
        if arguments.max_chambers is not None and processed >= arguments.max_chambers:
            break
        payload = inventory_chamber(
            atlas,
            arguments.discriminant,
            source_hash,
            chamber_index,
            ambient_gram,
            dual_gram,
            integral_scaled_dual_gram,
        )
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        processed += 1
        weighted = payload["oriented_weighted_wall_group"]
        if (
            arguments.progress_every
            and (
                processed % arguments.progress_every == 0
                or weighted["strictly_extends_archived_AutD"]
            )
        ):
            print(
                json.dumps(
                    {
                        "chamber": chamber_index,
                        "processed_this_run": processed,
                        "wall_count": payload["wall_count"],
                        "archived_order": payload["archived_AutD"]["order"],
                        "oriented_weighted_order": weighted["order"],
                        "strict_extension": weighted[
                            "strictly_extends_archived_AutD"
                        ],
                        "proper_high_classes": sum(
                            record["order"] > 1
                            and record["fixed_rank"] >= MIN_FIXED_RANK
                            for record in weighted["classes"]
                        ),
                    },
                    sort_keys=True,
                ),
                flush=True,
            )

    inventory = _aggregate(
        arguments.discriminant,
        source,
        source_hash,
        atlas,
        case_directory,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(inventory, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(inventory["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
