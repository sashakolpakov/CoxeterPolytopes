"""Probe an exact hyperideal divisor direction and certify its first missing wall.

Starting from the fixed ample class ``h`` and the explicit spacelike extreme
ray ``v`` of the 181-wall V4 outer cone, take the exact dyadic point
``h + (377595/32) v`` immediately before the positive light root on that
line.  Shimada's induced-nef atlas locates a genuine ``_minustwo`` boundary.
The script reconstructs its norm-two golden root, closes it under V4, and
tests the two previously certified hyperideal witness orbits against the
enlarged wall system.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from scipy.optimize import linprog
from sympy import Matrix

from coxeter_polytopes.exact.coxeter import classify_root_pair
from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.geometric.golden_halfspaces import (
    golden_quadratic_value,
    inverse_golden_matrix,
)
from coxeter_polytopes.verify.chamber import (
    _dot,
    _full_rank_kernel,
    _mat_vec,
    _pair,
    _row_rank,
    verify_hyperideal_extreme_ray,
)
from scripts.analyze_shimada_disc11_v4_root_orbit import (
    DEFAULT_EXTENSION,
    EXPECTED_TRANSFORMED_ROOTS_SHA256,
    NEW_HYPERIDEAL_RAY,
    _canonical_sha,
    _encoded_vector,
    _known_hyperideal_orbit,
    _positive_canonical,
    _transported_roots,
    _v4,
)
from scripts.coupled_discriminant11_candidate import field_gram
from scripts.diagnose_h11_vinberg_prefix import HYPERIDEAL_WITNESS
from scripts.fold_h11_seed_via_shimada import (
    _row_product,
    _walk_until_boundary_or_interior,
)
from scripts.frame_h11_vinberg_prefix_on_k3 import DEFAULT_TRANSCRIPT
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas
from scripts.shimada_disc11_conjugate import (
    AMPLE_FIELD_COEFFICIENTS,
    DISC11_SOURCE_SHA256,
)


TARGET_T = Fraction(377595, 32)
EXPECTED_REPRESENTATIVE = 1
EXPECTED_BOUNDARY_LABEL = (4, 1)
EXPECTED_BOUNDARY_FUNCTIONAL = (
    Golden(-217, -351),
    Golden(51, 83),
    Golden(5, 9),
    Golden(-5, -9),
    Golden(-83, -133),
    Golden(0, -1),
    Golden(5, 8),
    Golden(-5, -8),
    Golden(88, 143),
    Golden(129, 208),
)
EXPECTED_NEW_ROOT = (
    Golden(-40, -65),
    Golden(4, 7),
    Golden(9, 14),
    Golden(22, 33),
    Golden(31, 47),
    Golden(14, 19),
    Golden(-4, -9),
    Golden(-19, -32),
    Golden(-37, -60),
    Golden(-24, -41),
)
EXPECTED_AMPLE_ROOT_PAIRING = Golden(-427, -712)
EXPECTED_ESCAPE_ROOT_PAIRING = Golden(-5, 4) / 11
EXPECTED_CROSSING_T = Golden(3275, 5268)
EXPECTED_FACE_SQUARE = Golden(62614, -38946) / 11
ARCHIVED_HYPERIDEAL_RAY = (
    Golden(-1),
    Golden(Fraction(1, 2)),
    Golden(0),
    Golden(0),
    Golden(0),
    Golden(0),
    Golden(19, -23) / 28,
    Golden(19, -23) / 28,
    Golden(-11, -11) / 28,
    Golden(-11, -11) / 28,
)
SECOND_TARGET_T = Fraction(284021, 32)
SECOND_EXPECTED_REPRESENTATIVE = 72
SECOND_EXPECTED_BOUNDARY_LABEL = (11, 0)
SECOND_EXPECTED_BOUNDARY_FUNCTIONAL = (
    Golden(-51, -83),
    Golden(0),
    Golden(0),
    Golden(0),
    Golden(-14, -24),
    Golden(-9, -13),
    Golden(9, 14),
    Golden(-5, -8),
    Golden(19, 31),
    Golden(28, 46),
)
SECOND_NEW_ROOT = (
    Golden(-8, -14),
    Golden(4, 7),
    Golden(0),
    Golden(0),
    Golden(0),
    Golden(-1),
    Golden(-6, -9),
    Golden(-6, -9),
    Golden(-9, -14),
    Golden(-9, -14),
)
SECOND_GALLERY_SHA256 = (
    "19c3e61399d8fb19b8499b46f0a8babc6acca25f9e890f5d6777370f0dfd2bae"
)
SECOND_EXPANDED_ROOTS_SHA256 = (
    "f5626ff344f74839fa551c5e13acf1bb77017ec3c05edf88042cda7583502e18"
)
CHECKPOINT_SCHEMA_VERSION = 1
FROZEN_TERMINAL_WALL_COUNT = 309
FROZEN_TERMINAL_ROOTS_SHA256 = (
    "d5b9030d49c16fd55f2186ff115ad23eb6afe99f489a62682049691cb069d197"
)
FROZEN_TERMINAL_PAIRINGS_SHA256 = (
    "2588b7add212ce29cc7df2021c255568d7524a754330d0176868d9a89dbd642a"
)
FROZEN_CHECKPOINT_PAYLOAD_SHA256 = (
    "133759a4077a3c4fa543fb259228c8b9d4f7b20512ea82adb0610f6de521a346"
)
FROZEN_TERMINAL_RAY = (
    Golden(-1),
    Golden(-164, 415) / 1021,
    Golden(-941, 196) / 1021,
    Golden(-840, -16) / 1021,
    Golden(-1781, 180) / 1021,
    Golden(-1781, 180) / 1021,
    Golden(-1781, 180) / 1021,
    Golden(-1781, 180) / 1021,
    Golden(-1781, 180) / 1021,
    Golden(-1680, -32) / 1021,
)
FROZEN_TERMINAL_NORM = Golden(5514638, -3408230) / 1042441
FROZEN_TERMINAL_ACTIVE_FACETS = (4, 5, 6, 8, 9, 11, 13, 249, 253)
FROZEN_DIRECTED_OUTPUT_HASHES = (
    "5e148e8580195ce3ef472e97180931d62f55077e6ded85feb835765a7722cd62",
    "d20812f443db2ad092bc63d1d92e355bd5216b5e72e968b79ed23f94fe02458a",
    "79e7ffe3024d3ebb6a3db093a1a63ef8b698a2890d8c0a094e459600cf8a8615",
    "84f82d8aa596ea08e02f9a687fb1dda344ca039abdcf797df7dfef4f0590b9fb",
    "f544665b3b5c05ccc3697accda1ea1333103797542310d6acc5920b02c0301fc",
    "62da600a4eb0abbf6924edc89e98d1902d4928442d2ace0c2915419441a62c24",
    "5971b49624313942c77d8b5cbf90c708a8626a8420fc60a811cbd1ece9b1d038",
    "572c5a9045c71fdcf1eeb7c7a3fdb57a625852247107b3ed3816a4f523eeb26f",
    "7b543ffc9e063117def8016a2fde13b7cf6a38b31deae03660b3a67ea769b2a6",
    "a8fbda524b50aebfa9b49ce064975a906bb8b388e4277b9064f8e07ca80350a5",
    "0422b751d2f0321ce1a51e34b35ac94956713aa7164abd6188d3a3cb931fa9b4",
    "9159f7c23226f13cf048d5143ae50f995513b1942160992c51a49531817cc426",
    "c9b0a43a5c914bf2d8a897067fb124371d5d81ec75d81498cb5119c54cd6ce4f",
    "fc6dd0977ea641ea1475a4b6460510c14991f70910fb06cb050f19aa3854b4d4",
    "e71d0356c4412d1df643b49b0a7c1262682b02a52ed278bcc4b163b3760948b9",
    "a55bb1f3ca300115bf0ed06c98466cc855bb2337356103e7660f94ebb55a66da",
    "ce1e435aa80028f96e89aa035afcc519cb1e1504cb5ad92be58aa8ab9eb7d14a",
    "b7906d1af3bd44d606340a6676c03c550c4f5d2338b8da81a661a0d5c4fef1ba",
    "a5115baa9b0d14f48d2c979ca4456358977d8ac3ff34816e2ab74028cd12d717",
    "34ce25f34a796da997d560292a6e40ff9e0434092bd0b97fe1cdf715ffb8cd44",
    "5c6d56a92901ce3fd399cf599967951c570bcf929dd556db39cf5e54560954d8",
    "6a1eb94cdbdd202a6edfbf87c822b989a84837ee7c62291a3484baf1b19afa74",
    "23ebb8d355324461aca1b92b4639fd8661dea9cc93ba5dc52d69f46d07a20b51",
    "f3de6773b07586cf049be885010ff1ac0963c1a8f3cc7c40e94f25ecdcec542f",
    "2a97ae8acfd0f2487659a437231499a2dacf646029a1d9993b951d01d4cd8d1e",
    "beb8f377aa7f7b58fd8a8c61ab6814fde438485cb75fcb7ff03e2c6262187104",
    "4c24a6858bb1b61d8f00ad074a0c69e789932895ff865dfb22207a146a5ed039",
    "47f1007e322908b95b52da39a8e25cd06d2e0d31ac36f8666828f0f781278cfc",
    "81cbbf829f5ae60606fcdc6dd22a3fed511579a6e9a89d90d2f487de3d14d65b",
    FROZEN_TERMINAL_ROOTS_SHA256,
)


def _decoded_vector(value: Sequence[dict[str, str]]) -> tuple[Golden, ...]:
    return tuple(Golden.from_json(entry) for entry in value)


def _roots_sha256(roots: Sequence[Sequence[Golden]]) -> str:
    return _canonical_sha([_encoded_vector(root) for root in roots])


def _closed_181(atlas: ShimadaDisc11Atlas, transcript: Path, extension: Path):
    roots97, frame = _transported_roots(atlas, transcript, extension)
    group = _v4(atlas)
    provenance: dict[tuple[Golden, ...], list[tuple[int, str]]] = {}
    for source, root in enumerate(roots97):
        for name, action in group:
            image = _row_product(root, action)
            provenance.setdefault(image, []).append((source, name))
    roots181 = tuple(provenance)
    if len(roots181) != 181:
        raise ArithmeticError("the certified V4 closure no longer has 181 roots")
    return roots97, roots181, provenance, group, frame


def _frozen_roots189(
    roots181: Sequence[tuple[Golden, ...]],
    group,
) -> tuple[tuple[Golden, ...], ...]:
    """Reconstruct the two certified directed V4 wall orbits after D181."""

    current = tuple(roots181)
    for representative in (EXPECTED_NEW_ROOT, SECOND_NEW_ROOT):
        orbit: dict[tuple[Golden, ...], None] = {}
        for _, action in group:
            orbit.setdefault(_row_product(representative, action), None)
        current = tuple((*current, *(item for item in orbit if item not in current)))
    if len(current) != 189 or _roots_sha256(current) != SECOND_EXPANDED_ROOTS_SHA256:
        raise ArithmeticError("the frozen 189-wall base changed")
    return current


def _pairing_record(
    ray: Sequence[Golden],
    roots: Sequence[Sequence[Golden]],
    new_indices: Sequence[int],
) -> dict[str, Any]:
    gram = tuple(tuple(entry for entry in row) for row in field_gram())
    pairings = tuple(_pair(ray, gram, root) for root in roots)
    positive = [index for index, value in enumerate(pairings) if value.sign() > 0]
    return {
        "ray": _encoded_vector(ray),
        "norm": _pair(ray, gram, ray).to_json(),
        "feasible_for_enlarged_outer_cone": not positive,
        "positive_cutting_facets": positive,
        "new_wall_orbit_pairings": [
            {
                "facet": index,
                "pairing": pairings[index].to_json(),
            }
            for index in new_indices
        ],
        "all_pairings_sha256": _canonical_sha(
            [value.to_json() for value in pairings]
        ),
    }


def _witness_orbit_test(
    base_ray: Sequence[Golden],
    group,
    roots: Sequence[Sequence[Golden]],
    new_indices: Sequence[int],
) -> dict[str, Any]:
    orbit: dict[tuple[Golden, ...], dict[str, Any]] = {}
    for name, action in group:
        image = _positive_canonical(_row_product(tuple(base_ray), action))
        record = _pairing_record(image, roots, new_indices)
        record["V4_image_by"] = name
        orbit.setdefault(image, record)
    records = list(orbit.values())
    return {
        "orbit_size": len(records),
        "surviving_rays": sum(
            record["feasible_for_enlarged_outer_cone"] for record in records
        ),
        "records": records,
    }


def _dyadic_before_light(
    ample: Sequence[Golden],
    ray: Sequence[Golden],
    gram: Sequence[Sequence[Golden]],
    denominator: int = 32,
) -> tuple[Fraction, Golden, Golden]:
    if denominator <= 0 or denominator & (denominator - 1):
        raise ValueError("the light bracket denominator must be a power of two")

    def square(numerator: int) -> Golden:
        parameter = Golden(Fraction(numerator, denominator))
        point = tuple(ample[i] + parameter * ray[i] for i in range(10))
        return golden_quadratic_value(gram, point)

    low = 0
    high = 1
    while square(high).sign() < 0:
        low = high
        high *= 2
    while high - low > 1:
        middle = (low + high) // 2
        if square(middle).sign() < 0:
            low = middle
        else:
            high = middle
    return Fraction(low, denominator), square(low), square(high)


def _norm_two_root_from_functional(
    functional: Sequence[Golden],
    gram: Sequence[Sequence[Golden]],
    ample: Sequence[Golden],
) -> tuple[tuple[Golden, ...], int]:
    inverse = inverse_golden_matrix(gram)
    negative_normal = tuple(
        -sum((inverse[i][j] * functional[j] for j in range(10)), Golden(0))
        for i in range(10)
    )
    candidates = []
    for exponent in range(-12, 13):
        candidate = tuple(entry * (PHI**exponent) for entry in negative_normal)
        if (
            all(entry.is_integral for entry in candidate)
            and _pair(candidate, gram, candidate) == 2
        ):
            if _pair(ample, gram, candidate).sign() > 0:
                candidate = tuple(-entry for entry in candidate)
            candidates.append((candidate, exponent))
    unique = {candidate: exponent for candidate, exponent in candidates}
    if len(unique) != 1:
        raise ArithmeticError(f"functional has {len(unique)} norm-two normalizations")
    root, exponent = next(iter(unique.items()))
    return root, exponent


def _random_lp_exact_escape_search(
    roots: Sequence[tuple[Golden, ...]],
    gram: Sequence[Sequence[Golden]],
    ample: Sequence[Golden],
    *,
    random_trials: int = 512,
) -> dict[str, Any]:
    """Discover numerical vertices, accepting only exact Q(phi) witnesses."""

    functionals = [_mat_vec(gram, root) for root in roots]
    phi_float = (1.0 + 5.0**0.5) / 2.0
    floating = np.array(
        [
            [float(entry.a) + float(entry.b) * phi_float for entry in row]
            for row in functionals
        ],
        dtype=float,
    )
    lengths = np.linalg.norm(floating, axis=1)
    normalized = floating / lengths[:, None]
    section = -normalized.sum(axis=0)
    rng = np.random.default_rng(20260805)
    objectives = [
        *(np.eye(10)[i] for i in range(10)),
        *(-np.eye(10)[i] for i in range(10)),
        *(rng.normal(size=10) for _ in range(random_trials)),
    ]
    exact: dict[tuple[Golden, ...], dict[str, Any]] = {}
    successful_lp = 0
    exact_vertices = 0
    for objective in objectives:
        solution = linprog(
            objective,
            A_ub=normalized,
            b_ub=np.zeros(len(roots)),
            A_eq=section.reshape(1, -1),
            b_eq=np.ones(1),
            bounds=[(None, None)] * 10,
            method="highs",
        )
        if not solution.success:
            continue
        successful_lp += 1
        approximate_values = normalized @ solution.x
        nearest = np.argsort(np.abs(approximate_values))
        chosen: list[int] = []
        current_rank = 0
        for raw_index in nearest[:60]:
            index = int(raw_index)
            candidate = [functionals[item] for item in (*chosen, index)]
            candidate_rank = _row_rank(candidate, 10)
            if candidate_rank > current_rank:
                chosen.append(index)
                current_rank = candidate_rank
            if current_rank == 9:
                break
        if current_rank != 9:
            continue
        ray = _full_rank_kernel([functionals[index] for index in chosen])
        if ray is None:
            continue
        pairings = tuple(_dot(ray, functional) for functional in functionals)
        if all(value.sign() <= 0 for value in pairings):
            oriented = ray
        elif all(value.sign() >= 0 for value in pairings):
            oriented = [-entry for entry in ray]
        else:
            continue
        exact_vertices += 1
        canonical = _positive_canonical(oriented)
        norm = _pair(canonical, gram, canonical)
        sheet = _pair(canonical, gram, ample)
        if norm.sign() < 0 and sheet.sign() < 0:
            continue
        exact_pairings = tuple(_pair(canonical, gram, root) for root in roots)
        exact.setdefault(
            canonical,
            {
                "ray": _encoded_vector(canonical),
                "kind": (
                    "hyperideal"
                    if norm.sign() > 0
                    else "cusp"
                    if norm.sign() == 0 and sheet.sign() < 0
                    else "opposite_sheet"
                ),
                "norm": norm.to_json(),
                "pairing_with_ample": sheet.to_json(),
                "active_facets": [
                    index
                    for index, value in enumerate(exact_pairings)
                    if not value
                ],
                "independent_active_facets": chosen,
                "all_wall_pairings_nonpositive": all(
                    value.sign() <= 0 for value in exact_pairings
                ),
                "wall_pairings_sha256": _canonical_sha(
                    [value.to_json() for value in exact_pairings]
                ),
            },
        )
        if exact:
            break
    return {
        "method": "floating_random_LP_active_set_discovery_then_exact_Qphi_verification",
        "floating_output_used_as_proof": False,
        "objective_count": len(objectives),
        "successful_LP_count": successful_lp,
        "exact_vertex_reconstructions": exact_vertices,
        "exact_noncompact_witness_count": len(exact),
        "exact_noncompact_witnesses": list(exact.values()),
    }


def _second_directed_probe(
    atlas: ShimadaDisc11Atlas,
    roots185: Sequence[tuple[Golden, ...]],
    group,
    gram: Sequence[Sequence[Golden]],
    ample: Sequence[Golden],
) -> dict[str, Any]:
    parameter, square, next_square = _dyadic_before_light(
        ample, ARCHIVED_HYPERIDEAL_RAY, gram
    )
    point = tuple(
        ample[i] + Golden(parameter) * ARCHIVED_HYPERIDEAL_RAY[i]
        for i in range(10)
    )
    old_pairings = tuple(_pair(point, gram, root) for root in roots185)
    location = _walk_until_boundary_or_interior(atlas, point, max_gallery=20_000)
    record: dict[str, Any] = {
        "t": str(parameter),
        "point": _encoded_vector(point),
        "square": square.to_json(),
        "next_dyadic_square": next_square.to_json(),
        "strict_on_all_185_walls": all(value.sign() < 0 for value in old_pairings),
        "old_pairings_sha256": _canonical_sha(
            [value.to_json() for value in old_pairings]
        ),
        "location_kind": location.kind,
        "representative": location.state.index,
        "auxiliary_gallery_length": len(location.gallery),
        "auxiliary_gallery_sha256": _canonical_sha(
            [list(step) for step in location.gallery]
        ),
    }
    if location.kind != "root_boundary" or location.boundary_functional is None:
        return record
    functional = location.boundary_functional
    root, unit_exponent = _norm_two_root_from_functional(functional, gram, ample)
    ray_root_pairing = _pair(ARCHIVED_HYPERIDEAL_RAY, gram, root)
    ample_root_pairing = _pair(ample, gram, root)
    crossing_t = -ample_root_pairing / ray_root_pairing
    face = tuple(
        ample[i] + crossing_t * ARCHIVED_HYPERIDEAL_RAY[i] for i in range(10)
    )
    face_pairings = tuple(_pair(face, gram, old_root) for old_root in roots185)
    root_orbit: dict[tuple[Golden, ...], list[str]] = {}
    for name, action in group:
        image = _row_product(root, action)
        root_orbit.setdefault(image, []).append(name)
    expanded = tuple((*roots185, *(item for item in root_orbit if item not in roots185)))
    new_indices = tuple(range(len(roots185), len(expanded)))
    pairwise_acute = all(
        _pair(left, gram, right).sign() <= 0
        for left in root_orbit
        for right in expanded
        if left != right
    )
    archived_after = _witness_orbit_test(
        ARCHIVED_HYPERIDEAL_RAY, group, expanded, new_indices
    )
    lp_search = _random_lp_exact_escape_search(expanded, gram, ample)
    expanded_sha = _canonical_sha([_encoded_vector(item) for item in expanded])
    face_square = golden_quadratic_value(gram, face)
    verified = (
        parameter == SECOND_TARGET_T
        and square.sign() < 0
        and next_square.sign() >= 0
        and all(value.sign() < 0 for value in old_pairings)
        and location.kind == "root_boundary"
        and location.state.index == SECOND_EXPECTED_REPRESENTATIVE
        and location.boundary_label == SECOND_EXPECTED_BOUNDARY_LABEL
        and functional == SECOND_EXPECTED_BOUNDARY_FUNCTIONAL
        and len(location.gallery) == 25
        and _canonical_sha([list(step) for step in location.gallery])
        == SECOND_GALLERY_SHA256
        and root == SECOND_NEW_ROOT
        and _pair(root, gram, root) == 2
        and ray_root_pairing.sign() > 0
        and _pair(point, gram, root).sign() > 0
        and _pair(face, gram, root) == 0
        and face_square.sign() < 0
        and all(value.sign() < 0 for value in face_pairings)
        and len(root_orbit) == 4
        and len(expanded) == 189
        and expanded_sha == SECOND_EXPANDED_ROOTS_SHA256
        and pairwise_acute
        and archived_after["surviving_rays"] == 0
    )
    record.update(
        {
            "verified": verified,
            "boundary_label": list(location.boundary_label),
            "boundary_functional": _encoded_vector(functional),
            "root": _encoded_vector(root),
            "normalizing_phi_exponent": unit_exponent,
            "root_norm": _pair(root, gram, root).to_json(),
            "ample_root_pairing": ample_root_pairing.to_json(),
            "escape_root_pairing": ray_root_pairing.to_json(),
            "target_root_pairing": _pair(point, gram, root).to_json(),
            "crossing_t": crossing_t.to_json(),
            "face_square": face_square.to_json(),
            "face_strict_on_all_185_old_walls": all(
                value.sign() < 0 for value in face_pairings
            ),
            "root_V4_orbit_size": len(root_orbit),
            "expanded_system_pairwise_acute": pairwise_acute,
            "expanded_wall_count": len(expanded),
            "expanded_roots_sha256": expanded_sha,
            "archived_witness_orbit_after_second_wall": archived_after,
            "exact_escape_search_after_189_walls": lp_search,
        }
    )
    return record


def _directed_boundary_once(
    atlas: ShimadaDisc11Atlas,
    old_roots: Sequence[tuple[Golden, ...]],
    ray: tuple[Golden, ...],
    group,
    gram: Sequence[Sequence[Golden]],
    ample: Sequence[Golden],
) -> tuple[dict[str, Any], tuple[tuple[Golden, ...], ...]]:
    parameter, square, next_square = _dyadic_before_light(ample, ray, gram)
    point = tuple(ample[i] + Golden(parameter) * ray[i] for i in range(10))
    old_pairings = tuple(_pair(point, gram, root) for root in old_roots)
    ray_before = _witness_orbit_test(ray, group, old_roots, ())
    location = _walk_until_boundary_or_interior(atlas, point, max_gallery=20_000)
    gallery = [list(step) for step in location.gallery]
    record: dict[str, Any] = {
        "input_wall_count": len(old_roots),
        "target_ray": _encoded_vector(ray),
        "target_ray_norm": _pair(ray, gram, ray).to_json(),
        "target_ray_V4_orbit_before_probe": ray_before,
        "dyadic_t": str(parameter),
        "dyadic_point": _encoded_vector(point),
        "dyadic_square": square.to_json(),
        "next_dyadic_square": next_square.to_json(),
        "strict_on_all_input_walls": all(value.sign() < 0 for value in old_pairings),
        "input_wall_pairings_sha256": _canonical_sha(
            [value.to_json() for value in old_pairings]
        ),
        "location_kind": location.kind,
        "representative": location.state.index,
        "auxiliary_gallery_length": len(gallery),
        "auxiliary_gallery_sha256": _canonical_sha(gallery),
        "auxiliary_gallery": gallery,
    }
    if location.kind != "root_boundary" or location.boundary_functional is None:
        record["verified_new_wall"] = False
        return record, tuple(old_roots)

    functional = location.boundary_functional
    root, exponent = _norm_two_root_from_functional(functional, gram, ample)
    ample_root = _pair(ample, gram, root)
    ray_root = _pair(ray, gram, root)
    crossing_t = -ample_root / ray_root
    face = tuple(ample[i] + crossing_t * ray[i] for i in range(10))
    face_pairings = tuple(_pair(face, gram, item) for item in old_roots)
    face_square = golden_quadratic_value(gram, face)

    physical_wall = atlas.state_wall_orbits(location.state)[
        location.boundary_label[0]
    ][location.boundary_label[1]]
    physical_matrix = Matrix(1, 20, physical_wall) * atlas.gram_inverse
    physical_integral = all(entry.is_Integer for entry in physical_matrix)
    physical_root = tuple(int(entry) for entry in physical_matrix)
    physical_square = (
        Matrix(1, 20, physical_root)
        * Matrix(atlas.gram)
        * Matrix(20, 1, physical_root)
    )[0]

    root_orbit: dict[tuple[Golden, ...], list[str]] = {}
    for name, action in group:
        image = _row_product(root, action)
        root_orbit.setdefault(image, []).append(name)
    expanded = tuple(
        (*old_roots, *(item for item in root_orbit if item not in old_roots))
    )
    new_indices = tuple(range(len(old_roots), len(expanded)))
    pairwise_acute = all(
        _pair(left, gram, right).sign() <= 0
        for left in root_orbit
        for right in expanded
        if left != right
    )
    target_after = _witness_orbit_test(ray, group, expanded, new_indices)
    expanded_sha = _canonical_sha([_encoded_vector(item) for item in expanded])
    verified = (
        square.sign() < 0
        and next_square.sign() >= 0
        and all(value.sign() < 0 for value in old_pairings)
        and ray_before["surviving_rays"] == ray_before["orbit_size"]
        and all(entry.is_integral for entry in root)
        and _pair(root, gram, root) == 2
        and ample_root.sign() < 0
        and ray_root.sign() > 0
        and _pair(point, gram, root).sign() > 0
        and _pair(face, gram, root) == 0
        and face_square.sign() < 0
        and all(value.sign() < 0 for value in face_pairings)
        and physical_integral
        and physical_square == -2
        and len(root_orbit) in (1, 2, 4)
        and len(expanded) > len(old_roots)
        and pairwise_acute
        and target_after["surviving_rays"] == 0
    )
    record.update(
        {
            "verified_new_wall": verified,
            "boundary_label": list(location.boundary_label),
            "boundary_functional": _encoded_vector(functional),
            "normalizing_phi_exponent": exponent,
            "field_root": _encoded_vector(root),
            "field_root_norm": _pair(root, gram, root).to_json(),
            "physical_NS_root": list(physical_root),
            "physical_NS_square": str(physical_square),
            "ample_root_pairing": ample_root.to_json(),
            "target_ray_root_pairing": ray_root.to_json(),
            "target_point_root_pairing": _pair(point, gram, root).to_json(),
            "crossing_t": crossing_t.to_json(),
            "strict_timelike_face_witness": {
                "point": _encoded_vector(face),
                "square": face_square.to_json(),
                "new_root_pairing_zero": _pair(face, gram, root) == 0,
                "strict_on_all_prior_walls": all(
                    value.sign() < 0 for value in face_pairings
                ),
                "prior_wall_pairings_sha256": _canonical_sha(
                    [value.to_json() for value in face_pairings]
                ),
            },
            "new_root_V4_orbit_size": len(root_orbit),
            "new_root_V4_orbit": [
                {
                    "root": _encoded_vector(item),
                    "obtained_by": root_orbit[item],
                }
                for item in root_orbit
            ],
            "new_distinct_wall_count": len(new_indices),
            "expanded_wall_count": len(expanded),
            "expanded_roots_sha256": expanded_sha,
            "expanded_system_pairwise_acute": pairwise_acute,
            "target_ray_V4_orbit_after_probe": target_after,
        }
    )
    return record, expanded


def _iterated_directed_search(
    atlas: ShimadaDisc11Atlas,
    roots189: tuple[tuple[Golden, ...], ...],
    initial_search: dict[str, Any],
    group,
    gram: Sequence[Sequence[Golden]],
    ample: Sequence[Golden],
    *,
    max_new_orbits: int = 30,
    iteration_start: int = 3,
) -> dict[str, Any]:
    current = roots189
    search = initial_search
    iterations = []
    terminal_failed_attempts = []
    for index in range(max_new_orbits):
        witnesses = search["exact_noncompact_witnesses"]
        if not witnesses:
            break
        witness = witnesses[0]
        base_ray = tuple(Golden.from_json(entry) for entry in witness["ray"])
        ray_orbit: dict[tuple[Golden, ...], str] = {}
        for name, action in group:
            image = _positive_canonical(_row_product(base_ray, action))
            ray_orbit.setdefault(image, name)
        selected = None
        expanded = current
        failed_attempts = []
        for ray, name in ray_orbit.items():
            probe_record, candidate_expanded = _directed_boundary_once(
                atlas, current, ray, group, gram, ample
            )
            probe_record["attempted_V4_image_by"] = name
            if (
                probe_record.get("verified_new_wall")
                and len(candidate_expanded) > len(current)
            ):
                selected = probe_record
                expanded = candidate_expanded
                break
            failed_attempts.append(probe_record)
        if selected is None:
            terminal_failed_attempts = failed_attempts
            break
        selected["iteration"] = index + iteration_start
        selected["escape_discovery"] = search
        selected["failed_V4_image_attempts_before_success"] = failed_attempts
        iterations.append(selected)
        current = expanded
        search = _random_lp_exact_escape_search(current, gram, ample)
    return {
        "max_new_orbits": max_new_orbits,
        "completed_new_orbits": sum(
            bool(item.get("verified_new_wall")) for item in iterations
        ),
        "iterations": iterations,
        "terminal_failed_V4_image_attempts": terminal_failed_attempts,
        "terminal_wall_count": len(current),
        "terminal_roots_sha256": _canonical_sha(
            [_encoded_vector(item) for item in current]
        ),
        "terminal_escape_search": search,
        "stopping_reason": (
            "no_exact_noncompact_witness_found_by_fixed_LP_schedule"
            if not search["exact_noncompact_witnesses"]
            else "all_V4_image_locations_failed"
            if terminal_failed_attempts
            else "new_orbit_limit_reached"
        ),
    }


def _failed_attempt_summary(record: dict[str, Any]) -> dict[str, Any]:
    """Retain the deterministic failure trace without duplicating large pairing tables."""

    return {
        "attempted_V4_image_by": record.get("attempted_V4_image_by"),
        "target_ray": record.get("target_ray"),
        "location_kind": record.get("location_kind"),
        "representative": record.get("representative"),
        "auxiliary_gallery_length": record.get("auxiliary_gallery_length"),
        "auxiliary_gallery_sha256": record.get("auxiliary_gallery_sha256"),
        "boundary_label": record.get("boundary_label"),
        "field_root": record.get("field_root"),
        "verified_new_wall": record.get("verified_new_wall", False),
    }


def _checkpoint_step(
    record: dict[str, Any],
    current: Sequence[tuple[Golden, ...]],
    group,
) -> tuple[dict[str, Any], tuple[tuple[Golden, ...], ...]]:
    """Compress one successful search record into a replayable exact step."""

    if record.get("verified_new_wall") is not True:
        raise ArithmeticError("cannot checkpoint an unverified directed wall")
    root = _decoded_vector(record["field_root"])
    root_orbit: dict[tuple[Golden, ...], list[str]] = {}
    for name, action in group:
        image = _row_product(root, action)
        root_orbit.setdefault(image, []).append(name)
    new_roots = tuple(item for item in root_orbit if item not in current)
    expanded = tuple((*current, *new_roots))
    output_hash = _roots_sha256(expanded)
    if (
        len(current) != record["input_wall_count"]
        or len(new_roots) != record["new_distinct_wall_count"]
        or len(expanded) != record["expanded_wall_count"]
        or output_hash != record["expanded_roots_sha256"]
    ):
        raise ArithmeticError("directed record does not replay while checkpointing")
    escape_search = record["escape_discovery"]
    escape_witnesses = escape_search["exact_noncompact_witnesses"]
    if not escape_witnesses:
        raise ArithmeticError("directed record has no source escape witness")
    strict_face = {
        key: value
        for key, value in record["strict_timelike_face_witness"].items()
        if key != "point"
    }
    step = {
        "iteration": record["iteration"],
        "input_wall_count": len(current),
        "input_roots_sha256": _roots_sha256(current),
        "escape_discovery": {
            "method": escape_search["method"],
            "floating_output_used_as_proof": escape_search[
                "floating_output_used_as_proof"
            ],
            "source_witness_active_facets": escape_witnesses[0]["active_facets"],
            "source_witness_ray_sha256": _canonical_sha(
                escape_witnesses[0]["ray"]
            ),
            "source_witness_wall_pairings_sha256": escape_witnesses[0][
                "wall_pairings_sha256"
            ],
        },
        "attempted_V4_image_by": record["attempted_V4_image_by"],
        "target_ray": record["target_ray"],
        "target_ray_norm": record["target_ray_norm"],
        "directed_probe": {
            "dyadic_t": record["dyadic_t"],
            "dyadic_square": record["dyadic_square"],
            "next_dyadic_square": record["next_dyadic_square"],
            "strict_on_all_input_walls": record["strict_on_all_input_walls"],
            "input_wall_pairings_sha256": record[
                "input_wall_pairings_sha256"
            ],
        },
        "shimada_location": {
            "kind": record["location_kind"],
            "representative": record["representative"],
            "boundary_label": record["boundary_label"],
            "boundary_functional": record["boundary_functional"],
            "auxiliary_gallery_length": record["auxiliary_gallery_length"],
            "auxiliary_gallery_sha256": record["auxiliary_gallery_sha256"],
            "auxiliary_gallery": record["auxiliary_gallery"],
        },
        "new_wall": {
            "normalizing_phi_exponent": record["normalizing_phi_exponent"],
            "field_root": record["field_root"],
            "field_root_norm": record["field_root_norm"],
            "physical_NS_root": record["physical_NS_root"],
            "physical_NS_square": record["physical_NS_square"],
            "ample_root_pairing": record["ample_root_pairing"],
            "target_ray_root_pairing": record["target_ray_root_pairing"],
            "target_point_root_pairing": record["target_point_root_pairing"],
            "crossing_t": record["crossing_t"],
            "strict_timelike_face_witness": strict_face,
            "V4_orbit_size": record["new_root_V4_orbit_size"],
        },
        "output": {
            "new_distinct_wall_count": len(new_roots),
            "wall_count": len(expanded),
            "roots_sha256": output_hash,
            "pairwise_acute": record["expanded_system_pairwise_acute"],
            "target_ray_orbit_surviving_after_wall": record[
                "target_ray_V4_orbit_after_probe"
            ]["surviving_rays"],
        },
        "failed_V4_image_attempts_before_success": [
            _failed_attempt_summary(item)
            for item in record["failed_V4_image_attempts_before_success"]
        ],
    }
    return step, expanded


def _seal_checkpoint(payload: dict[str, Any]) -> dict[str, Any]:
    return {
        "payload_sha256": _canonical_sha(payload),
        "payload": payload,
    }


def _extract_checkpoint(document: dict[str, Any]) -> tuple[dict[str, Any], str]:
    """Accept a bare seal or an artifact/result containing one."""

    if "payload" in document and "payload_sha256" in document:
        payload = document["payload"]
        digest = document["payload_sha256"]
    elif "checkpoint" in document:
        return _extract_checkpoint(document["checkpoint"])
    elif "frozen_checkpoint" in document:
        return _extract_checkpoint(document["frozen_checkpoint"])
    else:
        raise ValueError("no sealed directed-wall checkpoint found")
    if _canonical_sha(payload) != digest:
        raise ValueError("directed-wall checkpoint payload hash mismatch")
    return payload, digest


def _build_checkpoint(
    iterated: dict[str, Any],
    roots_start: Sequence[tuple[Golden, ...]],
    group,
    *,
    source_pins: dict[str, Any] | None = None,
    prefix_payload: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a compact, exact, resumable transcript from raw probe records."""

    current = tuple(roots_start)
    if prefix_payload is None:
        base_wall_count = len(current)
        base_roots_sha256 = _roots_sha256(current)
        steps: list[dict[str, Any]] = []
        pins = dict(source_pins or {})
    else:
        base_wall_count = prefix_payload["base"]["wall_count"]
        base_roots_sha256 = prefix_payload["base"]["roots_sha256"]
        steps = list(prefix_payload["steps"])
        pins = dict(prefix_payload["source_pins"])
        if source_pins is not None and pins != source_pins:
            raise ValueError("resume source pins differ from the checkpoint")
    for record in iterated["iterations"]:
        step, current = _checkpoint_step(record, current, group)
        steps.append(step)
    terminal_search = iterated["terminal_escape_search"]
    terminal_witnesses = terminal_search["exact_noncompact_witnesses"]
    payload = {
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "analysis_id": "shimada_disc11_directed_v4_wall_orbits",
        "source_pins": pins,
        "base": {
            "wall_count": base_wall_count,
            "roots_sha256": base_roots_sha256,
        },
        "steps": steps,
        "terminal": {
            "wall_count": len(current),
            "roots_sha256": _roots_sha256(current),
            "escape_search": terminal_search,
            "escape_witness_count": len(terminal_witnesses),
            "stopping_reason": iterated["stopping_reason"],
            "failed_V4_image_attempts": [
                _failed_attempt_summary(item)
                for item in iterated["terminal_failed_V4_image_attempts"]
            ],
        },
        "proof_scope": (
            "exact_facet_lower_bound_and_exact_noncompactness_of_the_recorded_"
            "finite_outer_cone_not_global_wall_completeness"
        ),
    }
    if base_wall_count == 189 and len(steps) == 30:
        output_hashes = tuple(step["output"]["roots_sha256"] for step in steps)
        if output_hashes != FROZEN_DIRECTED_OUTPUT_HASHES:
            raise ArithmeticError("the frozen 30-step output hash sequence changed")
        if (
            payload["terminal"]["wall_count"] != FROZEN_TERMINAL_WALL_COUNT
            or payload["terminal"]["roots_sha256"]
            != FROZEN_TERMINAL_ROOTS_SHA256
        ):
            raise ArithmeticError("the frozen 309-wall terminal state changed")
    seal = _seal_checkpoint(payload)
    if (
        base_wall_count == 189
        and len(steps) == 30
        and payload["terminal"]["wall_count"] == FROZEN_TERMINAL_WALL_COUNT
        and seal["payload_sha256"] != FROZEN_CHECKPOINT_PAYLOAD_SHA256
    ):
        raise ArithmeticError("the frozen 30-step checkpoint payload changed")
    return seal


def _verify_checkpoint(
    document: dict[str, Any],
    atlas: ShimadaDisc11Atlas,
    roots_start: Sequence[tuple[Golden, ...]],
    group,
    gram: Sequence[Sequence[Golden]],
    ample: Sequence[Golden],
    *,
    replay_galleries: bool = True,
) -> tuple[tuple[tuple[Golden, ...], ...], dict[str, Any], dict[str, Any]]:
    """Replay a sealed transcript with exact Q(phi) and physical NS checks."""

    payload, digest = _extract_checkpoint(document)
    if payload.get("schema_version") != CHECKPOINT_SCHEMA_VERSION:
        raise ValueError("unsupported directed-wall checkpoint schema")
    current = tuple(roots_start)
    if payload["base"] != {
        "wall_count": len(current),
        "roots_sha256": _roots_sha256(current),
    }:
        raise ValueError("checkpoint base does not match the reconstructed roots")
    replayed_galleries = 0
    replayed_gallery_steps = 0
    for offset, step in enumerate(payload["steps"]):
        if step["input_wall_count"] != len(current):
            raise ArithmeticError(f"step {offset} input wall count changed")
        if step["input_roots_sha256"] != _roots_sha256(current):
            raise ArithmeticError(f"step {offset} input root hash changed")

        ray = _decoded_vector(step["target_ray"])
        ray_pairings = tuple(_pair(ray, gram, root) for root in current)
        if (
            _pair(ray, gram, ray) != Golden.from_json(step["target_ray_norm"])
            or _pair(ray, gram, ray).sign() <= 0
            or not all(value.sign() <= 0 for value in ray_pairings)
        ):
            raise ArithmeticError(f"step {offset} target is not hyperideal feasible")
        directed = step["directed_probe"]
        parameter = Fraction(directed["dyadic_t"])
        point = tuple(ample[i] + Golden(parameter) * ray[i] for i in range(10))
        old_pairings = tuple(_pair(point, gram, root) for root in current)
        dyadic_t, square, next_square = _dyadic_before_light(ample, ray, gram)
        if (
            dyadic_t != parameter
            or square != Golden.from_json(directed["dyadic_square"])
            or next_square != Golden.from_json(directed["next_dyadic_square"])
            or square.sign() >= 0
            or next_square.sign() < 0
            or not all(value.sign() < 0 for value in old_pairings)
            or _canonical_sha([value.to_json() for value in old_pairings])
            != directed["input_wall_pairings_sha256"]
        ):
            raise ArithmeticError(f"step {offset} directed point check failed")

        location_record = step["shimada_location"]
        if replay_galleries:
            location = _walk_until_boundary_or_interior(
                atlas, point, max_gallery=20_000
            )
            gallery = [list(item) for item in location.gallery]
            replayed_galleries += 1
            replayed_gallery_steps += len(gallery)
            if (
                location.kind != location_record["kind"]
                or location.state.index != location_record["representative"]
                or list(location.boundary_label or ())
                != location_record["boundary_label"]
                or _encoded_vector(location.boundary_functional or ())
                != location_record["boundary_functional"]
                or gallery != location_record["auxiliary_gallery"]
                or _canonical_sha(gallery)
                != location_record["auxiliary_gallery_sha256"]
            ):
                raise ArithmeticError(f"step {offset} Shimada gallery changed")
        else:
            gallery = location_record["auxiliary_gallery"]
            if (
                len(gallery) != location_record["auxiliary_gallery_length"]
                or _canonical_sha(gallery)
                != location_record["auxiliary_gallery_sha256"]
            ):
                raise ArithmeticError(f"step {offset} gallery hash changed")
            location = None

        wall = step["new_wall"]
        functional = _decoded_vector(location_record["boundary_functional"])
        root, exponent = _norm_two_root_from_functional(functional, gram, ample)
        if (
            root != _decoded_vector(wall["field_root"])
            or exponent != wall["normalizing_phi_exponent"]
            or not all(entry.is_integral for entry in root)
            or _pair(root, gram, root) != 2
        ):
            raise ArithmeticError(f"step {offset} field root reconstruction failed")

        if location is not None:
            label = tuple(location_record["boundary_label"])
            physical_wall = atlas.state_wall_orbits(location.state)[label[0]][label[1]]
            physical_matrix = Matrix(1, 20, physical_wall) * atlas.gram_inverse
            if not all(entry.is_Integer for entry in physical_matrix):
                raise ArithmeticError(f"step {offset} physical root is not integral")
            physical_root = tuple(int(entry) for entry in physical_matrix)
            physical_square = (
                Matrix(1, 20, physical_root)
                * Matrix(atlas.gram)
                * Matrix(20, 1, physical_root)
            )[0]
            if (
                list(physical_root) != wall["physical_NS_root"]
                or str(physical_square) != wall["physical_NS_square"]
                or physical_square != -2
            ):
                raise ArithmeticError(f"step {offset} physical (-2)-class changed")

        ample_root = _pair(ample, gram, root)
        ray_root = _pair(ray, gram, root)
        crossing_t = -ample_root / ray_root
        face = tuple(ample[i] + crossing_t * ray[i] for i in range(10))
        face_pairings = tuple(_pair(face, gram, item) for item in current)
        face_record = wall["strict_timelike_face_witness"]
        if (
            ample_root != Golden.from_json(wall["ample_root_pairing"])
            or ample_root.sign() >= 0
            or ray_root != Golden.from_json(wall["target_ray_root_pairing"])
            or ray_root.sign() <= 0
            or crossing_t != Golden.from_json(wall["crossing_t"])
            or _pair(face, gram, root) != 0
            or _pair(face, gram, face) != Golden.from_json(face_record["square"])
            or _pair(face, gram, face).sign() >= 0
            or not all(value.sign() < 0 for value in face_pairings)
            or _canonical_sha([value.to_json() for value in face_pairings])
            != face_record["prior_wall_pairings_sha256"]
        ):
            raise ArithmeticError(f"step {offset} strict face check failed")

        root_orbit: dict[tuple[Golden, ...], list[str]] = {}
        for name, action in group:
            image = _row_product(root, action)
            root_orbit.setdefault(image, []).append(name)
        new_roots = tuple(item for item in root_orbit if item not in current)
        expanded = tuple((*current, *new_roots))
        output = step["output"]
        pairwise_acute = all(
            _pair(left, gram, right).sign() <= 0
            for left in root_orbit
            for right in expanded
            if left != right
        )
        target_after = _witness_orbit_test(ray, group, expanded, ())
        if (
            len(root_orbit) != wall["V4_orbit_size"]
            or len(new_roots) != output["new_distinct_wall_count"]
            or len(expanded) != output["wall_count"]
            or _roots_sha256(expanded) != output["roots_sha256"]
            or not pairwise_acute
            or output["pairwise_acute"] is not True
            or target_after["surviving_rays"]
            != output["target_ray_orbit_surviving_after_wall"]
            or target_after["surviving_rays"] != 0
        ):
            raise ArithmeticError(f"step {offset} V4 orbit closure failed")
        current = expanded

    terminal = payload["terminal"]
    if (
        terminal["wall_count"] != len(current)
        or terminal["roots_sha256"] != _roots_sha256(current)
    ):
        raise ArithmeticError("terminal root system changed")
    terminal_witnesses = terminal["escape_search"]["exact_noncompact_witnesses"]
    if len(terminal_witnesses) != terminal["escape_witness_count"]:
        raise ArithmeticError("terminal escape witness count changed")
    witness_record = terminal_witnesses[0] if terminal_witnesses else None
    terminal_kind = "none_found"
    if witness_record is not None:
        witness = _decoded_vector(witness_record["ray"])
        pairings = tuple(_pair(witness, gram, root) for root in current)
        active = [index for index, value in enumerate(pairings) if not value]
        norm = _pair(witness, gram, witness)
        if (
            not all(value.sign() <= 0 for value in pairings)
            or norm != Golden.from_json(witness_record["norm"])
            or norm.sign() < 0
            or active != witness_record["active_facets"]
            or _canonical_sha([value.to_json() for value in pairings])
            != witness_record["wall_pairings_sha256"]
        ):
            raise ArithmeticError("terminal escape witness changed")
        terminal_kind = "hyperideal" if norm.sign() > 0 else "cusp"

    coxeter_relation_counts: dict[str, int] = {}
    for i, left in enumerate(current):
        for right in current[i + 1 :]:
            relation = classify_root_pair(
                Golden(2), _pair(left, gram, right), Golden(2)
            )
            if relation.kind == "finite_unrecognized":
                raise ArithmeticError("terminal acute system has a non-Coxeter angle")
            key = (
                f"m{relation.label}"
                if relation.kind == "finite"
                else relation.kind
            )
            coxeter_relation_counts[key] = coxeter_relation_counts.get(key, 0) + 1

    frozen_309_prefix = (
        payload["base"]["wall_count"] == 189
        and len(payload["steps"]) >= 30
        and tuple(
            step["output"]["roots_sha256"] for step in payload["steps"][:30]
        )
        == FROZEN_DIRECTED_OUTPUT_HASHES
        and payload["steps"][29]["output"]["wall_count"]
        == FROZEN_TERMINAL_WALL_COUNT
        and payload["steps"][29]["output"]["roots_sha256"]
        == FROZEN_TERMINAL_ROOTS_SHA256
    )
    terminal_equals_frozen_309 = (
        len(current) == FROZEN_TERMINAL_WALL_COUNT
        and _roots_sha256(current) == FROZEN_TERMINAL_ROOTS_SHA256
        and witness_record is not None
        and _decoded_vector(witness_record["ray"]) == FROZEN_TERMINAL_RAY
        and Golden.from_json(witness_record["norm"]) == FROZEN_TERMINAL_NORM
        and tuple(witness_record["active_facets"])
        == FROZEN_TERMINAL_ACTIVE_FACETS
        and witness_record["wall_pairings_sha256"]
        == FROZEN_TERMINAL_PAIRINGS_SHA256
    )
    audit = {
        "verified": True,
        "payload_sha256": digest,
        "base_wall_count": payload["base"]["wall_count"],
        "replayed_orbits": len(payload["steps"]),
        "replayed_galleries": replayed_galleries,
        "replayed_gallery_steps": replayed_gallery_steps,
        "terminal_wall_count": len(current),
        "terminal_roots_sha256": _roots_sha256(current),
        "terminal_escape_kind": terminal_kind,
        "terminal_coxeter_relation_counts": coxeter_relation_counts,
        "frozen_309_prefix_match": frozen_309_prefix,
        "terminal_equals_frozen_309": terminal_equals_frozen_309,
    }
    return current, audit, payload


def probe(
    disc11: Path,
    transcript: Path = DEFAULT_TRANSCRIPT,
    extension: Path = DEFAULT_EXTENSION,
    *,
    max_new_orbits: int = 30,
) -> dict[str, Any]:
    source_sha = hashlib.sha256(disc11.read_bytes()).hexdigest()
    if source_sha != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")
    atlas = ShimadaDisc11Atlas(disc11)
    roots97, roots181, provenance, group, frame = _closed_181(
        atlas, transcript, extension
    )
    gram = tuple(tuple(entry for entry in row) for row in field_gram())
    ample = tuple(AMPLE_FIELD_COEFFICIENTS)
    point = tuple(
        ample[i] + Golden(TARGET_T) * NEW_HYPERIDEAL_RAY[i]
        for i in range(10)
    )
    next_point = tuple(
        ample[i] + Golden(TARGET_T + Fraction(1, 32)) * NEW_HYPERIDEAL_RAY[i]
        for i in range(10)
    )
    point_square = golden_quadratic_value(gram, point)
    next_square = golden_quadratic_value(gram, next_point)
    old_pairings = tuple(_pair(point, gram, root) for root in roots181)
    old_strict = all(value.sign() < 0 for value in old_pairings)

    location = _walk_until_boundary_or_interior(atlas, point, max_gallery=20_000)
    if location.boundary_functional is None:
        raise RuntimeError(f"directed probe ended at {location.kind}, not a root boundary")
    functional = location.boundary_functional
    inverse = inverse_golden_matrix(gram)
    negative_normal = tuple(
        -sum((inverse[i][j] * functional[j] for j in range(10)), Golden(0))
        for i in range(10)
    )
    reconstructed_root = tuple(entry / PHI for entry in negative_normal)
    ample_root_pairing = _pair(ample, gram, reconstructed_root)
    escape_root_pairing = _pair(
        NEW_HYPERIDEAL_RAY, gram, reconstructed_root
    )
    crossing_t = -ample_root_pairing / escape_root_pairing
    face_point = tuple(
        ample[i] + crossing_t * NEW_HYPERIDEAL_RAY[i] for i in range(10)
    )
    face_square = golden_quadratic_value(gram, face_point)
    face_old_pairings = tuple(_pair(face_point, gram, root) for root in roots181)

    # The actual Shimada wall is an integral (-2)-root in the physical NS
    # lattice.  Recover it independently from the transported dual row.
    physical_wall = atlas.state_wall_orbits(location.state)[
        location.boundary_label[0]
    ][location.boundary_label[1]]
    physical_root_matrix = Matrix(1, 20, physical_wall) * atlas.gram_inverse
    physical_root_integral = all(entry.is_Integer for entry in physical_root_matrix)
    physical_root = tuple(int(entry) for entry in physical_root_matrix)
    physical_square = (
        Matrix(1, 20, physical_root)
        * Matrix(atlas.gram)
        * Matrix(20, 1, physical_root)
    )[0]

    new_orbit: dict[tuple[Golden, ...], list[str]] = {}
    for name, action in group:
        image = _row_product(reconstructed_root, action)
        new_orbit.setdefault(image, []).append(name)
    new_roots = tuple(new_orbit)
    expanded = tuple((*roots181, *(root for root in new_roots if root not in provenance)))
    new_indices = tuple(range(len(roots181), len(expanded)))
    pairwise_acute = all(
        _pair(left, gram, right).sign() <= 0
        for left in new_roots
        for right in expanded
        if left != right
    )

    # Reorient the archived witness once using its original 97 inequalities.
    archived = verify_hyperideal_extreme_ray(
        roots97, gram, HYPERIDEAL_WITNESS, dimension=9
    ).ray
    archived_pairings = tuple(_pair(archived, gram, root) for root in roots97)
    if all(value.sign() >= 0 for value in archived_pairings):
        archived = tuple(-entry for entry in archived)
    elif not all(value.sign() <= 0 for value in archived_pairings):
        raise ArithmeticError("archived witness has mixed original signs")

    archived_test = _witness_orbit_test(
        archived, group, expanded, new_indices
    )
    new_witness_test = _witness_orbit_test(
        NEW_HYPERIDEAL_RAY, group, expanded, new_indices
    )
    second_probe = _second_directed_probe(
        atlas, expanded, group, gram, ample
    )
    second_root_orbit: dict[tuple[Golden, ...], list[str]] = {}
    for name, action in group:
        image = _row_product(SECOND_NEW_ROOT, action)
        second_root_orbit.setdefault(image, []).append(name)
    roots189 = tuple(
        (*expanded, *(item for item in second_root_orbit if item not in expanded))
    )
    if (
        len(roots189) != 189
        or _canonical_sha([_encoded_vector(item) for item in roots189])
        != SECOND_EXPANDED_ROOTS_SHA256
    ):
        raise ArithmeticError("the frozen 189-wall system changed")
    iterated = _iterated_directed_search(
        atlas,
        roots189,
        second_probe["exact_escape_search_after_189_walls"],
        group,
        gram,
        ample,
        max_new_orbits=max_new_orbits,
    )
    checkpoint = _build_checkpoint(
        iterated,
        roots189,
        group,
        source_pins={
            "disc11_sha256": source_sha,
            "transcript_sha256": frame["transcript_sha256"],
            "extension_sha256": frame["extension_sha256"],
            "transported_97_roots_sha256": frame[
                "transported_roots_sha256"
            ],
            "base_189_roots_sha256": SECOND_EXPANDED_ROOTS_SHA256,
        },
    )
    gallery = [list(step) for step in location.gallery]
    verified = (
        frame["transported_roots_sha256"] == EXPECTED_TRANSFORMED_ROOTS_SHA256
        and point_square.sign() < 0
        and next_square.sign() >= 0
        and old_strict
        and location.kind == "root_boundary"
        and location.state.index == EXPECTED_REPRESENTATIVE
        and location.boundary_label == EXPECTED_BOUNDARY_LABEL
        and functional == EXPECTED_BOUNDARY_FUNCTIONAL
        and len(gallery) == 67
        and reconstructed_root == EXPECTED_NEW_ROOT
        and all(entry.is_integral for entry in reconstructed_root)
        and _pair(reconstructed_root, gram, reconstructed_root) == 2
        and ample_root_pairing == EXPECTED_AMPLE_ROOT_PAIRING
        and escape_root_pairing == EXPECTED_ESCAPE_ROOT_PAIRING
        and crossing_t == EXPECTED_CROSSING_T
        and _pair(face_point, gram, reconstructed_root) == 0
        and face_square == EXPECTED_FACE_SQUARE
        and face_square.sign() < 0
        and all(value.sign() < 0 for value in face_old_pairings)
        and _pair(point, gram, reconstructed_root).sign() > 0
        and physical_root_integral
        and physical_square == -2
        and len(new_roots) in (1, 2, 4)
        and len(expanded) == len(roots181) + len(new_indices)
        and pairwise_acute
        and second_probe.get("verified") is True
        and all(
            item.get("verified_new_wall") is True
            for item in iterated["iterations"]
        )
        and len(iterated["iterations"]) <= max_new_orbits
    )
    return {
        "verified": verified,
        "source_sha256": source_sha,
        **frame,
        "directed_point": {
            "t": str(TARGET_T),
            "coordinates": _encoded_vector(point),
            "square": point_square.to_json(),
            "next_dyadic_square": next_square.to_json(),
            "strict_on_all_181_old_walls": old_strict,
            "old_wall_pairings_sha256": _canonical_sha(
                [value.to_json() for value in old_pairings]
            ),
        },
        "shimada_location": {
            "kind": location.kind,
            "representative": location.state.index,
            "boundary_label": list(location.boundary_label),
            "boundary_functional": _encoded_vector(functional),
            "auxiliary_gallery_length": len(gallery),
            "auxiliary_gallery_sha256": _canonical_sha(gallery),
            "auxiliary_gallery": gallery,
            "physical_state_wall_sha256": _canonical_sha(
                atlas.state_wall_orbits(location.state)
            ),
        },
        "new_wall": {
            "field_root": _encoded_vector(reconstructed_root),
            "field_norm": _pair(reconstructed_root, gram, reconstructed_root).to_json(),
            "ample_pairing": ample_root_pairing.to_json(),
            "escape_ray_pairing": escape_root_pairing.to_json(),
            "target_pairing": _pair(point, gram, reconstructed_root).to_json(),
            "exact_crossing_parameter": crossing_t.to_json(),
            "strict_timelike_face_witness": {
                "coordinates": _encoded_vector(face_point),
                "square": face_square.to_json(),
                "new_wall_pairing_zero": _pair(
                    face_point, gram, reconstructed_root
                ) == 0,
                "strict_on_all_181_old_walls": all(
                    value.sign() < 0 for value in face_old_pairings
                ),
                "old_wall_pairings_sha256": _canonical_sha(
                    [value.to_json() for value in face_old_pairings]
                ),
            },
            "physical_NS_root": list(physical_root),
            "physical_NS_square": str(physical_square),
            "V4_orbit_size": len(new_roots),
            "V4_orbit": [
                {
                    "root": _encoded_vector(root),
                    "obtained_by": new_orbit[root],
                }
                for root in new_roots
            ],
        },
        "expanded_wall_system": {
            "old_count": len(roots181),
            "new_distinct_count": len(new_indices),
            "total_count": len(expanded),
            "roots_sha256": _canonical_sha(
                [_encoded_vector(root) for root in expanded]
            ),
            "pairwise_acute": pairwise_acute,
        },
        "archived_hyperideal_orbit_after_new_walls": archived_test,
        "new_hyperideal_orbit_after_new_walls": new_witness_test,
        "second_probe_on_surviving_archived_orbit": second_probe,
        "iterated_directed_escape_search": iterated,
        "frozen_checkpoint": checkpoint,
        "conclusion": (
            "new_missing_K3_curve_wall_found_but_an_exact_hyperideal_witness_survives"
            if archived_test["surviving_rays"] or new_witness_test["surviving_rays"]
            else "both_known_hyperideal_witness_orbits_are_cut"
        ),
        "scope": "one_directed_dyadic_probe_not_global_wall_completeness",
    }


def _checkpoint_source_pins(
    disc11: Path,
    transcript: Path,
    extension: Path,
    frame: dict[str, Any],
) -> dict[str, Any]:
    return {
        "disc11_sha256": hashlib.sha256(disc11.read_bytes()).hexdigest(),
        "transcript_sha256": frame["transcript_sha256"],
        "extension_sha256": frame["extension_sha256"],
        "transported_97_roots_sha256": frame["transported_roots_sha256"],
        "base_189_roots_sha256": SECOND_EXPANDED_ROOTS_SHA256,
    }


def _checkpoint_context(
    disc11: Path,
    transcript: Path,
    extension: Path,
):
    source_sha = hashlib.sha256(disc11.read_bytes()).hexdigest()
    if source_sha != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")
    atlas = ShimadaDisc11Atlas(disc11)
    _, roots181, _, group, frame = _closed_181(atlas, transcript, extension)
    roots189 = _frozen_roots189(roots181, group)
    gram = tuple(tuple(entry for entry in row) for row in field_gram())
    ample = tuple(AMPLE_FIELD_COEFFICIENTS)
    return atlas, roots189, group, gram, ample, frame


def verify_checkpoint_file(
    disc11: Path,
    checkpoint_path: Path,
    transcript: Path = DEFAULT_TRANSCRIPT,
    extension: Path = DEFAULT_EXTENSION,
    *,
    replay_galleries: bool = True,
) -> dict[str, Any]:
    document = json.loads(checkpoint_path.read_text())
    atlas, roots189, group, gram, ample, frame = _checkpoint_context(
        disc11, transcript, extension
    )
    _, audit, payload = _verify_checkpoint(
        document,
        atlas,
        roots189,
        group,
        gram,
        ample,
        replay_galleries=replay_galleries,
    )
    expected_pins = _checkpoint_source_pins(disc11, transcript, extension, frame)
    if payload["source_pins"] != expected_pins:
        raise ValueError("checkpoint source pins do not match the supplied inputs")
    return audit


def resume_checkpoint_file(
    disc11: Path,
    checkpoint_path: Path,
    transcript: Path = DEFAULT_TRANSCRIPT,
    extension: Path = DEFAULT_EXTENSION,
    *,
    additional_orbits: int,
    replay_galleries: bool = True,
) -> dict[str, Any]:
    if additional_orbits < 0:
        raise ValueError("additional_orbits must be nonnegative")
    document = json.loads(checkpoint_path.read_text())
    atlas, roots189, group, gram, ample, frame = _checkpoint_context(
        disc11, transcript, extension
    )
    current, audit, payload = _verify_checkpoint(
        document,
        atlas,
        roots189,
        group,
        gram,
        ample,
        replay_galleries=replay_galleries,
    )
    expected_pins = _checkpoint_source_pins(disc11, transcript, extension, frame)
    if payload["source_pins"] != expected_pins:
        raise ValueError("checkpoint source pins do not match the supplied inputs")
    terminal_search = payload["terminal"]["escape_search"]
    iteration_start = (
        payload["steps"][-1]["iteration"] + 1 if payload["steps"] else 3
    )
    iterated = _iterated_directed_search(
        atlas,
        current,
        terminal_search,
        group,
        gram,
        ample,
        max_new_orbits=additional_orbits,
        iteration_start=iteration_start,
    )
    checkpoint = _build_checkpoint(
        iterated,
        current,
        group,
        source_pins=expected_pins,
        prefix_payload=payload,
    )
    return {
        "resume_base_audit": audit,
        "additional_search": iterated,
        "checkpoint": checkpoint,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    parser.add_argument("--transcript", type=Path, default=DEFAULT_TRANSCRIPT)
    parser.add_argument("--extension", type=Path, default=DEFAULT_EXTENSION)
    parser.add_argument("--max-new-orbits", type=int, default=30)
    parser.add_argument("--checkpoint-only", action="store_true")
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--verify-checkpoint", type=Path)
    modes.add_argument("--resume-checkpoint", type=Path)
    parser.add_argument("--additional-orbits", type=int, default=0)
    parser.add_argument(
        "--output",
        type=Path,
        help="write the JSON result atomically to this path instead of stdout",
    )
    parser.add_argument(
        "--skip-gallery-replay",
        action="store_true",
        help="verify exact roots and inequalities but trust the hashed Shimada galleries",
    )
    args = parser.parse_args()
    if args.verify_checkpoint is not None:
        result = verify_checkpoint_file(
            args.disc11,
            args.verify_checkpoint,
            args.transcript,
            args.extension,
            replay_galleries=not args.skip_gallery_replay,
        )
    elif args.resume_checkpoint is not None:
        result = resume_checkpoint_file(
            args.disc11,
            args.resume_checkpoint,
            args.transcript,
            args.extension,
            additional_orbits=args.additional_orbits,
            replay_galleries=not args.skip_gallery_replay,
        )
    else:
        result = probe(
            args.disc11,
            args.transcript,
            args.extension,
            max_new_orbits=args.max_new_orbits,
        )
        if args.checkpoint_only:
            result = result["frozen_checkpoint"]
    if args.checkpoint_only or args.verify_checkpoint is not None:
        rendered = json.dumps(result, sort_keys=True, separators=(",", ":"))
    else:
        rendered = json.dumps(result, indent=2, sort_keys=True)
    if args.output is None:
        print(rendered)
    else:
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(rendered + "\n")
        temporary.replace(args.output)
        print(
            json.dumps(
                {
                    "output": str(args.output),
                    "verified": bool(result.get("verified", True)),
                    "payload_sha256": result.get("payload_sha256")
                    or result.get("checkpoint", {}).get("payload_sha256"),
                },
                sort_keys=True,
            )
        )


if __name__ == "__main__":
    main()
