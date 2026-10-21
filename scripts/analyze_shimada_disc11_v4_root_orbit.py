"""Close the 97 certified H11/K3 walls under the exact V4 and audit the cone.

The input walls are a proven facet lower bound, not a complete wall list.
Their V4 orbit is again a proven facet lower bound because every V4 element
is an enhanced K3 automorphism preserving the complete eigen-nef chamber.
If this finite outer cone were compact, the unknown complete chamber inside
it would be compact as well.  Otherwise the script records every remaining
noncompact projective extreme ray exactly.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
from typing import Any, Sequence

from coxeter_polytopes.exact.coxeter import classify_root_pair
from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.geometric.golden_halfspaces import golden_quadratic_value
from coxeter_polytopes.search.h11_vinberg import h11_basis_gram
from coxeter_polytopes.verify.chamber import (
    ChamberVerification,
    _canonical_ray,
    _dot,
    _full_rank_kernel,
    _mat_vec,
    _pair,
    _row_rank,
    diagnose_chamber_cone,
    verify_hyperideal_extreme_ray,
)
from scripts.certify_shimada_disc11_bugaenko_v4 import (
    _bugaenko_reversal_isometry,
    _field_involution_a,
    _fold_with_word,
    _identity,
    _lift_reversal,
    _matrix_product,
    _row_product,
)
from scripts.coupled_discriminant11_candidate import field_gram
from scripts.fold_h11_seed_via_shimada import fold_point_to_eigen_nef
from scripts.frame_h11_vinberg_prefix_on_k3 import (
    DEFAULT_TRANSCRIPT,
    _extension_roots,
    _transcript_roots,
    _vinberg_strict_point,
)
from scripts.diagnose_h11_vinberg_prefix import HYPERIDEAL_WITNESS
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas
from scripts.shimada_disc11_conjugate import (
    AMPLE_FIELD_COEFFICIENTS,
    DISC11_SOURCE_SHA256,
)


DEFAULT_EXTENSION = Path("data/candidates/h11_vinberg_extension_bound21.json")
EXPECTED_TRANSFORMED_ROOTS_SHA256 = (
    "84e102afa057691bdd067b4de068a65c60e5f554d1a539ccded845da77a350e7"
)
NEW_HYPERIDEAL_INDEPENDENT_ACTIVE = (5, 6, 9, 4, 13, 12, 14, 25, 35)
NEW_HYPERIDEAL_ACTIVE = (4, 5, 6, 9, 12, 13, 14, 21, 25, 35)
NEW_HYPERIDEAL_RAY = (
    Golden(-65, -80) / 121,
    Golden(27, -4) / 121,
    Golden(-1, 3) / 11,
    Golden(-7, 10) / 11,
    Golden(-8, 13) / 11,
    Golden(-6, 7) / 11,
    Golden(-4, 1) / 11,
    Golden(-7, -1) / 11,
    Golden(-5, -7) / 11,
    Golden(-1),
)
NEW_HYPERIDEAL_NORM = Golden(978, -602) / 1331


def _canonical_sha(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _encoded_vector(vector: Sequence[Golden]) -> list[dict[str, str]]:
    return [entry.to_json() for entry in vector]


def _encoded_matrix(matrix: Sequence[Sequence[Golden]]) -> list[list[dict[str, str]]]:
    return [_encoded_vector(row) for row in matrix]


def _transported_roots(
    atlas: ShimadaDisc11Atlas,
    transcript_path: Path,
    extension_path: Path,
) -> tuple[tuple[tuple[Golden, ...], ...], dict[str, Any]]:
    transcript_bytes = transcript_path.read_bytes()
    transcript = json.loads(transcript_bytes)
    extension_bytes = extension_path.read_bytes()
    extension = json.loads(extension_bytes)
    transcript_sha = hashlib.sha256(transcript_bytes).hexdigest()
    extension_sha = hashlib.sha256(extension_bytes).hexdigest()
    if extension.get("base_prefix", {}).get("sha256") != transcript_sha:
        raise ValueError("extension is not pinned to the base H11 prefix")
    roots = (*_transcript_roots(transcript), *_extension_roots(extension))
    if len(roots) != 97:
        raise ValueError("expected the certified 97-wall prefix")
    initial_roots = tuple(tuple(-entry for entry in root) for root in roots)
    fold = fold_point_to_eigen_nef(
        atlas,
        _vinberg_strict_point(),
        max_root_reflections=200,
        max_gallery=2000,
    )
    if fold.terminal.kind != "interior":
        raise RuntimeError("the 97-wall frame did not reach a strict K3 chamber interior")
    transformed = tuple(_row_product(root, fold.row_action) for root in initial_roots)
    encoded = [_encoded_vector(root) for root in transformed]
    digest = _canonical_sha(encoded)
    if digest != EXPECTED_TRANSFORMED_ROOTS_SHA256:
        raise ValueError("transported 97-root hash changed")
    gram = tuple(tuple(entry for entry in row) for row in field_gram())
    if not all(_pair(root, gram, root) == 2 for root in transformed):
        raise ArithmeticError("a transported wall is not a norm-two root")
    if not all(
        _pair(tuple(AMPLE_FIELD_COEFFICIENTS), gram, root).sign() < 0
        for root in transformed
    ):
        raise ArithmeticError("a transported root has the wrong ample orientation")
    return transformed, {
        "transcript_sha256": transcript_sha,
        "extension_sha256": extension_sha,
        "transported_roots_sha256": digest,
        "frame_root_reflections": sum(
            record["result"] == "root_boundary" for record in fold.records
        ),
        "frame_terminal_representative": fold.terminal.state.index,
    }


def _v4(atlas: ShimadaDisc11Atlas):
    _, reversal = _bugaenko_reversal_isometry()
    lift_plus = _lift_reversal(reversal, 1)
    lift_minus = _lift_reversal(reversal, -1)
    plus, _, plus_terminal, plus_word = _fold_with_word(atlas, lift_plus)
    minus, _, minus_terminal, minus_word = _fold_with_word(atlas, lift_minus)
    a = _field_involution_a()
    if (
        plus_terminal.kind != "interior"
        or minus_terminal.kind != "interior"
        or len(plus_word) != 88
        or len(minus_word) != 88
        or _matrix_product(plus, minus) != a
        or _matrix_product(minus, plus) != a
    ):
        raise ArithmeticError("the certified V4 reconstruction failed")
    return (
        ("identity", _identity(10)),
        ("A", a),
        ("B_plus", plus),
        ("B_minus", minus),
    )


def _positive_canonical(vector: Sequence[Golden]) -> tuple[Golden, ...]:
    pivot = next(entry for entry in vector if entry)
    scale = pivot if pivot.sign() > 0 else -pivot
    return tuple(entry / scale for entry in vector)


def _extreme_ray_record(
    active: Sequence[int],
    roots: Sequence[Sequence[Golden]],
    gram: Sequence[Sequence[Golden]],
    ample: Sequence[Golden],
) -> dict[str, Any]:
    functionals = [_mat_vec(gram, root) for root in roots]
    chosen: list[int] = []
    current_rank = 0
    for index in active:
        candidate = [functionals[item] for item in (*chosen, index)]
        candidate_rank = _row_rank(candidate, 10)
        if candidate_rank > current_rank:
            chosen.append(index)
            current_rank = candidate_rank
        if current_rank == 9:
            break
    if current_rank != 9:
        raise ArithmeticError("an extreme-ray active set has rank below nine")
    ray = _full_rank_kernel([functionals[index] for index in chosen])
    if ray is None:
        raise ArithmeticError("failed to reconstruct an extreme ray")
    values = [_dot(ray, functional) for functional in functionals]
    if all(value.sign() <= 0 for value in values):
        oriented = ray
    elif all(value.sign() >= 0 for value in values):
        oriented = [-entry for entry in ray]
    else:
        raise ArithmeticError("reconstructed extreme ray violates a wall")
    canonical = _positive_canonical(oriented)
    pairings = tuple(_pair(canonical, gram, root) for root in roots)
    actual_active = tuple(index for index, value in enumerate(pairings) if not value)
    if actual_active != tuple(active):
        raise ArithmeticError("reconstructed active set changed")
    if not all(value.sign() <= 0 for value in pairings):
        raise ArithmeticError("canonical extreme ray lost its orientation")
    return {
        "active_facets": list(active),
        "independent_active_facets": chosen,
        "ray": _encoded_vector(canonical),
        "norm": _pair(canonical, gram, canonical).to_json(),
        "pairing_with_ample": _pair(canonical, gram, ample).to_json(),
        "all_wall_pairings_nonpositive": True,
        "wall_pairings_sha256": _canonical_sha(
            [value.to_json() for value in pairings]
        ),
    }


def _known_hyperideal_orbit(
    roots97: Sequence[Sequence[Golden]],
    closed_roots: Sequence[Sequence[Golden]],
    group,
    gram: Sequence[Sequence[Golden]],
    ample: Sequence[Golden],
) -> tuple[dict[str, Any], ...]:
    witness = verify_hyperideal_extreme_ray(
        roots97, gram, HYPERIDEAL_WITNESS, dimension=9
    )
    base_pairings = tuple(_pair(witness.ray, gram, root) for root in roots97)
    if all(value.sign() <= 0 for value in base_pairings):
        base_ray = witness.ray
    elif all(value.sign() >= 0 for value in base_pairings):
        base_ray = tuple(-entry for entry in witness.ray)
    else:
        raise ArithmeticError("the archived witness has mixed original-wall signs")
    rays: dict[tuple[Golden, ...], tuple[str, tuple[Golden, ...]]] = {}
    for name, action in group:
        image = _row_product(base_ray, action)
        rays.setdefault(_positive_canonical(image), (name, image))
    records = []
    for canonical, (name, _) in rays.items():
        pairings = tuple(_pair(canonical, gram, root) for root in closed_roots)
        records.append(
            {
                "V4_image_by": name,
                "ray": _encoded_vector(canonical),
                "norm": _pair(canonical, gram, canonical).to_json(),
                "pairing_with_ample": _pair(canonical, gram, ample).to_json(),
                "active_facets": [
                    index for index, value in enumerate(pairings) if not value
                ],
                "all_V4_closed_wall_pairings_nonpositive": all(
                    value.sign() <= 0 for value in pairings
                ),
                "wall_pairings_sha256": _canonical_sha(
                    [value.to_json() for value in pairings]
                ),
            }
        )
    return tuple(records)


def _active_diagram(
    active: Sequence[int],
    roots: Sequence[Sequence[Golden]],
    gram: Sequence[Sequence[Golden]],
) -> list[dict[str, Any]]:
    edges = []
    for local_i, global_i in enumerate(active):
        for local_j, global_j in enumerate(active[local_i + 1 :], local_i + 1):
            product = _pair(roots[global_i], gram, roots[global_j])
            relation = classify_root_pair(Golden(2), product, Golden(2))
            if relation.kind == "finite" and relation.label == 2:
                continue
            edges.append(
                {
                    "local_facets": [local_i, local_j],
                    "closed_root_facets": [global_i, global_j],
                    "product": product.to_json(),
                    "kind": relation.kind,
                    "label": relation.label,
                    "q": relation.q.to_json(),
                }
            )
    return edges


def _new_hyperideal_orbit_certificate(
    closed_roots: Sequence[tuple[Golden, ...]],
    provenance: dict[tuple[Golden, ...], list[tuple[int, str]]],
    group,
    gram: Sequence[Sequence[Golden]],
    ample: Sequence[Golden],
) -> dict[str, Any]:
    raw_pairings = tuple(
        _pair(NEW_HYPERIDEAL_RAY, gram, root) for root in closed_roots
    )
    raw_active = tuple(i for i, value in enumerate(raw_pairings) if not value)
    independent_rows = [
        _mat_vec(gram, closed_roots[index])
        for index in NEW_HYPERIDEAL_INDEPENDENT_ACTIVE
    ]
    raw_verified = (
        len(NEW_HYPERIDEAL_RAY) == 10
        and _pair(NEW_HYPERIDEAL_RAY, gram, NEW_HYPERIDEAL_RAY)
        == NEW_HYPERIDEAL_NORM
        and NEW_HYPERIDEAL_NORM.sign() > 0
        and all(value.sign() <= 0 for value in raw_pairings)
        and raw_active == NEW_HYPERIDEAL_ACTIVE
        and _row_rank(independent_rows, 10) == 9
    )

    orbit: dict[tuple[Golden, ...], dict[str, Any]] = {}
    for name, action in group:
        image = _positive_canonical(_row_product(NEW_HYPERIDEAL_RAY, action))
        pairings = tuple(_pair(image, gram, root) for root in closed_roots)
        active = tuple(i for i, value in enumerate(pairings) if not value)
        active_rank = _row_rank(
            [_mat_vec(gram, closed_roots[index]) for index in active], 10
        )
        orbit.setdefault(
            image,
            {
                "V4_image_by": name,
                "ray": _encoded_vector(image),
                "norm": _pair(image, gram, image).to_json(),
                "pairing_with_ample": _pair(image, gram, ample).to_json(),
                "active_facets": list(active),
                "active_rank": active_rank,
                "all_181_wall_pairings_nonpositive": all(
                    value.sign() <= 0 for value in pairings
                ),
                "wall_pairings_sha256": _canonical_sha(
                    [value.to_json() for value in pairings]
                ),
                "active_facet_provenance": [
                    {
                        "closed_root_facet": index,
                        "sources": [
                            {
                                "source_97_facet": source,
                                "V4_image_by": group_name,
                            }
                            for source, group_name in provenance[closed_roots[index]]
                        ],
                    }
                    for index in active
                ],
            },
        )
    orbit_records = list(orbit.values())
    orbit_verified = all(
        Golden.from_json(record["norm"]).sign() > 0
        and record["active_rank"] == 9
        and record["all_181_wall_pairings_nonpositive"]
        for record in orbit_records
    )
    return {
        "verified": raw_verified and orbit_verified,
        "raw_ray": _encoded_vector(NEW_HYPERIDEAL_RAY),
        "raw_norm": NEW_HYPERIDEAL_NORM.to_json(),
        "independent_active_facets": list(NEW_HYPERIDEAL_INDEPENDENT_ACTIVE),
        "actual_active_facets": list(NEW_HYPERIDEAL_ACTIVE),
        "active_diagram_edges": _active_diagram(
            NEW_HYPERIDEAL_ACTIVE, closed_roots, gram
        ),
        "V4_orbit_size": len(orbit_records),
        "V4_orbit": orbit_records,
    }


def _independent_seed_first(
    roots: Sequence[tuple[Golden, ...]],
    gram: Sequence[Sequence[Golden]],
) -> tuple[tuple[Golden, ...], ...]:
    """Avoid a combinatorial seed search in incremental double description."""

    functionals = [_mat_vec(gram, root) for root in roots]
    chosen: list[int] = []
    current_rank = 0
    for index, functional in enumerate(functionals):
        candidate = [functionals[item] for item in chosen] + [functional]
        candidate_rank = _row_rank(candidate, 10)
        if candidate_rank > current_rank:
            chosen.append(index)
            current_rank = candidate_rank
        if current_rank == 10:
            break
    if current_rank != 10:
        raise ArithmeticError("V4-closed walls do not span rank ten")
    selected = set(chosen)
    return tuple(roots[index] for index in (*chosen, *(i for i in range(len(roots)) if i not in selected)))


def _diagnose_with_progress(
    roots: Sequence[tuple[Golden, ...]],
    gram: Sequence[Sequence[Golden]],
    ample: Sequence[Golden],
) -> ChamberVerification:
    """The audited incremental DD algorithm with a per-facet progress stream."""

    ambient_rank = 10
    functionals = [_mat_vec(gram, root) for root in roots]
    seed = tuple(range(ambient_rank))
    if _row_rank([functionals[index] for index in seed], ambient_rank) != ambient_rank:
        raise ArithmeticError("the scheduled first ten walls are not independent")
    processed = list(seed)
    rays: dict[tuple[Golden, ...], tuple[list[Golden], frozenset[int]]] = {}
    for omitted in seed:
        active = tuple(index for index in seed if index != omitted)
        ray = _full_rank_kernel([functionals[index] for index in active])
        if ray is None:
            raise ArithmeticError("independent seed has a dependent proper subset")
        values = [_dot(ray, functionals[index]) for index in seed]
        if all(value.sign() <= 0 for value in values):
            oriented = ray
        elif all(value.sign() >= 0 for value in values):
            oriented = [-entry for entry in ray]
        else:
            raise ArithmeticError("seed wall cone has inconsistent orientation")
        rays[_canonical_ray(oriented)] = (
            oriented,
            frozenset(index for index in seed if not _dot(oriented, functionals[index])),
        )

    started = time.monotonic()
    rank_cache: dict[frozenset[int], int] = {}
    for new_facet in range(ambient_rank, len(roots)):
        functional = functionals[new_facet]
        negative = []
        zero = []
        positive = []
        for ray, active in rays.values():
            value = _dot(ray, functional)
            target = negative if value.sign() < 0 else positive if value.sign() > 0 else zero
            target.append((ray, active, value))
        updated: dict[tuple[Golden, ...], tuple[list[Golden], frozenset[int]]] = {}
        for ray, active, _ in (*negative, *zero):
            new_active = active | ({new_facet} if not _dot(ray, functional) else set())
            updated[_canonical_ray(ray)] = (ray, frozenset(new_active))
        adjacent_pairs = 0
        for negative_ray, negative_active, negative_value in negative:
            for positive_ray, positive_active, positive_value in positive:
                common = negative_active & positive_active
                common_rank = rank_cache.get(common)
                if common_rank is None:
                    common_rank = _row_rank(
                        [functionals[index] for index in common], ambient_rank
                    )
                    rank_cache[common] = common_rank
                if common_rank != ambient_rank - 2:
                    continue
                adjacent_pairs += 1
                ray = [
                    positive_value * negative_ray[i]
                    - negative_value * positive_ray[i]
                    for i in range(ambient_rank)
                ]
                key = _canonical_ray(ray)
                actual = frozenset(
                    index
                    for index in (*processed, new_facet)
                    if not _dot(ray, functionals[index])
                )
                updated[key] = (ray, actual)
        rays = updated
        processed.append(new_facet)
        print(
            json.dumps(
                {
                    "processed_facets": len(processed),
                    "rays": len(rays),
                    "kept": len(negative),
                    "cut": len(positive),
                    "boundary": len(zero),
                    "new_adjacent_pairs": adjacent_pairs,
                    "elapsed_seconds": round(time.monotonic() - started, 3),
                },
                sort_keys=True,
            ),
            file=sys.stderr,
            flush=True,
        )

    finite: set[tuple[int, ...]] = set()
    cusps: set[tuple[int, ...]] = set()
    hyperideal: set[tuple[int, ...]] = set()
    opposite: set[tuple[int, ...]] = set()
    for ray, active_set in rays.values():
        active = tuple(sorted(active_set))
        norm_sign = _pair(ray, gram, ray).sign()
        sheet_sign = _pair(ray, gram, ample).sign()
        if norm_sign < 0:
            (finite if sheet_sign < 0 else opposite).add(active)
        elif norm_sign == 0:
            (cusps if sheet_sign < 0 else opposite).add(active)
        else:
            hyperideal.add(active)
    return ChamberVerification(
        dimension=9,
        facet_count=len(roots),
        finite_vertices=tuple(sorted(finite)),
        cusps=tuple(sorted(cusps)),
        hyperideal_rays=tuple(sorted(hyperideal)),
        opposite_sheet_rays=tuple(sorted(opposite)),
    )


def _cdd_discover_exact_escape_rays(
    roots: Sequence[tuple[Golden, ...]],
    gram: Sequence[Sequence[Golden]],
    ample: Sequence[Golden],
) -> dict[str, Any]:
    """Use floating cdd only to propose active sets, then certify them exactly."""

    functionals = [_mat_vec(gram, root) for root in roots]
    # A deep convergent to phi; cdd on this installation emits floating
    # output, so neither this substitution nor cdd is trusted in the proof.
    phi_approximation = Fraction(267914296 + 165580141, 267914296)

    def rational_approximation(value: Golden) -> Fraction:
        return value.a + value.b * phi_approximation

    rows = [
        "0 "
        + " ".join(
            str(-rational_approximation(entry)) for entry in functional
        )
        for functional in functionals
    ]
    cdd_input = (
        "H-representation\nbegin\n"
        f"{len(rows)} 11 rational\n"
        + "\n".join(rows)
        + "\nend\n"
    )
    completed = subprocess.run(
        ["cddexec", "--rep"],
        input=cdd_input,
        text=True,
        capture_output=True,
        check=True,
    )
    lines = completed.stdout.splitlines()
    begin = next(i for i, line in enumerate(lines) if line.strip() == "begin")
    header = lines[begin + 1].split()
    ray_count = int(header[0])
    numeric_rays = []
    for line in lines[begin + 2 : begin + 2 + ray_count]:
        values = [float(value) for value in line.split()]
        if len(values) != 11 or values[0] != 0.0:
            continue
        numeric_rays.append(values[1:])

    exact: dict[tuple[Golden, ...], dict[str, Any]] = {}
    phi_float = (1.0 + 5.0**0.5) / 2.0
    float_functionals = [
        [float(entry.a) + float(entry.b) * phi_float for entry in functional]
        for functional in functionals
    ]
    for numeric in numeric_rays:
        approximate_values = [
            sum(numeric[j] * functional[j] for j in range(10))
            for functional in float_functionals
        ]
        nearest = sorted(range(len(roots)), key=lambda i: abs(approximate_values[i]))
        chosen: list[int] = []
        current_rank = 0
        for index in nearest[:40]:
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
        canonical = _positive_canonical(oriented)
        norm = _pair(canonical, gram, canonical)
        sheet = _pair(canonical, gram, ample)
        if norm.sign() < 0 and sheet.sign() < 0:
            continue
        exact_pairings = tuple(_pair(canonical, gram, root) for root in roots)
        record = {
            "active_facets": [
                index for index, value in enumerate(exact_pairings) if not value
            ],
            "independent_active_facets": chosen,
            "ray": _encoded_vector(canonical),
            "norm": norm.to_json(),
            "pairing_with_ample": sheet.to_json(),
            "kind": (
                "hyperideal"
                if norm.sign() > 0
                else "cusp"
                if norm.sign() == 0 and sheet.sign() < 0
                else "opposite_sheet"
            ),
            "all_wall_pairings_nonpositive": all(
                value.sign() <= 0 for value in exact_pairings
            ),
            "wall_pairings_sha256": _canonical_sha(
                [value.to_json() for value in exact_pairings]
            ),
        }
        exact.setdefault(canonical, record)
    return {
        "method": "floating_cdd_active_sets_followed_by_exact_Qphi_verification",
        "floating_output_used_as_proof": False,
        "floating_cdd_ray_count": ray_count,
        "parsed_floating_ray_count": len(numeric_rays),
        "exact_noncompact_witness_count": len(exact),
        "exact_noncompact_witnesses": list(exact.values()),
    }


def analyze(
    disc11: Path,
    transcript_path: Path = DEFAULT_TRANSCRIPT,
    extension_path: Path = DEFAULT_EXTENSION,
    *,
    full_census: bool = False,
    progress: bool = False,
    cdd_ray_search: bool = False,
) -> dict[str, Any]:
    source_sha = hashlib.sha256(disc11.read_bytes()).hexdigest()
    if source_sha != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")
    atlas = ShimadaDisc11Atlas(disc11)
    roots97, frame = _transported_roots(atlas, transcript_path, extension_path)
    group = _v4(atlas)
    gram = tuple(tuple(entry for entry in row) for row in field_gram())
    ample = tuple(AMPLE_FIELD_COEFFICIENTS)

    closed: dict[tuple[Golden, ...], list[tuple[int, str]]] = {}
    for source_index, root in enumerate(roots97):
        for name, action in group:
            image = _row_product(root, action)
            if _pair(image, gram, image) != 2:
                raise ArithmeticError("a V4 root image lost norm two")
            if _pair(ample, gram, image).sign() >= 0:
                raise ArithmeticError("a V4 root image lost its chamber orientation")
            closed.setdefault(image, []).append((source_index, name))
    closed_roots = tuple(closed)
    index = {root: i for i, root in enumerate(closed_roots)}
    permutations = {
        name: tuple(index[_row_product(root, action)] for root in closed_roots)
        for name, action in group
    }
    action_verified = (
        all(
            permutations[name][permutations[name][i]] == i
            for name, _ in group
            for i in range(len(closed_roots))
        )
        and all(
            permutations["B_plus"][permutations["B_minus"][i]]
            == permutations["A"][i]
            for i in range(len(closed_roots))
        )
    )

    pairwise_acute = all(
        _pair(left, gram, right).sign() <= 0
        for i, left in enumerate(closed_roots)
        for right in closed_roots[i + 1 :]
    )
    if not pairwise_acute:
        raise ArithmeticError("V4 closure is not an acute simple-root system")
    if golden_quadratic_value(gram, ample).sign() >= 0:
        raise ArithmeticError("fixed ample point is not timelike")

    known_hyperideal_orbit = _known_hyperideal_orbit(
        roots97, closed_roots, group, gram, ample
    )
    known_witness_survives = any(
        record["all_V4_closed_wall_pairings_nonpositive"]
        for record in known_hyperideal_orbit
    )
    witness_orbit_consistent = len(
        {
            bool(record["all_V4_closed_wall_pairings_nonpositive"])
            for record in known_hyperideal_orbit
        }
    ) == 1
    new_hyperideal = _new_hyperideal_orbit_certificate(
        closed_roots, closed, group, gram, ample
    )

    geometry = None
    noncompact_records = None
    cdd_search = (
        _cdd_discover_exact_escape_rays(closed_roots, gram, ample)
        if cdd_ray_search
        else None
    )
    if full_census:
        distance_ordered = tuple(
            sorted(closed_roots, key=lambda root: -_pair(ample, gram, root))
        )
        census_roots = _independent_seed_first(distance_ordered, gram)
        geometry = (
            _diagnose_with_progress(census_roots, gram, ample)
            if progress
            else diagnose_chamber_cone(census_roots, gram, ample, dimension=9)
        )
        noncompact_records = {
            "cusps": [
                _extreme_ray_record(active, census_roots, gram, ample)
                for active in geometry.cusps
            ],
            "hyperideal_rays": [
                _extreme_ray_record(active, census_roots, gram, ample)
                for active in geometry.hyperideal_rays
            ],
            "opposite_sheet_rays": [
                _extreme_ray_record(active, census_roots, gram, ample)
                for active in geometry.opposite_sheet_rays
            ],
        }
    encoded_closed = [_encoded_vector(root) for root in closed_roots]
    orbit_sizes: dict[int, int] = {}
    seen: set[int] = set()
    for start in range(len(closed_roots)):
        if start in seen:
            continue
        orbit = {permutation[start] for permutation in permutations.values()}
        seen.update(orbit)
        orbit_sizes[len(orbit)] = orbit_sizes.get(len(orbit), 0) + 1
    compact = (
        not any(noncompact_records.values())
        if noncompact_records is not None
        else False
    )
    return {
        "verified": (
            action_verified
            and pairwise_acute
            and (geometry is None or geometry.facet_count == len(closed_roots))
            and len(seen) == len(closed_roots)
            and witness_orbit_consistent
            and new_hyperideal["verified"]
        ),
        "source_sha256": source_sha,
        **frame,
        "input_facet_count": len(roots97),
        "V4_closed_distinct_facet_count": len(closed_roots),
        "closed_roots_sha256": _canonical_sha(encoded_closed),
        "V4_action_verified": action_verified,
        "V4_orbit_size_counts": {
            str(size): count for size, count in sorted(orbit_sizes.items())
        },
        "pairwise_acute": pairwise_acute,
        "archived_hyperideal_witness_survives_V4_closure": known_witness_survives,
        "archived_hyperideal_witness_orbit_consistent": witness_orbit_consistent,
        "archived_hyperideal_witness_V4_orbit": list(known_hyperideal_orbit),
        "new_exact_hyperideal_witness": new_hyperideal,
        "cdd_proposed_exact_escape_search": cdd_search,
        "outer_cone": {
            "full_census_run": geometry is not None,
            "finite_ray_count": len(geometry.finite_vertices) if geometry else None,
            "cusp_count": len(geometry.cusps) if geometry else None,
            "hyperideal_ray_count": len(geometry.hyperideal_rays) if geometry else None,
            "certified_hyperideal_ray_count_lower_bound": new_hyperideal[
                "V4_orbit_size"
            ],
            "opposite_sheet_ray_count": len(geometry.opposite_sheet_rays) if geometry else None,
            "compact": compact,
            "finite_vertex_active_facets": [
                list(active) for active in geometry.finite_vertices
            ] if geometry else None,
            "noncompact_witnesses": noncompact_records,
        },
        "conclusion": (
            "V4_orbit_outer_cone_is_compact_so_complete_eigen_nef_chamber_is_compact"
            if compact is True
            else "V4_orbit_outer_cone_still_has_exact_hyperideal_escape_rays"
        ),
        "scope": "exact_V4_closure_of_97_certified_facets_not_global_wall_completeness",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    parser.add_argument("--transcript", type=Path, default=DEFAULT_TRANSCRIPT)
    parser.add_argument("--extension", type=Path, default=DEFAULT_EXTENSION)
    parser.add_argument("--full-census", action="store_true")
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--cdd-ray-search", action="store_true")
    args = parser.parse_args()
    print(
        json.dumps(
            analyze(
                args.disc11,
                args.transcript,
                args.extension,
                full_census=args.full_census,
                progress=args.progress,
                cdd_ray_search=args.cdd_ray_search,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
