"""Exact global facet reduction for the two cyclic Fermat-char-3 survivors."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
import sympy as sp

from scripts.search_apery_fermi_chamber_folds import _exact_nonnegative_dependence
from scripts.search_ss_fermat_quartic_char3_folds import (
    _coxeter_label,
    _primitive_positive as primitive_positive,
    parse_source,
)


ROOT = Path(__file__).resolve().parents[1]
CENSUS = ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_element_class_screen.json"
OUTPUT = ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_survivor_screen.json"


def wall_cast(matrix: np.ndarray) -> dict[str, object]:
    source = parse_source()
    walls = np.asarray(source["walls"], dtype=np.int64)
    wall_types = source["wall_types"]
    ample = np.asarray(source["ample"], dtype=np.int64).reshape(-1)
    gram = sp.Matrix(source["gram"])
    lookup = {tuple(map(int, wall)): index for index, wall in enumerate(walls)}
    permutation = tuple(
        lookup[tuple(map(int, image))] for image in walls @ matrix.T
    )
    unseen = set(range(len(walls)))
    orbits: list[tuple[int, ...]] = []
    while unseen:
        index = min(unseen)
        image = permutation[index]
        orbit = (index,) if image == index else (index, image)
        if any(permutation[wall] not in orbit for wall in orbit):
            raise ArithmeticError("not an involution wall orbit")
        unseen.difference_update(orbit)
        orbits.append(orbit)
    ray_sources: dict[tuple[int, ...], list[int]] = {}
    for orbit_index, orbit in enumerate(orbits):
        ray = primitive_positive(np.sum(walls[list(orbit)], axis=0), ample)
        ray_sources.setdefault(ray, []).append(orbit_index)
    rays = tuple(ray_sources)
    restricted = sp.Matrix(rays)
    _, coordinate_columns = restricted.rref()
    coordinates = restricted[:, list(coordinate_columns)]
    fixed_rank = gram.rows - (sp.Matrix(matrix) - sp.eye(gram.rows)).rank()
    if coordinates.cols != fixed_rank or coordinates.rank() != fixed_rank:
        raise ArithmeticError("coordinate projection is not injective on the fixed span")
    return {
        "gram": gram,
        "ample": ample,
        "walls": walls,
        "wall_types": wall_types,
        "orbits": orbits,
        "rays": rays,
        "ray_sources": ray_sources,
        "restricted": restricted,
        "coordinates": coordinates,
    }


def fraction_text(value: sp.Rational | Fraction | int) -> str:
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def strict_facets(cast: dict[str, object]) -> tuple[int, ...]:
    gram: sp.Matrix = cast["gram"]  # type: ignore[assignment]
    restricted: sp.Matrix = cast["restricted"]  # type: ignore[assignment]
    ample = np.asarray(cast["ample"], dtype=np.int64)
    inverse = gram.inv()
    denominator = int(sp.ilcm(*[sp.Rational(value).q for value in inverse]))
    metric = np.asarray(inverse * denominator, dtype=np.int64)
    rays = np.asarray(restricted.tolist(), dtype=np.int64)
    ray_gram = (rays @ metric) @ rays.T
    norms = np.diag(ray_gram)
    pairings = rays @ ample
    facets = []
    for index, norm_value in enumerate(norms):
        norm = int(norm_value)
        if norm >= 0:
            continue
        margins = int(pairings[index]) * ray_gram[index] - norm * pairings
        if np.all(np.delete(margins, index) > 0):
            facets.append(index)
    return tuple(facets)


def exact_global_reduction(
    cast: dict[str, object], facets: tuple[int, ...]
) -> tuple[list[dict[str, object]], tuple[int, ...]]:
    coordinates: sp.Matrix = cast["coordinates"]  # type: ignore[assignment]
    facet_set = set(facets)
    certificates: list[dict[str, object]] = []
    unresolved: list[int] = []
    for target in range(coordinates.rows):
        if target in facet_set:
            continue
        local_indices = [*facets, target]
        dependence = _exact_nonnegative_dependence(
            coordinates[local_indices, :], len(facets)
        )
        if dependence is None:
            unresolved.append(target)
            continue
        certificates.append(
            {
                "restricted_ray_index": target + 1,
                "terms": [
                    {
                        "restricted_ray_index": facets[support] + 1,
                        "coefficient": fraction_text(coefficient),
                    }
                    for support, coefficient in dependence
                ],
            }
        )
    return certificates, tuple(unresolved)


def positive_facet_gram(
    cast: dict[str, object], facets: tuple[int, ...]
) -> sp.Matrix:
    gram: sp.Matrix = cast["gram"]  # type: ignore[assignment]
    restricted: sp.Matrix = cast["restricted"]  # type: ignore[assignment]
    normals = restricted[list(facets), :]
    return -(normals * gram.inv() * normals.T)


def strict_witness_certificates(
    cast: dict[str, object], facets: tuple[int, ...]
) -> list[dict[str, object]]:
    gram: sp.Matrix = cast["gram"]  # type: ignore[assignment]
    restricted: sp.Matrix = cast["restricted"]  # type: ignore[assignment]
    ample = np.asarray(cast["ample"], dtype=np.int64)
    inverse = gram.inv()
    denominator = int(sp.ilcm(*[sp.Rational(value).q for value in inverse]))
    metric = np.asarray(inverse * denominator, dtype=np.int64)
    rays = np.asarray(restricted.tolist(), dtype=np.int64)
    ray_gram = (rays @ metric) @ rays.T
    norms = np.diag(ray_gram)
    pairings = rays @ ample
    ample_square = int((sp.Matrix(1, gram.rows, ample.tolist()) * gram * sp.Matrix(ample.tolist()))[0])
    certificates = []
    for index in facets:
        norm = int(norms[index])
        margins = int(pairings[index]) * ray_gram[index] - norm * pairings
        other = np.delete(margins, index)
        witness_square = Fraction(
            ample_square * norm - int(pairings[index]) ** 2 * denominator,
            norm,
        )
        if norm >= 0 or witness_square <= 0 or int(other.min()) <= 0:
            raise ArithmeticError("strict facet witness failed exact replay")
        certificates.append(
            {
                "restricted_ray_index": index + 1,
                "normal_square": fraction_text(Fraction(norm, denominator)),
                "ample_pairing": int(pairings[index]),
                "witness_square": fraction_text(witness_square),
                "minimum_integer_strict_margin": int(other.min()),
            }
        )
    return certificates


def coxeter_pair_gate(positive: sp.Matrix) -> dict[str, object]:
    denominator = int(sp.ilcm(*[sp.Rational(value).q for value in positive]))
    integral = np.asarray(positive * denominator, dtype=np.int64)
    labels: Counter[str] = Counter()
    bad: list[dict[str, object]] = []
    for left in range(positive.rows):
        for right in range(left):
            # The helper uses negative-square inward normals.  Negating the
            # positive-wall Gram restores that convention.
            product = -int(integral[left, right])
            left_norm = -int(integral[left, left])
            right_norm = -int(integral[right, right])
            label = _coxeter_label(product, left_norm, right_norm)
            if label is None:
                bad.append(
                    {
                        "facet_positions": [right + 1, left + 1],
                        "normalized_squared_product": fraction_text(
                            Fraction(product * product, left_norm * right_norm)
                        ),
                    }
                )
            else:
                labels[label] += 1
    return {
        "pair_count": positive.rows * (positive.rows - 1) // 2,
        "label_counts": dict(sorted(labels.items())),
        "bad_pair_count": len(bad),
        "bad_pairs": bad,
        "all_pairs_coxeter": not bad,
    }


def main() -> None:
    census = json.loads(CENSUS.read_text(encoding="utf-8"))
    survivors = [
        record
        for record in census["high_fixed_rank_casts"]
        if int(record["class_size"]) in {4536, 9072}
    ]
    if len(survivors) != 2:
        raise ArithmeticError(f"expected two survivors, got {len(survivors)}")
    archived_results = []
    for record in sorted(survivors, key=lambda item: item["fixed_rank"]):
        identity = (record["class_size"], record["fixed_rank"])
        print("CAST", identity, flush=True)
        representative = np.asarray(record["representative_matrix"], dtype=np.int64)
        cast = wall_cast(representative)
        strict = strict_facets(cast)
        witnesses = strict_witness_certificates(cast, strict)
        certificates, unresolved = exact_global_reduction(cast, strict)
        exact = strict if not unresolved else ()
        positive = positive_facet_gram(cast, strict)
        pair_gate = coxeter_pair_gate(positive)
        if not pair_gate["all_pairs_coxeter"]:
            raise ArithmeticError("a frozen survivor failed the exact Coxeter pair gate")
        print(
            "REDUCTION",
            identity,
            "rays",
            len(cast["rays"]),
            "strict",
            len(strict),
            "exact",
            len(exact),
            "redundant",
            len(certificates),
            "unresolved",
            len(unresolved),
            flush=True,
        )
        if unresolved:
            continue
        volume_input = {
            "status": {"all_facet_pairs_are_coxeter": True},
            "group_class": {
                "class_size": record["class_size"],
                "fixed_rank": record["fixed_rank"],
                "hyperbolic_dimension": record["hyperbolic_dimension"],
            },
            "facet_count": len(strict),
            "facet_restricted_ray_indices": [index + 1 for index in strict],
            "positive_wall_gram": [
                [fraction_text(value) for value in positive.row(row)]
                for row in range(positive.rows)
            ],
        }
        volume_path = (
            ROOT
            / "data/enumerations"
            / f"ss_fermat_quartic_char3_c2_h{record['hyperbolic_dimension']}_coxeter_screen.json"
        )
        volume_path.write_text(
            json.dumps(volume_input, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        archived_results.append(
            {
                "class_size": record["class_size"],
                "fixed_rank": record["fixed_rank"],
                "hyperbolic_dimension": record["hyperbolic_dimension"],
                "representative_archived_ns_matrix": representative.tolist(),
                "wall_orbit_count": len(cast["orbits"]),
                "candidate_ray_count": len(cast["rays"]),
                "facet_count": len(strict),
                "facet_restricted_ray_indices": [index + 1 for index in strict],
                "strict_facet_witnesses": witnesses,
                "redundant_ray_count": len(certificates),
                "redundancy_certificates": certificates,
                "unresolved_count": len(unresolved),
                "exact_complete": not unresolved,
                "coxeter_pair_gate": pair_gate,
                "coxeter_screen": str(volume_path.relative_to(ROOT)),
                "decision": "complete_irredundant_coxeter_polytope_requires_volume_audit",
            }
        )
    OUTPUT.write_text(
        json.dumps(
            {
                "schema": "ss_fermat_quartic_char3_c2_survivor_screen_v1",
                "source_full_stabilizer_census": str(CENSUS.relative_to(ROOT)),
                "surface_automorphism_required": False,
                "all_5944_published_walls_screened": True,
                "artifact_stops_before_volume_gate": True,
                "casts": archived_results,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print("WROTE", OUTPUT)


if __name__ == "__main__":
    main()
