"""Certify the common (-4)-reflection mechanism in three Shimada H18 stocks.

The fresh weighted-wall searches for discriminants 11, 15, and 16 produced
large families of rational involutions with fixed rank 19.  This script
replays every such H18 case and proves that its nonidentity action is exactly
reflection in a primitive divisor vector ``v`` satisfying

    v^2 = -4,   div(v) = 1.

Consequently ``NS intersect v^perp`` has absolute discriminant ``4*d`` when
``NS`` has absolute discriminant ``d``.  The action itself has denominator 2,
so this is a divisor-lattice construction, not an integral surface
automorphism.  The script also ties every member to its existing exact
non-Coxeter certificate; it makes no finite-volume claim.
"""

from __future__ import annotations

from collections import Counter
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.inventory_shimada_oriented_weighted_wall_symmetries import CONFIGS
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_three_discriminant_h18_minus4_reflection_family.json"
)
SCHEMA = "coxeter.shimada_three_discriminant_h18_minus4_reflection_family.v1"
AMBIENT_RANK = 20
EXPECTED_COUNTS = {11: 216, 15: 144, 16: 574}
EXPECTED_SMITH = {
    11: [1] * 18 + [44],
    15: [1] * 18 + [60],
    16: [1] * 16 + [2, 4, 8],
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _matrix(payload: list[list[str]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _matrix_denominator(matrix: sp.Matrix) -> int:
    return int(sp.ilcm(*(sp.Rational(value).q for value in matrix)))


def _primitive_row(row: sp.Matrix) -> sp.Matrix:
    if row.rows != 1 or not any(value != 0 for value in row):
        raise ArithmeticError("a nonzero rational row is required")
    denominator = int(sp.ilcm(*(sp.Rational(value).q for value in row)))
    entries = [int(denominator * sp.Rational(value)) for value in row]
    divisor = math.gcd(*(abs(value) for value in entries))
    entries = [value // divisor for value in entries]
    first = next(value for value in entries if value)
    if first < 0:
        entries = [-value for value in entries]
    return sp.Matrix(1, len(entries), entries)


def _action_for_generic_case(
    chamber_record: dict[str, Any], case_payload: dict[str, Any]
) -> sp.Matrix:
    subgroup = {
        tuple(map(int, permutation))
        for permutation in case_payload["case"]["subgroup_wall_permutations"]
    }
    if len(subgroup) != 2:
        raise ArithmeticError("an H18 family member is not an involution subgroup")
    identity = tuple(range(int(case_payload["case"]["published_wall_count"])))
    nonidentity = next(permutation for permutation in subgroup if permutation != identity)
    matches = [
        record
        for record in chamber_record["oriented_weighted_wall_group"]["all_actions"]
        if tuple(map(int, record["permutation"])) == nonidentity
    ]
    if len(matches) != 1:
        raise ArithmeticError("the involution has no unique inventory lift")
    return _matrix(matches[0]["primal_row_action"])


def _certify_case(
    *,
    discriminant: int,
    manifest_record: dict[str, Any],
    case_payload: dict[str, Any],
    ambient_gram: sp.Matrix,
    action: sp.Matrix,
) -> dict[str, Any]:
    case_id = str(manifest_record["id"])
    identity = sp.eye(AMBIENT_RANK)
    difference = action - identity
    if difference.rank() != 1:
        raise ArithmeticError(f"{case_id}: the involution does not have rank-one defect")
    vector = _primitive_row(next(difference.row(i) for i in range(AMBIENT_RANK) if any(difference.row(i))))
    square = sp.Rational((vector * ambient_gram * vector.T)[0])
    covector = vector * ambient_gram
    divisibility = math.gcd(*(abs(int(value)) for value in covector))
    reflection = identity - sp.Rational(2, square) * ambient_gram * vector.T * vector

    fixed = case_payload["saturated_fixed_divisor_lattice"]
    basis = _matrix(fixed["basis_primal_rows"])
    fixed_gram = _matrix(fixed["gram"])
    smith = list(map(int, fixed["smith_diagonal_absolute"]))
    fixed_determinant = abs(int(fixed["determinant"]))
    predicted_determinant = (
        discriminant * abs(int(square)) // (divisibility * divisibility)
    )
    non_coxeter_pairs = int(case_payload["coxeter_gate"]["non_coxeter_pair_count"])

    checks = {
        "ambient_determinant": int(ambient_gram.det()) == -discriminant,
        "subgroup_is_C2": int(manifest_record["subgroup_order"]) == 2,
        "fixed_rank_is_19": int(manifest_record["fixed_rank"]) == 19,
        "hyperbolic_dimension_is_18": int(manifest_record["hyperbolic_dimension"]) == 18,
        "action_is_involution": action * action == identity,
        "action_is_Lorentz_isometry": action * ambient_gram * action.T == ambient_gram,
        "action_has_rank_one_defect": difference.rank() == 1,
        "primitive_vector_is_anti_invariant": vector * action == -vector,
        "primitive_vector_square_is_minus_4": square == -4,
        "primitive_vector_divisibility_is_1": divisibility == 1,
        "action_equals_exact_vector_reflection": action == reflection,
        "action_denominator_is_2": _matrix_denominator(action) == 2,
        "action_is_not_integral": any(not value.is_Integer for value in action),
        "fixed_basis_rank_is_19": basis.rank() == 19,
        "fixed_basis_is_pointwise_fixed": basis * action == basis,
        "fixed_basis_is_vector_orthogonal": basis * ambient_gram * vector.T == sp.zeros(19, 1),
        "fixed_gram_replays": basis * ambient_gram * basis.T == fixed_gram,
        "fixed_determinant_matches_orthogonal_formula": fixed_determinant == predicted_determinant,
        "fixed_determinant_is_4d": fixed_determinant == 4 * discriminant,
        "fixed_smith_form_matches_family": smith == EXPECTED_SMITH[discriminant],
        "existing_case_is_not_Coxeter": case_payload["decision"]["classification"] == "not_coxeter",
        "strict_forbidden_Coxeter_pair_exists": non_coxeter_pairs >= 1,
    }
    failed = [name for name, value in checks.items() if not value]
    if failed:
        raise ArithmeticError(f"{case_id}: failed checks {failed}")

    return {
        "id": case_id,
        "discriminant": discriminant,
        "chamber_index": int(manifest_record["chamber_index"]),
        "class_index": int(manifest_record["class_index"]),
        "primitive_minus4_vector_primal_row": [int(value) for value in vector],
        "vector_square": int(square),
        "vector_divisibility": divisibility,
        "reflection_action_denominator": _matrix_denominator(action),
        "fixed_lattice_determinant_absolute": fixed_determinant,
        "fixed_lattice_smith_diagonal_absolute": smith,
        "certified_non_Coxeter_pair_count": non_coxeter_pairs,
        "checks": checks,
    }


def _generic_records(discriminant: int) -> tuple[list[dict[str, Any]], dict[str, str]]:
    inventory_path = (
        ROOT
        / "data"
        / "enumerations"
        / f"shimada_disc{discriminant}_oriented_weighted_wall_symmetry_inventory.json"
    )
    manifest_path = (
        ROOT
        / "data"
        / "enumerations"
        / f"shimada_disc{discriminant}_oriented_weighted_high_section_census.json"
    )
    case_directory = (
        ROOT
        / "data"
        / "enumerations"
        / f"shimada_disc{discriminant}_oriented_weighted_high_sections"
    )
    inventory = json.loads(inventory_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    chambers = {
        int(record["chamber_index"]): record for record in inventory["chambers"]
    }
    atlas = ShimadaDisc11Atlas(CONFIGS[discriminant]["source"])
    ambient_gram = sp.Matrix(atlas.gram)
    records = []
    for item in manifest["cases"]:
        if int(item["hyperbolic_dimension"]) != 18:
            continue
        case_path = ROOT / item["case_certificate"]
        payload = json.loads(case_path.read_text())
        if not (
            payload["case"]["weighted_group_strictly_extends_archived_AutD"]
            and not payload["case"]["representative_contained_in_archived_AutD"]
        ):
            raise ArithmeticError(f"{item['id']}: H18 member is not a fresh class")
        action = _action_for_generic_case(chambers[int(item["chamber_index"])], payload)
        records.append(
            _certify_case(
                discriminant=discriminant,
                manifest_record=item,
                case_payload=payload,
                ambient_gram=ambient_gram,
                action=action,
            )
        )
    return records, {
        "primary_atlas": str(CONFIGS[discriminant]["source"]),
        "primary_atlas_sha256": _sha256(CONFIGS[discriminant]["source"]),
        "inventory": str(inventory_path.relative_to(ROOT)),
        "inventory_sha256": _sha256(inventory_path),
        "manifest": str(manifest_path.relative_to(ROOT)),
        "manifest_sha256": _sha256(manifest_path),
        "case_directory": str(case_directory.relative_to(ROOT)),
    }


def _disc16_records() -> tuple[list[dict[str, Any]], dict[str, str]]:
    discriminant = 16
    manifest_path = (
        ROOT
        / "data"
        / "enumerations"
        / "shimada_disc16_fresh_weighted_high_section_census.json"
    )
    case_directory = (
        ROOT
        / "data"
        / "enumerations"
        / "shimada_disc16_fresh_weighted_high_sections"
    )
    inventory_path = (
        ROOT
        / "data"
        / "enumerations"
        / "shimada_disc16_weighted_wall_symmetry_inventory.json"
    )
    manifest = json.loads(manifest_path.read_text())
    atlas = ShimadaDisc11Atlas(CONFIGS[discriminant]["source"])
    ambient_gram = sp.Matrix(atlas.gram)
    records = []
    for item in manifest["cases"]:
        if int(item["hyperbolic_dimension"]) != 18:
            continue
        case_path = ROOT / item["case_certificate"]
        payload = json.loads(case_path.read_text())
        actions = payload["case"]["nonidentity_rational_actions"]
        if len(actions) != 1:
            raise ArithmeticError(f"{item['id']}: H18 member has no unique involution")
        action = _matrix(actions[0]["primal_row_action"])
        records.append(
            _certify_case(
                discriminant=discriminant,
                manifest_record=item,
                case_payload=payload,
                ambient_gram=ambient_gram,
                action=action,
            )
        )
    return records, {
        "primary_atlas": str(CONFIGS[discriminant]["source"]),
        "primary_atlas_sha256": _sha256(CONFIGS[discriminant]["source"]),
        "inventory": str(inventory_path.relative_to(ROOT)),
        "inventory_sha256": _sha256(inventory_path),
        "manifest": str(manifest_path.relative_to(ROOT)),
        "manifest_sha256": _sha256(manifest_path),
        "case_directory": str(case_directory.relative_to(ROOT)),
    }


def main() -> None:
    all_records: list[dict[str, Any]] = []
    sources: dict[str, dict[str, str]] = {}
    for discriminant in (11, 15):
        records, source = _generic_records(discriminant)
        all_records.extend(records)
        sources[str(discriminant)] = source
        print(json.dumps({"discriminant": discriminant, "certified": len(records)}))
    records, source = _disc16_records()
    all_records.extend(records)
    sources["16"] = source
    print(json.dumps({"discriminant": 16, "certified": len(records)}))

    count_by_discriminant = Counter(record["discriminant"] for record in all_records)
    if dict(count_by_discriminant) != EXPECTED_COUNTS:
        raise ArithmeticError(
            f"family count changed: {dict(count_by_discriminant)} != {EXPECTED_COUNTS}"
        )
    vector_count_by_discriminant = {
        str(discriminant): len(
            {
                tuple(record["primitive_minus4_vector_primal_row"])
                for record in all_records
                if record["discriminant"] == discriminant
            }
        )
        for discriminant in EXPECTED_COUNTS
    }
    bad_pair_histogram = Counter(
        record["certified_non_Coxeter_pair_count"] for record in all_records
    )
    payload = {
        "schema": SCHEMA,
        "scope": {
            "description": "all fresh H18 classes in the completed discriminant 11, 15, and 16 weighted-wall searches",
            "surface_automorphism_required": False,
            "finite_volume_claimed": False,
            "ordinary_hyperbolic_Coxeter_polytope_claimed": False,
        },
        "sources": sources,
        "theorem": {
            "statement": "Every one of the 934 fresh H18 fixed sections is the fixed hyperplane of a rational reflection in a primitive divisor vector v with v^2=-4 and div(v)=1.",
            "fixed_lattice_discriminant_formula": "|disc(NS intersect v^perp)| = |disc(NS)|*|v^2|/div(v)^2 = 4d",
            "integrality_obstruction": "The reflection action has denominator 2 in every case and therefore is not an integral Neron-Severi isometry.",
            "Coxeter_outcome": "Every section has at least one exact forbidden facet-pair ratio and hence is not an ordinary Coxeter polytope.",
        },
        "summary": {
            "case_count": len(all_records),
            "case_count_by_discriminant": {
                str(key): value for key, value in sorted(count_by_discriminant.items())
            },
            "distinct_primitive_vector_count_by_discriminant": vector_count_by_discriminant,
            "reflection_denominator_histogram": dict(
                sorted(Counter(record["reflection_action_denominator"] for record in all_records).items())
            ),
            "vector_square_histogram": dict(
                sorted(Counter(record["vector_square"] for record in all_records).items())
            ),
            "vector_divisibility_histogram": dict(
                sorted(Counter(record["vector_divisibility"] for record in all_records).items())
            ),
            "fixed_lattice_determinant_histogram": dict(
                sorted(Counter(record["fixed_lattice_determinant_absolute"] for record in all_records).items())
            ),
            "certified_non_Coxeter_pair_count_histogram": {
                str(key): value for key, value in sorted(bad_pair_histogram.items())
            },
            "minimum_certified_non_Coxeter_pair_count": min(bad_pair_histogram),
            "all_checks_pass": all(
                all(record["checks"].values()) for record in all_records
            ),
        },
        "records": all_records,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["summary"], sort_keys=True))
    print(OUTPUT)


if __name__ == "__main__":
    main()
