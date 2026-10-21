#!/usr/bin/env python3
"""Certify the fixed Enriques K3-cover lattice and a universal lift obstruction.

The mathematical theorem recorded in
``results/enriques_variable_lift_package_theorem.md`` has both discrete and
period-theoretic hypotheses.  This script certifies the completely finite
part which is easy to get subtly wrong:

* the invariant/anti-invariant decomposition of the K3 lattice and its glue;
* the exact equations for a pair of exchanged ``(-2)`` lifts of a downstairs
  Enriques root; and
* the failure of the tempting canonical simultaneous lift for every one of
  the seven Kondō--Martin full-nef wall systems.

The last item is an obstruction, not a nonexistence theorem.  It proves that
the downstairs wall Gram does not by itself determine a viable K3-cover
period: nontrivial two-sheet intersection data are indispensable.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

from coxeter_polytopes.exact.matrix import inertia


ROOT = Path(__file__).resolve().parents[1]
RAW_DATA = (
    ROOT
    / "data"
    / "external"
    / "shimada_l10l26note"
    / "L10L26compdata2.txt"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "enriques_variable_lift_package.json"
)

RAW_DATA_SHA256 = "218c8188e467bde29369cdeda3771d767c3f70105d8fb56def8d698b90e75a6c"

TYPE_TO_RECORD = {
    "I": "12A",
    "II": "12B",
    "III": "20B",
    "IV": "20F",
    "V": "20A",
    "VI": "20E",
    "VII": "20C",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _extract_integer_array(record_text: str, field: str) -> list[Any]:
    matches = tuple(re.finditer(rf"\b{re.escape(field)}\s*:=", record_text))
    if len(matches) != 1:
        raise ArithmeticError(f"record has {len(matches)} occurrences of {field}")
    opening = record_text.find("[", matches[0].end())
    if opening < 0:
        raise ArithmeticError(f"missing array after {field}")
    level = 0
    for position in range(opening, len(record_text)):
        character = record_text[position]
        if character == "[":
            level += 1
        elif character == "]":
            level -= 1
            if level == 0:
                value = ast.literal_eval(record_text[opening : position + 1])
                if not isinstance(value, list):
                    raise ArithmeticError(f"{field} is not a list")
                return value
    raise ArithmeticError(f"truncated array for {field}")


def parse_gap_records(path: Path = RAW_DATA) -> dict[str, dict[str, list[Any]]]:
    """Parse Shimada's records without importing either seven-type certifier."""

    chunks = re.split(r"\)\s*,\s*rec\(", path.read_text(encoding="utf-8"))
    if len(chunks) != 16:
        raise ArithmeticError(f"expected 16 GAP records, found {len(chunks)}")
    records: dict[str, dict[str, list[Any]]] = {}
    for chunk in chunks:
        names = re.findall(r'\bname\s*:=\s*"([^"]+)"', chunk)
        if len(names) != 1 or names[0] in records:
            raise ArithmeticError(f"bad record name census {names}")
        records[names[0]] = {
            field: _extract_integer_array(chunk, field)
            for field in ("GramL10", "configuration", "walls")
        }
    return records


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[i, j]) for j in range(matrix.cols)]
        for i in range(matrix.rows)
    ]


def _column(values: Iterable[int]) -> sp.Matrix:
    return sp.Matrix([int(value) for value in values])


def _gcd(values: Iterable[int]) -> int:
    result = 0
    for value in values:
        result = math.gcd(result, abs(int(value)))
    return result


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=sp.ZZ)
    return [
        abs(int(smith[i, i]))
        for i in range(min(smith.rows, smith.cols))
        if smith[i, i] != 0
    ]


def _pair_histogram(configuration: sp.Matrix) -> dict[str, int]:
    counts = Counter(
        int(configuration[i, j])
        for i in range(configuration.rows)
        for j in range(i + 1, configuration.cols)
    )
    return {str(value): counts[value] for value in sorted(counts)}


def _fixed_cover_model(gram_m: sp.Matrix) -> tuple[dict[str, Any], sp.Matrix, sp.Matrix]:
    """Return the exact cover certificate together with Lambda Gram and tau."""

    u = sp.Matrix([[0, 1], [1, 0]])
    gram_lambda = sp.diag(u, gram_m, gram_m)
    identity_2 = sp.eye(2)
    identity_10 = sp.eye(10)
    zero_2_10 = sp.zeros(2, 10)
    zero_10_2 = sp.zeros(10, 2)

    tau = sp.Matrix.vstack(
        sp.Matrix.hstack(-identity_2, zero_2_10, zero_2_10),
        sp.Matrix.hstack(zero_10_2, sp.zeros(10), identity_10),
        sp.Matrix.hstack(zero_10_2, identity_10, sp.zeros(10)),
    )
    invariant = sp.Matrix.vstack(zero_2_10, identity_10, identity_10)
    anti_invariant = sp.Matrix.vstack(
        sp.Matrix.hstack(identity_2, zero_2_10),
        sp.Matrix.hstack(zero_10_2, identity_10),
        sp.Matrix.hstack(zero_10_2, -identity_10),
    )

    gram_s = invariant.T * gram_lambda * invariant
    gram_t = anti_invariant.T * gram_lambda * anti_invariant
    combined = sp.Matrix.hstack(invariant, anti_invariant)
    cross = invariant.T * gram_lambda * anti_invariant

    if tau * tau != sp.eye(22) or tau.T * gram_lambda * tau != gram_lambda:
        raise ArithmeticError("tau is not an involutive K3-lattice isometry")
    if tau * invariant != invariant or tau * anti_invariant != -anti_invariant:
        raise ArithmeticError("wrong tau eigenspaces")
    if gram_s != 2 * gram_m or gram_t != sp.diag(u, 2 * gram_m):
        raise ArithmeticError("wrong invariant/anti-invariant Gram")
    if cross != sp.zeros(10, 12):
        raise ArithmeticError("eigenspaces are not orthogonal")

    half_sum_checks = []
    for index in range(10):
        half_sum = (invariant[:, index] + anti_invariant[:, index + 2]) / 2
        expected = sp.zeros(22, 1)
        expected[index + 2, 0] = 1
        half_sum_checks.append(half_sum == expected)
    if not all(half_sum_checks):
        raise ArithmeticError("the ten half-sum glue generators changed")

    model = {
        "ambient": {
            "description": "Lambda = U + M + M, M = L10",
            "rank": 22,
            "determinant": int(gram_lambda.det()),
            "inertia": asdict(inertia(_rows(gram_lambda))),
            "tau_formula": "tau(u,x,y)=(-u,y,x)",
            "tau_is_involutive_isometry": True,
        },
        "invariant_lattice_S": {
            "description": "S = {(0,x,x)} isometric to M(2)",
            "rank": 10,
            "inertia": asdict(inertia(_rows(gram_s))),
            "absolute_discriminant": abs(int(gram_s.det())),
            "primitive_embedding_smith_diagonal": _smith_diagonal(invariant),
        },
        "anti_invariant_lattice_T": {
            "description": "T = {(u,x,-x)} isometric to U + M(2)",
            "rank": 12,
            "inertia": asdict(inertia(_rows(gram_t))),
            "absolute_discriminant": abs(int(gram_t.det())),
            "primitive_embedding_smith_diagonal": _smith_diagonal(anti_invariant),
        },
        "glue": {
            "S_perp_T": True,
            "index_of_S_plus_T_in_Lambda": abs(int(combined.det())),
            "quotient": "(Z/2Z)^10",
            "ten_half_sum_generators_integral": all(half_sum_checks),
            "half_sum_formula": "(s_i+t_i)/2=(0,e_i,0)",
        },
    }
    return model, gram_lambda, tau


def analyze() -> dict[str, Any]:
    if sha256(RAW_DATA) != RAW_DATA_SHA256:
        raise AssertionError("Shimada raw data hash changed")
    records = parse_gap_records()
    gram_m = sp.Matrix(records["12A"]["GramL10"])
    if abs(int(gram_m.det())) != 1:
        raise ArithmeticError("invalid L10 Gram")
    if asdict(inertia(_rows(gram_m))) != {"positive": 1, "negative": 9, "zero": 0}:
        raise ArithmeticError("L10 inertia changed")
    for raw_name in TYPE_TO_RECORD.values():
        if sp.Matrix(records[raw_name]["GramL10"]) != gram_m:
            raise ArithmeticError(f"{raw_name} uses a different L10 Gram")

    cover_model, gram_lambda, tau = _fixed_cover_model(gram_m)
    type_certificates = []
    for type_name, raw_name in TYPE_TO_RECORD.items():
        walls = sp.Matrix(records[raw_name]["walls"])
        configuration = sp.Matrix(records[raw_name]["configuration"])
        if walls * gram_m * walls.T != configuration:
            raise ArithmeticError(f"{raw_name}: configuration reconstruction failed")
        if walls.rank() != 10:
            raise ArithmeticError(f"{raw_name}: walls do not span M over Q")
        if any(_gcd(walls.row(i)) != 1 for i in range(walls.rows)):
            raise ArithmeticError(f"{raw_name}: imprimitive wall")

        canonical = sp.zeros(22, walls.rows)
        for index in range(walls.rows):
            canonical[2:12, index] = _column(walls.row(index))
        conjugate = tau * canonical
        invariant_sums = canonical + conjugate
        anti_differences = canonical - conjugate
        canonical_gram = canonical.T * gram_lambda * canonical
        cross_gram = canonical.T * gram_lambda * conjugate
        delta_gram = anti_differences.T * gram_lambda * anti_differences

        if canonical_gram != configuration or cross_gram != sp.zeros(walls.rows):
            raise ArithmeticError(f"{raw_name}: bad canonical lift intersections")
        if delta_gram != 2 * configuration:
            raise ArithmeticError(f"{raw_name}: delta Gram is not twice downstairs Gram")
        if any(
            (canonical[:, i].T * gram_lambda * canonical[:, i])[0] != -2
            or (conjugate[:, i].T * gram_lambda * conjugate[:, i])[0] != -2
            or (invariant_sums[:, i].T * gram_lambda * invariant_sums[:, i])[0] != -4
            or (anti_differences[:, i].T * gram_lambda * anti_differences[:, i])[0] != -4
            or 2 * canonical[:, i] != invariant_sums[:, i] + anti_differences[:, i]
            for i in range(walls.rows)
        ):
            raise ArithmeticError(f"{raw_name}: an individual lift identity failed")

        configuration_inertia = inertia(_rows(configuration))
        delta_inertia = inertia(_rows(delta_gram))
        expected = {
            "positive": 1,
            "negative": 9,
            "zero": walls.rows - 10,
        }
        if asdict(configuration_inertia) != expected or asdict(delta_inertia) != expected:
            raise ArithmeticError(f"{raw_name}: unexpected canonical delta inertia")

        type_certificates.append(
            {
                "nikulin_kondo_type": type_name,
                "raw_shimada_record": raw_name,
                "facet_count": walls.rows,
                "wall_rank": int(walls.rank()),
                "all_downstairs_walls_primitive_minus_two_roots": True,
                "downstairs_intersection_histogram": _pair_histogram(configuration),
                "canonical_lift": {
                    "formula": "b_i=(0,r_i,0), tau(b_i)=(0,0,r_i)",
                    "all_individual_pairs_are_disjoint_exchanged_minus_two_classes": True,
                    "b_i_dot_b_j_equals_downstairs_configuration": True,
                    "b_i_dot_tau_b_j_is_zero": True,
                    "delta_gram_equals_twice_downstairs_configuration": True,
                    "delta_span_rank": int(anti_differences.rank()),
                    "delta_gram_inertia": asdict(delta_inertia),
                    "delta_span_signature": {"positive": 1, "negative": 9},
                    "orthogonal_complement_in_T_inertia": {
                        "positive": 1,
                        "negative": 1,
                        "zero": 0,
                    },
                    "positive_directions_left_for_period_plane": 1,
                    "required_positive_directions_for_K3_period_plane": 2,
                    "simultaneous_period_exists": False,
                },
            }
        )

    return {
        "schema": "enriques_variable_lift_package_v1",
        "scope": {
            "finite_certificate": (
                "fixed K3-cover lattice, lift identities, and canonical simultaneous-lift obstruction"
            ),
            "not_certified_by_this_json": [
                "existence of a noncanonical anti-invariant lift package",
                "period-map surjectivity",
                "effectivity or irreducibility of lift roots",
                "full-nef cone equality",
            ],
        },
        "raw_input": {
            "path": str(RAW_DATA.relative_to(ROOT)),
            "sha256": RAW_DATA_SHA256,
        },
        "downstairs_lattice_M": {
            "description": "L10 = U + E8(-1)",
            "rank": 10,
            "inertia": asdict(inertia(_rows(gram_m))),
            "determinant": int(gram_m.det()),
        },
        "fixed_K3_cover_model": cover_model,
        "general_lift_identities": {
            "invariant_pullback": "v_i=(0,r_i,r_i) in S, v_i^2=-4",
            "anti_invariant_choice": "delta_i in T, delta_i^2=-4",
            "glue_gate": "b_i=(v_i+delta_i)/2 lies in Lambda",
            "exchanged_pair": "tau(b_i)=(v_i-delta_i)/2",
            "consequences": "b_i^2=tau(b_i)^2=-2 and b_i.tau(b_i)=0",
            "two_wall_formula": (
                "if m=r_i.r_j, a=b_i.b_j, c=b_i.tau(b_j), then "
                "a+c=m and delta_i.delta_j=2(a-c)=2(2a-m)"
            ),
            "nonnegative_curve_distribution": {
                "m=0": [0],
                "m=1": [-2, 2],
                "m=2": [-4, 0, 4],
            },
        },
        "canonical_simultaneous_lift_obstruction": {
            "statement": (
                "If the downstairs roots span L10 over Q, putting every lift on one "
                "fixed sheet gives delta Gram 2*Gram(R), of signature (1,9); its "
                "orthogonal complement in T has signature (1,1), so it cannot contain "
                "the positive two-plane required by a K3 period."
            ),
            "all_seven_wall_systems_span_L10": True,
            "all_seven_canonical_simultaneous_lifts_obstructed": True,
            "interpretation": (
                "Individual liftability is not simultaneous liftability; viable Enriques "
                "realizations require nontrivial two-sheet intersection/glue data not "
                "contained in the downstairs Gram matrix."
            ),
        },
        "types": type_certificates,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(analyze(), indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.is_file() or args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived certificate differs: {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(hashlib.sha256(rendered.encode("utf-8")).hexdigest())


if __name__ == "__main__":
    main()
