#!/usr/bin/env python3
"""Certify the determinant obstruction to an exact-fixed pure order-four lift.

The archived Roulleau ``U+E8(-1)+A1(-1)^4`` transposition has a
rank-13 fixed lattice F of determinant 16.  The order-four Lefschetz
relations leave two possible cohomological characters for a purely
non-symplectic lift.  The exhaustive Brandhorst--Hofmann classification
contains one action of each character; their fixed determinants are 64
and 256.  Since an exact H12 fixed-section lift requires the global fixed
lattice to have rank 13, it would give a full-rank containment F <= S(g);
the index-determinant identity rules out both actions.

The optional ``--brandhorst-hofmann-zip`` argument replays the complete
79-row external order-four scan from the CC-BY K3Groups archive.  The
default certificate remains offline and verifies the pinned target Gram
matrices extracted by that scan.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import io
import json
import math
from pathlib import Path
import re
from typing import Any
import zipfile

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ as POLY_ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/roulleau_r14_ue8a1x4_all_subgroup_folds.json"
OUTPUT = ROOT / "data/enumerations/roulleau_r14_c2t_pure_order4_obstruction.json"
EXPECTED_SOURCE_SHA256 = (
    "5c39aa26fc4e2c7f2af5bb08815532ec9cce7b2222279c65f257e383f356f4ef"
)

BH_DOI = "10.5281/zenodo.10410278"
BH_VERSION = "2.1"
BH_ZIP_MEMBER = "K3Groups/results/purely_ns/order04.txt"
BH_ZIP_MD5 = "143047b9cb76f4ad4cb2aadd5959d879"
BH_ZIP_SHA256 = "c1272db2c0fe2ae0a788aed7745201fccd2a1337b6c2d85555880c434cab4d7c"
BH_ORDER4_SHA256 = "d5a61603c1539278e8411e188656c9417662336ff4282b714e939c40d38c94e8"
BH_ORDER4_ROW_COUNT = 79


BH_TARGETS: dict[tuple[int, int, int], dict[str, Any]] = {
    (13, 3, 3): {
        "id": "0.4.2.9",
        "row_one_based": 24,
        "row_sha256_with_newline": (
            "17f47b1c67d9fd724b8af3fd15b29a71571213b8fda90b953d8c21dbb06994d0"
        ),
        "multiplicity_in_complete_order4_file": 1,
        "fixed_gram": [
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
            [0, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, -2, -1, 1, -1, -1, 0, 1, 1, 0],
            [0, 0, 0, 0, -1, -2, 1, -1, -1, 0, 1, 1, 0],
            [0, 0, 0, 0, 1, 1, -2, 1, 1, 0, -1, -1, 0],
            [0, 0, 0, 0, -1, -1, 1, -2, -1, -1, 1, 1, 0],
            [0, 0, 0, 0, -1, -1, 1, -1, -2, 0, 1, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, -1, 0, -2, 1, 1, 0],
            [0, 0, 0, 0, 1, 1, -1, 1, 1, 1, -2, -1, 0],
            [0, 0, 0, 0, 1, 1, -1, 1, 1, 1, -1, -2, 0],
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ],
        "fixed_determinant": 64,
        "fixed_smith": [*[1] * 10, 4, 4, 4],
    },
    (13, 7, 1): {
        "id": "0.4.0.2",
        "row_one_based": 2,
        "row_sha256_with_newline": (
            "b4f4f40fb706258e76568063140df8824a5686f63f386d359eebee46d1951e18"
        ),
        "multiplicity_in_complete_order4_file": 1,
        "fixed_gram": [
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2],
            [0, -2, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
            [0, 1, -2, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
            [0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, -1, -1, -2, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, -2, 1, 0, 0, 0, -1, 0],
            [0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, -1, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0],
            [0, 1, -1, 0, 0, 0, 0, 0, 0, -1, -2, 0, 0],
            [0, 0, 0, 0, 0, 0, -1, 1, -1, 0, 0, -2, 0],
            [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ],
        "fixed_determinant": 256,
        "fixed_smith": [*[1] * 6, *[2] * 6, 4],
    },
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def _domain_matrix(matrix: sp.Matrix) -> DomainMatrix:
    if any(sp.Rational(value).q != 1 for value in matrix):
        raise AssertionError("an integral Smith computation received fractions")
    return DomainMatrix.from_Matrix(sp.Matrix(matrix)).convert_to(POLY_ZZ)


def saturated_row_kernel(equations: sp.Matrix) -> sp.Matrix:
    """Return a saturated Z-basis for row vectors x with x*equations=0."""

    transposed = equations.T
    denominators = [sp.Rational(value).q for value in transposed]
    denominator = int(sp.ilcm(*denominators)) if denominators else 1
    integral = sp.Matrix(denominator * transposed)
    diagonal, _, right = smith_normal_decomp(_domain_matrix(integral))
    diagonal_matrix = diagonal.to_Matrix()
    right_matrix = right.to_Matrix()
    zero_columns = [
        column
        for column in range(integral.cols)
        if all(
            diagonal_matrix[row, column] == 0
            for row in range(diagonal_matrix.rows)
        )
    ]
    basis = right_matrix[:, zero_columns].T
    if basis * equations != sp.zeros(basis.rows, equations.cols):
        raise AssertionError("Smith-kernel replay failed")
    return basis


def smith_diagonal(matrix: sp.Matrix) -> list[int]:
    normal = smith_normal_form(matrix, domain=ZZ)
    return [abs(int(normal[index, index])) for index in range(min(normal.shape))]


def same_row_lattice(left: sp.Matrix, right: sp.Matrix) -> bool:
    if left.shape != right.shape or left.rank() != left.rows:
        return False
    pivots = right.rref()[1]
    change = left[:, list(pivots)] * right[:, list(pivots)].inv()
    return (
        all(sp.Rational(value).q == 1 for value in change)
        and abs(int(change.det())) == 1
        and change * right == left
    )


def _parse_bh_row(line: str) -> dict[str, Any]:
    exact = re.sub(r"(?<![\w.])(\d+)/(\d+)", r"Fraction(\1,\2)", line)
    pairs = eval(exact, {"__builtins__": {}, "Fraction": Fraction})
    return dict(pairs)


def _rational(value: Any) -> sp.Rational:
    if isinstance(value, Fraction):
        return sp.Rational(value.numerator, value.denominator)
    return sp.Rational(value)


def _bh_target_from_row(data: dict[str, Any]) -> tuple[tuple[int, int, int], sp.Matrix]:
    ambient_action = sp.Matrix(22, 22, data["distinguished_generator"])
    fixed_rank = 22 - (ambient_action - sp.eye(22)).rank()
    anti_rank = 22 - (ambient_action + sp.eye(22)).rank()
    gaussian_rank = (22 - fixed_rank - anti_rank) // 2
    character = (fixed_rank, anti_rank, gaussian_rank)

    basis = sp.Matrix(22, 22, list(map(_rational, data["basis"])))
    ambient_gram = sp.Matrix(
        22, 22, list(map(_rational, data["inner_product_matrix"]))
    )
    gram = sp.simplify(basis * ambient_gram * basis.T)
    action = sp.simplify(basis * ambient_action * basis.inv())
    if not (
        gram.det() == -1
        and exact_inertia(gram) == (3, 19, 0)
        and all(sp.Rational(value).q == 1 for value in gram)
        and all(sp.Rational(value).q == 1 for value in action)
        and action * gram * action.T == gram
        and action**4 == sp.eye(22)
    ):
        raise AssertionError("Brandhorst--Hofmann K3 row replay failed")
    fixed_basis = saturated_row_kernel(action - sp.eye(22))
    fixed_gram = sp.simplify(fixed_basis * gram * fixed_basis.T)
    return character, fixed_gram


def replay_brandhorst_hofmann_zip(path: Path) -> dict[str, Any]:
    """Replay all 79 pure order-four rows from the published archive."""

    zip_bytes = path.read_bytes()
    if not (
        hashlib.md5(zip_bytes).hexdigest() == BH_ZIP_MD5  # noqa: S324; provenance only
        and sha256_bytes(zip_bytes) == BH_ZIP_SHA256
    ):
        raise AssertionError("the pinned Brandhorst--Hofmann ZIP changed")
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as archive:
        order4_bytes = archive.read(BH_ZIP_MEMBER)
    if sha256_bytes(order4_bytes) != BH_ORDER4_SHA256:
        raise AssertionError("the pinned pure order-four data file changed")
    lines = order4_bytes.decode("utf-8").splitlines()
    if len(lines) != BH_ORDER4_ROW_COUNT:
        raise AssertionError("the pure order-four row count changed")

    hits: dict[tuple[int, int, int], list[dict[str, Any]]] = {
        character: [] for character in BH_TARGETS
    }
    for row_number, line in enumerate(lines, 1):
        data = _parse_bh_row(line)
        ambient_action = sp.Matrix(22, 22, data["distinguished_generator"])
        fixed_rank = 22 - (ambient_action - sp.eye(22)).rank()
        anti_rank = 22 - (ambient_action + sp.eye(22)).rank()
        character = (
            fixed_rank,
            anti_rank,
            (22 - fixed_rank - anti_rank) // 2,
        )
        if character not in hits:
            continue
        replayed_character, fixed_gram = _bh_target_from_row(data)
        if replayed_character != character:
            raise AssertionError("character replay changed")
        hits[character].append(
            {
                "id": data["id"],
                "row_one_based": row_number,
                "row_sha256_with_newline": sha256_bytes((line + "\n").encode()),
                "fixed_gram": fixed_gram,
            }
        )

    for character, expected in BH_TARGETS.items():
        observed = hits[character]
        if len(observed) != expected["multiplicity_in_complete_order4_file"]:
            raise AssertionError("the classified target multiplicity changed")
        row = observed[0]
        if not (
            row["id"] == expected["id"]
            and row["row_one_based"] == expected["row_one_based"]
            and row["row_sha256_with_newline"]
            == expected["row_sha256_with_newline"]
            and row["fixed_gram"] == sp.Matrix(expected["fixed_gram"])
        ):
            raise AssertionError("the classified target row changed")
    return {
        "zip_sha256": BH_ZIP_SHA256,
        "order4_file_sha256": BH_ORDER4_SHA256,
        "row_count": len(lines),
        "target_multiplicities": {
            str(character): len(rows) for character, rows in hits.items()
        },
    }


def _verify_pinned_classification_targets() -> None:
    for character, target in BH_TARGETS.items():
        gram = sp.Matrix(target["fixed_gram"])
        if not (
            gram.shape == (13, 13)
            and exact_inertia(gram) == (1, 12, 0)
            and abs(int(gram.det())) == target["fixed_determinant"]
            and smith_diagonal(gram) == target["fixed_smith"]
            and target["multiplicity_in_complete_order4_file"] == 1
            and character[0] == 13
        ):
            raise AssertionError("a pinned classification target is inconsistent")


def _index_square_gate(subdeterminant: int, superdeterminant: int) -> dict[str, Any]:
    ratio = sp.Rational(subdeterminant, superdeterminant)
    integer_square = (
        ratio.q == 1 and math.isqrt(int(ratio.p)) ** 2 == int(ratio.p)
    )
    return {
        "det_F_over_det_S": str(ratio),
        "required_identity": "det(F)=[S:F]^2*det(S)",
        "is_nonnegative_integer_square": integer_square,
        "containment_possible": integer_square,
    }


def build_certificate(external_zip: Path | None = None) -> dict[str, Any]:
    source_hash = sha256(SOURCE)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise AssertionError("the pinned R14/27 source certificate changed")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = next(
        row for row in source["cases"] if row["id"] == "R14-UE8A1x4-C2T-H12-23"
    )
    vinberg = case["geometry"]["vinberg"]
    if not (
        source["all_checks_pass"]
        and case["all_checks_pass"]
        and all(case["checks"].values())
        and vinberg["classification"]
        == "finite_volume_noncompact_Coxeter_polytope"
        and vinberg["dimension"] == 12
        and vinberg["facet_count"] == 23
        and vinberg["finite_vertex_count"] == 320
        and vinberg["ideal_vertex_count"] == 17
        and vinberg["elliptic_edge_count"] == 2308
        and vinberg["bad_edge_count"] == 0
        and vinberg["bad_parabolic_neighborhood_count"] == 0
    ):
        raise AssertionError("the upstream H12/23 Vinberg theorem changed")

    n_gram = sp.Matrix(
        source["source_reconstruction"]["authoritative_lane"]["basis_gram"]
    )
    wall_columns = sp.Matrix(
        source["source_reconstruction"]["authoritative_lane"][
            "wall_coordinate_columns"
        ]
    )
    permutation = [
        value - 1 for value in case["subgroup"]["elements_one_based"][1]
    ]
    column_action = sp.Matrix.hstack(
        *(wall_columns[:, permutation[index]] for index in range(14))
    )
    n_action = column_action.T
    archived_fixed_basis = sp.Matrix(
        case["fixed_lattice"]["saturated_basis_columns_in_NS_basis"]
    ).T
    fixed_gram = sp.simplify(archived_fixed_basis * n_gram * archived_fixed_basis.T)
    computed_fixed_basis = saturated_row_kernel(n_action - sp.eye(14))
    if not (
        n_gram.det() == -16
        and exact_inertia(n_gram) == (1, 13, 0)
        and n_action**2 == sp.eye(14)
        and n_action * n_gram * n_action.T == n_gram
        and sp.factor(n_action.charpoly().as_expr())
        == (sp.Symbol("lambda") - 1) ** 13 * (sp.Symbol("lambda") + 1)
        and archived_fixed_basis.rows == 13
        and same_row_lattice(archived_fixed_basis, computed_fixed_basis)
        and fixed_gram.det() == 16
        and smith_diagonal(fixed_gram) == [*[1] * 10, 2, 2, 4]
    ):
        raise AssertionError("the literal archived fixed lattice changed")

    # Artebani--Sarti's pure order-four Lefschetz relations are
    # r+l+2m=22 and 4*alpha=r-l-2.  Here r=rank(F)=13.
    lefschetz_possibilities = []
    for anti_rank in range(23):
        remaining = 22 - 13 - anti_rank
        numerator = 13 - anti_rank - 2
        if remaining < 0 or remaining % 2 or numerator < 0 or numerator % 4:
            continue
        gaussian_rank = remaining // 2
        alpha = numerator // 4
        if gaussian_rank >= 1:
            lefschetz_possibilities.append((13, anti_rank, gaussian_rank, alpha))
    if lefschetz_possibilities != [(13, 3, 3, 2), (13, 7, 1, 1)]:
        raise AssertionError("the pure order-four Lefschetz gate changed")

    _verify_pinned_classification_targets()
    if external_zip is not None:
        replay_brandhorst_hofmann_zip(external_zip)

    classification_rows = []
    for fixed_rank, anti_rank, gaussian_rank, alpha in lefschetz_possibilities:
        character = (fixed_rank, anti_rank, gaussian_rank)
        target = BH_TARGETS[character]
        gate = _index_square_gate(16, target["fixed_determinant"])
        if gate["containment_possible"]:
            raise AssertionError("a classified fixed lattice passed the index gate")
        classification_rows.append(
            {
                "character_r_l_m": list(character),
                "lefschetz_alpha": alpha,
                "classification_id": target["id"],
                "multiplicity_in_complete_order4_file": target[
                    "multiplicity_in_complete_order4_file"
                ],
                "row_one_based": target["row_one_based"],
                "row_sha256_with_newline": target["row_sha256_with_newline"],
                "fixed_determinant": target["fixed_determinant"],
                "fixed_smith": target["fixed_smith"],
                "full_rank_containment_gate": gate,
            }
        )

    certificate: dict[str, Any] = {
        "schema": "coxeter.roulleau_r14_c2t.pure_order4_obstruction.v1",
        "all_checks_pass": True,
        "source_file_sha256": source_hash,
        "case_id": "R14-UE8A1x4-C2T-H12-23",
        "theorem": (
            "The archived R14 transposition defining the finite-volume H12/23 "
            "section has no purely non-symplectic order-four K3 automorphism "
            "lift whose global fixed lattice has rank 13"
        ),
        "archived_divisor_action": {
            "lattice": "U+E8(-1)+A1(-1)^4",
            "rank_signature_determinant": [14, [1, 13], -16],
            "action_characteristic_polynomial": "(lambda-1)^13*(lambda+1)",
            "fixed_rank_determinant": [13, 16],
            "fixed_smith": [*[1] * 10, 2, 2, 4],
        },
        "pure_order4_lefschetz_gate": {
            "relations": ["r+l+2m=22", "4*alpha=r-l-2"],
            "fixed_rank_r": 13,
            "surviving_r_l_m_alpha": [
                list(possibility) for possibility in lefschetz_possibilities
            ],
        },
        "brandhorst_hofmann_classification": {
            "citation": "Brandhorst--Hofmann, Finite subgroups of automorphisms of K3 surfaces",
            "doi": BH_DOI,
            "archive_version": BH_VERSION,
            "license": "CC-BY-4.0",
            "zip_md5": BH_ZIP_MD5,
            "zip_sha256": BH_ZIP_SHA256,
            "member": BH_ZIP_MEMBER,
            "member_sha256": BH_ORDER4_SHA256,
            "complete_row_count": BH_ORDER4_ROW_COUNT,
            "full_external_scan_replay_available": True,
            "full_external_scan_command": (
                "python3 scripts/certify_roulleau_r14_c2t_pure_order4_obstruction.py "
                "--brandhorst-hofmann-zip /path/to/K3Groups.zip --check-archive"
            ),
            "surviving_rows": classification_rows,
        },
        "determinant_obstruction": {
            "logic": (
                "An exact-fixed lift fixes every vector of F and has global "
                "fixed rank 13. Thus F is a finite-index sublattice of S(g), "
                "so det(F)=[S(g):F]^2 det(S(g))."
            ),
            "classified_fixed_determinants": [64, 256],
            "archived_fixed_determinant": 16,
            "both_index_square_tests_fail": True,
        },
        "hyperbolic_polytope_handoff": {
            "dimension": 12,
            "facets": 23,
            "classification": "finite_volume_noncompact",
            "elliptic_edges": 2308,
            "finite_vertices": 320,
            "cusps": 17,
            "bad_edges": 0,
            "bad_parabolic_neighborhoods": 0,
            "section_to_odd_unimodular_full_group_index": 8320,
        },
        "scope": {
            "exact_rank13_fixed_pure_order4_lift_ruled_out": True,
            "pure_order4_lift_with_enlarged_global_fixed_lattice_ruled_out": False,
            "symplectic_square_order4_lift_ruled_out": False,
            "higher_order_lift_ruled_out": False,
            "abstract_hyperbolic_polytope_invalidated": False,
            "special_period_search_needed": False,
            "route_status": "retired_by_complete_classification_and_determinant_gate",
        },
    }
    canonical = json.dumps(certificate, sort_keys=True, separators=(",", ":"))
    certificate["payload_sha256"] = hashlib.sha256(canonical.encode()).hexdigest()
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--brandhorst-hofmann-zip", type=Path)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    certificate = build_certificate(args.brandhorst_hofmann_zip)
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            raise AssertionError("the archived pure order-four obstruction changed")
        print(f"archive OK: {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(args.output)
    print(certificate["payload_sha256"])


if __name__ == "__main__":
    main()
