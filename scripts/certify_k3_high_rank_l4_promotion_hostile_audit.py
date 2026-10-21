#!/usr/bin/env python3
"""Hostile logical audit of the proposed high-rank K3Root L4 promotion.

The upstream package proves a useful marked K3 carrier/return theorem and a
correct 19-row finite-volume classification.  This audit asks the narrower
question imposed by the refined inverse mandate: does the *classification*
come from AG information independent of the reflective-lattice input, and
does it therefore create a new theorem back on the polytope class?

It also freezes a concrete warning against the shorthand
"carrier nonempty iff finite volume".  The even lattice

    U + E7(-1) + A2(-1)

has an explicit primitive K3-lattice embedding and hence full-nef K3
carriers, but it is outside the 19 finite-volume rows.  The correct statement
uses the *finite-Aut/polyhedral-effective-cone/Cox-finite* carrier locus.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
UPSTREAM_PATH = (
    ROOT / "data/enumerations/k3_high_rank_full_root_classification.json"
)
UPSTREAM_REPORT = (
    ROOT / "results/k3_high_rank_full_root_classification_transfer.md"
)
OUTPUT_PATH = (
    ROOT / "data/enumerations/k3_high_rank_l4_promotion_hostile_audit.json"
)

ROULLEAU_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
AHL_TEXT = Path("/private/tmp/ahl_k3.txt")

PINNED_EXTERNAL_HASHES = {
    "roulleau_tex": (
        "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
    ),
    "ahl_text": (
        "3278ddaf1215222b60206d012575c6b6fb08486128a344168453314af439a91c"
    ),
}

# One-based inclusive source windows.  Their hashes pin the exact dependency
# witnesses without copying the excerpts into the machine artifact.
PINNED_EXCERPT_HASHES = {
    "roulleau_introduction_lines_197_218": (
        "8bc2148efbe02d72bd3b153045ebf5a9bd49daf859e41cd288601e2d6c3cf476"
    ),
    "ahl_theorems_2_11_2_12_lines_500_530": (
        "36e6a0bd37abde088da6aaf62f426d5f439fd2a9679c4823cd11b5736e686a5e"
    ),
}


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _window_hash(path: Path, first: int, last: int) -> str:
    lines = path.read_text(errors="replace").splitlines()
    payload = "\n".join(lines[first - 1 : last]).encode()
    return _sha256_bytes(payload)


def _bareiss_det(matrix: list[list[int]]) -> int:
    """Exact integer determinant by fraction-free Bareiss elimination."""

    a = [row[:] for row in matrix]
    n = len(a)
    if n == 0:
        return 1
    sign = 1
    previous = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k + 1, n) if a[i][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign *= -1
        pivot_value = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot_value - a[i][k] * a[k][j]
                if numerator % previous:
                    raise ArithmeticError("Bareiss division ceased to be exact")
                a[i][j] = numerator // previous
        previous = pivot_value
        for i in range(k + 1, n):
            a[i][k] = 0
    return sign * a[-1][-1]


def _cartan_a2() -> list[list[int]]:
    return [[2, -1], [-1, 2]]


def _cartan_e7() -> list[list[int]]:
    # Chain 0--1--2--3--4--5 with node 6 attached to node 2.
    matrix = [[2 if i == j else 0 for j in range(7)] for i in range(7)]
    for i in range(5):
        matrix[i][i + 1] = matrix[i + 1][i] = -1
    matrix[2][6] = matrix[6][2] = -1
    return matrix


def _cartan_e8() -> list[list[int]]:
    # Chain 0--...--6 with node 7 attached to node 2.
    matrix = [[2 if i == j else 0 for j in range(8)] for i in range(8)]
    for i in range(6):
        matrix[i][i + 1] = matrix[i + 1][i] = -1
    matrix[2][7] = matrix[7][2] = -1
    return matrix


def _principal(matrix: list[list[int]], indices: list[int]) -> list[list[int]]:
    return [[matrix[i][j] for j in indices] for i in indices]


def _block_diag(*blocks: list[list[int]]) -> list[list[int]]:
    size = sum(len(block) for block in blocks)
    result = [[0 for _ in range(size)] for _ in range(size)]
    offset = 0
    for block in blocks:
        for i, row in enumerate(block):
            for j, value in enumerate(row):
                result[offset + i][offset + j] = value
        offset += len(block)
    return result


def _negative(matrix: list[list[int]]) -> list[list[int]]:
    return [[-entry for entry in row] for row in matrix]


def check_external_sources() -> dict[str, Any]:
    """Recheck the exact source passages used to diagnose circularity."""

    for path in (ROULLEAU_TEX, AHL_TEXT):
        if not path.exists():
            raise FileNotFoundError(path)
    if _sha256(ROULLEAU_TEX) != PINNED_EXTERNAL_HASHES["roulleau_tex"]:
        raise AssertionError("Roulleau source hash changed")
    if _sha256(AHL_TEXT) != PINNED_EXTERNAL_HASHES["ahl_text"]:
        raise AssertionError("AHL source hash changed")

    roulleau_window = _window_hash(ROULLEAU_TEX, 197, 218)
    ahl_window = _window_hash(AHL_TEXT, 500, 530)
    assert (
        roulleau_window
        == PINNED_EXCERPT_HASHES["roulleau_introduction_lines_197_218"]
    )
    assert (
        ahl_window
        == PINNED_EXCERPT_HASHES["ahl_theorems_2_11_2_12_lines_500_530"]
    )

    roulleau_text = ROULLEAU_TEX.read_text(errors="replace")
    ahl_text = AHL_TEXT.read_text(errors="replace")
    assert "The result of their classification is a list of $118$" in roulleau_text
    assert "In the present paper, we exhibit for each of these surfaces" in roulleau_text
    assert "The hyperbolic lattices satisfying (ii) have been classi" in ahl_text
    assert "together with Theorem 2.11, give the following" in ahl_text
    assert "Theorem 2.12" in ahl_text

    return {
        "checked": True,
        "all_checks_pass": True,
        "full_source_sha256": PINNED_EXTERNAL_HASHES,
        "excerpt_sha256": PINNED_EXCERPT_HASHES,
        "logical_witness": {
            "roulleau": (
                "the atlas starts from the Nikulin--Vinberg Picard-lattice "
                "classification and then supplies geometric constructions"
            ),
            "ahl": (
                "Theorem 2.12 is explicitly obtained by combining Theorem "
                "2.11 with the prior classification of O(N)/W(N)-finite "
                "hyperbolic lattices"
            ),
        },
    }


def _counterexample_lattice_audit(upstream: dict[str, Any]) -> dict[str, Any]:
    """Audit U+E7(-1)+A2(-1) as carrier-positive but not FV."""

    u = [[0, 1], [1, 0]]
    e7 = _cartan_e7()
    e8 = _cartan_e8()
    a2 = _cartan_a2()
    gram = _block_diag(u, _negative(e7), _negative(a2))

    # Positive definiteness of E7 is exact by Sylvester's criterion.
    e7_leading_minors = [
        _bareiss_det([row[:k] for row in e7[:k]]) for k in range(1, 8)
    ]
    assert e7_leading_minors == [2, 3, 4, 5, 6, 7, 2]
    assert len(gram) == 11
    assert all(len(row) == 11 for row in gram)
    assert all(gram[i][i] % 2 == 0 for i in range(11))
    assert _bareiss_det(gram) == 6

    # A literal primitive embedding into 3U+2E8(-1): use one U; in the first
    # E8 simple-root Z-basis use nodes [0,...,5,7], an E7 subdiagram; in the
    # second use the adjacent simple roots 0 and 1.  Coordinate subsets of
    # a Z-basis are primitive.
    e7_inside_e8 = _principal(e8, [0, 1, 2, 3, 4, 5, 7])
    a2_inside_e8 = _principal(e8, [0, 1])
    assert e7_inside_e8 == e7
    assert a2_inside_e8 == a2
    assert _bareiss_det(e8) == 1

    label = "U+E7(-1)+A2(-1)"
    listed_labels = {case["lattice"] for case in upstream["cases"]}
    assert label not in listed_labels

    return {
        "lattice": label,
        "rank": 11,
        "hyperbolic_dimension": 10,
        "signature": [1, 10],
        "even": True,
        "determinant": 6,
        "e7_leading_principal_minors": e7_leading_minors,
        "primitive_k3_embedding": {
            "ambient": "3U+2E8(-1)",
            "witness": (
                "one U summand; E7 on first-E8 simple-root nodes "
                "[0,1,2,3,4,5,7]; an A2 on second-E8 nodes [0,1]"
            ),
            "coordinate_sublattices_primitive": True,
        },
        "exact_picard_k3_carrier_nonempty": True,
        "full_root_chamber_returns_as_nef": True,
        "in_nineteen_finite_volume_rows": False,
        "finite_volume": False,
        "finite_volume_reason": (
            "the accepted Nikulin--Vinberg exhaustion says every finite-volume "
            "rank-11 member is one of its four rows; this lattice is not one"
        ),
        "logical_use": (
            "counterexample to raw carrier-nonemptiness iff finite volume; "
            "not a counterexample to the decorated finite-Aut carrier iff"
        ),
    }


def build_artifact(source_check: dict[str, Any] | None = None) -> dict[str, Any]:
    upstream = json.loads(UPSTREAM_PATH.read_text())
    upstream_report = UPSTREAM_REPORT.read_text()
    assert upstream["all_checks_pass"] is True
    assert upstream["surface_free_class"]["name"] == "K3Root[11,20]"
    assert upstream["classification"][
        "rank_11_through_19_finite_volume_type_count"
    ] == 19
    assert upstream["classification"]["rank_20_finite_volume_type_count"] == 0
    assert upstream["classification"]["worldwide_novelty"] is False
    assert upstream["transfer_theorem"]["U"]["status"] == "closed"
    assert upstream["transfer_theorem"]["R"]["status"].startswith("closed")
    assert upstream["transfer_theorem"]["C"]["status"] == "closed"
    assert upstream["transfer_theorem"]["V"]["status"].startswith("closed")
    assert upstream["transfer_theorem"]["V"]["noncompact_all_19"] is True
    assert upstream["transfer_theorem"]["strict_m3c"]["status"] == "open"
    assert "Nikulin--Vinberg" in upstream_report
    assert "classification output" in upstream_report

    source_hashes = upstream["source_check"]["sha256"]
    assert source_hashes["roulleau_tex"] == PINNED_EXTERNAL_HASHES["roulleau_tex"]
    assert source_hashes["ahl_text"] == PINNED_EXTERNAL_HASHES["ahl_text"]
    assert upstream["source_check"]["ahl_witness"]["result"] == "Theorem 2.12"

    counterexample = _counterexample_lattice_audit(upstream)
    cases = upstream["cases"]
    assert len(cases) == 19
    assert all(case["contains_U"] for case in cases)
    assert all(
        case["literal_vinberg_status"]["decision"]
        == "finite_volume_noncompact"
        for case in cases
    )
    assert sum(
        case["literal_vinberg_status"]["ideal_vertex_count"] for case in cases
    ) == 1156

    checks = {
        "non_singleton_surface_free_arithmetic_class": True,
        "essential_surjectivity_to_declared_full_nef_k3_species": True,
        "oriented_embedding_and_period_for_all_19": True,
        "same_marked_lattice_and_complete_nef_chamber_return": True,
        "no_hidden_minus_two_roots_by_complete_arrangement": True,
        "literal_finite_volume_for_all_positive_rows": True,
        "all_positive_rows_noncompact": True,
        "raw_carrier_nonemptiness_does_not_detect_finite_volume": True,
        "nineteen_zero_exhaustion_is_independent_ag_native": False,
        "new_polytope_theorem_from_ag_pushback": False,
        "strict_l4": False,
    }

    return {
        "schema": "coxeter.k3_high_rank_l4_promotion_hostile_audit.v1",
        "artifact_id": "k3_high_rank_l4_promotion_hostile_audit",
        "all_checks_pass": True,
        "upstream": {
            "path": str(UPSTREAM_PATH.relative_to(ROOT)),
            "sha256": _sha256(UPSTREAM_PATH),
            "report_path": str(UPSTREAM_REPORT.relative_to(ROOT)),
            "report_sha256": _sha256(UPSTREAM_REPORT),
            "mathematical_classification_valid": True,
        },
        "strict_score": {
            "C": (
                "PASS: non-singleton arithmetic lattice/chamber class, though "
                "primitive K3-lattice embeddability is built into membership"
            ),
            "Car": (
                "PASS for the declared full-nef K3 carrier species: primitive "
                "embedding plus periods realizes every ambient member; the "
                "nineteen finite-volume rows have one oriented embedding orbit"
            ),
            "Return": (
                "PASS for integrally marked data: K3 Weyl transitivity returns "
                "the same complete all-minus-two chamber with no hidden roots"
            ),
            "AG": (
                "FAIL_INDEPENDENT_EXHAUSTION: period/Torelli, finite-Aut and "
                "Cox-ring equivalences are AG, but the 19/0 list is imported "
                "from the Nikulin--Vinberg O(N)/W_2(N)-finite lattice list"
            ),
            "Pushback": (
                "FAIL_NEW_THEOREM: the finite list and no-H19 conclusion restate "
                "the same reflective-lattice classification; the literal "
                "Vinberg replays verify positives but cannot exhaust rank 20"
            ),
        },
        "verdict": {
            "strict_l4": False,
            "promotion": "rejected",
            "best_honest_status": (
                "class-wide marked K3 carrier/return and AG volume dictionary, "
                "plus a replayed published high-rank classification; not L4"
            ),
            "not_a_refutation_of_upstream_mathematics": True,
            "worldwide_new_abstract_polytopes": 0,
        },
        "classification_dependency": {
            "ahl_theorem_2_11_role": (
                "identifies polyhedral Eff, finite O(N)/W_2(N), and finite Aut"
            ),
            "ahl_theorem_2_12_role": (
                "imports the previously classified hyperbolic lattices and gives "
                "the rank counts, including zero at rank 20"
            ),
            "roulleau_role": (
                "starts from the Nikulin--Vinberg 118-lattice list and supplies "
                "geometric models, curve configurations, and moduli information"
            ),
            "literal_vinberg_role": (
                "independently certifies finite volume/noncompactness of the 19 "
                "listed chambers, but does not prove the list exhaustive"
            ),
            "circularity_for_strict_l4": (
                "the property being classified on the surface side is equivalent "
                "to the reflective-lattice property, and the actual exhaustion "
                "is sourced from the reflective-lattice classification"
            ),
        },
        "correct_finite_volume_iff": {
            "false_shorthand": (
                "finite volume iff the unqualified exact-Picard K3 carrier locus "
                "is nonempty"
            ),
            "correct_statement": (
                "finite volume iff some/every exact-Picard full-nef K3 carrier "
                "has finite Aut iff Eff is rational polyhedral iff Cox is finitely "
                "generated"
            ),
            "counterexample": counterexample,
        },
        "compact_noncompact_separation": {
            "compact_finite_volume_locus_rank_11_20": "empty",
            "compact_reason": (
                "Meyer makes every rank-at-least-five indefinite rational lattice "
                "isotropic; a finite-volume complete K3 root chamber therefore "
                "has a nef isotropic cusp"
            ),
            "compact_reason_is_ag_native": False,
            "finite_volume_noncompact_locus": {
                "type_count": 19,
                "hyperbolic_dimensions": [10, 18],
                "total_ideal_vertices_in_literal_replays": 1156,
                "no_H19_type": True,
            },
        },
        "missing_l4_gate": {
            "required": (
                "an AG-native classification/boundedness theorem for the marked "
                "K3 carrier loci that does not invoke the Nikulin--Vinberg "
                "2-reflective lattice exhaustion, followed by a polytope "
                "conclusion not already contained in that exhaustion"
            ),
            "acceptable_forms": [
                "an independent geometric classification via surface models, genus-one fibrations, curve configurations, or Cox algebras that derives the finite list",
                "a broader raw-P class on which AG effectivity/moduli removes or constructs cases not fixed by the input Coxeter classification",
                "a genuinely new facet/dimension restriction or positive terminal realization obtained only after the AG stage",
            ],
            "not_sufficient": [
                "period existence",
                "Weyl re-marking",
                "finite Aut iff finite reflection index",
                "literal Vinberg verification of the already listed positives",
                "Roulleau geometric models indexed by the inherited 118-lattice list",
            ],
        },
        "source_check": source_check
        or {
            "checked": False,
            "reason": (
                "upstream pinned hashes and witness metadata used; run with "
                "--check-sources to re-read the external source copies"
            ),
            "full_source_sha256": PINNED_EXTERNAL_HASHES,
            "excerpt_sha256": PINNED_EXCERPT_HASHES,
        },
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-sources", action="store_true")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()

    source_check = check_external_sources() if args.check_sources else None
    artifact = build_artifact(source_check)
    if args.write:
        OUTPUT_PATH.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    if args.check_archive:
        frozen = json.loads(OUTPUT_PATH.read_text())
        artifact = build_artifact(frozen["source_check"])
        if artifact != frozen:
            raise ArithmeticError("frozen hostile-audit archive changed")
    print(json.dumps(artifact, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
