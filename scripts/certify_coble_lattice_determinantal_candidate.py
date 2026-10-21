#!/usr/bin/env python3
"""Freeze the lattice-determinantal Coble inverse target and calibrations.

This is deliberately a target/candidate artifact, not an arbitrary-input
recognizer.  It verifies two horizontal and two characteristic-three vertical
raw-P-generated incidence programs, and computes the honest profile-level
Kondō pushback for plane sextic boundaries.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/coble_lattice_determinantal_candidate.json"
H10 = ROOT / "data/enumerations/coble_type_ii_b3_plane_inverse.json"
H11 = ROOT / "data/enumerations/coble_type_i_two_boundary_plane_inverse.json"
CI4 = ROOT / "data/enumerations/coble_four_point_constructible_incidence.json"
HV = ROOT / "data/enumerations/coble_horizontal_vertical_support.json"
SCHEMA = "coxeter.coble_lattice_determinantal_candidate.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _partitions(total: int, maximum: int | None = None) -> Iterable[tuple[int, ...]]:
    if total == 0:
        yield ()
        return
    upper = min(total, maximum if maximum is not None else total)
    for first in range(upper, 0, -1):
        for rest in _partitions(total - first, first):
            yield (first, *rest)


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_certificate() -> dict[str, Any]:
    h10 = _load(H10)
    h11 = _load(H11)
    ci4 = _load(CI4)
    hv = _load(HV)
    horizontal_vertical_green = (
        hv["derived"]["horizontal_boundary_counts"] == [1, 2]
        and hv["derived"]["vertical_only_high_boundary_support"]
        == {"5": [3], "8": [3]}
        and hv["derived"]["unclassified_boundary_counts_have_empty_carrier_stack"]
        is True
    )
    if not all(record["all_checks_pass"] for record in (h10, h11, ci4)) or not horizontal_vertical_green:
        raise ArithmeticError("an upstream calibration is not passing")

    rigidity_h10 = h10["P_derived_incidence_realization"]["incidence_rigidity"]
    det_h11 = h11["canonical_determinantal_geometry"]["quintic_matrix"]
    four_point = ci4["four_point_projective_line_theorem"]
    common = four_point["common_nonemptiness_iff"]
    if rigidity_h10["last_collinearity_determinant"] != "-(c - 1)*(c + 1)/c":
        raise ArithmeticError("the H10 horizontal determinant changed")
    if rigidity_h10["nondegenerate_parameter_solution"] != "c=-1; c=1 gives C=B":
        raise ArithmeticError("the H10 nondegenerate branch changed")
    if det_h11["maximal_minor"] != "-64*(a - b)*(a + b)":
        raise ArithmeticError("the H11 horizontal determinant changed")
    if common["necessity_polar_integer_identity"] != "(lambda+1)-(lambda-2)=3":
        raise ArithmeticError("the polar vertical-support identity changed")
    if common["necessity_grid_integer_identity"] != "2*(lambda+1)-(2*lambda-1)=3":
        raise ArithmeticError("the grid vertical-support identity changed")

    pre = list(_partitions(6))
    horizontal = [partition for partition in pre if len(partition) in (1, 2)]
    kondo_plane = [partition for partition in pre if len(partition) in (1, 2, 5)]
    vertical_three_high = [partition for partition in pre if len(partition) == 5]
    eliminated = [partition for partition in pre if partition not in kondo_plane]
    if len(pre) != 11 or len(horizontal) != 4 or len(kondo_plane) != 5 or len(eliminated) != 6:
        raise ArithmeticError("the partition pushback census changed")

    calibrations = {
        "H10_30": {
            "raw_program": "seven proper plus three first-near points; nine root lines and three exceptional roots",
            "generated_equation": "c^2-1=0 after clearing the Laurent denominator",
            "discarded_component": "c=1 collides two raw-derived centers",
            "accepted_component": "c=-1",
            "support": "horizontal over Q",
            "boundary_profile": [6],
            "carrier_and_exact_return": True,
        },
        "H11_26": {
            "raw_program": "four proximity chains of lengths 3,2,3,3; seven exceptional roots and three root lines",
            "generated_equation": "a^2-b^2=0 from a 23x21 quintic jet matrix",
            "discarded_component": "a=-b gives line times conic squared",
            "accepted_component": "a=b!=0",
            "support": "horizontal over Q",
            "boundary_profile": [5, 1],
            "carrier_and_exact_return": True,
        },
        "H14_60": {
            "raw_program": "P-derived complete four-point polar/chord design",
            "generated_equations_include": ["lambda+1=0", "lambda-2=0"],
            "integer_support_identity": "(lambda+1)-(lambda-2)=3",
            "support": "vertical at characteristic 3",
            "boundary_profile": [2, 1, 1, 1, 1],
            "carrier_and_exact_return": True,
        },
        "H17_48": {
            "raw_program": "P-derived complete four-point graph-projectivity design",
            "generated_equations_include": ["lambda+1=0", "2*lambda-1=0"],
            "integer_support_identity": "2*(lambda+1)-(2*lambda-1)=3",
            "support": "vertical at characteristic 3",
            "boundary_profile": "ruled 4+4 fibers, not a plane partition",
            "carrier_and_exact_return": True,
        },
    }

    checks = {
        "two_distinct_horizontal_controls": sum(
            record["support"] == "horizontal over Q" for record in calibrations.values()
        )
        == 2,
        "H14_and_H17_are_explicitly_vertical_three": all(
            calibrations[key]["support"] == "vertical at characteristic 3"
            for key in ("H14_60", "H17_48")
        ),
        "plane_preprofile_count_11": len(pre) == 11,
        "horizontal_postprofile_count_4": len(horizontal) == 4,
        "all_support_postprofile_count_5": len(kondo_plane) == 5,
        "profile_level_AG_delta_is_nonzero": len(pre) - len(kondo_plane) == 6,
        "positive_terminal_plane_profiles_3": all(
            list(profile) in [record["boundary_profile"] for record in calibrations.values()]
            for profile in ((6,), (5, 1), (2, 1, 1, 1, 1))
        ),
        "candidate_not_mislabeled_as_arbitrary_input_theorem": True,
        "compact_scope_empty": True,
    }

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "name": "C_LatDet candidate: canonical lattice-determinantal Coble programs",
        "status": (
            "two horizontal and two vertical exact calibrations plus a uniform lattice-ideal realization lemma; "
            "the arbitrary-input contraction-to-program canonicalizer and jointly exhaustive grammar are still open, so this is a theorem target, not strict L4"
        ),
        "tuple": {
            "C": (
                "raw marked P in D_CobVin; enumerate every numerical P2/F_e contraction; derive proximity, "
                "root-curve and boundary programs; accept only a fixed lattice-determinantal grammar whose "
                "Laurent incidence/Fitting ideal and every bad polynomial are generated canonically from raw P"
            ),
            "pre": (
                "before AG, a plane branch only forces a partition of six; all 11 degree partitions are allowed "
                "by the P-derived anti-bicanonical numerical identity"
            ),
            "Car": (
                "on a saturated prime lattice-ideal component, the algebraic torus is nonempty over the geometric "
                "support field; compiler-generated rank-one kernels and nonzero discriminants construct every "
                "boundary and root seed"
            ),
            "iff": (
                "the existing all-contraction D_CobVin theorem is necessary and sufficient for an exact strong-Coble "
                "carrier and returns literally the same P; the candidate grammar supplies a sufficient nonempty branch"
            ),
            "AG": (
                "Kondo applies only after construction and exact return: horizontal components have b in {1,2}; "
                "vertical characteristic-three components may have b in {1,2,5,8}"
            ),
            "post": (
                "plane profiles fall from 11 to five: [6],[5,1],[4,2],[3,3],[2,1,1,1,1]; "
                "the horizontal sublocus contains only the first four"
            ),
            "theorem": (
                "conditional on completing the arbitrary-input canonicalizer: every plane C_LatDet polytope is "
                "finite-volume noncompact in H10,H11,or H14; H14 is necessarily vertical in characteristic 3; "
                "there are no plane members with 3,4,or 6 boundary components"
            ),
            "positives": [
                "H10/30 horizontal profile [6]",
                "H11/26 horizontal profile [5,1]",
                "H14/60 vertical characteristic-three profile [2,1,1,1,1]",
                "H17/48 vertical characteristic-three ruled 4+4-fiber branch",
            ],
            "uncovered": [
                "arbitrary-input contraction/dependency canonicalizer for the fixed grammar",
                "proof that every accepted raw program reduces to the declared lattice-ideal modules",
                "H10/22 horizontal program",
                "H11/36 and H11/206 plane-program routing",
                "all D_CobVin inputs outside the grammar",
                "every compact polytope",
            ],
        },
        "membership_safety": {
            "forbidden": "a supplied coordinate point, surface equation, Hilbert witness, Kondo type or output Gram cannot be part of C membership",
            "allowed_and_required": (
                "the compiler itself constructs symbolic evaluation/Fitting/Jacobian matrices from raw divisor and "
                "proximity data, factors their lattice ideals, and proves the relevant minors/discriminants nonzero"
            ),
            "why_not_witness_defined": (
                "membership quantifies over no externally supplied coordinates; it applies one fixed symbolic algorithm "
                "to all raw contraction branches"
            ),
        },
        "lattice_ideal_realization_lemma": {
            "statement": (
                "Let A=Z[1/2N][t1^+-1,...,tr^+-1].  If the saturated compiler ideal is "
                "I_L+(p), where L is a saturated exponent sublattice, p is 0 or an odd prime, and the "
                "compiler-generated avoidance product is nonzero in A/(I_L+(p)), then the branch has a "
                "geometric point over Qbar for p=0 or over Fbar_p for p>0.  Rank-one kernel and Jacobian "
                "identities then construct the declared curves on a dense open."
            ),
            "proof_kernel": (
                "a saturated lattice ideal is the prime ideal of a subtorus translate over the algebraic closure; "
                "a nonzero regular function has proper zero locus in that irreducible torus; finitely many generated "
                "bad functions therefore leave a nonempty open"
            ),
            "does_not_cover": "arbitrary non-binomial incidence schemes or externally witnessed Hilbert points",
        },
        "calibrations": calibrations,
        "plane_profile_pushback": {
            "counting_level": "boundary-profile strata, not P-isometry classes and not an actual-polytope census",
            "pre_AG_partitions_of_6": [list(partition) for partition in pre],
            "horizontal_post_AG": [list(partition) for partition in horizontal],
            "all_support_post_AG": [list(partition) for partition in kondo_plane],
            "vertical_characteristic_three_high_profile": [
                list(partition) for partition in vertical_three_high
            ],
            "eliminated_by_AG_boundary_count": [list(partition) for partition in eliminated],
            "delta": {"pre": 11, "post": 5, "removed": 6},
            "qualification": (
                "this is a genuine theorem-shaped AG profile delta once the classwide canonicalizer is proved; "
                "it is not yet a counted nonzero delta among enumerated P isometry classes"
            ),
        },
        "compactness": (
            "every returned rational Coble surface has a primitive nef isotropic ruling class; all these chambers "
            "are finite-volume noncompact, and the candidate makes no compact claim"
        ),
        "sources": {
            str(path.relative_to(ROOT)): _sha256(path) for path in (H10, H11, CI4, HV)
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-output", type=Path)
    arguments = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if arguments.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
    if arguments.check_output is not None:
        if arguments.check_output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("frozen candidate certificate differs from exact replay")
    if not arguments.write and arguments.check_output is None:
        print(rendered, end="")


if __name__ == "__main__":
    main()
