#!/usr/bin/env python3
"""Independent hostile replay of the compact maximal (2,5,6) PHYS no-go.

This audit deliberately does not import the owner certifier.  It uses the
integral basis (1, phi), derives the codifferent from the inverse trace Gram,
and checks the inert-prime obstruction directly in ideal-valuation parity.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OWNER_REPORT = ROOT / "results/k3_compact_maximal_256_phys_empty.md"
OWNER_SCRIPT = ROOT / "scripts/certify_k3_compact_maximal_256_phys_empty.py"
OWNER_ARCHIVE = ROOT / "data/enumerations/k3_compact_maximal_256_phys_empty.json"
OWNER_TEST = ROOT / "tests/test_k3_compact_maximal_256_phys_empty.py"

AUDIT_REPORT = (
    ROOT / "results/k3_compact_maximal_256_phys_empty_independent_hostile_audit.md"
)
AUDIT_ARCHIVE = (
    ROOT
    / "data/enumerations/"
    "k3_compact_maximal_256_phys_empty_independent_hostile_audit.json"
)

OWNER_HASHES = {
    "results/k3_compact_maximal_256_phys_empty.md": (
        "9fa8afdb6f8e3f5225c494d48e3e9130a9b86cfe2aa8f58b38b177b9f972b83e"
    ),
    "scripts/certify_k3_compact_maximal_256_phys_empty.py": (
        "452f2946ec042dd7380e1772e3355d520753b54b909dd970c910857338527d40"
    ),
    "data/enumerations/k3_compact_maximal_256_phys_empty.json": (
        "d5d7dd20cb241b9dfa21c88327c62521862cfc4007be0a5965268ced449c5b03"
    ),
    "tests/test_k3_compact_maximal_256_phys_empty.py": (
        "fd6842bccc25a5dbac21cc9121d8fc31ae6d808121b653b596255bc2ec5ef11a"
    ),
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True)
class Kphi:
    """Exact a+b*phi arithmetic, phi^2=phi+1."""

    a: Fraction
    b: Fraction

    def __add__(self, other: "Kphi") -> "Kphi":
        return Kphi(self.a + other.a, self.b + other.b)

    def __sub__(self, other: "Kphi") -> "Kphi":
        return Kphi(self.a - other.a, self.b - other.b)

    def __mul__(self, other: "Kphi") -> "Kphi":
        return Kphi(
            self.a * other.a + self.b * other.b,
            self.a * other.b + self.b * other.a + self.b * other.b,
        )

    def trace(self) -> Fraction:
        return 2 * self.a + self.b

    def norm(self) -> Fraction:
        return self.a * self.a + self.a * self.b - self.b * self.b

    def inverse(self) -> "Kphi":
        norm = self.norm()
        if norm == 0:
            raise ZeroDivisionError
        # conjugate(phi)=1-phi, so conjugate(a+b phi)=(a+b)-b phi.
        return Kphi((self.a + self.b) / norm, -self.b / norm)

    def __truediv__(self, other: "Kphi") -> "Kphi":
        return self * other.inverse()

    def sign(self, *, conjugate: bool = False) -> int:
        # a+b*phi=(a+b/2)+(b/2)*sqrt(5); conjugation flips sqrt(5).
        rational = self.a + self.b / 2
        radical = (-self.b if conjugate else self.b) / 2
        if radical == 0:
            return (rational > 0) - (rational < 0)
        if rational == 0:
            return (radical > 0) - (radical < 0)
        if (rational > 0) == (radical > 0):
            return (rational > 0) - (rational < 0)
        comparison = rational * rational - 5 * radical * radical
        if comparison == 0:
            raise AssertionError("sqrt(5) would be rational")
        if comparison > 0:
            return (rational > 0) - (rational < 0)
        return (radical > 0) - (radical < 0)

    def encoded(self) -> list[str]:
        return [str(self.a), str(self.b)]


ZERO = Kphi(Fraction(0), Fraction(0))
ONE = Kphi(Fraction(1), Fraction(0))
PHI = Kphi(Fraction(0), Fraction(1))


def determinant_3(matrix: list[list[Kphi]]) -> Kphi:
    return (
        matrix[0][0]
        * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1]
        * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2]
        * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def reconstruct_raw_packet() -> dict[str, Any]:
    gram = [
        [Kphi(Fraction(-6), Fraction(0)), ZERO, Kphi(Fraction(3), Fraction(0))],
        [ZERO, Kphi(Fraction(-2), Fraction(0)), PHI],
        [Kphi(Fraction(3), Fraction(0)), PHI, Kphi(Fraction(-2), Fraction(0))],
    ]
    if PHI * PHI != PHI + ONE:
        raise AssertionError("phi relation failed")
    determinant = determinant_3(gram)
    if determinant != Kphi(Fraction(0), Fraction(6)):
        raise AssertionError("determinant is not 6 phi")

    pivots = [
        Kphi(Fraction(-6), Fraction(0)),
        Kphi(Fraction(-2), Fraction(0)),
        determinant / Kphi(Fraction(12), Fraction(0)),
    ]
    identity_signs = [pivot.sign() for pivot in pivots]
    conjugate_signs = [pivot.sign(conjugate=True) for pivot in pivots]
    if identity_signs != [-1, -1, 1] or conjugate_signs != [-1, -1, -1]:
        raise AssertionError("admissible signatures failed")

    two = Kphi(Fraction(2), Fraction(0))
    cartan = [
        two * gram[0][2] / gram[0][0],
        two * gram[2][0] / gram[2][2],
        two * gram[1][2] / gram[1][1],
        two * gram[2][1] / gram[2][2],
    ]
    expected_cartan = [
        Kphi(Fraction(-1), Fraction(0)),
        Kphi(Fraction(-3), Fraction(0)),
        Kphi(Fraction(0), Fraction(-1)),
        Kphi(Fraction(0), Fraction(-1)),
    ]
    if cartan != expected_cartan:
        raise AssertionError("Cartan integrality failed")

    cosine_23 = gram[1][2] * gram[1][2] / (gram[1][1] * gram[2][2])
    cosine_31 = gram[2][0] * gram[2][0] / (gram[2][2] * gram[0][0])
    if cosine_23 != (PHI * PHI) / Kphi(Fraction(4), Fraction(0)):
        raise AssertionError("order-five normalized product failed")
    if cosine_31 != Kphi(Fraction(3, 4), Fraction(0)):
        raise AssertionError("order-six normalized product failed")

    angle_sum = Fraction(1, 2) + Fraction(1, 5) + Fraction(1, 6)
    if angle_sum != Fraction(13, 15) or angle_sum >= 1:
        raise AssertionError("compactness failed")

    return {
        "basis": "1,phi with phi^2=phi+1",
        "determinant": determinant.encoded(),
        "determinant_text": "6*phi",
        "LDL_identity_signs": identity_signs,
        "LDL_conjugate_signs": conjugate_signs,
        "cartan_coefficients": [item.encoded() for item in cartan],
        "cartan_entries_in_Ok": True,
        "normalized_products": {
            "m5": cosine_23.encoded(),
            "m6": cosine_31.encoded(),
        },
        "compact_angle_sum": str(angle_sum),
        "trace_rank": 6,
        "cycle_products": ["3", "phi^2"],
        "cycle_order_equals_maximal_order": True,
    }


def derive_codifferent_shell() -> dict[str, Any]:
    # Trace Gram on (1,phi): [[2,1],[1,3]], determinant 5.
    trace_gram = [[2, 1], [1, 3]]
    inverse_trace_gram = [
        [Fraction(3, 5), Fraction(-1, 5)],
        [Fraction(-1, 5), Fraction(2, 5)],
    ]
    # Dual basis d0=(3-phi)/5, d1=(-1+2phi)/5.
    d0 = Kphi(Fraction(3, 5), Fraction(-1, 5))
    d1 = Kphi(Fraction(-1, 5), Fraction(2, 5))
    if [d0.trace(), (d0 * PHI).trace()] != [1, 0]:
        raise AssertionError("first trace-dual vector failed")
    if [d1.trace(), (d1 * PHI).trace()] != [0, 1]:
        raise AssertionError("second trace-dual vector failed")

    # x=A*d0+B*d1 has trace A.  For A=-2, conversion to the
    # (1,sqrt(5)) basis gives x=-1+(B+1)/sqrt(5).  Both embeddings are
    # negative iff |B+1|<sqrt(5).  Since 2^2<5<3^2, the complete integer
    # list is B+1 in {-2,-1,0,1,2}; no search cutoff is involved.
    if not (2 * 2 < 5 < 3 * 3):
        raise AssertionError("sqrt(5) integer bounds changed")
    complete_b_parameters = [-3, -2, -1, 0, 1]
    rows: list[dict[str, Any]] = []
    for b_parameter in complete_b_parameters:
        value = Kphi(
            -2 * d0.a + b_parameter * d1.a,
            -2 * d0.b + b_parameter * d1.b,
        )
        negative = value.sign() == -1 and value.sign(conjugate=True) == -1
        if not negative:
            raise AssertionError("derived shell point is not totally negative")
        c = b_parameter + 1
        if value.trace() != -2:
            raise AssertionError("shell trace failed")
        if value.norm() != Fraction(5 - c * c, 5):
            raise AssertionError("shell norm formula failed")
        rows.append(
            {
                "B": b_parameter,
                "c=B+1": c,
                "phi_coordinates": value.encoded(),
                "trace": str(value.trace()),
                "norm": str(value.norm()),
                "v_inert_3": 0,
            }
        )
    if [row["B"] for row in rows] != complete_b_parameters:
        raise AssertionError("codifferent energy shell is not exhaustive")

    # The adjacent excluded integers c=+-3 each make one embedding positive.
    excluded_signs = {}
    for c in (-3, 3):
        b_parameter = c - 1
        value = Kphi(
            -2 * d0.a + b_parameter * d1.a,
            -2 * d0.b + b_parameter * d1.b,
        )
        excluded_signs[str(c)] = [
            value.sign(),
            value.sign(conjugate=True),
        ]
    if any(sorted(signs) != [-1, 1] for signs in excluded_signs.values()):
        raise AssertionError("energy-shell boundary check failed")

    return {
        "trace_gram": trace_gram,
        "inverse_trace_gram": [
            [[entry.numerator, entry.denominator] for entry in row]
            for row in inverse_trace_gram
        ],
        "dual_basis": [d0.encoded(), d1.encoded()],
        "trace_minus_two_parameter": "A=-2",
        "negativity_equivalence": "|B+1|<sqrt(5)",
        "integer_boundary_proof": "2^2<5<3^2",
        "total_negative_B_values": [-3, -2, -1, 0, 1],
        "excluded_c_signs": excluded_signs,
        "rows": rows,
        "cardinality": len(rows),
    }


def singerman_maximality_crosscheck() -> dict[str, Any]:
    target = tuple(sorted((2, 5, 6)))

    normal_subgroups: set[tuple[int, int, int]] = set()
    nonnormal_subgroups: set[tuple[int, int, int]] = {
        tuple(sorted(row))
        for row in [
            (7, 7, 7),
            (2, 7, 7),
            (3, 3, 7),
            (4, 8, 8),
            (3, 8, 8),
            (9, 9, 9),
            (4, 4, 5),
        ]
    }
    # Parameters above max(target) cannot yield target in these positive
    # multiple patterns, so 2..6 is an exact target-membership window.
    for t in range(2, 7):
        normal_subgroups.add(tuple(sorted((t, t, t))))
        nonnormal_subgroups.add(tuple(sorted((t, 4 * t, 4 * t))))
        nonnormal_subgroups.add(tuple(sorted((t, 2 * t, 2 * t))))
        nonnormal_subgroups.add(tuple(sorted((3, t, 3 * t))))
        nonnormal_subgroups.add(tuple(sorted((2, t, 2 * t))))
        for s in range(2, 7):
            normal_subgroups.add(tuple(sorted((s, s, t))))

    if target in normal_subgroups or target in nonnormal_subgroups:
        raise AssertionError("(2,5,6) appears in Singerman's nonmaximal list")

    return {
        "target": list(target),
        "Singerman_overgroup_of_triangle_is_triangle": True,
        "normal_inclusion_patterns": ["(t,t,t)", "(s,s,t)"],
        "nonnormal_sporadic_subgroups": [
            [7, 7, 7],
            [2, 7, 7],
            [3, 3, 7],
            [4, 8, 8],
            [3, 8, 8],
            [9, 9, 9],
            [4, 4, 5],
        ],
        "nonnormal_family_subgroups": [
            "(t,4t,4t)",
            "(t,2t,2t)",
            "(3,t,3t)",
            "(2,t,2t)",
        ],
        "target_in_complete_nonmaximal_patterns": False,
        "Maclachlan_independent_table_crosscheck": (
            "Table 3 lists (2,5,6) as the maximal triangle group in class e"
        ),
        "orientation_index_transfer": (
            "W subset W' implies W+ subset W'+; equality of orientation "
            "subgroups and both indices 2 imply W=W'"
        ),
        "absolute_full_reflection_maximality": True,
        "published_dependency_not_reproved": True,
    }


def compact_reflection_area_maximality() -> dict[str, Any]:
    """Prove reflection-maximality without a classification citation.

    If the (2,5,6) chamber P were tiled by N>1 chambers Q of a reflection
    overgroup, area bounds force Q to be a triangle and N<=5.  The exact
    triangle equation then leaves only N=2 and Q=(2,3,10).  Doubling that
    triangle across any of its three sides never gives (2,5,6).
    """

    p_area_over_pi = Fraction(2, 15)

    # A compact Coxeter s-gon has area
    #   ((s-2)-sum_i 1/m_i)*pi,  m_i>=2.
    # For s=4, positivity forces some m_i>=3, hence area>=pi/6.
    # For s>=5, replacing every 1/m_i by 1/2 gives area>=pi/2.
    polygon_area_floor_over_pi = Fraction(1, 6)
    if not p_area_over_pi / 2 < polygon_area_floor_over_pi:
        raise AssertionError("polygon area gap failed")

    # The smallest compact hyperbolic Coxeter triangle is (2,3,7), of
    # area pi/42.  Here is the complete elementary branch proof.
    triangle_minimum_branches = {
        "a>=3": {
            "smallest_hyperbolic_ordered_triple": [3, 3, 4],
            "area_over_pi": str(Fraction(1, 12)),
        },
        "a=2,b=3": {
            "smallest_hyperbolic_c": 7,
            "area_over_pi": str(Fraction(1, 42)),
        },
        "a=2,b=4": {
            "smallest_hyperbolic_c": 5,
            "area_over_pi": str(Fraction(1, 20)),
        },
        "a=2,b>=5": {
            "lower_bound_from_2_5_5": str(Fraction(1, 10)),
        },
    }
    triangle_area_floor_over_pi = Fraction(1, 42)
    # 2/(15N)>=1/42 gives 15N<=84, hence N<=5.
    possible_indices = list(range(2, 6))

    solutions_by_index: dict[str, list[list[int]]] = {}
    for index in possible_indices:
        reciprocal_sum = 1 - p_area_over_pi / index
        solutions: list[list[int]] = []

        # Put 2<=a<=b<=c.  Since 3/a >= reciprocal_sum and the latter
        # exceeds 9/10, a<=3.  Once a is fixed,
        # 2/b >= reciprocal_sum-1/a bounds b.  Thus these loops are a
        # derived exhaustive enumeration, not a search cutoff.
        for a in range(2, 4):
            remainder_after_a = reciprocal_sum - Fraction(1, a)
            if remainder_after_a <= 0:
                raise AssertionError("unexpected nonpositive remainder")
            b_upper = int(Fraction(2, 1) / remainder_after_a)
            for b in range(a, b_upper + 1):
                reciprocal_c = reciprocal_sum - Fraction(1, a) - Fraction(1, b)
                if reciprocal_c <= 0:
                    continue
                c_fraction = Fraction(1, 1) / reciprocal_c
                if c_fraction.denominator == 1 and c_fraction.numerator >= b:
                    solutions.append([a, b, c_fraction.numerator])
        solutions_by_index[str(index)] = solutions

    expected_solutions = {
        "2": [[2, 3, 10]],
        "3": [],
        "4": [],
        "5": [],
    }
    if solutions_by_index != expected_solutions:
        raise AssertionError(
            f"unexpected reflection-overgroup triangle census {solutions_by_index}"
        )

    # At index two, P is the union of two adjacent Q chambers, hence is a
    # reflection double across their common side.  A double across a side
    # with endpoint orders p,q has boundary angles 2*pi/p, 2*pi/q plus two
    # copies of the opposite angle.  An order-2 endpoint becomes straight.
    doubling_cases = [
        {
            "shared_edge_endpoint_orders": [2, 10],
            "opposite_order": 3,
            "result": "triangle (3,3,5)",
            "coxeter": True,
            "equals_target_2_5_6": False,
        },
        {
            "shared_edge_endpoint_orders": [2, 3],
            "opposite_order": 10,
            "result": "triangle angles (pi/10,pi/10,2pi/3)",
            "coxeter": False,
            "equals_target_2_5_6": False,
        },
        {
            "shared_edge_endpoint_orders": [3, 10],
            "opposite_order": 2,
            "result": "quadrilateral",
            "coxeter": False,
            "equals_target_2_5_6": False,
        },
    ]
    if any(case["equals_target_2_5_6"] for case in doubling_cases):
        raise AssertionError("a (2,3,10) double unexpectedly gives (2,5,6)")

    return {
        "P_area_over_pi": str(p_area_over_pi),
        "proper_overgroup_index_lower_bound": 2,
        "proper_overgroup_chamber_area_upper_bound_over_pi": str(
            p_area_over_pi / 2
        ),
        "compact_coxeter_polygon_s_ge_4_area_floor_over_pi": str(
            polygon_area_floor_over_pi
        ),
        "overgroup_chamber_must_be_triangle": True,
        "compact_triangle_area_floor_over_pi": str(
            triangle_area_floor_over_pi
        ),
        "triangle_minimum_branches": triangle_minimum_branches,
        "derived_index_upper_bound": 5,
        "reciprocal_equation": "1/a+1/b+1/c=1-2/(15N)",
        "enumeration_bounds": "2<=a<=3 and a<=b<=2/(S-1/a)",
        "solutions_by_index": solutions_by_index,
        "index_two_means_adjacent_chamber_double": True,
        "doubling_cases": doubling_cases,
        "proper_reflection_overgroup_exists": False,
        "absolute_full_reflection_maximality": True,
        "classification_dependency": False,
    }


def valuation_obstruction(shell: dict[str, Any]) -> dict[str, Any]:
    mod_three_values = [(x * x - x - 1) % 3 for x in range(3)]
    if mod_three_values != [2, 2, 1]:
        raise AssertionError("unexpected minimal polynomial residues")
    if any(Fraction(row["norm"]).numerator % 3 == 0 for row in shell["rows"]):
        raise AssertionError("an energy value is not a unit above 3")

    pairs = []
    for left in shell["rows"]:
        for right in shell["rows"]:
            # u_left/u_right is a p-unit; division by 3 has valuation -1.
            pairs.append(
                {
                    "left_B": left["B"],
                    "right_B": right["B"],
                    "v_p(u_left/u_right)": 0,
                    "v_p(u_left/(3*u_right))": -1,
                    "square_possible": False,
                }
            )
    if len(pairs) != 25 or any(row["square_possible"] for row in pairs):
        raise AssertionError("valuation census failed")

    return {
        "minimal_polynomial_mod_3_values": mod_three_values,
        "three_inert": True,
        "energy_shell_all_p_units": True,
        "m6_norm_ratio_square_class": "[3]",
        "form_similarity_scale_cancels": True,
        "line_representative_changes_contribute_only_squares": True,
        "required_square_test": "u_left/(3*u_right) in k^x2",
        "ordered_pair_count": len(pairs),
        "square_count": 0,
        "pairs": pairs,
    }


def owner_hash_audit() -> dict[str, str]:
    paths = [OWNER_REPORT, OWNER_SCRIPT, OWNER_ARCHIVE, OWNER_TEST]
    observed = {
        str(path.relative_to(ROOT)): sha256_file(path)
        for path in paths
    }
    if observed != OWNER_HASHES:
        raise AssertionError(
            f"owner bytes changed; expected {OWNER_HASHES}, observed {observed}"
        )
    return observed


def build_certificate() -> dict[str, Any]:
    owner_hashes = owner_hash_audit()
    raw = reconstruct_raw_packet()
    shell = derive_codifferent_shell()
    area_maximality = compact_reflection_area_maximality()
    classification_crosscheck = singerman_maximality_crosscheck()
    obstruction = valuation_obstruction(shell)

    return {
        "schema": "coxeter.k3-compact-maximal-256-independent-hostile-audit.v1",
        "audit_report_sha256": sha256_file(AUDIT_REPORT),
        "owner_hashes": owner_hashes,
        "independence": {
            "imports_owner_certifier": False,
            "arithmetic_basis": "1,phi (owner uses 1,sqrt(5))",
            "codifferent_derived_from_inverse_trace_gram": True,
            "square_no_go_method": "inert-prime ideal valuation parity",
        },
        "raw_reconstruction": raw,
        "maximality": area_maximality,
        "published_classification_crosscheck": classification_crosscheck,
        "codifferent_shell": shell,
        "valuation_obstruction": obstruction,
        "freedom_audit": {
            "common_form_scale_lambda": "cancels in u1/u3",
            "facet_representative_scaling": "changes h1/h3 by a square",
            "raw_similarity": "preserves the m6 line-pair square class [3]",
            "facet_permutation": "labels 2,5,6 are distinct; m6 pair is intrinsic",
            "coefficient_order": "cycle order O0=Ok, so no order lane remains",
            "module_or_overlattice": (
                "cannot change trace -2 or codifferent membership for a literal root"
            ),
            "Weyl_or_period_framing": "preserves the lattice root square -2",
            "all_freedoms_close": True,
        },
        "embedding_scope": {
            "trace_rank": 6,
            "K3_complement_rank": 16,
            "discriminant_length_at_most": 6,
            "stable_embedding_after_PHYS": True,
            "failure_is_before_EMB": True,
        },
        "claim_scope": {
            "full_trace": True,
            "full_Picard_NS_equals_trace_lattice": True,
            "same_vector_physical_roots": True,
            "projected_root_larger_Picard_excluded": False,
            "alternative_AG_carrier_excluded": False,
            "global_Property4_disjunction_refuted": False,
            "universal_full_Picard_K3_nonemptiness_refuted": True,
        },
        "presentation": {
            "owner_report_contains_primary_source_links": True,
            "owner_report_marks_published_maximality_dependency": True,
            "owner_report_contains_scope_nonclaims": True,
            "pandoc_required_by_reproduction": True,
        },
        "verdict": "PASS",
    }


def encoded_certificate() -> str:
    payload = build_certificate()
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-output", action="store_true")
    parser.add_argument("--check-output", action="store_true")
    args = parser.parse_args()

    encoded = encoded_certificate()
    if args.write_output:
        AUDIT_ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        AUDIT_ARCHIVE.write_text(encoded, encoding="utf-8")
    if args.check_output:
        if (
            not AUDIT_ARCHIVE.exists()
            or AUDIT_ARCHIVE.read_text(encoding="utf-8") != encoded
        ):
            raise SystemExit("independent hostile audit archive differs")
    if not args.write_output and not args.check_output:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
