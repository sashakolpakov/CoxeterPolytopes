"""Independent hostile checks for the canonical-chamber transfer synthesis."""

from __future__ import annotations

from fractions import Fraction
import hashlib
from math import comb, isqrt
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/canonical_chamber_rigidity_absw_transfer_theorem.md"
AUDIT = (
    ROOT
    / "results/canonical_chamber_rigidity_absw_transfer_theorem_hostile_audit.md"
)
PFIRST = ROOT / "results/rational_pfirst_effectivity_realizability_no_go.md"
TIMELIKE = ROOT / "results/rational_pfirst_timelike_anticanonical_vertex_closure.md"
TIMELIKE_AUDIT = (
    ROOT / "results/rational_pfirst_timelike_anticanonical_vertex_closure_hostile_audit.md"
)
CUSP = ROOT / "results/rational_pfirst_anticanonical_cusp_saturation_closure.md"
CUSP_AUDIT = (
    ROOT / "results/rational_pfirst_anticanonical_cusp_saturation_closure_hostile_audit.md"
)
DIRECT = ROOT / "results/rational_direct_full_ns_n11_dimension_theorem.md"
DIRECT_AUDIT = ROOT / "results/rational_direct_full_ns_n11_dimension_theorem_hostile_audit.md"
SEMIAMPLENESS = ROOT / "results/pfirst_weak_log_semiampleness_certificate.md"
WEAK_LOG_TRANSFER = ROOT / "results/pfirst_weak_log_complete_carrier_transfer.md"
WEAK_LOG_AUDIT = (
    ROOT / "results/pfirst_weak_log_complete_carrier_transfer_independent_hostile_audit.md"
)
K3_INTERVAL = ROOT / "results/rank4_k3_exact_maximal_eigen_interval_infinite_fields.md"
K3_INTERVAL_AUDIT = (
    ROOT / "results/rank4_k3_exact_maximal_eigen_interval_infinite_fields_hostile_scope_audit.md"
)
PINS = {
    REPORT: "850d2b5832c3f7bd9fa85c39a6e84449c8f1fb5ede6c98de931f966218da3767",
    PFIRST: "9c1b123caaee41baee9be7344aff28e9b17b8f8a12aeda2c443eece6d48ffd8a",
    TIMELIKE: "6450ab0629c539fcb3fff00e6aa2bdd0ddf194f17bdaaa2c1f2265d2b4aebb29",
    TIMELIKE_AUDIT: "2d7eddae058f81f06fded910b4d8ed45219b30e011af5be6a232f20655df0499",
    CUSP: "9d4142c8c4647bfced451ac1ff29184a6a278d9a51cdb1a0638276f3eb813032",
    CUSP_AUDIT: "1cc36249914c82ba42145e64cb0e1a43a888d4ac0c02d6b02084cade3ea498d8",
    DIRECT: "a5dbd6e2fbbd7cd55025f5af324de6b778a196cf4f72cf06903a1dce2e4ae27f",
    DIRECT_AUDIT: "9293c924bca806212ccac9a776ff1f477d9f9009c9ba1ce328dd81b1bc193193",
    SEMIAMPLENESS: "a96d6fd0332b3acf9e15afeda559d489436084321f8d8e523683c2e06fb73b58",
    WEAK_LOG_TRANSFER: "bff338b5e2f4dde2c918742a853656661c2ac074684938031ebbc95df8fd14b2",
    WEAK_LOG_AUDIT: "ce7cc977b2172418a7d6667c4a4a6e8ac8c9a3f3ec150d151ee91f0f7109cb3f",
    K3_INTERVAL: "9351875964fbc73da81b55ec3ff68db54af3dbe1ec90272555cff1938ab13c95",
    K3_INTERVAL_AUDIT: "3204a7bcfc580e24a0a015dd46f74e11d284a31ab33f42acdb6fa70c27220264",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compact(value: str) -> str:
    return " ".join(value.split())


def pair(left: tuple[int, ...], right: tuple[int, ...]) -> Fraction:
    """The rational Lorentz form diag(1,-1,-1)."""

    return Fraction(left[0] * right[0]) - sum(
        Fraction(a * b) for a, b in zip(left[1:], right[1:])
    )


def cartan(normals: tuple[tuple[int, ...], ...]) -> list[list[Fraction]]:
    return [
        [2 * pair(left, right) / pair(left, left) for right in normals]
        for left in normals
    ]


def cycle_product(matrix: list[list[Fraction]], cycle: tuple[int, ...]):
    result = Fraction(1)
    for index, source in enumerate(cycle):
        target = cycle[(index + 1) % len(cycle)]
        result *= matrix[source][target]
    return result


def test_frozen_inputs_are_exact() -> None:
    assert {path: sha256(path) for path in PINS} == PINS


def test_cartan_cycles_survive_independent_wall_scaling_exactly() -> None:
    normals = ((0, 1, 0), (1, 2, 0), (1, 1, 2))
    matrix = cartan(normals)
    scales = (Fraction(3, 5), Fraction(-7, 2), Fraction(11, 13))
    scaled = [
        [scales[j] * matrix[i][j] / scales[i] for j in range(3)]
        for i in range(3)
    ]
    for cycle in ((0,), (0, 1), (1, 2), (0, 1, 2), (0, 2, 1)):
        assert cycle_product(scaled, cycle) == cycle_product(matrix, cycle)
        assert cycle_product(matrix, cycle).denominator > 0

    # A common form scale cancels in every Cartan quotient.
    form_scale = Fraction(-17, 19)
    common_scaled = [
        [
            2 * form_scale * pair(left, right)
            / (form_scale * pair(left, left))
            for right in normals
        ]
        for left in normals
    ]
    assert common_scaled == matrix

    # Length-two normalization introduces square roots edgewise, but around
    # a closed cycle they cancel to this rational expression.
    cycle = (0, 1, 2)
    numerator = Fraction(1)
    denominator = Fraction(1)
    for index, source in enumerate(cycle):
        target = cycle[(index + 1) % len(cycle)]
        numerator *= pair(normals[source], normals[target])
        denominator *= -pair(normals[source], normals[source])
    assert numerator / denominator == Fraction(-1, 6)


def test_minimal_field_scope_and_scaling_formula_are_frozen() -> None:
    report = compact(REPORT.read_text(encoding="utf-8"))
    required = (
        "minimal Vinberg ground field",
        r"\((\lambda_j/\lambda_i)a_{ij}\)",
        "invariant under arbitrary nonzero wall rescaling",
        "common scaling of the intersection form",
        "irrationally normalized length-two wall vectors",
        "distinguished eigen-nef",
    )
    for phrase in required:
        assert phrase in report


def test_three_branch_transfer_has_the_actual_carrier_hypotheses() -> None:
    report = compact(REPORT.read_text(encoding="utf-8"))
    required = (
        r"Let \(n\ge2\)",
        "rational-polyhedral complete Mori cone",
        "Square-zero tangent Mori rays are allowed",
        "smooth projective **complex** rational surface",
        "possibly empty disjoint reduced boundary",
        r"\((-4)\)-curves",
        "adapted exponent-two host",
        "finite affine cover",
        "compatible global sections",
        "exact localized unit identities",
        "No complete Mori list is assumed",
        "entire integral trace lattice",
        "integral pairing isometries",
        "fixing the canonical class",
        "finite numerical Coxeter orbit",
        "not a property inferred from the existence of a weak NS embedding",
    )
    for phrase in required:
        assert phrase in report


def test_weak_log_branch_derives_nefness_completeness_and_rank_in_order() -> None:
    report = compact(REPORT.read_text(encoding="utf-8"))
    required = (
        r"make \(mD\) base-point-free and hence make \(D\) nef",
        "classifies every hidden negative curve",
        r"boundary/\((-2)\)/\((-1)\) adapted integral reflection alphabet",
        "proposed chamber is the entire nef chamber, without a complete Mori list",
        r"\(1/4(1,1)\) points",
        r"D^2=10-\rho+b",
        r"\widehat c_2=\rho+2-{7b\over4}",
        r"\(\rho\le26\), hence \(n\le25\)",
    )
    for phrase in required:
        assert phrase in report

    # Replay the elimination independently, including the empty-boundary
    # specialization.  The theorem keeps the uniform n<=25 conclusion.
    for rho in range(1, 80):
        for b in range(0, 80):
            if 10 - rho + b >= 0 and 4 * rho + 8 - 7 * b >= 0:
                assert rho <= 26
            if b == 0 and 10 - rho >= 0:
                assert rho <= 10


def test_weak_log_source_is_the_stronger_complete_transfer() -> None:
    report = compact(REPORT.read_text(encoding="utf-8"))
    assert "complete weak-log transfer proves items 2--3 for the stronger certified branch" in report
    transfer = compact(WEAK_LOG_TRANSFER.read_text(encoding="utf-8"))
    for phrase in (
        "strictly weaker on the geometric input than the complete-Mori weak-incidence theorem",
        "possibly empty disjoint union",
        "does not require a prior hidden-curve census",
        "arbitrary PB-Atlas B1/B2a-to-cover adapter is currently a precise contract",
    ):
        assert phrase in transfer
    semiampleness = compact(SEMIAMPLENESS.read_text(encoding="utf-8"))
    assert "exact adapter contract, not an implemented end-to-end adapter" in semiampleness
    assert "complete toric \\(\\mathbb F_4\\) provider" in semiampleness


def test_numeric_branch_maximum_and_characteristic_scope() -> None:
    # Direct <=11, weak <=25, and complex full trace <=19 because d>=1.
    ceilings = {"direct": 11, "weak": 25, "trace_complex": 20 // 1 - 1}
    assert max(ceilings.values()) == 25
    assert ceilings["trace_complex"] == 19
    assert 22 // 1 - 1 == 21

    report = compact(REPORT.read_text(encoding="utf-8"))
    assert r"d(n+1)\le20" in report
    assert r"last-row **numerical** ceiling changes to \(22\)" in report
    assert "does not promote the complex weak-incidence branch" in report
    assert "exhaustive complex compiler" in report
    assert "global positive-characteristic compiler" not in report


def test_bounded_degree_does_not_make_the_field_set_finite() -> None:
    # Distinct squarefree radicands give distinct quadratic fields.  This
    # finite prefix can be extended through all primes while degree stays 2.
    primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
    assert len(set(primes)) == len(primes)
    for prime in primes:
        assert all(prime % divisor for divisor in range(2, isqrt(prime) + 1))
    for i, left in enumerate(primes):
        for right in primes[i + 1 :]:
            # left/right cannot be a rational square for distinct primes.
            assert left * right != isqrt(left * right) ** 2

    # In contrast, bounded degree plus a uniform conjugate bound traps every
    # integral minimal polynomial in a finite coefficient box.
    for degree in range(1, 24):
        box_size = 1
        for j in range(1, degree + 1):
            coefficient_bound = comb(degree, j) * 14**j
            box_size *= 2 * coefficient_bound + 1
        assert isinstance(box_size, int) and box_size > 0
    assert 23**4 < 10**6


def test_exact_rank_one_k3_family_is_used_with_all_scope_limits() -> None:
    report = compact(REPORT.read_text(encoding="utf-8"))
    required = (
        "explicit rank-four K3 family",
        "exact finite-length distinguished eigen-nef intervals",
        "complete restricted-root saturation",
        r"declared \(\mathbf Z[t]\)-trace host",
        r"primitive ample class of fixed square \(14\)",
        "does not cover arbitrary overorders or commensurable hosts",
        r"does not settle the required \(n\ge2\) case",
        "fixed polarization square do not by themselves bound field height",
    )
    for phrase in required:
        assert phrase in report
    k3 = compact(K3_INTERVAL.read_text(encoding="utf-8"))
    assert "hyperbolic-dimension-one counterfamily" in k3
    assert "No maximality across arbitrary overorders" in k3


def test_finite_gram_and_absw_premises_are_not_weakened() -> None:
    report = compact(REPORT.read_text(encoding="utf-8"))
    finite_gram = (
        "bounded finite-type",
        "total** degrees",
        "global finite return grammar",
        "smooth proper",
        "Signed cycles must have both positive and negative parts bounded",
        "A polarization degree by itself",
    )
    absw = (
        "isometry classes** of integral reflection hosts",
        "complete normalized Gram packets",
        "all** integral hyperbolic reflections",
        "locally finite complete mirror arrangement",
        "one chamber orbit",
        r"\(L(2)\subseteq M\subseteq L(2)^*=\frac12L\)",
        "if it merely produces finitely many trace hosts",
    )
    for phrase in finite_gram + absw:
        assert phrase in report


def test_rank_three_counterexample_is_used_only_as_a_no_go() -> None:
    report = compact(REPORT.read_text(encoding="utf-8"))
    no_go = compact(PFIRST.read_text(encoding="utf-8"))
    assert "rank-three canonical-orientation counterexample" in report
    assert "Bare maximality cannot be substituted" in report
    assert "even absolute reflection" in no_go
    assert "maximality do not force a rational full-nef realization" in no_go
    assert r"\(K_S\) would be nef" in no_go
    assert "unbounded family" not in report


def test_property_four_progress_has_exact_ranges_and_next_boundary() -> None:
    report = compact(REPORT.read_text(encoding="utf-8"))
    for phrase in (
        r"ordinary anticanonical sequence \(I_{1,n}\), \(3\le n\le9\)",
        r"For \(3\le n\le8\)",
        "timelike vertex \\(-K\\)",
        r"For \(n=9\)",
        "primitive anticanonical cusp",
        r"twice-primitive \((-4)\) facets",
        "nonunimodular adapted hosts",
    ):
        assert phrase in report
    timelike_audit = compact(TIMELIKE_AUDIT.read_text(encoding="utf-8"))
    cusp_audit = compact(CUSP_AUDIT.read_text(encoding="utf-8"))
    assert r"PASS for \(I_{1,n}\), \(3\le n\le8\)" in timelike_audit
    assert "PASS" in cusp_audit and r"\(I_{1,9}\)" in cusp_audit


def test_every_local_markdown_dependency_resolves() -> None:
    report = REPORT.read_text(encoding="utf-8")
    targets = re.findall(r"\(([^)]+\.md)\)", report)
    assert targets
    missing = [target for target in targets if not (REPORT.parent / target).exists()]
    assert missing == []


def test_audit_records_conditional_scope_and_no_placeholders() -> None:
    audit = compact(AUDIT.read_text(encoding="utf-8"))
    for phrase in (
        "PASS at the stated carrier-relative scope",
        "global conclusion remains conditional on an exhaustive exact carrier compiler",
        "minimal Vinberg field",
        "The certified weak-log proof does not hide a Mori census",
        "Field degree is not field finiteness",
        "ABSW-type finiteness and property-4 scope",
        PINS[REPORT],
        PINS[PFIRST],
        PINS[SEMIAMPLENESS],
        PINS[WEAK_LOG_TRANSFER],
        PINS[K3_INTERVAL],
    ):
        assert phrase in audit
    combined = REPORT.read_text(encoding="utf-8") + audit
    assert "PLACEHOLDER" not in combined
    assert "d(n+1)le" not in combined
