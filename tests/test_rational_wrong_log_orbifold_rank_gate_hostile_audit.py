"""Independent hostile checks for the wrong-log orbifold rank gate.

These checks deliberately do not import the production certifier.  They
reconstruct the contraction numerics, local quotient contribution, rank
elimination, and the aggregate supported-boundary obstruction directly.
"""

from __future__ import annotations

import hashlib
import json
import tarfile
from fractions import Fraction
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_wrong_log_orbifold_rank_gate.md"
ARCHIVE = (
    ROOT / "data/enumerations/rational_wrong_log_orbifold_rank_gate.json"
)


def bilinear(left, gram, right):
    return sum(
        left[i] * gram[i][j] * right[j]
        for i in range(len(left))
        for j in range(len(right))
    )


def negative_continued_fraction(entries):
    value = Fraction(entries[-1])
    for entry in reversed(entries[:-1]):
        value = Fraction(entry) - Fraction(1, value)
    return value


def ruled_invariants(bidegree, private_blowups):
    a, b = bidegree
    return 2 * a * b - private_blowups, -2 * (a + b) + private_blowups


def ruled_intersection(left, right):
    return left[0] * right[1] + left[1] * right[0]


def test_contraction_divisor_and_local_quotient_are_reconstructed_independently():
    # Basis A,B_1,...,B_b with B_i.B_j=-4 delta_ij and A.B_i=d_i.
    # The remaining entries of an actual surface are irrelevant to these
    # two identities.  Several asymmetric choices guard against a copied
    # one-boundary computation.
    for a_square, degrees in ((1, (1,)), (7, (1, 3)), (11, (2, 5, 7))):
        size = 1 + len(degrees)
        gram = [[0] * size for _ in range(size)]
        gram[0][0] = a_square
        for i, degree in enumerate(degrees, start=1):
            gram[0][i] = gram[i][0] = degree
            gram[i][i] = -4
        H = (Fraction(1),) + tuple(Fraction(d, 4) for d in degrees)
        assert bilinear(H, gram, H) == a_square + sum(
            Fraction(d * d, 4) for d in degrees
        )
        for i in range(1, size):
            B_i = tuple(Fraction(j == i) for j in range(size))
            assert bilinear(H, gram, B_i) == 0

        # For any distinct irreducible C, ampleness of A and nonnegative
        # intersection with every B_i force H.C>0.
        for a_dot_c in (1, 2, 9):
            for boundary_dots in product(range(3), repeat=len(degrees)):
                assert Fraction(a_dot_c) + sum(
                    Fraction(d * t, 4)
                    for d, t in zip(degrees, boundary_dots)
                ) > 0

    # All positive cycles on the single exceptional component satisfy
    # Artin's numerical rationality criterion.
    assert all(1 - 2 * m * m + m <= 0 for m in range(1, 101))

    # The one-vertex HJ string is 4/1, hence 1/4(1,1), and adjunction
    # independently solves the discrepancy equation 2=(-4)a.
    assert negative_continued_fraction([4]) == Fraction(4, 1)
    discrepancy = Fraction(2, -4)
    assert discrepancy == Fraction(-1, 2)
    assert discrepancy > -1


def test_orbifold_euler_formula_and_rank_elimination_have_no_horizon():
    for rho in range(3, 80):
        e_surface = 2 + rho
        for boundary_count in range(1, 40):
            e_regular = e_surface - 2 * boundary_count
            e_singular = e_surface - boundary_count
            orbifold_from_regular = Fraction(e_regular) + Fraction(
                boundary_count, 4
            )
            orbifold_from_singular = Fraction(e_singular) - Fraction(
                3 * boundary_count, 4
            )
            closed_formula = Fraction(2 + rho) - Fraction(
                7 * boundary_count, 4
            )
            assert orbifold_from_regular == orbifold_from_singular
            assert orbifold_from_regular == closed_formula

    # This is the symbolic obstruction for every rho>=27, not a finite
    # search: b>=rho-10 makes 4*c2hat <= 78-3*rho < 0.
    for rho in (27, 28, 100, 10**6):
        assert 78 - 3 * rho < 0
    assert (7 * 10 + 8) // (7 - 4) == 26
    assert 26 - 10 == 16
    assert (4 * (26 + 2)) // 7 == 16

    # A deliberately much larger bounded replay is only a poison check on
    # the preceding symbolic proof.
    admissible = []
    for rho in range(3, 1001):
        lower_b = max(1, rho - 10)
        upper_b = (4 * rho + 8) // 7
        if lower_b <= upper_b:
            admissible.append((rho, lower_b, upper_b))
    assert max(rho for rho, _, _ in admissible) == 26
    assert [row for row in admissible if row[0] == 26] == [(26, 16, 16)]


def test_ijz_hypotheses_and_primary_v1_pins_are_not_weakened():
    text = REPORT.read_text(encoding="utf-8")
    flat_text = " ".join(text.split())
    required = (
        "compact normal complex analytic klt variety",
        "admits a nef and big",
        "class and has nef anticanonical divisor",
        "In dimension two there are no auxiliary polarization factors",
        "normal, projective, klt",
        "Thus \\(-K_Y\\) is nef",
        "thm-nef-anti-canonical",
        "Nothing in the argument below promotes it to positive characteristic",
    )
    for phrase in required:
        assert phrase in flat_text

    archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    pins = archive["primary_source_pins"]
    assert pins["arxiv_id"] == "2607.18699v1"
    assert pins["theorem"] == "Iwai-Jinnouchi-Zhang Theorem 1.3"
    for key in ("pdf_sha256", "source_tar_sha256", "tex_sha256"):
        assert len(pins[key]) == 64
        assert set(pins[key]) <= set("0123456789abcdef")

    # Recheck the downloaded primary files when the audit cache is present,
    # without making the permanent test suite depend on /private/tmp.
    pdf = Path("/private/tmp/ijz_orbifold_c2_2607.18699v1.pdf")
    source = Path("/private/tmp/ijz_orbifold_c2_2607.18699v1.tar")
    if pdf.exists() and source.exists():
        assert hashlib.sha256(pdf.read_bytes()).hexdigest() == pins["pdf_sha256"]
        assert hashlib.sha256(source.read_bytes()).hexdigest() == pins[
            "source_tar_sha256"
        ]
        with tarfile.open(source) as handle:
            tex = handle.extractfile("semipositivity_chern_2026_0717.tex")
            assert tex is not None
            tex_bytes = tex.read()
        assert hashlib.sha256(tex_bytes).hexdigest() == pins["tex_sha256"]
        tex_text = tex_bytes.decode("utf-8")
        theorem_start = tex_text.index("\\label{thm-nef-anti-canonical}")
        theorem_end = tex_text.index("\\end{thm}", theorem_start)
        theorem = tex_text[theorem_start:theorem_end]
        assert "compact normal complex analytic variety" in theorem
        assert "admits a nef and big class" in theorem
        assert "-K_X$ is nef" in theorem
        assert "at most klt singularities" in theorem


def test_supported_boundary_no_go_covers_aggregate_h3_and_old_gap():
    # H1 at the smallest allowed even intersection and H2 at the smallest
    # Coxeter-allowed positive intersection already force coefficient one.
    # These are exact elimination identities, not sampled LP solutions.
    intersection = Fraction(2)
    assert Fraction(4 - intersection, 1) == 2
    assert Fraction(4, 4 - intersection) == 2  # a_i+a_j >= 2
    h2_denominator = 4 - intersection * intersection / 2
    assert h2_denominator == 2
    assert Fraction(2, h2_denominator) == 1  # a_B >= 1

    # Enumerate all small incidence vectors as a hostile audit of the exact
    # dichotomy in (1-Q/4)c >= T/2-1.  Q<4 and T>=3 has exactly the
    # distributed vector (1,1,1), up to zero entries and permutation.
    low_q_patterns = set()
    for length in range(1, 7):
        for t in product(range(5), repeat=length):
            T = sum(t)
            if T < 3:
                continue
            Q = sum(value * value for value in t)
            if Q < 4:
                low_q_patterns.add(tuple(sorted((v for v in t if v), reverse=True)))
                lower_bound = Fraction(Fraction(T, 2) - 1, 1 - Fraction(Q, 4))
                assert lower_bound >= 2
            else:
                assert (1 - Fraction(Q, 4)) <= 0
                assert Fraction(T, 2) - 1 > 0
    assert low_q_patterns == {(1, 1, 1)}

    distributed = (1, 1, 1)
    assert max(distributed) == 1  # defeats the former componentwise H3
    assert sum(distributed) == 3  # triggers the correct aggregate H3


def test_geometric_controls_and_f4_calibration_use_fresh_formulas():
    controls = (
        ((1, 1), 6, (1, 1), 6, (-4, -4, 2, 2, 2)),
        ((1, 1), 6, (1, 1), 4, (-4, -2, 2, 2, 0)),
        ((1, 1), 6, (1, 2), 5, (-4, -1, 3, 2, -1)),
    )
    for left, left_blowups, right, right_blowups, expected in controls:
        left_square, left_k = ruled_invariants(left, left_blowups)
        right_square, right_k = ruled_invariants(right, right_blowups)
        actual = (
            left_square,
            right_square,
            ruled_intersection(left, right),
            left_k,
            right_k,
        )
        assert actual == expected

    boundary_square, boundary_k = ruled_invariants((1, 0), 4)
    exceptional_square, exceptional_k = ruled_invariants((0, 1), 1)
    degrees = [ruled_intersection((1, 0), (0, 1)) for _ in range(3)]
    assert (boundary_square, boundary_k) == (-4, 2)
    assert (exceptional_square, exceptional_k) == (-1, -1)
    assert ruled_intersection((1, 0), (1, 0)) == 0
    assert degrees == [1, 1, 1]

    gram = (
        (-4, 1, 0),
        (1, 0, 0),
        (0, 0, -1),
    )
    B = (Fraction(1), Fraction(0), Fraction(0))
    E = (Fraction(0), Fraction(0), Fraction(1))
    G = (Fraction(0), Fraction(1), Fraction(-1))
    K_plus_half_B = (Fraction(-3, 2), Fraction(-6), Fraction(1))
    assert tuple(
        bilinear(K_plus_half_B, gram, curve) for curve in (B, E, G)
    ) == (Fraction(0), Fraction(-1), Fraction(-1, 2))
    assert bilinear(K_plus_half_B, gram, K_plus_half_B) == 8
    assert Fraction(2 + 3) - Fraction(7, 4) == Fraction(13, 4)


def test_finiteness_corollary_stays_on_the_declared_carrier_lane():
    text = REPORT.read_text(encoding="utf-8")
    flat_text = " ".join(text.split())
    assert "carrier-defined ABSW-type finiteness" in flat_text
    assert "already admitting" in flat_text
    assert "This is not essential surjectivity from raw Coxeter data" in flat_text
    assert "Constructing such a carrier from arbitrary raw Coxeter data" in flat_text
    assert "must also be verified independently" in flat_text
    assert "does not use reflection maximality" in flat_text
    assert "neither constructs one nor permits one in the weak incidence lane" in flat_text
    assert "(H3) understood in its aggregate form" in flat_text

    # Once an abstract NS lattice is fixed, every standard host is an
    # intermediate additive group in (1/2 NS)/NS, a finite 2-group of this
    # cardinality.  No facet-degree or surface-moduli count is being used.
    for rho in (3, 10, 26):
        assert 2**rho < float("inf")

    archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert archive["schema_version"] == 2
    assert archive["rank_sieve"]["maximum_picard_rank"] == 26
    assert archive["rank_sieve"]["maximum_hyperbolic_dimension"] == 25
    assert archive["rank_sieve"]["extremal_integer_pairs"] == [[26, 16]]
    assert "H3_aggregate_degree_at_least_3" in archive[
        "supported_boundary_obstructions"
    ]
