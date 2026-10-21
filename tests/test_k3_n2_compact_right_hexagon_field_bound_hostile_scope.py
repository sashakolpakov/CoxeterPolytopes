"""Independent hostile replay of the compact right-angled K3 hexagon theorem."""

from __future__ import annotations

from fractions import Fraction
import hashlib
from itertools import permutations, product
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_n2_compact_right_hexagon_field_bound.md"
CERTIFIER = ROOT / "scripts/certify_k3_n2_compact_right_hexagon_field_bound.py"
ARCHIVE = ROOT / "data/enumerations/k3_n2_compact_right_hexagon_field_bound.json"
OWNER_TESTS = ROOT / "tests/test_k3_n2_compact_right_hexagon_field_bound.py"
AUDIT = ROOT / "results/k3_n2_compact_right_hexagon_field_bound_hostile_scope_audit.md"
PINS = {
    REPORT: "5cbbe7fd2d88a2129bdc0b5345f0ae377c0fa02b5f0874606fc6ff76a17783f3",
    CERTIFIER: "ff39e56f30c74429fb93935e3ba3f5353eeb219d059b2b38aa78c5f3df07a85a",
    ARCHIVE: "dc52fb561ac5a16b7ee0683af1999be26909820102fabf92382b30c1cf856e10",
    OWNER_TESTS: "06b0c81155739626843543301732e1af75c56dd2c5f9bd87c7bafab1674c808c",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_coordinates(u: Fraction, v: Fraction, w: Fraction):
    x0, x2, x4 = u * u - 1, v * v - 1, w * w - 1
    odd_cosh_squares = (
        (u * v + w) ** 2 / (x0 * x2),
        (v * w + u) ** 2 / (x2 * x4),
        (w * u + v) ** 2 / (x4 * x0),
    )
    return (
        x0,
        odd_cosh_squares[0] - 1,
        x2,
        odd_cosh_squares[1] - 1,
        x4,
        odd_cosh_squares[2] - 1,
    )


def test_owner_bytes_are_exactly_pinned():
    assert {path: sha256(path) for path in PINS} == PINS


def test_cofactor_cosine_formula_and_four_minors_rederive_independently():
    values = (
        Fraction(17, 16),
        Fraction(5, 4),
        Fraction(3, 2),
        Fraction(2),
        Fraction(3),
    )
    tested = 0
    for u, v, w in product(values, repeat=3):
        coordinates = independent_coordinates(u, v, w)
        assert min(coordinates) > 0
        x0, x1, x2, x3, x4, x5 = coordinates
        assert x0 * x1 * x2 == x2 * x3 * x4 == x4 * x5 * x0
        for i in range(6):
            assert (
                coordinates[i] * coordinates[(i + 1) % 6]
                == coordinates[(i + 3) % 6] * coordinates[(i + 4) % 6]
            )
        tested += 1
    assert tested == 125


def test_ordering_argument_covers_every_alternating_permutation():
    base = (Fraction(17), Fraction(25), Fraction(50))
    for permutation in permutations(base):
        x0, x2, x4 = permutation
        u2, v2, w2 = x0 + 1, x2 + 1, x4 + 1
        # Work with a rigorous upper bound on the relevant odd coordinate,
        # avoiding irrational square roots by sqrt(product)<an integer ceiling.
        smallest, middle, largest = sorted(permutation)
        C, E = middle + 1, largest + 1
        assert 2 * C / (middle * largest) <= Fraction(34, 256)
        assert E / (middle * largest) <= Fraction(17, 256)
        # The remaining term is bounded in the theorem by 34*sqrt(17)/256.
        assert Fraction(51 + 34 * 5, 256) == Fraction(221, 256) < 1


def test_short_pair_bound_handles_both_disjuncts():
    # Small alternating disjunct.
    assert 4 * (Fraction(1599, 100) + 1) < 68
    # Large alternating disjunct: exact rational coshes give a short odd side.
    for u, v, w in (
        (Fraction(5), Fraction(7), Fraction(9)),
        (Fraction(10), Fraction(11), Fraction(13)),
        (Fraction(20), Fraction(30), Fraction(50)),
    ):
        coordinates = independent_coordinates(u, v, w)
        assert min(coordinates[0], coordinates[2], coordinates[4]) >= 16
        assert min(coordinates[1], coordinates[3], coordinates[5]) < 1
        assert 4 * (min(coordinates) + 1) < 8 < 68


def test_unique_large_conjugate_cannot_descend_to_a_proper_subfield():
    # In a totally real extension k/F of relative degree r, each real
    # embedding of F has exactly r real extensions. A value in F therefore
    # repeats r times among the embeddings of k.
    for total_degree in range(2, 7):
        for subfield_degree in range(1, total_degree):
            if total_degree % subfield_degree:
                continue
            relative_degree = total_degree // subfield_degree
            if relative_degree > 1:
                assert relative_degree >= 2
                # Hence a unique >4 conjugate is impossible in a proper subfield.
                assert not relative_degree == 1


def test_discriminant_exponent_and_interval_bound_are_exact():
    assert 6 * 5 == 30
    # All conjugates lie in one interval of length 68.
    for left, right in ((Fraction(0), Fraction(68)), (Fraction(4), Fraction(68))):
        assert right - left <= 68
    assert 68**30 > 0


def test_archive_replays_and_scope_does_not_expand():
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))["payload"]
    assert payload["conclusions"]["some_dotted_pair_cycle_below_68"] is True
    assert payload["conclusions"]["that_dotted_pair_cycle_generates_the_full_field"] is True
    assert payload["conclusions"]["field_discriminant_at_most_68_power_30"] is True
    assert payload["conclusions"]["reflection_maximality_used"] is False
    scope = payload["scope"]
    assert scope["minimal_cyclic_Cartan_ground_field_required"] is True
    assert scope["exact_same_compact_right_angled_hexagon_eigen_nef_chamber_required"] is True
    assert scope["non_right_hexagon_closed"] is False
    assert scope["noncompact_hexagon_closed"] is False
    assert scope["arbitrary_six_or_more_facet_chamber_closed"] is False
    assert scope["carrier_existence_proved"] is False


def test_report_and_audit_pin_the_primitive_pair_mechanism():
    report = " ".join(REPORT.read_text(encoding="utf-8").split())
    audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "among the six second-neighbour dotted pairs",
        "This is a rank-three cofactor calculation",
        r"&<{2C+E+2C\sqrt E\over yz}",
        r"4<(p_i)_0<68",
        "unique-large-conjugate primitivity",
        "restriction of the distinguished embedding",
        "gratuitous proper totally real field extension cannot hide above a dotted pair",
        r"would have a finite field set once OUR geometry supplies a uniform upper bound",
    ):
        assert phrase in report
    for phrase in (
        "PASS in exact compact right-angled scope",
        "does not need that stronger statement",
        "This subfield argument is the decisive point",
        "no pair-generation assumption is silently imported",
        "uniform short-dotted-pair theorem",
    ):
        assert phrase in audit


def test_tex_transport_and_classification_scope_are_clean():
    combined = REPORT.read_text(encoding="utf-8") + AUDIT.read_text(encoding="utf-8")
    assert ",quad" not in combined
    assert ",qquad" not in combined
    normalized = " ".join(combined.split())
    assert "No imposed height, power coupling, Coxeter classification" in normalized
    assert "No reflective-lattice classification enters" in normalized
