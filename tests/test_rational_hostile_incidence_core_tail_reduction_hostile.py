import hashlib
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SYNTHESIS = ROOT / "results/rational_hostile_incidence_core_tail_reduction.md"
AUDIT = (
    ROOT
    / "results/rational_hostile_incidence_core_tail_reduction_hostile_audit.md"
)

EXPECTED_HASHES = {
    "results/rational_hostile_incidence_core_tail_reduction.md": (
        "092e2258867d243c748d2995f3243b4e4886c0b54b3effe3e492b0326b29175d"
    ),
    "results/rational_hostile_incidence_core_tail_reduction_hostile_audit.md": (
        "df1d7406894177be147db2b0881bc149de598e01e2704970af68ba7e844dff5d"
    ),
    "results/rational_finite_vertex_discrepancy_tail_gate.md": (
        "817d30726dad66f7512c9eaa115823f34c2f38d7e86224876b4d50613b9a4548"
    ),
    "results/rational_finite_vertex_discrepancy_tail_gate_hostile_audit.md": (
        "8a20e688292a2e5abca3bbc492c27e430d18153ee967d0cf74eb30ec37afc507"
    ),
    "results/rational_h1_intersecting_minus4_log_core.md": (
        "5c296bbaf3893e7f9f0021cc20ea101cc71618c07b2dad1255a41f2beea6b927"
    ),
    "results/rational_h1_intersecting_minus4_log_core_hostile_audit.md": (
        "b13b452d98908fdf47e45c3faf90d0c6761881688ff7d1e144711e745aa288b2"
    ),
    "results/rational_h2_minus4_minus2_genus_one_core.md": (
        "1464b327d588f9a5e8c5b2a6809667194d27c09a69aad5805250df4053c9634b"
    ),
    "results/rational_h2_minus4_minus2_genus_one_core_hostile_audit.md": (
        "d90b45deb81f4305bc417642328f5077dc5c1c49657022e30364864e4428af1c"
    ),
    "results/rational_h3_hostile_incidence_trichotomy.md": (
        "c0bae8e9447f0793d21a059fd43c13b23bd4b26ff75d5947925965fcab2e5d44"
    ),
    "results/rational_h3_hostile_incidence_trichotomy_hostile_audit.md": (
        "c9096bae163c6339b02014a2c4c71c5e4cbba5df40c45ae45bf0a8aecdfd5029"
    ),
    "results/rational_wrong_log_orbifold_rank_gate.md": (
        "36ef3a14968300119f4f0178e9dad11003220d4ee2879cc49c9f0494dc961afa"
    ),
    "results/rational_wrong_log_orbifold_rank_gate_hostile_audit.md": (
        "cce94b8013d853e2e13212fa20c09ecb3d042a15a0b335be32546ba4d6cc41a2"
    ),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _dot(left: list[int], gram: list[list[int]], right: list[int]) -> int:
    return sum(
        left[i] * gram[i][j] * right[j]
        for i in range(len(left))
        for j in range(len(right))
    )


def test_frozen_dependency_hashes_and_links():
    for relative, expected in EXPECTED_HASHES.items():
        path = ROOT / relative
        assert path.is_file(), relative
        assert _sha256(path) == expected, relative

    synthesis = SYNTHESIS.read_text(encoding="utf-8")
    for dependency in (
        "rational_finite_vertex_discrepancy_tail_gate.md",
        "rational_finite_vertex_discrepancy_tail_gate_hostile_audit.md",
        "rational_h1_intersecting_minus4_log_core.md",
        "rational_h1_intersecting_minus4_log_core_hostile_audit.md",
        "rational_h2_minus4_minus2_genus_one_core.md",
        "rational_h2_minus4_minus2_genus_one_core_hostile_audit.md",
        "rational_h3_hostile_incidence_trichotomy.md",
        "rational_h3_hostile_incidence_trichotomy_hostile_audit.md",
        "rational_wrong_log_orbifold_rank_gate.md",
        "rational_wrong_log_orbifold_rank_gate_hostile_audit.md",
    ):
        assert f"]({dependency})" in synthesis


def test_H1_H2_are_only_minimal_finite_angle_cores():
    # H1 roots B_i/2 have norm -2 and pairing I/2.
    h1_ratio_I2 = Fraction(1 * 1, 2 * 2)
    h1_ratio_I4 = Fraction(2 * 2, 2 * 2)
    assert h1_ratio_I2 == Fraction(1, 4)  # cos(pi/3)^2
    assert h1_ratio_I4 == 1  # parallel

    # H2 roots B/2,R have norms -2,-4 and pairing I.
    h2_ratio_I2 = Fraction(2 * 2, 2 * 4)
    h2_ratio_I3 = Fraction(3 * 3, 2 * 4)
    assert h2_ratio_I2 == Fraction(1, 2)  # cos(pi/4)^2
    assert h2_ratio_I3 > 1  # ultraparallel/dotted

    # K^2=10-rho; the remaining constants come from 2K.D+D^2.
    assert 10 + 2 * 4 - 4 == 14
    assert 10 + 2 * 2 - 2 == 12


def test_H3_star_and_host_factor_two_conventions():
    for incidences in ([1, 1, 1], [1, 1, 1, 1], [2, 1, 1]):
        t = list(incidences)
        size = 1 + len(t)
        physical = [[0] * size for _ in range(size)]
        physical[0][0] = -1
        for i, value in enumerate(t, start=1):
            physical[i][i] = -4
            physical[0][i] = physical[i][0] = value
        F = [4] + t
        Q = sum(value * value for value in t)
        T = sum(t)
        K_degrees = [-1] + [2] * len(t)
        assert [_dot(F, physical, [int(j == i) for j in range(size)]) for i in range(size)] == [Q - 4] + [0] * len(t)
        assert _dot(F, physical, F) == 4 * (Q - 4)
        assert sum(a * b for a, b in zip(K_degrees, F, strict=True)) == -4 + 2 * T

        host = [[0] * size for _ in range(size)]
        for i in range(size):
            host[i][i] = -2
        for i, value in enumerate(t, start=1):
            host[0][i] = host[i][0] = value
        r = [2] + t  # F/2=2E+sum t_i(B_i/2)
        assert [_dot(r, host, [int(j == i) for j in range(size)]) for i in range(size)] == [Q - 4] + [0] * len(t)
        assert _dot(r, host, r) == 2 * (Q - 4)

    # q_M=2q_NS, with roots R,E,B/2 respectively.
    physical_rows = {
        "R": (-2, 0, Fraction(1)),
        "E": (-1, -1, Fraction(1)),
        "B": (-4, 2, Fraction(1, 2)),
    }
    colors = {}
    for name, (square, canonical_degree, scale) in physical_rows.items():
        colors[name] = (
            2 * square * scale * scale,
            2 * canonical_degree * scale,
        )
    assert colors == {
        "R": (Fraction(-4), Fraction(0)),
        "E": (Fraction(-2), Fraction(-2)),
        "B": (Fraction(-2), Fraction(2)),
    }


def test_Bm_one_mixed_orbit_spans_but_does_not_bound_tail():
    for length in (1, 2, 5, 17):
        projection = list(range(1, length + 1))
        # A sign flip gives twice one coordinate; permutations move it to all
        # coordinates. These axis vectors are visibly independent.
        axes = []
        for index in range(length):
            flipped = projection[:]
            flipped[index] *= -1
            difference = [
                left - right
                for left, right in zip(projection, flipped, strict=True)
            ]
            assert difference[index] != 0
            assert sum(value != 0 for value in difference) == 1
            axes.append(difference)
        assert len(axes) == length
        assert all(axes[i][i] != 0 for i in range(length))


def test_synthesis_scope_phrases_and_clean_bytes():
    raw = SYNTHESIS.read_bytes()
    audit_raw = AUDIT.read_bytes()
    for blob in (raw, audit_raw):
        assert b"\r" not in blob
        assert not any(byte < 32 and byte not in (9, 10) for byte in blob)
    text = raw.decode("utf-8")
    for needle in (
        "negative-definite packet admitting the declared normal projective",
        "separate semiampleness/abundance gate",
        "minimal finite-angle hostile incidences",
        "if it is also \\(\\mathbf Q\\)-Gorenstein",
        "supplied \\(\\mathbf Q\\)-Gorenstein contraction",
        "orbit of any one nonzero exceptional projection",
        "does **not** bound their number",
        "dimension-bound existence theorem",
    ):
        assert needle in text
    assert "entire minimal hostile boundary has been reduced to one" not in text
