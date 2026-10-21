from pathlib import Path
import shutil
import subprocess

from scripts.certify_rational_i1n_inverse_design import (
    COXITER_SHA256,
    DIMENSIONS,
    EXPECTED_FACETS,
    arithmetic_genus,
    build_certificate,
    canonical_pairing,
    coxiter_filename,
    enumerate_simple_roots,
    flex_obstruction_record,
    locate_named_roots,
    lorentz_pair,
    named_classes,
    parse_coxiter,
    p12_p13_decomposition,
    p16_cone_certificate,
    p16_decomposition,
    p16_extreme_rays,
    p16_inequality_rows,
    q_tail_cone_certificate,
    q_tail_decomposition,
    q_tail_extreme_rays,
    q_tail_inequality_rows,
    relation_edges,
)


COXITER_ROOT = Path("/private/tmp/coxiter-p118/graphs/group units qf")
SINGULAR_CERTIFICATE = Path("scripts/certify_rational_i1n_curves.sing")


def test_streamed_roots_have_the_expected_norms_and_adjunction_data():
    for n in DIMENSIONS:
        roots = enumerate_simple_roots(n)
        assert len(roots) == EXPECTED_FACETS[n]
        assert all(lorentz_pair(root, root) in (1, 2) for root in roots)
        assert all(arithmetic_genus(root) >= 0 for root in roots)
        for root in roots[n:]:
            multiplicities = tuple(-entry for entry in root[1:])
            assert all(
                multiplicities[i] >= multiplicities[i + 1]
                for i in range(n - 1)
            )
            expected = -3 * root[0] + sum(multiplicities)
            assert canonical_pairing(root) == expected


def test_every_pinned_coxiter_relation_graph_matches_exactly_when_available():
    if not COXITER_ROOT.exists():
        return
    import hashlib

    for n in DIMENSIONS:
        path = COXITER_ROOT / coxiter_filename(n)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == COXITER_SHA256[n]
        facets, dimension, published = parse_coxiter(path)
        reconstructed = relation_edges(enumerate_simple_roots(n))
        assert (facets, dimension) == (EXPECTED_FACETS[n], n)
        assert len(published) == len(reconstructed)
        assert set(published) == set(reconstructed)


def test_explicit_cluster_realizes_all_named_walls_through_n16():
    expected = {
        12: {"L", "C"},
        13: {"L", "C"},
        14: {"L", "C", "Q"},
        15: {"L", "C", "Q"},
        16: {"L", "C", "Q", "S"},
    }
    for n, names in expected.items():
        labels = locate_named_roots(enumerate_simple_roots(n), n)
        assert {name for name, label in labels.items() if label is not None} == names

    classes = named_classes(16)
    assert -lorentz_pair(classes["L"], classes["C"]) == 0
    assert -lorentz_pair(classes["C"], classes["Q"]) == 0
    assert -lorentz_pair(classes["C"], classes["S"]) == 0
    assert -lorentz_pair(classes["Q"], classes["S"]) == 0
    assert [lorentz_pair(classes[name], classes[name]) for name in ("L", "C", "Q", "S")] == [2, 2, 2, 1]


def test_explicit_forms_are_absolutely_irreducible_when_singular_is_available():
    executable = shutil.which("singular")
    if executable is None:
        return
    completed = subprocess.run(
        [executable, "-q", str(SINGULAR_CERTIFICATE)],
        check=True,
        capture_output=True,
        text=True,
    )
    for name in ("C", "Q", "S"):
        assert f"{name}_NONCONSTANT_FACTORS=1" in completed.stdout
        assert f"{name}_FACTOR_EXPONENT=1" in completed.stdout


def test_flex_cubic_obstruction_starts_exactly_at_n17():
    for n in range(12, 17):
        assert flex_obstruction_record(n, enumerate_simple_roots(n)) is None
    for n in range(17, 20):
        roots = enumerate_simple_roots(n)
        obstruction = flex_obstruction_record(n, roots)
        assert obstruction is not None
        assert all(obstruction["wall_labels"].values())
        assert obstruction["intersection_numbers"] == {
            "L.C": 0,
            "C.Q": 0,
            "C.A": 1,
            "Q.A": 0,
        }


def test_p12_p13_effective_cone_formula_is_exact():
    samples = {
        12: [(4, (2, *([1] * 11))), (6, (2, 2, *([1] * 10)))],
        13: [(4, (2, *([1] * 12))), (9, (3, 3, 3, *([1] * 10)))],
    }
    for n, classes in samples.items():
        for degree, multiplicities in classes:
            x, y, exceptional = p12_p13_decomposition(n, degree, multiplicities)
            assert min((x, y, *exceptional)) >= 0
            assert x + 3 * y == degree
            base = [
                (x if i < 3 else 0) + (y if i < 11 else 0)
                for i in range(n)
            ]
            previous = 0
            reconstructed = []
            for i, coefficient in enumerate(exceptional):
                reconstructed.append(base[i] - (coefficient - previous))
                previous = coefficient
            assert tuple(reconstructed) == multiplicities


def test_p14_p15_exact_extreme_rays_all_decompose():
    for n, expected in ((14, 94), (15, 145)):
        rows = q_tail_inequality_rows(n)
        rays = q_tail_extreme_rays(n)
        assert len(rays) == expected
        for ray in rays:
            assert all(sum(a * b for a, b in zip(row, ray)) >= 0 for row in rows)
            x, y, w, exceptional = q_tail_decomposition(n, ray[0], ray[1:])
            assert min((x, y, w, *exceptional)) >= 0
            assert x + 3 * y + 4 * w == ray[0]
            base = [
                (x if i < 3 else 0)
                + (y if i < 11 else 0)
                + w
                + (w if i == 0 else 0)
                for i in range(n)
            ]
            previous = 0
            reconstructed = []
            for i, coefficient in enumerate(exceptional):
                reconstructed.append(base[i] - (coefficient - previous))
                previous = coefficient
            assert tuple(reconstructed) == ray[1:]
        certificate = q_tail_cone_certificate(n)
        assert certificate["primitive_extreme_ray_count"] == expected
        assert certificate[
            "every_ray_has_exact_nonnegative_L_C_Q_exceptional_decomposition"
        ]


def test_p16_exact_extreme_rays_all_decompose():
    rows = p16_inequality_rows()
    rays = p16_extreme_rays()
    assert len(rays) == 325
    for ray in rays:
        assert all(sum(a * b for a, b in zip(row, ray)) >= 0 for row in rows)
        x, y, w, v, exceptional = p16_decomposition(ray[0], ray[1:])
        assert min((x, y, w, v, *exceptional)) >= 0
        assert x + 3 * y + 4 * w + 6 * v == ray[0]
    certificate = p16_cone_certificate()
    assert certificate["primitive_extreme_ray_count"] == 325
    assert certificate[
        "every_ray_has_exact_nonnegative_L_C_Q_S_exceptional_decomposition"
    ]


def test_certificate_separates_wall_realization_from_nef_completeness():
    certificate = build_certificate(COXITER_ROOT if COXITER_ROOT.exists() else None)
    for n in range(12, 17):
        status = certificate["dimensions"][str(n)]["status"]
        assert status["wall_realization"].startswith("explicit")
        assert status["nef_cone_completeness"].startswith("proved")
        volume = certificate["dimensions"][str(n)]["vinberg_volume"]
        assert volume["finite_volume"] and not volume["compact"]
    for n in range(17, 20):
        status = certificate["dimensions"][str(n)]["status"]
        assert status["wall_realization"].startswith("impossible")
        obstruction = certificate["dimensions"][str(n)][
            "flex_cubic_obstruction"
        ]
        upgrade = obstruction["marking_independent_mmp_upgrade"]
        expected_orbit = {17: 8, 18: 24, 19: 120}[n]
        assert upgrade["canonical_class_count"] == expected_orbit
        assert upgrade["rank_one_P2_endpoints"] == expected_orbit
        assert upgrade["rank_two_endpoints"] == {
            "F1_nonminimal": expected_orbit,
            "F2_minimal": expected_orbit,
        }
        assert upgrade["weighted_diagram_automorphism_order"] == expected_orbit
        assert obstruction["scope"]["unconditional"].startswith(
            "nonexistence of a complete facet-cone carrier"
        )
    assert "abstract Vinberg finite-volume" in certificate["dimensions"]["17"][
        "status"
    ]["coxeter_polytope"]
    for n in (18, 19):
        assert "Vinberg--Kaplinskaja" in certificate["dimensions"][str(n)]["status"][
            "coxeter_polytope"
        ]
    assert certificate["global_conclusion"]["new_hyperbolic_polytope"] is False
    assert certificate["global_conclusion"]["positive_algebro_geometric_realizations"] == [12, 13, 14, 15, 16]
    assert certificate["global_conclusion"][
        "rational_carrier_nonexistence_dimensions"
    ] == [17, 18, 19]
