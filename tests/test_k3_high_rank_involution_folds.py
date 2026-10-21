import hashlib
import json

from scripts.search_k3_high_rank_involution_folds import (
    ARCHIVE_PATH,
    SOURCE_PATH,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(ARCHIVE_PATH.read_text())


def test_exact_weighted_graph_involution_inventory_replays_quickly() -> None:
    archived = _archive()
    replay = build_certificate(inventory_only=True)
    assert archived["source_sha256"] == hashlib.sha256(SOURCE_PATH.read_bytes()).hexdigest()
    assert replay["inventory_sha256"] == archived["inventory_sha256"]
    assert replay["inventory_sha256"] == (
        "671358ce1466142dccbad743b76735c3a9ed1558c3fca47cfa123a06a9b722af"
    )
    assert replay["summary"]["weighted_graph_automorphism_counts"] == {
        "15": 12,
        "16": 4,
        "17": 2,
    }
    assert replay["summary"]["nonidentity_involution_counts"] == {
        "15": 7,
        "16": 3,
        "17": 1,
    }
    assert replay["summary"]["dimension_greater_than_nine_involution_counts"] == {
        "15": 3,
        "16": 1,
        "17": 0,
    }
    assert replay["all_checks_pass"] is True


def test_frozen_positive_folds_pass_every_exact_geometry_gate() -> None:
    archived = _archive()
    assert archived["artifact_status"] == "exact_search_and_geometry_certified"
    assert archived["all_checks_pass"] is True
    geometries = [
        geometry
        for case in archived["cases"]
        for geometry in case["geometry_certificates"]
    ]
    assert [
        (
            geometry["ambient_rho"],
            geometry["hyperbolic_dimension"],
            geometry["facet_count"],
            geometry["finite_vertex_count"],
            geometry["ideal_vertex_count"],
        )
        for geometry in geometries
    ] == [
        (15, 13, 19, 178, 7),
        (16, 14, 17, 89, 5),
    ]
    assert all(geometry["all_checks_pass"] for geometry in geometries)
    assert all(geometry["finite_volume_certified"] for geometry in geometries)
    assert all(geometry["compact"] is False for geometry in geometries)
    assert all(all(geometry["checks"].values()) for geometry in geometries)
    assert all(
        geometry["neron_severi_action"]["all_integral_lattice_checks_pass"]
        for geometry in geometries
    )
    assert all(
        not geometry["neron_severi_action"][
            "discriminant_action_is_scalar_plus_or_minus_one"
        ]
        for geometry in geometries
    )

    rank15, rank16 = geometries
    assert rank15["coxeter"]["edge_counts"] == {
        "2": 146,
        "3": 14,
        "4": 2,
        "6": 0,
        "infinity": 6,
        "ultraparallel": 3,
    }
    assert rank16["coxeter"]["edge_counts"] == {
        "2": 119,
        "3": 14,
        "4": 2,
        "6": 0,
        "infinity": 1,
        "ultraparallel": 0,
    }
    assert rank15["vinberg_volume_criterion"]["edge_extension_histogram"] == {
        "2": 1323
    }
    assert rank16["vinberg_volume_criterion"]["edge_extension_histogram"] == {
        "2": 704
    }
    assert rank15["vinberg_volume_criterion"]["canonical_full_audit_sha256"] == (
        "f08875e243a45d110620af4b4d86cb33c5b417a2b286e524ecdc1d9e1015bcff"
    )
    assert rank16["vinberg_volume_criterion"]["canonical_full_audit_sha256"] == (
        "8a9ae49a6c2fd562b336fc0203726fd02436332eaf25e2a74a094f9571543619"
    )


def test_positive_folds_are_graph_symmetries_not_surface_automorphisms() -> None:
    archived = _archive()
    geometries = [
        geometry
        for case in archived["cases"]
        for geometry in case["geometry_certificates"]
    ]
    assert all(
        geometry["surface_automorphism_realization"]["status"]
        == "not_realized_on_general_K3_special_period_unassessed"
        for geometry in geometries
    )
    assert archived["summary"][
        "all_certified_positive_folds_are_surface_automorphism_fixed"
    ] is False
    assert archived["summary"]["novelty_status"] == (
        "unresolved; these exact folded polytopes are not claimed new"
    )


def test_every_conjugate_fold_has_an_exact_orbit_gram_transport() -> None:
    archived = _archive()
    assert all(
        transport["orbit_gram_congruence_verified"]
        for case in archived["cases"]
        for conjugacy_class in case["conjugacy_classes"]
        for transport in conjugacy_class["member_transports"]
    )
    # Three rank-15 involutions are conjugate and give the same H^13 fold;
    # the rank-16 H^14 involution is a singleton conjugacy class.
    positive_classes = [
        conjugacy_class
        for case in archived["cases"]
        for conjugacy_class in case["conjugacy_classes"]
        if conjugacy_class["eligible_for_full_geometry_gate"]
    ]
    assert [conjugacy_class["size"] for conjugacy_class in positive_classes] == [3, 1]
    assert archived["summary"]["certified_finite_volume_involution_count"] == 4
