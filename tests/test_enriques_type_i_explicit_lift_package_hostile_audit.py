from __future__ import annotations

from functools import lru_cache
import hashlib
import json

from scripts.audit_enriques_type_i_explicit_lift_package_hostile import (
    OUTPUT,
    analyze,
)


ARCHIVE_SHA256 = "9a79f9036c8383c4bee67ce71cae36008c34d1507ee94a49ee3adfe5d32e7ee6"


@lru_cache(maxsize=1)
def _payload() -> dict:
    return analyze()


def test_hostile_archive_replays_and_freezes() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert hashlib.sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    assert payload["verdict"] == "PASS/FREEZE"
    assert payload["all_checks_pass"] is True


def test_row_coordinate_involution_and_saturated_rootless_anti_lattice() -> None:
    convention = _payload()["row_column_convention"]
    assert convention["row_isometry_identity_tau_G_tauT"] is True
    assert convention["column_isometry_identity_tauT_G_tau"] is False
    assert convention["column_convention_maximum_error"] == 152
    anti = _payload()["anti_lattice"]
    assert anti["rank"] == 9
    assert anti["determinant"] == -1024
    assert anti["smith_diagonal"] == [2] * 8 + [4]
    assert anti["basis_primitive_in_SX_smith_diagonal"] == [1] * 9
    assert anti["selected_difference_span_smith_in_anti_basis"] == [1] * 9
    assert anti["all_norms_divisible_by_4"] is True
    assert anti["minus_two_root_count"] == 0


def test_independent_reflection_and_LDL_exhaustions_match_production_sets() -> None:
    fiber = _payload()["independent_fiber_replay"]
    assert fiber["engine"] == "simple-reflection closure, not qfminim"
    assert fiber["root_count"] == 242
    assert fiber["all_lift_intersection_histogram"] == {"0": 246, "1": 28, "2": 2}
    l4 = _payload()["independent_L4_replay"]
    assert l4["engine"] == "independent exact rational LDL ellipsoid recursion"
    assert l4["bound"] == 378
    assert l4["including_zero_leaf_count"] == 473
    assert l4["nonzero_vector_count"] == 472
    assert l4["degree_zero_root_count"] == 0
    assert l4["negative_component_witness_count"] == 0
    assert l4["positive_root_degree_histogram"] == {
        "1": 20,
        "2": 20,
        "3": 24,
        "4": 28,
        "5": 36,
        "6": 44,
        "7": 64,
    }


def test_surface_existence_is_literature_backed_and_scope_is_not_overclaimed() -> None:
    scope = _payload()["geometric_provenance_and_scope"]
    assert scope["serialized_datum_kind"] == "lattice/Borcherds context, not surface equations"
    assert scope["oscar_1_4_loader_tuple_order"] == "SY, SX, L26, w, u"
    assert scope["oscar_fourth_entry_is_weyl_vector"] is True
    assert scope["machine_certificate_alone_constructs_surface_equations"] is False
    assert scope["actual_surface_identification_supplies_the_bridge_from_SX_in_II_1_25_to_NS_of_a_K3_cover"] is True
    assert scope["new_abstract_polytope_claimed"] is False
    assert scope["new_surface_claimed"] is False
