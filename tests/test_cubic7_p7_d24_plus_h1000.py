import json
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_bna1 import _functional_value, _row_rank
from coxeter_polytopes.search.cubic7_d24_plus_coefficients import (
    DELTA7,
    P7_D24_SEED,
    _reflection_integral,
    d24_coefficient_pair,
    d24_height_zero_simple_roots_orthonormal,
)
from coxeter_polytopes.search.cubic7_dn import d24_plus_membership


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/cubic7_p7_d24_plus_h1000.json"


def _element(payload: dict[str, str]) -> Cubic7:
    return Cubic7.from_json(payload)


def _vector(payload: list[dict[str, str]]) -> tuple[Cubic7, ...]:
    return tuple(_element(value) for value in payload)


def _functional(wall: tuple[Cubic7, ...]) -> tuple[Cubic7, ...]:
    return (-P7_D24_SEED.timelike_coefficient * wall[0], *wall[1:])


def _payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_h1000_census_and_prefix_scope_are_frozen() -> None:
    payload = _payload()
    assert payload["seed_id"] == P7_D24_SEED.seed_id
    assert payload["height_bound"] == 1000
    assert payload["control_level_count"] == payload["height_packet_count"] == 94
    assert payload["nonempty_packet_count"] == 3
    assert payload["accepted_root_count"] == 155
    assert payload["spinor_accepted_count"] == 10
    assert [
        (packet["height"]["display"], packet["candidate_count"], packet["accepted_count"])
        for packet in payload["height_packets"]
        if packet["accepted_count"]
    ] == [
        ("7 + 22*theta + 10*theta^2", 1, 1),
        ("78 + 220*theta + 98*theta^2", 77, 77),
        ("122 + 342*theta + 152*theta^2", 77, 77),
    ]
    nonempty = [packet for packet in payload["height_packets"] if packet["accepted_count"]]
    assert [packet["accepted_cosets"] for packet in nonempty] == [
        {"D": 1},
        {"D": 72, "spinor-parity-1": 4, "spinor-parity-2": 1},
        {"D": 72, "spinor-parity-1": 4, "spinor-parity-7": 1},
    ]
    assert payload["frontier"]["all_controls_at_height_at_most_1000_exhausted"]
    assert payload["discovery"] is None
    assert "does not imply global nonreflectivity" in payload["scope"]
    assert not payload["geometry_conclusion"]["full_double_description_claimed"]


def test_h1000_archived_controls_and_roots_are_exact() -> None:
    payload = _payload()
    coefficient = P7_D24_SEED.timelike_coefficient
    canonical_norms = set(P7_D24_SEED.canonical_mirror_norms)
    controls = set()
    for row in payload["control_levels"]:
        height = _element(row["height"]["exact"])
        norm = _element(row["norm"]["exact"])
        control = _element(row["control"]["exact"])
        target = _element(row["spatial_square_target"]["exact"])
        assert norm in canonical_norms
        assert height == coefficient * control * control / norm
        assert target == norm + coefficient * control * control
        assert (2 * coefficient * control / norm).is_integral
        assert height.sign(0) > 0 and (Cubic7(1000) - height).sign(0) >= 0
        controls.add((height, norm, control, target))
    assert len(controls) == 94

    roots = []
    for row in payload["accepted_roots"]:
        vector = _vector(row["vector_orthonormal"])
        norm = _element(row["norm"]["exact"])
        control = _element(row["control"]["exact"])
        target = _element(row["spatial_square_target"]["exact"])
        assert len(vector) == 25 and vector[0] == control
        assert d24_coefficient_pair(P7_D24_SEED, vector, vector) == norm
        assert sum((value * value for value in vector[1:]), start=Cubic7()) == target
        membership = d24_plus_membership(vector[1:])
        assert membership == (
            "spinor" if row["coset"].startswith("spinor-parity-") else row["coset"]
        )
        assert _reflection_integral(P7_D24_SEED, norm, control, vector[1:])
        roots.append(vector)
    assert len(roots) == 155

    earlier = list(d24_height_zero_simple_roots_orthonormal())
    offset = 0
    for packet in payload["same_height_audit"]:
        count = packet["accepted_count"]
        current = roots[offset : offset + count]
        assert all(
            d24_coefficient_pair(P7_D24_SEED, root, wall).sign(0) <= 0
            for root in current
            for wall in earlier
        )
        positive = sum(
            d24_coefficient_pair(P7_D24_SEED, left, right).sign(0) > 0
            for index, left in enumerate(current)
            for right in current[index + 1 :]
        )
        assert positive == packet["positive_pairing_count"] == 0
        earlier_values = tuple(
            d24_coefficient_pair(P7_D24_SEED, root, wall)
            for root in current
            for wall in earlier
        )
        assert packet["earlier_wall_pairing_census"] == {
            "negative": sum(value.sign(0) < 0 for value in earlier_values),
            "zero": sum(not value for value in earlier_values),
            "positive": 0,
        }
        assert packet["passed"]
        earlier.extend(current)
        offset += count
    assert offset == 155


def test_h300_rays_are_checked_against_both_new_packets_exactly() -> None:
    payload = _payload()
    roots = [_vector(row["vector_orthonormal"]) for row in payload["accepted_roots"]]
    packet2 = roots[1:78]
    packet3 = roots[78:155]
    audit = payload["old_h300_hyperideal_ray_audit"]
    assert audit["old_ray_count"] == len(audit["rays"]) == 23
    survives2 = survives3 = 0
    for row in audit["rays"]:
        ray = _vector(row["vector"])
        values2 = tuple(d24_coefficient_pair(P7_D24_SEED, ray, wall) for wall in packet2)
        values3 = tuple(d24_coefficient_pair(P7_D24_SEED, ray, wall) for wall in packet3)
        assert sum(value.sign(0) > 0 for value in values2) == row[
            "packet_2_positive_pairing_count"
        ]
        assert sum(not value for value in values2) == row["packet_2_zero_pairing_count"]
        assert max(values2) == _element(row["packet_2_max_pairing"]["exact"])
        assert sum(value.sign(0) > 0 for value in values3) == row[
            "packet_3_positive_pairing_count"
        ]
        assert sum(not value for value in values3) == row["packet_3_zero_pairing_count"]
        assert max(values3) == _element(row["packet_3_max_pairing"]["exact"])
        survives_packet2 = all(value.sign(0) <= 0 for value in values2)
        survives_both = survives_packet2 and all(value.sign(0) <= 0 for value in values3)
        assert survives_packet2 == row["survives_packet_2"]
        assert survives_both == row["survives_both_new_packets"]
        survives2 += survives_packet2
        survives3 += survives_both
    assert survives2 == audit["survives_packet_2_count"]
    assert survives3 == audit["survives_both_new_packets_count"]


def test_h1000_hyperideal_ray_certificate_is_entirely_exact() -> None:
    payload = _payload()
    certificate = payload["certified_hyperideal_extreme_ray"]
    witness = _vector(certificate["vector"])
    roots = [_vector(row["vector_orthonormal"]) for row in payload["accepted_roots"]]
    walls = (*d24_height_zero_simple_roots_orthonormal(), *roots)
    functionals = tuple(_functional(wall) for wall in walls)
    pairings = tuple(_functional_value(functional, witness) for functional in functionals)
    active = tuple(index for index, value in enumerate(pairings) if not value)
    norm = d24_coefficient_pair(P7_D24_SEED, witness, witness)
    assert len(walls) == certificate["wall_count"] == 179
    assert all(value.sign(0) <= 0 for value in pairings)
    assert active == tuple(certificate["active_facets"])
    assert _row_rank(tuple(functionals[index] for index in active), 25) == 24
    assert certificate["active_rank"] == 24
    assert norm == _element(certificate["norm"]["exact"])
    assert norm.sign(0) > 0
    assert certificate["wall_pairing_census"] == {
        "negative": sum(value.sign(0) < 0 for value in pairings),
        "zero": sum(not value for value in pairings),
        "positive": 0,
    }
    assert certificate["certificate_is_entirely_exact"]
    assert not payload["geometry_conclusion"]["finite_volume"]
    assert not payload["geometry_conclusion"]["compact"]


def test_complete_norm_gate_is_unchanged() -> None:
    payload = _payload()
    assert P7_D24_SEED.canonical_mirror_norms == (
        Cubic7(1),
        Cubic7(2),
        DELTA7,
        2 * DELTA7,
    )
    assert payload["canonical_mirror_norms"] == [
        {
            "display": str(norm),
            "exact": norm.to_json(),
            "signs": [norm.sign(place) for place in range(3)],
        }
        for norm in P7_D24_SEED.canonical_mirror_norms
    ]
    assert "complete canonical norms" in payload["complete_mirror_norm_proof"]
