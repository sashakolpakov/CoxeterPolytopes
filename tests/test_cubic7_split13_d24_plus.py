import json
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_bna1 import _functional_value, _row_rank
from coxeter_polytopes.search.cubic7_d24_plus_coefficients import (
    _reflection_integral,
    d24_coefficient_pair,
)
from coxeter_polytopes.search.cubic7_d24_plus_split13 import (
    P13_D24_SEEDS,
    P13_R7,
    P13_R10,
    P13_R8,
    enumerate_split13_control_levels,
    enumerate_split13_level,
    split13_height_zero_patterns,
    split13_prime,
)
from coxeter_polytopes.search.cubic7_d24_plus_unit import (
    d24_height_zero_simple_roots_orthonormal,
)
from coxeter_polytopes.search.cubic7_dn import d24_plus_membership
from coxeter_polytopes.search.cubic7_seed_manifest import (
    initial_cubic7_diagonal_seeds,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/cubic7_split13_d24_plus_h300.json"


def _element(payload: dict[str, str]) -> Cubic7:
    return Cubic7.from_json(payload)


def _vector(payload: list[dict[str, str]]) -> tuple[Cubic7, ...]:
    return tuple(_element(value) for value in payload)


def _payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def _functional(seed, wall):
    return (-seed.timelike_coefficient * wall[0], *wall[1:])


def test_split13_d24_seeds_match_manifest_and_remain_distinguished() -> None:
    manifest = initial_cubic7_diagonal_seeds()[3:6]
    primes = (P13_R7, P13_R10, P13_R8)
    assert tuple(seed.coefficient for seed in manifest) == primes
    assert tuple(split13_prime(seed) for seed in P13_D24_SEEDS) == primes
    assert [seed.support for seed in P13_D24_SEEDS] == [
        "p13(theta=7)",
        "p13(theta=10)",
        "p13(theta=8)",
    ]
    assert len(set(P13_D24_SEEDS)) == 3
    for seed, prime in zip(P13_D24_SEEDS, primes, strict=True):
        assert prime.norm() == 13
        assert all(prime.sign(place) > 0 for place in range(3))
        assert seed.admissible
        assert seed.canonical_mirror_norms == (
            Cubic7(1),
            Cubic7(2),
            prime,
            2 * prime,
        )
        assert "complete norm ideals" in seed.complete_norm_class_proof


def test_height_zero_is_exactly_the_classical_d24_orbit() -> None:
    for seed in P13_D24_SEEDS:
        patterns = split13_height_zero_patterns(seed)
        assert len(patterns) == 1
        norm, coset, spatial = patterns[0]
        assert norm == 2 and coset == "D"
        assert spatial[:2] == (Cubic7(1), Cubic7(1))
        assert not any(spatial[2:])


def test_h300_control_census_and_frontier_are_exact() -> None:
    payload = _payload()
    assert payload["height_bound"] == 300
    assert payload["seed_order"] == [seed.seed_id for seed in P13_D24_SEEDS]
    records = {record["seed_id"]: record for record in payload["seeds"]}
    expected_controls = {
        P13_D24_SEEDS[0].seed_id: (46, [2, 7, 9, 28]),
        P13_D24_SEEDS[1].seed_id: (46, [3, 6, 9, 28]),
        P13_D24_SEEDS[2].seed_id: (47, [3, 7, 9, 28]),
    }
    for seed in P13_D24_SEEDS:
        record = records[seed.seed_id]
        levels = enumerate_split13_control_levels(seed, 300)
        count, norm_counts = expected_controls[seed.seed_id]
        assert len(levels) == record["control_level_count"] == count
        assert len({level.height for level in levels}) == len(levels)
        assert [row["count"] for row in record["control_counts_by_norm"]] == norm_counts
        archived = {
            (
                _element(row["height"]["exact"]),
                _element(row["norm"]["exact"]),
                _element(row["control"]["exact"]),
                _element(row["spatial_square_target"]["exact"]),
            )
            for row in record["control_levels"]
        }
        assert archived == {
            (level.height, level.norm, level.control, level.target) for level in levels
        }
        assert record["frontier"]["all_controls_at_height_at_most_300_exhausted"]
    assert payload["comparison"]["finite_volume_survivors"] == []
    assert payload["comparison"]["compact_survivors"] == []
    assert "No global nonreflectivity conclusion" in payload["scope"]


def test_nonempty_packets_and_p_divisible_gate_replay_exactly() -> None:
    payload = _payload()
    records = {record["seed_id"]: record for record in payload["seeds"]}
    expected_packets = {
        P13_D24_SEEDS[0].seed_id: [
            ("19 + 54*theta + 24*theta^2", 1, 1),
            ("38 + 108*theta + 48*theta^2", 5, 4),
        ],
        P13_D24_SEEDS[1].seed_id: [
            ("33 + 94*theta + 42*theta^2", 1, 1),
        ],
        P13_D24_SEEDS[2].seed_id: [
            ("6 + 18*theta + 8*theta^2", 1, 1),
        ],
    }
    expected_p_roots = {
        P13_D24_SEEDS[0].seed_id: 0,
        P13_D24_SEEDS[1].seed_id: 1,
        P13_D24_SEEDS[2].seed_id: 1,
    }
    for seed in P13_D24_SEEDS:
        record = records[seed.seed_id]
        packets = [packet for packet in record["height_packets"] if packet["accepted_count"]]
        assert [
            (packet["height"]["display"], packet["candidate_count"], packet["accepted_count"])
            for packet in packets
        ] == expected_packets[seed.seed_id]
        levels = enumerate_split13_control_levels(seed, 300)
        for packet in packets:
            height = _element(packet["height"]["exact"])
            candidates = tuple(
                root
                for level in levels
                if level.height == height
                for root in enumerate_split13_level(seed, level)
            )
            assert len(candidates) == packet["candidate_count"]

        prime = split13_prime(seed)
        roots = [_vector(row["vector_orthonormal"]) for row in record["accepted_roots"]]
        norms = [_element(row["norm"]["exact"]) for row in record["accepted_roots"]]
        assert sum(norm in (prime, 2 * prime) for norm in norms) == expected_p_roots[
            seed.seed_id
        ]
        for row, root, norm in zip(record["accepted_roots"], roots, norms, strict=True):
            assert d24_coefficient_pair(seed, root, root) == norm
            assert _reflection_integral(seed, norm, root[0], root[1:])
            membership = d24_plus_membership(root[1:])
            assert membership == (
                "spinor" if row["coset"].startswith("spinor-parity-") else row["coset"]
            )
            if norm in (prime, 2 * prime):
                reduced = tuple(value / prime for value in root[1:])
                assert all(value.is_integral for value in reduced)
                assert d24_plus_membership(reduced) in ("D", "spinor")


def test_simultaneous_height_selection_is_acute() -> None:
    payload = _payload()
    records = {record["seed_id"]: record for record in payload["seeds"]}
    for seed in P13_D24_SEEDS:
        record = records[seed.seed_id]
        roots = [_vector(row["vector_orthonormal"]) for row in record["accepted_roots"]]
        earlier = list(d24_height_zero_simple_roots_orthonormal())
        cursor = 0
        for packet, audit in zip(
            (row for row in record["height_packets"] if row["accepted_count"]),
            record["same_height_audit"],
            strict=True,
        ):
            current = roots[cursor : cursor + packet["accepted_count"]]
            within = tuple(
                d24_coefficient_pair(seed, left, right)
                for index, left in enumerate(current)
                for right in current[index + 1 :]
            )
            prior = tuple(
                d24_coefficient_pair(seed, root, wall)
                for root in current
                for wall in earlier
            )
            assert audit["within_packet"] == {
                "negative": sum(value.sign(0) < 0 for value in within),
                "zero": sum(not value for value in within),
                "positive": 0,
            }
            assert audit["against_earlier_walls"] == {
                "negative": sum(value.sign(0) < 0 for value in prior),
                "zero": sum(not value for value in prior),
                "positive": 0,
            }
            assert all(value.sign(0) <= 0 for value in (*within, *prior))
            earlier.extend(current)
            cursor += len(current)
        assert cursor == len(roots)


def test_exact_cones_and_hyperideal_witnesses_are_frozen() -> None:
    payload = _payload()
    records = {record["seed_id"]: record for record in payload["seeds"]}
    expected = {
        P13_D24_SEEDS[0].seed_id: [(25, 25, 1, 0, 24), (29, 895, 3, 0, 892)],
        P13_D24_SEEDS[1].seed_id: [(25, 25, 2, 0, 23)],
        P13_D24_SEEDS[2].seed_id: [(25, 25, 2, 0, 23)],
    }
    for seed in P13_D24_SEEDS:
        record = records[seed.seed_id]
        snapshots = record["geometry"]["snapshots"]
        assert [
            (
                row["wall_count"],
                row["extreme_ray_count"],
                row["finite_vertex_count"],
                row["cusp_count"],
                row["hyperideal_ray_count"],
            )
            for row in snapshots
        ] == expected[seed.seed_id]
        final = record["geometry"]["final"]
        assert not final["finite_volume"] and not final["compact"]

        roots = [_vector(row["vector_orthonormal"]) for row in record["accepted_roots"]]
        walls = (*d24_height_zero_simple_roots_orthonormal(), *roots)
        functionals = tuple(_functional(seed, wall) for wall in walls)
        witness = record["geometry"]["hyperideal_witness"]
        vector = _vector(witness["vector"])
        pairings = tuple(_functional_value(functional, vector) for functional in functionals)
        active = tuple(index for index, value in enumerate(pairings) if not value)
        norm = d24_coefficient_pair(seed, vector, vector)
        assert all(value.sign(0) <= 0 for value in pairings)
        assert active == tuple(witness["active_facets"])
        assert _row_rank(tuple(functionals[index] for index in active), 25) == 24
        assert witness["active_rank"] == 24
        assert norm == _element(witness["norm"]["exact"])
        assert norm.sign(0) > 0
        assert witness["wall_pairing_census"] == {
            "negative": sum(value.sign(0) < 0 for value in pairings),
            "zero": sum(not value for value in pairings),
            "positive": 0,
        }
