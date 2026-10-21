import json
from hashlib import sha256
from pathlib import Path

from coxeter_polytopes.exact.silver import SILVER_ONE, SQRT2, Silver
from coxeter_polytopes.search.silver_d24_plus import (
    _integer_determinant,
    d24_plus_gram,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_p7_d24_plus_h100.json"
H300_ARTIFACT = ROOT / "data/enumerations/silver_p7_d24_plus_h300.json"
H500_ARTIFACT = ROOT / "data/enumerations/silver_p7_d24_plus_h500.json"
from coxeter_polytopes.search.silver_p7_d24_plus import (
    EPSILON,
    P7_DELTA,
    P7_TIMELIKE,
    H500_REPLACEMENT_WITNESS_VECTOR,
    SILVER_P7_D24_SEED,
    analyze_silver_p7_d24_prefix,
    audit_first_packet_hyperideal_survival,
    audit_frozen_h100_witness,
    certify_full_prefix_hyperideal_ray,
    certify_h300_replacement_hyperideal_ray,
    certify_h500_replacement_hyperideal_ray,
    enumerate_p7_d24_control_levels,
    p7_factorization_audit,
    run_silver_p7_d24_prefix,
    silver_d24_height_zero_simple_roots,
    silver_p7_d24_pair,
    silver_p7_d24_vector_divisible_by_sqrt2,
)


def test_prime_seven_factorization_norm_classes_and_glue_basis() -> None:
    assert p7_factorization_audit()
    assert EPSILON * P7_TIMELIKE == P7_DELTA
    assert P7_TIMELIKE.norm() == -7 and P7_DELTA.norm() == 7
    assert SILVER_P7_D24_SEED.canonical_mirror_norms == (
        SILVER_ONE,
        Silver(2),
        P7_DELTA,
        2 * P7_DELTA,
    )
    gram = d24_plus_gram()
    assert _integer_determinant(gram) == 1
    assert all(entry.denominator == 1 for row in gram for entry in row)
    assert all(gram[index][index].numerator % 2 == 0 for index in range(24))


def test_h100_control_and_glue_census_is_complete() -> None:
    levels = enumerate_p7_d24_control_levels(100)
    search = run_silver_p7_d24_prefix(100)
    assert len(levels) == len(search.level_enumerations) == 25
    assert len(search.packets) == 17

    # A delta-norm mirror has spatial part in delta*D24+.  Every such target
    # is therefore divided by delta^2 before glue enumeration.  Through H=100
    # the integral reduced targets contain no D24+ pattern.
    prime_levels = [
        enumeration
        for enumeration in search.level_enumerations
        if enumeration.level.norm in (P7_DELTA, 2 * P7_DELTA)
    ]
    assert len(prime_levels) == 19
    assert all(not enumeration.candidates for enumeration in prime_levels)
    assert all(
        enumeration.reduced_target is None
        or enumeration.reduced_target.is_integral
        for enumeration in prime_levels
    )

    nonempty = [packet for packet in search.packets if packet.accepted_roots]
    assert [packet.height for packet in nonempty] == [
        Silver(5, 4),
        Silver(31, 22),
    ]
    assert [
        (
            packet.raw_reflective_count,
            packet.imprimitive_sqrt2_count,
            packet.candidate_count,
            len(packet.accepted_roots),
        )
        for packet in nonempty
    ] == [(4, 1, 3, 3), (30, 1, 29, 28)]
    assert len(search.accepted_roots) == 31
    assert sum(root.coset != "D24" for root in search.accepted_roots) == 5


def test_norm_two_imprimitive_mirrors_reduce_to_norm_one() -> None:
    search = run_silver_p7_d24_prefix(100)
    for height in (Silver(5, 4), Silver(31, 22)):
        enumerations = [
            enumeration
            for enumeration in search.level_enumerations
            if enumeration.level.height == height
        ]
        norm_one = next(
            enumeration.candidates[0]
            for enumeration in enumerations
            if enumeration.level.norm == 1
        )
        norm_two = next(
            enumeration
            for enumeration in enumerations
            if enumeration.level.norm == 2
        )
        assert len(norm_two.imprimitive_vectors) == 1
        raw = norm_two.imprimitive_vectors[0]
        assert silver_p7_d24_vector_divisible_by_sqrt2(raw)
        quotient = tuple(value / SQRT2 for value in raw)
        assert quotient == norm_one.vector
        assert silver_p7_d24_pair(raw, raw) == 2
        assert silver_p7_d24_pair(quotient, quotient) == 1
    assert all(
        not silver_p7_d24_vector_divisible_by_sqrt2(root.vector)
        for root in search.accepted_roots
    )


def test_first_complete_packet_exact_cone() -> None:
    search = run_silver_p7_d24_prefix(100)
    first = next(packet for packet in search.packets if packet.accepted_roots)
    assert [(root.norm, root.control, root.coset) for root in first.accepted_roots] == [
        (Silver(1), Silver(1, 1), "D24"),
        (Silver(2), Silver(2, 1), "D24"),
        (Silver(2), Silver(2, 1), "D24"),
    ]
    geometry = analyze_silver_p7_d24_prefix(search, 1)
    assert [
        (
            snapshot.wall_count,
            snapshot.extreme_ray_count,
            snapshot.finite_vertex_count,
            snapshot.cusp_count,
            snapshot.hyperideal_ray_count,
        )
        for snapshot in geometry.snapshots
    ] == [(27, 113, 4, 0, 109)]
    assert [
        (
            insertion.facet_index,
            insertion.negative_ray_count,
            insertion.zero_ray_count,
            insertion.positive_ray_count,
            insertion.output_ray_count,
        )
        for insertion in geometry.insertions
    ] == [(25, 3, 0, 22, 69), (26, 47, 0, 22, 113)]
    assert (
        geometry.ray_stream_sha256
        == "b5066aaa308ccfc7a48fb8c5897f4fe4613e4683364c5d0fd2def8d81158d24d"
    )
    assert not geometry.finite_volume and not geometry.compact


def test_h100_prefix_has_a_surviving_exact_hyperideal_extreme_ray() -> None:
    search = run_silver_p7_d24_prefix(100)
    witness = certify_full_prefix_hyperideal_ray(search)
    assert witness.later_wall_count == 28
    assert witness.surviving_hyperideal_ray_count == 2
    assert witness.final_active_rank == 24
    assert witness.new_zero_facets == (41, 46, 47, 48, 52)
    assert witness.ray.norm == Silver(48, -16) / 7
    assert witness.ray.vector == (
        Silver(20, 12) / 7,
        Silver(2, 2),
        Silver(2, 2),
        Silver(2),
        Silver(2),
        *(Silver(0) for _ in range(20)),
    )
    assert all(value.sign(0) <= 0 for value in witness.later_pairings)
    assert sum(not value for value in witness.later_pairings) == 5
    assert witness.ray.norm.sign(0) > 0


def test_frozen_h100_artifact_records_the_full_certificate() -> None:
    payload = json.loads(ARTIFACT.read_text())
    assert payload["height_bound"] == 100
    assert payload["enumeration"]["control_level_count"] == 25
    assert payload["enumeration"]["accepted_root_count"] == 31
    assert payload["enumeration"]["final_wall_count"] == 55
    snapshot = payload["first_complete_packet_geometry"]["snapshots"][0]
    assert (
        snapshot["wall_count"],
        snapshot["extreme_ray_count"],
        snapshot["finite_vertex_count"],
        snapshot["cusp_count"],
        snapshot["hyperideal_ray_count"],
    ) == (27, 113, 4, 0, 109)
    assert (
        payload["first_complete_packet_geometry"]["ray_stream_sha256"]
        == "b5066aaa308ccfc7a48fb8c5897f4fe4613e4683364c5d0fd2def8d81158d24d"
    )
    certificate = payload["full_h100_nonfinite_volume_certificate"]
    assert certificate["later_wall_count"] == 28
    assert certificate["all_later_pairings_nonpositive"]
    assert certificate["chosen_ray"]["final_active_rank"] == 24
    assert not certificate["finite_volume"] and not certificate["compact"]


def test_h300_wall_census_and_first_new_nonempty_height() -> None:
    search = run_silver_p7_d24_prefix(300)
    assert len(search.control_levels) == 48
    assert len(search.packets) == 33
    assert len(search.accepted_roots) == 324
    assert len(search.walls) == 348
    nonempty = [packet for packet in search.packets if packet.accepted_roots]
    assert [
        (
            packet.height,
            packet.raw_reflective_count,
            packet.imprimitive_sqrt2_count,
            packet.candidate_count,
            len(packet.accepted_roots),
        )
        for packet in nonempty
    ] == [
        (Silver(5, 4), 4, 1, 3, 3),
        (Silver(31, 22), 30, 1, 29, 28),
        (Silver(56, 40), 1, 0, 1, 1),
        (Silver(62, 44), 100, 0, 100, 87),
        (Silver(82, 58), 5, 0, 5, 5),
        (Silver(124, 88), 232, 0, 232, 200),
    ]
    first_new = next(
        packet
        for packet in nonempty
        if packet.height.sign(0) > 0
        and (packet.height - Silver(100)).sign(0) > 0
    )
    assert first_new.height == Silver(56, 40)
    assert len(first_new.accepted_roots) == 1
    root = first_new.accepted_roots[0]
    assert root.norm == 2 * P7_DELTA
    assert root.control == Silver(12, 8)
    assert root.coset == "D24"


def test_h300_old_witness_cut_and_initial_ray_survival_scope() -> None:
    search = run_silver_p7_d24_prefix(300)
    audit = audit_frozen_h100_witness(search)
    assert audit.first_cut_height == Silver(62, 44)
    assert len(audit.first_cut_packet_positive_indices) == 4
    assert [
        (
            packet.height,
            packet.wall_count,
            packet.negative_pairing_count,
            packet.zero_pairing_count,
            packet.positive_pairing_count,
            packet.minimum_pairing,
            packet.maximum_pairing,
        )
        for packet in audit.packets
    ] == [
        (Silver(31, 22), 28, 23, 5, 0, Silver(-12, -8), Silver(0)),
        (Silver(56, 40), 1, 1, 0, 0, Silver(-4, -2), Silver(-4, -2)),
        (Silver(62, 44), 87, 79, 4, 4, Silver(-28, -20), Silver(4)),
        (Silver(82, 58), 5, 4, 0, 1, Silver(-28, -18), Silver(4, 2)),
        (Silver(124, 88), 200, 183, 14, 3, Silver(-40, -28), 2 * SQRT2),
    ]


def test_h300_has_a_direct_exact_replacement_hyperideal_ray() -> None:
    search = run_silver_p7_d24_prefix(300)
    witness = certify_h300_replacement_hyperideal_ray(search)
    coordinate = Silver(172, -120) / 7
    assert witness.vector == (
        Silver(1),
        *(coordinate for _ in range(19)),
        *(Silver(0) for _ in range(5)),
    )
    assert witness.norm == Silver(1109345, -784418) / 49
    assert witness.norm.sign(0) > 0
    assert witness.wall_count == 348
    assert all(value.sign(0) <= 0 for value in witness.pairings)
    assert len(witness.active_facets) == 24
    assert witness.active_rank == 24
    assert witness.active_facets == (
        *range(18),
        *range(19, 24),
        59,
    )
    survival = audit_first_packet_hyperideal_survival(search)
    assert [
        (
            item.height,
            item.cumulative_later_wall_count,
            item.surviving_hyperideal_ray_count,
        )
        for item in survival
    ] == [
        (Silver(31, 22), 28, 2),
        (Silver(56, 40), 29, 2),
        (Silver(62, 44), 116, 0),
        (Silver(82, 58), 121, 0),
        (Silver(124, 88), 321, 0),
    ]


def test_frozen_h300_artifact_certifies_prefix_nonfinite_volume() -> None:
    payload = json.loads(H300_ARTIFACT.read_text())
    assert payload["height_bound"] == 300
    assert payload["enumeration"]["control_level_count"] == 48
    assert payload["enumeration"]["height_packet_count"] == 33
    assert payload["enumeration"]["accepted_root_count"] == 324
    assert payload["enumeration"]["final_wall_count"] == 348
    assert [
        item["accepted_count"]
        for item in payload["enumeration"]["nonempty_packet_summary"]
    ] == [3, 28, 1, 87, 5, 200]
    assert payload["h100_witness_extension_audit"]["first_cut_height"][
        "display"
    ] == "62 + 44*sqrt(2)"
    assert [
        item["surviving_hyperideal_ray_count"]
        for item in payload["first_packet_hyperideal_ray_survival"]
    ] == [2, 2, 0, 0, 0]
    geometry = payload["geometry_status"]
    assert geometry["h300_finite_volume"] is False
    assert geometry["h300_compact"] is False
    assert geometry["eventual_chamber_status"] == "unresolved"
    replacement = payload["h300_replacement_hyperideal_ray"]
    assert replacement["wall_count"] == 348
    assert replacement["all_wall_pairings_nonpositive"]
    assert replacement["zero_pairing_count"] == 24
    assert replacement["active_rank"] == 24
    assert replacement["norm"]["display"] == (
        "1109345/49 - 784418/49*sqrt(2)"
    )
    assert replacement["finite_volume"] is False
    assert replacement["compact"] is False
    assert payload["v1_h100_lineage"]["sha256"]["artifact"] == (
        "ec01a747418b4fb0ecbfa25040cbae6d3f6369746f533c5b5447ff9ac11c5fbd"
    )


def test_frozen_h500_extension_and_replacement_hyperideal_ray() -> None:
    payload = json.loads(H500_ARTIFACT.read_text())
    enumeration = payload["enumeration"]
    assert payload["height_bound"] == 500
    assert enumeration["control_level_count"] == 64
    assert enumeration["height_packet_count"] == 43
    assert enumeration["accepted_root_count"] == 4350
    assert enumeration["final_wall_count"] == 4374
    assert enumeration["extension_control_level_count"] == 16
    assert enumeration["extension_accepted_root_count"] == 4026
    assert [
        (
            item["height"]["display"],
            item["primitive_candidate_count"],
            item["accepted_count"],
        )
        for item in enumeration["new_nonempty_packet_summary"]
    ] == [
        ("181 + 128*sqrt(2)", 4090, 4003),
        ("206 + 146*sqrt(2)", 7, 6),
        ("248 + 176*sqrt(2)", 24, 17),
    ]

    roots = [
        root
        for packet in enumeration["extension_packets"]
        for root in packet["accepted_roots"]
    ]
    digest = sha256()
    for root in roots:
        digest.update(
            json.dumps(root, sort_keys=True, separators=(",", ":")).encode()
        )
        digest.update(b"\n")
    assert digest.hexdigest() == enumeration["extension_accepted_root_stream_sha256"]

    accelerator = payload["exact_dense_packet_accelerator"]
    assert accelerator["simultaneous_height_semantics"]
    assert accelerator["same_height_unordered_pair_count"] == 8_010_154
    assert (
        accelerator["maximum_proved_intermediate_absolute_bound"]
        < accelerator["int64_max"]
    )
    assert (
        accelerator["maximum_abs_pair_coefficient"]
        <= accelerator["exact_sign_square_safety_limit"]
    )
    assert payload["h300_replacement_ray_cut_audit"]["positive_pairing_count"] == 74

    h300 = json.loads(H300_ARTIFACT.read_text())
    old_walls = [
        tuple(Silver.from_json(value) for value in root["vector_orthonormal"])
        for packet in h300["enumeration"]["packets"]
        for root in packet["accepted_roots"]
    ]
    new_walls = [
        tuple(Silver.from_json(value) for value in root["vector_orthonormal"])
        for root in roots
    ]
    walls = (*silver_d24_height_zero_simple_roots(), *old_walls, *new_walls)
    certificate = certify_h500_replacement_hyperideal_ray(walls)
    assert certificate.vector == H500_REPLACEMENT_WITNESS_VECTOR
    assert certificate.norm == Silver(693, -490) / 5
    assert certificate.norm.sign(0) > 0
    assert certificate.wall_count == 4374
    assert certificate.active_facets == (
        *range(4),
        *range(5, 24),
        4361,
    )
    assert certificate.active_rank == 24
    assert all(value.sign(0) <= 0 for value in certificate.pairings)
    assert payload["geometry_status"]["h500_finite_volume"] is False
    assert payload["geometry_status"]["h500_compact"] is False
    assert payload["geometry_status"]["eventual_chamber_status"] == "unresolved"
