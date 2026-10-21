import json
from pathlib import Path

from sympy import Matrix

from coxeter_polytopes.exact.cubic7 import THETA7, Cubic7
from coxeter_polytopes.search.cubic7_odd_leech_section import (
    CANONICAL_MIRROR_NORMS,
    COMPLETE_MIRROR_NORM_PROOF,
    FIRST_PACKET_CONTROL,
    FIRST_PACKET_HEIGHT,
    FIRST_PACKET_NORM,
    FIRST_PACKET_SCALE,
    FIRST_PACKET_TARGET,
    RADICAL_NORM_EXCLUSION,
    SECTION_BASIS_SCALED,
    SECTION_GRAM,
    SECTION_WITNESS_SCALED,
    iter_minimal_section_vectors_scaled,
    odd_leech_section_first_packet_audit,
    odd_leech_section_model_audit,
    section_controls_through_first_packet,
    section_rank_three_gram_types,
    section_smith_invariants,
)
from coxeter_polytopes.search.silver_odd_leech import (
    odd_leech_pair_scaled,
    odd_leech_raw_member,
)
ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "data/enumerations/cubic7_odd_leech_section_h23_first_packet.json"
)


def _payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_exact_rank23_section_basis_determinant_and_smith_form() -> None:
    audit = odd_leech_section_model_audit()
    assert len(SECTION_BASIS_SCALED) == len(SECTION_GRAM) == 23
    assert all(sum(row) == 0 and odd_leech_raw_member(row) for row in SECTION_BASIS_SCALED)
    assert audit.determinant == 3
    assert audit.smith_invariants == section_smith_invariants() == (1,) * 22 + (3,)
    assert all(value > 0 for value in audit.leading_principal_minors)


def test_complete_minimal_shell_is_weight12_golay_section() -> None:
    audit = odd_leech_section_model_audit()
    shell = tuple(iter_minimal_section_vectors_scaled())
    assert len(shell) == audit.minimal_shell_count == 2576
    assert audit.minimum == 3
    assert audit.minimal_shell_span_rank == 23
    assert all(sum(vector) == 0 for vector in shell)
    assert all(odd_leech_pair_scaled(vector, vector) == 3 for vector in shell)
    assert audit.minimal_shell_digest == (
        "2860f05b59483cf95a28433e76915d6862cf05fd1348303e938233033e2923b9"
    )


def test_mod3_pairing_radical_excludes_height_zero_norms() -> None:
    audit = odd_leech_section_model_audit()
    generator = Matrix(audit.radical_generator_coordinates)
    assert audit.radical_mod_three_rank == 22
    assert all(int(value) % 3 == 0 for value in Matrix(SECTION_GRAM) * generator)
    assert audit.radical_generator_raw_residue == (1,) * 24
    assert audit.radical_minimum_lower_bound == 7
    assert "no norm 3 or 6 vector" in RADICAL_NORM_EXCLUSION


def test_complete_mirror_norm_alphabet_and_height_zero_census() -> None:
    audit = odd_leech_section_first_packet_audit()
    assert CANONICAL_MIRROR_NORMS == tuple(Cubic7(value) for value in (1, 2, 3, 6))
    assert "complete canonical norms" in COMPLETE_MIRROR_NORM_PROOF
    assert [len(types) for _norm, types in audit.height_zero_ambient_types] == [0, 0, 1, 1]
    assert not any(types for _norm, types in audit.height_zero_gated_types)
    assert audit.height_zero_wall_count == 0
    assert section_rank_three_gram_types(Cubic7(3), 3)[0].u_square == 3
    assert section_rank_three_gram_types(Cubic7(6), 3)[0].u_square == 6


def test_complete_first_control_census_and_unique_realized_packet() -> None:
    audit = odd_leech_section_first_packet_audit()
    controls = section_controls_through_first_packet()
    assert controls == audit.controls
    assert len(controls) == 5
    nonempty = [level for level in controls if level.necessary_gram_types]
    assert len(nonempty) == 1
    first = nonempty[0]
    assert (
        first.height,
        first.norm,
        first.control,
        first.target,
    ) == (
        FIRST_PACKET_HEIGHT,
        FIRST_PACKET_NORM,
        FIRST_PACKET_CONTROL,
        FIRST_PACKET_TARGET,
    )
    assert first.necessary_gram_types[0].matrix == (
        (3, 3, 0),
        (3, 3, 0),
        (0, 0, 0),
    )
    assert first.necessary_gram_types[0].rank == 1


def test_all_2576_first_roots_are_primitive_integral_and_mutually_acute() -> None:
    audit = odd_leech_section_first_packet_audit()
    assert abs(FIRST_PACKET_CONTROL.norm()) == abs(FIRST_PACKET_SCALE.norm()) == 1
    assert (
        -THETA7 * FIRST_PACKET_CONTROL**2 + 3 * FIRST_PACKET_SCALE**2
        == FIRST_PACKET_NORM
        == 2
    )
    assert audit.first_packet_root_count == 2576
    assert audit.distinct_spatial_pairings == (-3, -1, 0, 1)
    assert all(pairing.sign(0) < 0 for pairing in audit.distinct_root_pairings)


def test_strict_hyperideal_witness_and_frozen_artifact_replay() -> None:
    audit = odd_leech_section_first_packet_audit()
    assert odd_leech_raw_member(SECTION_WITNESS_SCALED)
    assert sum(SECTION_WITNESS_SCALED) == 0
    assert odd_leech_pair_scaled(SECTION_WITNESS_SCALED, SECTION_WITNESS_SCALED) == 4
    assert audit.model.witness_shell_pairing_histogram == ((-1, 672), (0, 1232), (1, 672))
    assert audit.witness_minimum_wall_margin == Cubic7(-1, 1, 1)
    assert audit.witness_minimum_wall_margin.sign(0) > 0
    assert audit.witness_hyperideal_excess == 4 - THETA7
    assert audit.witness_hyperideal_excess.sign(0) > 0
    assert not audit.finite_volume and not audit.compact
    payload = _payload()
    assert payload["status"] == (
        "FIRST_PACKET_INFINITE_VOLUME; EVENTUAL_CHAMBER_UNRESOLVED"
    )
    assert payload["spatial_lattice"]["determinant"] == audit.model.determinant
    assert payload["spatial_lattice"]["minimal_shell_digest"] == (
        audit.model.minimal_shell_digest
    )
    assert payload["complete_control_census"]["control_count"] == len(audit.controls)
    assert payload["first_packet"]["root_count"] == audit.first_packet_root_count
    assert payload["strict_hyperideal_certificate"]["finite_volume"] is False
    assert payload["strict_hyperideal_certificate"]["compact"] is False
