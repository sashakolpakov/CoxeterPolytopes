from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

from coxeter_polytopes.exact.golden import PHI, Golden
from coxeter_polytopes.search.golden_odd_leech_section import (
    CANONICAL_MIRROR_NORMS,
    section_norm_five_shell_audit,
    section_norm_four_shell_audit,
    twice_unit_section_first_packet_audit,
    unit_section_first_packet_audit,
)
from coxeter_polytopes.search.silver_odd_leech import (
    OddLeechGramTriple,
    catalog_golay_codewords,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/golden_odd_leech_section_h23_first_packets.json"


def test_golden_section_local_alphabet_and_height_zero():
    audit = unit_section_first_packet_audit()
    assert audit.model.determinant == 3
    assert audit.model.smith_invariants == (1,) * 22 + (3,)
    assert audit.model.minimum == 3
    assert audit.model.minimal_shell_count == 2576
    assert CANONICAL_MIRROR_NORMS == tuple(Golden(n) for n in (1, 2, 3, 6))
    assert [len(types) for _norm, types in audit.height_zero_ambient_types] == [
        0, 0, 1, 1
    ]
    assert all(not types for _norm, types in audit.height_zero_gated_types)


def test_norm_four_and_five_shell_counts_from_golay_weight_enumerator():
    # Clean-room zero-sum counts directly from every catalogue word.  This
    # does not use the shell component totals frozen by the search module.
    words = catalog_golay_codewords()
    weights = Counter(word.bit_count() for word in words)
    balanced_octad_signs = sum(
        1 for signs in range(1 << 8) if signs.bit_count() == 4
    )
    norm_four = weights[8] * balanced_octad_signs + 24 * 23
    odd_norm_four = 0
    norm_five = 0
    for word in words:
        signs = tuple(-1 if word & (1 << index) else 1 for index in range(24))
        baseline_sum = sum(signs)
        odd_norm_four += sum(
            baseline_sum - 4 * sign == 0 for sign in signs
        )
        norm_five += sum(
            baseline_sum - 4 * (signs[left] + signs[right]) == 0
            for left, right in combinations(range(24), 2)
        )
    assert odd_norm_four == 0
    assert norm_four == 53_682
    assert norm_five == 553_104
    assert section_norm_four_shell_audit().total_count == norm_four
    assert section_norm_five_shell_audit().total_count == norm_five


def test_phi_first_packet_and_strict_hyperideal_witness():
    audit = unit_section_first_packet_audit()
    nonempty = [level for level in audit.controls if level.necessary_gram_triples]
    assert len(audit.controls) == 7
    assert len(nonempty) == 1
    assert audit.first_height == Golden(Fraction(3, 2), Fraction(5, 2))
    assert audit.first_norm == 2
    assert audit.first_control == 1 + PHI
    assert audit.first_target == Golden(5, 5)
    assert audit.first_gram_type == OddLeechGramTriple(0, 5, 0)
    assert audit.first_packet_root_count == 553_104
    assert audit.packet_pairing_upper_bound == -2 * PHI
    assert audit.witness_spatial_scale == Fraction(2, 3)
    assert audit.witness_square == Golden(Fraction(16, 9), -1)
    assert audit.closest_wall_pairing == -1
    assert audit.strict_hyperideal and not audit.finite_volume and not audit.compact


def test_twice_phi_first_packet_and_strict_hyperideal_witness():
    audit = twice_unit_section_first_packet_audit()
    nonempty = [level for level in audit.controls if level.necessary_gram_triples]
    assert len(audit.controls) == 2
    assert len(nonempty) == 1
    assert audit.first_height == Golden(1, 2)
    assert audit.first_norm == 2
    assert audit.first_control == PHI
    assert audit.first_target == Golden(4, 4)
    assert audit.first_gram_type == OddLeechGramTriple(0, 4, 0)
    assert audit.first_packet_root_count == 53_682
    assert audit.packet_pairing_upper_bound == -2 * PHI
    assert audit.witness_spatial_scale == Fraction(4, 3)
    assert audit.witness_square == Golden(Fraction(16, 3), -2)
    assert audit.closest_wall_pairing == Golden(-2, Fraction(2, 3))
    assert audit.strict_hyperideal and not audit.finite_volume and not audit.compact


def test_frozen_section_schema_distinguishes_gram_rank_from_radical_dimension():
    model = json.loads(ARTIFACT.read_text())["section_model"]
    assert model["gram_rank_mod_three"] == 22
    assert model["pairing_radical_dimension_mod_three"] == 1
    assert "radical_mod_three_rank" not in model
