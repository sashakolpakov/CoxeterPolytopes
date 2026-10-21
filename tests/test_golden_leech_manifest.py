import json
from pathlib import Path

from coxeter_polytopes.exact.golden import ONE, PHI, Golden
from coxeter_polytopes.search.silver_leech import LeechGramTriple
from coxeter_polytopes.search.golden_leech_manifest import (
    CLASS_AND_UNIT_PROOF,
    HEIGHT_ZERO_MINIMUM_TRACE,
    HEIGHT_ZERO_PROOF,
    MIRROR_CLASS_PROOF,
    candidate_golden_leech_gram_triples,
    first_necessary_packet,
    full_manifest_first_packet_obstructions,
    golden_ideal_atoms,
    golden_ideal_squareclasses,
    golden_leech_ideal_squareclass_manifest,
    has_short_binary_combination,
    mirror_norm_classes,
)


REPOSITORY = Path(__file__).resolve().parents[1]


def test_golden_prime_ideal_and_squarefree_lineage_census() -> None:
    assert "norm-Euclidean" in CLASS_AND_UNIT_PROOF
    assert HEIGHT_ZERO_MINIMUM_TRACE == 8
    assert "there are no height-zero mirrors" in HEIGHT_ZERO_PROOF
    atoms = golden_ideal_atoms(128)
    assert len(atoms) == 28
    lookup = {atom.label: atom for atom in atoms}
    assert lookup["p2"].ideal_norm == 4
    assert lookup["p2"].splitting_type == "inert"
    assert lookup["p5"].ideal_norm == 5
    assert lookup["p5"].splitting_type == "ramified"
    assert lookup["p7"].ideal_norm == 49
    assert lookup["p11+"].generator.conjugate() == lookup["p11-"].generator

    rows = golden_ideal_squareclasses(128)
    assert len(rows) == 47
    by_lineage = {row.atom_labels: row for row in rows}
    assert by_lineage[()].coefficient == -ONE + PHI
    assert by_lineage[("p2",)].coefficient == -2 + 2 * PHI
    assert by_lineage[("p5",)].coefficient == -1 + 2 * PHI
    assert by_lineage[("p11+",)].coefficient == -1 + 3 * PHI
    assert by_lineage[("p11-",)].coefficient == -2 + 3 * PHI
    assert by_lineage[("p11+", "p11-")].coefficient == -11 + 11 * PHI
    assert all(row.coefficient.norm() == -row.ideal_norm for row in rows)
    assert all(row.coefficient.sign(0) > 0 > row.coefficient.sign(1) for row in rows)
    for row in rows:
        odd_prime_count = sum(label != "p2" for label in row.atom_labels)
        expected_class_count = 2**odd_prime_count
        if "p2" not in row.atom_labels:
            expected_class_count *= 2
        assert len(mirror_norm_classes(row.coefficient)) == expected_class_count


def test_inert_dyadic_and_odd_local_mirror_classes_are_complete() -> None:
    assert "prime (2) is inert" in MIRROR_CLASS_PROOF
    assert tuple(row.canonical_norm for row in mirror_norm_classes(-ONE + PHI)) == (
        Golden(1),
        Golden(2),
    )
    # When (2) divides squarefree (a), q has dyadic valuation exactly one.
    assert tuple(row.canonical_norm for row in mirror_norm_classes(-2 + 2 * PHI)) == (
        Golden(2),
    )
    ramified = mirror_norm_classes(-1 + 2 * PHI)
    assert tuple(row.canonical_norm for row in ramified) == (
        Golden(1),
        Golden(2),
        2 + PHI,
        4 + 2 * PHI,
    )
    assert tuple(row.odd_spatial_divisor for row in ramified) == (
        Golden(1),
        Golden(1),
        2 + PHI,
        2 + PHI,
    )
    # The two primes over 11 vary independently, giving 2*2^2 classes.
    assert len(mirror_norm_classes(-11 + 11 * PHI)) == 8


def test_golden_leech_gram_equations_and_short_vector_filter() -> None:
    triples = candidate_golden_leech_gram_triples(14 + 10 * PHI)
    assert triples == (
        LeechGramTriple(10, 4, 3),
        LeechGramTriple(8, 6, 2),
        LeechGramTriple(6, 8, 1),
        LeechGramTriple(4, 10, 0),
    )
    assert has_short_binary_combination(LeechGramTriple(32, 14, 21))
    assert not has_short_binary_combination(LeechGramTriple(4, 10, 5))


def test_all_norm128_directions_have_complete_first_necessary_packets() -> None:
    rows = golden_leech_ideal_squareclass_manifest(128)
    assert len(rows) == 47
    lookup = {row.coefficient: row for row in rows}
    assert lookup[-ONE + PHI].first_necessary_height == 6 + 10 * PHI
    assert lookup[-2 + 2 * PHI].first_necessary_height == 1 + 2 * PHI
    assert lookup[-2 + 2 * PHI].first_packet_shell_norms == (4,)
    assert lookup[-1 + 2 * PHI].first_necessary_height == 11 + 18 * PHI
    assert lookup[-1 + 2 * PHI].first_packet_shell_norms == (6,)
    assert lookup[-5 + 6 * PHI].first_necessary_height == 122 + 198 * PHI
    assert lookup[-5 + 6 * PHI].first_packet_triple_count == 12
    assert first_necessary_packet(-11 + 11 * PHI)[0] == 100 + 162 * PHI


def test_every_norm128_first_packet_has_a_strict_hyperideal_witness() -> None:
    rows = full_manifest_first_packet_obstructions(128)
    assert len(rows) == 47
    assert all(row.direction_square in (4, 6, 8, 10) for row in rows)
    assert all(row.minimum_wall_margin.sign(0) > 0 for row in rows)
    assert all(row.witness_hyperideal_excess.sign(0) > 0 for row in rows)
    tightest = min(rows, key=lambda row: row.limiting_hyperideal_excess)
    assert tightest.coefficient == -5 + 6 * PHI
    assert tightest.atom_labels == ("p41-",)
    assert tightest.direction_square == 4
    assert tightest.limiting_coordinate == Golden(-92, 94) / 55
    assert tightest.limiting_hyperideal_excess == Golden(3373) / 121 - Golden(10398, 0) * PHI / 605
    assert tightest.shrink_denominator == 142


def test_artifact_states_direction_and_nonterminal_scope() -> None:
    artifact = json.loads(
        (REPOSITORY / "data/enumerations/golden_leech_low_norm_manifest.json").read_text()
    )
    assert "K-square" in artifact["direction_interpretation"]
    assert "not an integral-host classification" in artifact["direction_interpretation"]
    assert "later mirrors can cut" in artifact["scope"]
    assert "does not prove global nonreflectivity" in artifact["scope"]
