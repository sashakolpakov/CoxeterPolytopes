from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_leech import LeechGramTriple
from coxeter_polytopes.search.silver_leech_manifest import (
    first_necessary_packet,
    full_manifest_first_packet_obstructions,
    has_short_binary_combination,
    mirror_norm_classes,
    norm89_priority_obstructions,
    primitive_admissible_coefficients,
    scalar_shell_priority_obstructions,
    silver_ideal_atoms,
    silver_ideal_squareclasses,
    silver_leech_ideal_squareclass_manifest,
    silver_leech_low_norm_manifest,
)


def test_primitive_norm128_squareclass_manifest_is_finite_and_canonical() -> None:
    coefficients = primitive_admissible_coefficients(128)
    assert len(coefficients) == 50
    assert coefficients[0].coefficient == -1 + SQRT2
    assert coefficients[-1].coefficient == 1 + 8 * SQRT2
    assert all(row.coefficient.norm() == -row.absolute_norm for row in coefficients)
    assert all(row.coefficient.sign(0) > 0 > row.coefficient.sign(1) for row in coefficients)


def test_full_squarefree_ideal_manifest_includes_content_and_split_choices() -> None:
    atoms = silver_ideal_atoms(128)
    assert len(atoms) == 32
    lookup = {atom.label: atom for atom in atoms}
    assert lookup["p2"].ideal_norm == 2
    assert lookup["p3"].ideal_norm == 9
    assert lookup["p3"].splitting_type == "inert"
    assert lookup["p7+"].ideal_norm == lookup["p7-"].ideal_norm == 7
    assert lookup["p7+"].generator.conjugate() == lookup["p7-"].generator

    squareclasses = silver_ideal_squareclasses(128)
    assert len(squareclasses) == 57
    by_lineage = {row.atom_labels: row for row in squareclasses}
    assert by_lineage[("p3",)].coefficient == -3 + 3 * SQRT2
    assert by_lineage[("p2", "p3")].coefficient == 3 * SQRT2
    assert by_lineage[("p7+",)].coefficient == 1 + 2 * SQRT2
    assert by_lineage[("p7-",)].coefficient == -1 + 2 * SQRT2
    assert all(row.coefficient.norm() == -row.ideal_norm for row in squareclasses)

    manifest = silver_leech_ideal_squareclass_manifest(128)
    assert len(manifest) == 57
    lookup_rows = {row.coefficient: row for row in manifest}
    assert lookup_rows[-3 + 3 * SQRT2].first_necessary_height == 21 + 15 * SQRT2
    assert lookup_rows[3 * SQRT2].first_necessary_height == 72 + 51 * SQRT2
    assert lookup_rows[3 + 6 * SQRT2].first_packet_shell_norms == (4,)
    assert lookup_rows[7 * SQRT2].first_necessary_height == 31 + 22 * SQRT2
    assert lookup_rows[7 * SQRT2].first_packet_shell_norms == (4,)


def test_p7_local_classes_and_complete_binary_minimum_filter() -> None:
    classes = mirror_norm_classes(-1 + 2 * SQRT2)
    assert tuple(row.canonical_norm for row in classes) == (
        Silver(1),
        Silver(2),
        3 + SQRT2,
        6 + 2 * SQRT2,
    )
    assert tuple(row.odd_spatial_divisor for row in classes) == (
        Silver(1),
        Silver(1),
        3 + SQRT2,
        3 + SQRT2,
    )
    assert has_short_binary_combination(LeechGramTriple(32, 14, 21))
    assert not has_short_binary_combination(LeechGramTriple(36, 14, 22))


def test_low_norm_manifest_first_packets_and_scalar_shortlist() -> None:
    manifest = silver_leech_low_norm_manifest(128)
    assert len(manifest) == 50
    lookup = {row.coefficient: row for row in manifest}
    assert lookup[-1 + SQRT2].first_necessary_height == 14 + 10 * SQRT2
    assert lookup[SQRT2].first_packet_shell_norms == (6,)
    assert lookup[-1 + 2 * SQRT2].first_packet_shell_norms == (4,)
    assert lookup[-6 + 7 * SQRT2].first_packet_shell_norms == (4,)
    assert lookup[-4 + 7 * SQRT2].first_packet_shell_norms == (6,)
    assert lookup[-3 + 7 * SQRT2].first_packet_triple_count == 5
    assert lookup[3 + 7 * SQRT2].first_packet_triple_count == 8
    assert first_necessary_packet(5 + 8 * SQRT2)[0] == 277 + 196 * SQRT2
    assert first_necessary_packet(7 + 9 * SQRT2)[0] == 335 + 237 * SQRT2


def test_norm89_priority_first_packets_have_exact_hyperideal_obstructions() -> None:
    first, second = norm89_priority_obstructions()
    assert first.coefficient == -3 + 7 * SQRT2
    assert first.first_height == 117 + 83 * SQRT2
    assert first.maximum_pairing == 21 + 14 * SQRT2
    assert first.wall_bound == 30 + 19 * SQRT2
    assert first.witness_coordinate == Silver(4) / 3
    assert first.wall_margin == 2 + SQRT2 / 3
    assert first.hyperideal_excess == Silver(91) / 9 - 7 * SQRT2
    assert first.hyperideal_excess.sign(0) > 0

    assert second.coefficient == 3 + 7 * SQRT2
    assert second.first_height == 219 + 155 * SQRT2
    assert second.maximum_pairing == 29 + 21 * SQRT2
    assert second.wall_bound == 54 + 37 * SQRT2
    assert second.witness_coordinate == Silver(9) / 5
    assert second.wall_margin == Silver(9, -4) / 5
    assert second.hyperideal_excess == Silver(249) / 25 - 7 * SQRT2
    assert second.hyperideal_excess.sign(0) > 0


def test_new_scalar_shell_rows_are_all_strictly_hyperideal() -> None:
    rows = scalar_shell_priority_obstructions()
    assert tuple(row.coefficient for row in rows) == (
        -6 + 7 * SQRT2,
        3 + 6 * SQRT2,
        -4 + 7 * SQRT2,
        7 * SQRT2,
    )
    assert tuple(row.first_height for row in rows) == (
        19 + 14 * SQRT2,
        33 + 24 * SQRT2,
        29 + 21 * SQRT2,
        31 + 22 * SQRT2,
    )
    assert all(row.maximum_pairing == 9 + 6 * SQRT2 for row in rows)
    assert all(row.wall_margin.sign(0) > 0 for row in rows)
    assert all(row.hyperideal_excess.sign(0) > 0 for row in rows)
    assert rows[0].hyperideal_excess == Silver(17457) / 196 - Silver(827, 0) * SQRT2 / 14
    assert rows[2].hyperideal_excess == Silver(5009) / 98 - Silver(211, 0) * SQRT2 / 7


def test_every_norm128_direction_prefix_has_a_strict_first_packet_patch() -> None:
    rows = full_manifest_first_packet_obstructions(128)
    assert len(rows) == 57
    assert all(row.direction_square in (4, 6, 8, 10) for row in rows)
    assert all(row.minimum_wall_margin.sign(0) > 0 for row in rows)
    assert all(row.witness_hyperideal_excess.sign(0) > 0 for row in rows)
    tightest = min(rows, key=lambda row: row.limiting_hyperideal_excess)
    assert tightest.coefficient == 3 + 8 * SQRT2
    assert tightest.ideal_norm == 119
    assert tightest.direction_square == 10
    assert tightest.limiting_coordinate == Silver(53) / 47 + Silver(5, 0) * SQRT2 / 94
    assert tightest.limiting_hyperideal_excess == (
        Silver(21588) / 2209 - Silver(15022, 0) * SQRT2 / 2209
    )
