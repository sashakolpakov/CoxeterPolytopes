import json
from pathlib import Path

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_leech_manifest_cleanroom import (
    CleanTriple,
    EXPLICIT_DIRECTION_VECTORS_SCALED,
    clean_atoms,
    clean_directions,
    clean_has_short_combination,
    clean_manifest,
    clean_mirror_classes,
    clean_prefix_witnesses,
    mog_scaled_square,
)


REPOSITORY = Path(__file__).resolve().parents[1]


def test_cleanroom_prime_ideal_and_direction_census() -> None:
    atoms = clean_atoms()
    assert len(atoms) == 32
    labels = {atom.label for atom in atoms}
    assert {"p2", "p3", "p5", "p11", "p7+", "p7-"} <= labels
    assert sum(atom.splitting_type == "ramified" for atom in atoms) == 1
    assert sum(atom.splitting_type == "inert" for atom in atoms) == 3
    assert sum(atom.splitting_type == "split" for atom in atoms) == 28

    directions = clean_directions()
    assert len(directions) == 57
    assert len({row.atom_labels for row in directions}) == 57
    assert all(row.coefficient.norm() == -row.ideal_norm for row in directions)
    assert all(row.coefficient.sign(0) > 0 > row.coefficient.sign(1) for row in directions)
    by_lineage = {row.atom_labels: row for row in directions}
    assert by_lineage[("p3",)].coefficient == -3 + 3 * SQRT2
    assert by_lineage[("p2", "p3")].coefficient == 3 * SQRT2
    assert by_lineage[("p7+",)].coefficient == 1 + 2 * SQRT2
    assert by_lineage[("p7-",)].coefficient == -1 + 2 * SQRT2
    assert by_lineage[("p7+", "p7-")].coefficient == -7 + 7 * SQRT2


def test_cleanroom_local_classes_keep_split_atoms_independent() -> None:
    directions = {row.atom_labels: row for row in clean_directions()}
    p7 = clean_mirror_classes(directions[("p7-",)])
    assert tuple(row.canonical_norm for row in p7) == (
        Silver(1),
        Silver(2),
        3 + SQRT2,
        6 + 2 * SQRT2,
    )
    assert tuple(row.odd_spatial_divisor for row in p7) == (
        Silver(1),
        Silver(1),
        3 + SQRT2,
        3 + SQRT2,
    )

    norm49 = clean_mirror_classes(directions[("p7+", "p7-")])
    assert len(norm49) == 8
    assert {row.odd_spatial_divisor for row in norm49} == {
        Silver(1),
        Silver(3, -1),
        Silver(3, 1),
        Silver(7),
    }
    inert = clean_mirror_classes(directions[("p3",)])
    assert {row.odd_spatial_divisor for row in inert} == {Silver(1), Silver(3)}


def test_cleanroom_binary_filter_exhausts_all_short_combinations() -> None:
    assert clean_has_short_combination(CleanTriple(32, 14, 21))
    assert not clean_has_short_combination(CleanTriple(36, 14, 22))
    rows = clean_manifest()
    assert all(
        not clean_has_short_combination(triple)
        for row in rows
        for component in row.first_packet_components
        for triple in component.surviving_triples
    )


def test_cleanroom_reconstructs_all_first_packets_and_gates() -> None:
    rows = clean_manifest()
    assert len(rows) == 57
    lookup = {row.coefficient: row for row in rows}
    assert lookup[-1 + SQRT2].first_necessary_height == 14 + 10 * SQRT2
    assert lookup[-1 + 2 * SQRT2].first_necessary_height == 5 + 4 * SQRT2
    assert lookup[-7 + 7 * SQRT2].first_necessary_height == 49 + 35 * SQRT2
    assert lookup[7 + 9 * SQRT2].first_necessary_height == 335 + 237 * SQRT2
    for row in rows:
        for component in row.first_packet_components:
            assert component.target_after_gate == (
                component.target_before_gate
                / (component.odd_spatial_divisor * component.odd_spatial_divisor)
            )


def test_cleanroom_all_57_witnesses_are_strict_prefix_only_certificates() -> None:
    assert tuple(map(mog_scaled_square, EXPLICIT_DIRECTION_VECTORS_SCALED)) == (4, 6, 8, 10)
    rows = clean_manifest()
    witnesses = clean_prefix_witnesses()
    assert len(witnesses) == 57
    by_coefficient = {row.coefficient: row for row in rows}
    for witness in witnesses:
        row = by_coefficient[witness.coefficient]
        assert witness.component_wall_bounds == tuple(
            row.coefficient * component.control / component.odd_spatial_divisor
            for component in row.first_packet_components
        )
        assert witness.minimum_wall_margin.sign(0) > 0
        assert witness.witness_hyperideal_excess.sign(0) > 0
    tightest = min(witnesses, key=lambda row: row.limiting_hyperideal_excess)
    assert tightest.coefficient == 3 + 8 * SQRT2
    assert tightest.direction_square == 10
    assert tightest.limiting_coordinate == Silver(53, 0) / 47 + SQRT2 * Silver(5, 0) / 94


def test_cleanroom_artifact_states_bounded_prefix_scope() -> None:
    artifact = json.loads(
        (REPOSITORY / "data/enumerations/silver_leech_low_norm_manifest_cleanroom.json").read_text()
    )
    assert artifact["counts"] == {
        "coefficient_directions": 57,
        "first_packets": 57,
        "prime_ideal_atoms": 32,
        "strict_first_prefix_witnesses": 57,
    }
    assert artifact["owner_backend_imported"] is False
    assert "later mirrors may cut" in artifact["volume_scope"]
    assert "not a final finite-volume" in artifact["volume_scope"]


def test_cleanroom_artifact_matches_owner_exactly_field_by_field() -> None:
    directory = REPOSITORY / "data/enumerations"
    owner = json.loads((directory / "silver_leech_low_norm_manifest.json").read_text())
    clean = json.loads(
        (directory / "silver_leech_low_norm_manifest_cleanroom.json").read_text()
    )
    assert owner["prime_ideal_atoms"] == clean["prime_ideal_atoms"]
    assert owner["ideal_squareclasses"] == clean["directions"]
    for owner_row, clean_row in zip(owner["rows"], clean["rows"], strict=True):
        assert owner_row["coefficient"] == clean_row["coefficient"]
        assert owner_row["absolute_norm"] == clean_row["ideal_norm"]
        assert owner_row["mirror_norm_classes"] == clean_row["mirror_norm_classes"]
        assert owner_row["first_necessary_height"] == clean_row["first_necessary_height"]
        assert owner_row["first_packet_components"] == clean_row["first_packet_components"]
    assert owner["all_first_packet_obstructions"] == clean["first_prefix_witnesses"]
