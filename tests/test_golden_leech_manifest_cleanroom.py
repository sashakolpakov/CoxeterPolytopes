import json
from dataclasses import asdict
from pathlib import Path

from coxeter_polytopes.exact.golden import ONE, PHI, Golden
from coxeter_polytopes.search.golden_leech_manifest import (
    full_manifest_first_packet_obstructions,
    golden_ideal_atoms,
    golden_ideal_squareclasses,
    golden_leech_ideal_squareclass_manifest,
)
from coxeter_polytopes.search.golden_leech_manifest_cleanroom import (
    EXPLICIT_DIRECTION_VECTORS_SCALED,
    clean_golden_atoms,
    clean_golden_candidate_triples,
    clean_golden_directions,
    clean_golden_has_short_combination,
    clean_golden_manifest,
    clean_golden_mirror_classes,
    clean_golden_prefix_witnesses,
    mog_scaled_square,
)


REPOSITORY = Path(__file__).resolve().parents[1]


def test_cleanroom_prime_ideal_and_direction_census() -> None:
    atoms = clean_golden_atoms()
    assert len(atoms) == 28
    assert sum(atom.splitting_type == "inert" for atom in atoms) == 3
    assert sum(atom.splitting_type == "ramified" for atom in atoms) == 1
    assert sum(atom.splitting_type == "split" for atom in atoms) == 24
    lookup = {atom.label: atom for atom in atoms}
    assert (lookup["p2"].ideal_norm, lookup["p2"].splitting_type) == (4, "inert")
    assert (lookup["p5"].ideal_norm, lookup["p5"].splitting_type) == (5, "ramified")
    assert lookup["p7"].ideal_norm == 49
    assert lookup["p11+"].generator.conjugate() == lookup["p11-"].generator

    directions = clean_golden_directions()
    assert len(directions) == 47
    assert len({row.atom_labels for row in directions}) == 47
    assert all(row.coefficient.norm() == -row.ideal_norm for row in directions)
    assert all(row.coefficient.sign(0) > 0 > row.coefficient.sign(1) for row in directions)
    by_lineage = {row.atom_labels: row for row in directions}
    assert by_lineage[()].coefficient == -ONE + PHI
    assert by_lineage[("p2",)].coefficient == -2 + 2 * PHI
    assert by_lineage[("p5",)].coefficient == -1 + 2 * PHI
    assert by_lineage[("p11+",)].coefficient == -1 + 3 * PHI
    assert by_lineage[("p11-",)].coefficient == -2 + 3 * PHI
    assert by_lineage[("p11+", "p11-")].coefficient == -11 + 11 * PHI


def test_cleanroom_inert_dyadic_and_independent_odd_classes() -> None:
    directions = {row.atom_labels: row for row in clean_golden_directions()}
    unit = clean_golden_mirror_classes(directions[()])
    assert tuple(row.canonical_norm for row in unit) == (Golden(1), Golden(2))

    # Since (2) is inert and divides the squarefree coefficient, the norm
    # valuation is exactly one and there is no spatial gate.
    dyadic = clean_golden_mirror_classes(directions[("p2",)])
    assert tuple(row.canonical_norm for row in dyadic) == (Golden(2),)
    assert tuple(row.odd_spatial_divisor for row in dyadic) == (Golden(1),)

    ramified = clean_golden_mirror_classes(directions[("p5",)])
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
    both_eleven = clean_golden_mirror_classes(directions[("p11+", "p11-")])
    assert len(both_eleven) == 8
    assert {row.odd_spatial_divisor for row in both_eleven} == {
        Golden(1),
        Golden(4, -1),
        Golden(3, 1),
        Golden(11),
    }


def test_cleanroom_binary_gram_filter_is_exhaustive() -> None:
    # Unit normalization of a height-zero mirror would require y^2=2.
    assert clean_golden_candidate_triples(Golden(2)) == ()
    triples = clean_golden_candidate_triples(14 + 10 * PHI)
    assert tuple((row.u_square, row.v_square, row.u_dot_v) for row in triples) == (
        (10, 4, 3),
        (8, 6, 2),
        (6, 8, 1),
        (4, 10, 0),
    )
    assert clean_golden_has_short_combination(type(triples[0])(32, 14, 21))
    assert not clean_golden_has_short_combination(type(triples[0])(4, 10, 5))
    assert all(
        not clean_golden_has_short_combination(triple)
        for row in clean_golden_manifest()
        for component in row.first_packet_components
        for triple in component.surviving_triples
    )


def test_cleanroom_reconstructs_all_first_packets() -> None:
    rows = clean_golden_manifest()
    assert len(rows) == 47
    assert all(len(row.first_packet_components) == 1 for row in rows)
    lookup = {row.coefficient: row for row in rows}
    assert lookup[-ONE + PHI].first_necessary_height == 6 + 10 * PHI
    assert lookup[-2 + 2 * PHI].first_necessary_height == 1 + 2 * PHI
    assert lookup[-1 + 2 * PHI].first_necessary_height == 11 + 18 * PHI
    assert lookup[-5 + 6 * PHI].first_necessary_height == 122 + 198 * PHI
    assert len(lookup[-5 + 6 * PHI].first_packet_components[0].surviving_triples) == 12
    assert lookup[-11 + 11 * PHI].first_necessary_height == 100 + 162 * PHI
    for row in rows:
        for component in row.first_packet_components:
            assert component.target_after_gate == (
                component.target_before_gate
                / (component.odd_spatial_divisor * component.odd_spatial_divisor)
            )


def test_cleanroom_all_prefix_witnesses_are_strict_and_p41_is_tightest() -> None:
    assert tuple(map(mog_scaled_square, EXPLICIT_DIRECTION_VECTORS_SCALED)) == (4, 6, 8, 10)
    rows = {row.coefficient: row for row in clean_golden_manifest()}
    witnesses = clean_golden_prefix_witnesses()
    assert len(witnesses) == 47
    for witness in witnesses:
        row = rows[witness.coefficient]
        assert witness.component_wall_bounds == tuple(
            row.coefficient * component.control / component.odd_spatial_divisor
            for component in row.first_packet_components
        )
        assert witness.minimum_wall_margin.sign(0) > 0
        assert witness.witness_hyperideal_excess.sign(0) > 0
    tightest = min(witnesses, key=lambda row: row.limiting_hyperideal_excess)
    assert tightest.coefficient == -5 + 6 * PHI
    assert tightest.atom_labels == ("p41-",)
    assert tightest.direction_square == 4
    assert tightest.component_pairing_bounds == (17 + 26 * PHI,)
    assert tightest.component_wall_bounds == (16 + 30 * PHI,)
    assert tightest.limiting_coordinate == Golden(-92, 94) / 55
    assert tightest.limiting_hyperideal_excess == (
        Golden(3373) / 121 - Golden(10398) * PHI / 605
    )
    assert tightest.shrink_denominator == 142
    assert tightest.minimum_wall_margin == Golden(8, 15) / 71


def test_cleanroom_matches_owner_backend_field_by_field() -> None:
    assert tuple(map(asdict, golden_ideal_atoms())) == tuple(map(asdict, clean_golden_atoms()))
    assert tuple(map(asdict, golden_ideal_squareclasses())) == tuple(
        map(asdict, clean_golden_directions())
    )
    owners = golden_leech_ideal_squareclass_manifest()
    cleans = clean_golden_manifest()
    for owner, clean in zip(owners, cleans, strict=True):
        assert owner.coefficient == clean.coefficient
        assert owner.ideal_norm == clean.ideal_norm
        assert owner.atom_labels == clean.atom_labels
        assert tuple(map(asdict, owner.mirror_norm_classes)) == tuple(
            map(asdict, clean.mirror_norm_classes)
        )
        assert owner.first_necessary_height == clean.first_necessary_height
        assert tuple(map(asdict, owner.first_packet_components)) == tuple(
            map(asdict, clean.first_packet_components)
        )
    assert tuple(map(asdict, full_manifest_first_packet_obstructions())) == tuple(
        map(asdict, clean_golden_prefix_witnesses())
    )


def test_cleanroom_artifact_records_bounded_nonterminal_scope() -> None:
    artifact = json.loads(
        (
            REPOSITORY
            / "data/enumerations/golden_leech_low_norm_manifest_cleanroom.json"
        ).read_text()
    )
    assert artifact["counts"] == {
        "coefficient_directions": 47,
        "first_packets": 47,
        "prime_ideal_atoms": 28,
        "strict_first_prefix_witnesses": 47,
    }
    assert artifact["owner_backend_imported"] is False
    assert "K-square" in artifact["coefficient_direction_scope"]
    assert "not an integral-host classification" in artifact["coefficient_direction_scope"]
    assert "later mirrors may cut" in artifact["volume_scope"]
    assert "does not decide final volume" in artifact["volume_scope"]


def test_cleanroom_artifact_matches_owner_artifact_exactly_on_audited_fields() -> None:
    directory = REPOSITORY / "data/enumerations"
    owner = json.loads((directory / "golden_leech_low_norm_manifest.json").read_text())
    clean = json.loads(
        (directory / "golden_leech_low_norm_manifest_cleanroom.json").read_text()
    )
    assert owner["prime_ideal_atoms"] == clean["prime_ideal_atoms"]
    assert owner["ideal_squareclasses"] == clean["directions"]
    for owner_row, clean_row in zip(owner["rows"], clean["rows"], strict=True):
        assert owner_row["coefficient"] == clean_row["coefficient"]
        assert owner_row["ideal_norm"] == clean_row["ideal_norm"]
        assert owner_row["atom_labels"] == clean_row["atom_labels"]
        assert owner_row["mirror_norm_classes"] == clean_row["mirror_norm_classes"]
        assert owner_row["first_necessary_height"] == clean_row["first_necessary_height"]
        assert owner_row["first_packet_components"] == clean_row["first_packet_components"]
    assert owner["all_first_packet_obstructions"] == clean["first_prefix_witnesses"]
