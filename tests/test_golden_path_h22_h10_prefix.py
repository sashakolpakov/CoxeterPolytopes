import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO
from coxeter_polytopes.exact.matrix import Inertia, determinant, inertia
from coxeter_polytopes.search.golden_path_host import (
    COEFFICIENT,
    DETERMINANT,
    P101,
    P5,
    WEYL_SQUARE,
    basis_vector,
    canonical_mirror_norms,
    cartan_matrix,
    golden_path_h10_audit,
    lattice_gram,
    pair,
    strict_interior_point,
    weyl_vector,
)
from scripts.certify_golden_path_h22_h10_prefix import encoded_payload


ARTIFACT = ROOT / "data/enumerations/golden_path_h22_h10_prefix.json"


def test_exact_rank23_cartan_symmetrization_and_two_real_signatures() -> None:
    cartan = cartan_matrix()
    gram = lattice_gram()
    assert len(cartan) == len(gram) == 23
    assert cartan[0][1] == -1
    assert cartan[1][0] == -COEFFICIENT
    assert all(cartan[i][i + 1] == cartan[i + 1][i] == -1 for i in range(1, 22))
    symmetrizer = (COEFFICIENT, *(ONE for _ in range(22)))
    assert gram == tuple(
        tuple(symmetrizer[i] * cartan[i][j] for j in range(23))
        for i in range(23)
    )
    assert determinant(gram) == DETERMINANT == Golden(-84, -172)
    assert DETERMINANT == -4 * P5 * Golden(4, 13)
    assert PHI * Golden(4, 13) == (ONE + PHI) * P101
    assert inertia(gram, embedding=0) == Inertia(22, 1, 0)
    assert inertia(gram, embedding=1) == Inertia(23, 0, 0)


def test_weyl_direction_and_complete_local_norm_candidate_alphabet() -> None:
    rho = weyl_vector()
    assert pair(rho, rho) == WEYL_SQUARE == Golden(-84, -172) / 23
    assert [pair(rho, basis_vector(i)) for i in range(1, 23)] == [ZERO] * 22
    assert pair(rho, basis_vector(0)) == WEYL_SQUARE
    assert canonical_mirror_norms() == tuple(
        Golden(*pair_value)
        for pair_value in (
            (1, 0),
            (2, 0),
            (2, 1),
            (4, 0),
            (4, 2),
            (8, 0),
            (8, 4),
            (9, 4),
            (23, -1),
            (16, 8),
            (18, 8),
            (46, -2),
            (36, 16),
            (92, -4),
            (72, 32),
            (184, -8),
        )
    )


def test_complete_h10_control_replay_and_chronological_wall_selection() -> None:
    audit = golden_path_h10_audit()
    assert len(audit.controls) == 383
    assert sum(row.recursion_nodes for row in audit.controls) == 19_760
    assert sum(row.energy_solution_count for row in audit.controls) == 6
    assert len(audit.candidates) == 3
    assert audit.control_stream_sha256 == (
        "1db06f9e673ea0f11e59ad220e9ca65b0408d98eb1c77e2f4701c5586bc4b1f6"
    )
    assert audit.candidate_stream_sha256 == (
        "f758c1b6a105df9b7817f0c58147d7270b4e467130ce9521c866a78923b02714"
    )

    seed, rejected, new = audit.decisions
    assert seed.accepted
    assert seed.candidate.height == Golden(-1, 22) / 23
    assert seed.candidate.norm == 4
    assert seed.candidate.vector == (PHI - 1, *(ZERO for _ in range(22)))
    assert rejected.accepted is False
    assert rejected.first_positive_wall == 22
    assert rejected.first_positive_pairing == 2 * PHI
    assert rejected.candidate.vector == (ONE, ONE, *(ZERO for _ in range(21)))
    assert new.accepted
    assert new.candidate.height == Golden(42, 86) / 23
    assert new.candidate.norm == 2
    assert new.candidate.vector == (
        ONE,
        ONE + PHI,
        ZERO,
        *(-PHI for _ in range(20)),
    )
    assert len(audit.walls) == 24


def test_true_24_wall_prefix_has_exact_hyperideal_extreme_ray() -> None:
    audit = golden_path_h10_audit()
    assert (
        len(audit.seed_geometry.finite_vertices),
        len(audit.seed_geometry.cusps),
        len(audit.seed_geometry.hyperideal_rays),
        len(audit.seed_geometry.opposite_sheet_rays),
    ) == (2, 0, 21, 0)
    assert (
        len(audit.geometry.finite_vertices),
        len(audit.geometry.cusps),
        len(audit.geometry.hyperideal_rays),
        len(audit.geometry.opposite_sheet_rays),
    ) == (2, 0, 41, 0)
    witness = audit.hyperideal_witness
    assert witness.active_facets == (*range(21), 23)
    assert witness.norm == Golden(5_708_212, 7_337_356) / (571 * 571)
    assert witness.norm.sign(0) > 0
    assert all(value.sign(0) <= 0 for value in witness.wall_pairings)
    interior = strict_interior_point()
    assert pair(interior, interior).sign(0) < 0
    assert all(pair(interior, wall).sign(0) < 0 for wall in audit.walls)


def test_frozen_json_is_a_byte_exact_replay() -> None:
    encoded = encoded_payload()
    assert ARTIFACT.read_bytes() == encoded
    payload = json.loads(encoded)
    assert payload["verdict"] == {
        "complete_through_height_10": True,
        "golden_H22_H29_nonexistence": "OPEN",
        "new_compact_example": False,
        "new_exact_Vinberg_wall": True,
        "reflectivity": "OPEN",
    }
