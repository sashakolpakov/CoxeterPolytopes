import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.geometric.k3 import K3_TRACE_SCALAR
from coxeter_polytopes.search.bugaenko_vinberg import (
    _fraction_inverse,
    _mat_vec,
    _pair,
)
from coxeter_polytopes.search.h11_vinberg import (
    _standard_e8_gram,
    standard_e8_roots_from_simple_reflections,
)
from coxeter_polytopes.search.h19_vinberg import (
    D,
    admissible_control_coordinates,
    build_h19_vinberg_transcript,
    h19_basis_gram,
    height_zero_roots,
    run_h19_vinberg_search,
)
from coxeter_polytopes.verify.chamber import verify_hyperideal_extreme_ray
from coxeter_polytopes.verify.vinberg import verify_h19_vinberg_transcript


TRANSCRIPT = (
    Path(__file__).parents[1] / "data" / "candidates" / "h19_vinberg_prefix.json"
)
HYPERIDEAL_WITNESS = (0, 1, 2, 3, 4, 5, 8, 9, 13)


def test_control_scan_is_exhaustive_through_bound_five() -> None:
    assert admissible_control_coordinates(5) == (
        Golden(1),
        PHI,
        1 + PHI,
        2 * PHI,
        2 + PHI,
        1 + 2 * PHI,
    )


def test_control_diagonalization_height_zero_and_bounded_counts() -> None:
    assert len(standard_e8_roots_from_simple_reflections()) == 240
    assert len(height_zero_roots()) == 242
    search = run_h19_vinberg_search(5)
    assert search.control_vector == (Golden(2), -PHI, *(Golden(0) for _ in range(8)))
    assert search.cone_root_count == 242
    assert len(search.cone_roots) == 9
    assert [len(level.roots) for level in search.levels] == [2, 1, 2, 1, 0, 7]
    assert len(search.roots) == 22
    assert [level.height for level in search.levels] == sorted(
        level.height for level in search.levels
    )
    assert search.levels[0].height == D / 4
    identity_eigenspace_square = K3_TRACE_SCALAR * (-2 * D)
    assert identity_eigenspace_square == Golden(14, 32) / 5
    assert identity_eigenspace_square.sign(0) > 0
    assert identity_eigenspace_square.trace() == 12


def test_archived_bounded_transcript_is_reproducible() -> None:
    archived = json.loads(TRANSCRIPT.read_text())
    assert archived == build_h19_vinberg_transcript(5)
    assert len(verify_h19_vinberg_transcript(archived).roots) == 22


def test_exact_hyperideal_extreme_ray_blocks_vinberg_stopping() -> None:
    search = run_h19_vinberg_search(5)
    rho = _mat_vec(_fraction_inverse(_standard_e8_gram()), [Golden(1)] * 8)
    inward = (Golden(0), Golden(1), *rho)
    interior = tuple(
        search.control_vector[i] + inward[i] / 10 for i in range(10)
    )
    assert _pair(interior, h19_basis_gram(), interior).sign(0) < 0
    assert all(
        _pair(interior, h19_basis_gram(), root).sign(0) < 0
        for root in search.roots
    )
    witness = verify_hyperideal_extreme_ray(
        search.roots, h19_basis_gram(), HYPERIDEAL_WITNESS, dimension=9
    )
    assert witness.norm == Golden(115, 0) / 2 - Golden(71, 0) * PHI / 2
    assert witness.norm.sign(0) > 0
    assert all(pairing.sign(0) <= 0 for pairing in witness.wall_pairings)
