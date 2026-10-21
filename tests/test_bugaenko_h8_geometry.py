import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.matrix import Inertia, determinant, inertia
from coxeter_polytopes.known_examples.bugaenko import (
    BUGAENKO_VERTICES,
    bugaenko_basis_gram,
    bugaenko_diagram_basis_gram,
    bugaenko_f_vector,
    bugaenko_h8_certificate,
    bugaenko_primitive_roots,
    bugaenko_simple_roots,
    bugaenko_wall_gram,
)
from coxeter_polytopes.verify.chamber import verify_chamber_geometry


def test_reconstructed_bugaenko_diagram_exactly() -> None:
    basis_gram = bugaenko_basis_gram()
    diagram_basis_gram = bugaenko_diagram_basis_gram()
    wall_gram = bugaenko_wall_gram()
    assert determinant(basis_gram) == -2 * PHI
    assert determinant(diagram_basis_gram) == 16 - 10 * PHI
    assert inertia(basis_gram) == Inertia(8, 1, 0)
    assert inertia(basis_gram, embedding=1) == Inertia(9, 0, 0)
    primitive_roots = bugaenko_primitive_roots()
    assert determinant([*primitive_roots[:8], primitive_roots[9]]) == 1
    assert determinant(bugaenko_simple_roots()[:9]) == -3 + 2 * PHI
    assert all(wall_gram[i][i] == 2 for i in range(11))
    assert wall_gram[0][1] == wall_gram[7][8] == -PHI
    assert wall_gram[9][10] == -2 * PHI
    expected_simple_edges = {
        (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (3, 9), (5, 10)
    }
    for i in range(11):
        for j in range(i):
            if (j, i) in {(0, 1), (7, 8)}:
                assert wall_gram[j][i] == -PHI
            elif (j, i) == (9, 10):
                assert wall_gram[j][i] == -2 * PHI
            elif (j, i) in expected_simple_edges:
                assert wall_gram[j][i] == Golden(-1)
            else:
                assert wall_gram[j][i] == Golden(0)


def test_independent_bugaenko_geometry_verification() -> None:
    result = verify_chamber_geometry(bugaenko_h8_certificate())
    assert result.dimension == 8
    assert result.facet_count == 11
    assert result.finite_vertices == BUGAENKO_VERTICES
    assert len(result.finite_vertices) == 41
    assert result.cusps == ()
    assert result.hyperideal_rays == ()
    assert result.opposite_sheet_rays == ()
    assert bugaenko_f_vector() == (41, 164, 316, 374, 294, 156, 54, 11, 1)


def test_archived_bugaenko_certificate_matches_generator() -> None:
    path = Path(__file__).resolve().parents[1] / "data/known_examples/bugaenko_compact_h8.json"
    with path.open(encoding="utf-8") as handle:
        assert json.load(handle) == bugaenko_h8_certificate()
