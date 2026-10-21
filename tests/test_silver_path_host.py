from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_path_host import (
    first_cut_roots,
    pair,
    silver_path_first_cut_audit,
)


def test_silver_path_first_cut_geometry() -> None:
    audit = silver_path_first_cut_audit()
    assert audit.determinant == Silver(-40, -42)
    assert audit.selected_inertia == Inertia(22, 1, 0)
    assert audit.conjugate_inertia == Inertia(23, 0, 0)
    assert (audit.finite_vertex_count, audit.cusp_count, audit.hyperideal_ray_count) == (6, 0, 97)


def test_silver_path_first_cut_roots() -> None:
    roots = first_cut_roots()
    assert [pair(root, root) for root in roots] == [Silver(2), Silver(2)]
    assert pair(*roots) == Silver(-2, -1)
