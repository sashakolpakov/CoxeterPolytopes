from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_a24_plus import (
    a24_basis_audit,
    a24_niemeier_membership,
    first_a24_patterns,
)


def test_a24_niemeier_glue_first_shell_is_exact() -> None:
    assert a24_basis_audit()
    invalid_mixed = (Cubic7(1) / 5,) * 5 + (Cubic7(0),) * 20
    assert a24_niemeier_membership(invalid_mixed) is None
    patterns = first_a24_patterns()
    assert all(a24_niemeier_membership(pattern.coordinates) == pattern.coset for pattern in patterns)
