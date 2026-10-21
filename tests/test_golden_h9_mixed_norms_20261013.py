"""Check the newly admitted integral reflections and the old default API."""
import json
from pathlib import Path
import pytest
from coxeter_polytopes.exact.golden import Golden as G, PHI, ZERO
from coxeter_polytopes.search.bugaenko_vinberg import _pair
from coxeter_polytopes.search.og10_golden_h9_vinberg import (
    _enumerate_height_level, _positive_definite_enumeration_data,
    h9_field_gram,
)


def test_first_norm_four_wall_is_the_known_pell_reflection():
    inverse, schur = _positive_definite_enumeration_data()
    level = _enumerate_height_level(PHI-1, (), inverse, schur, norm=4)
    root = (PHI-1, *(ZERO for _ in range(8)), PHI)
    assert level.roots == (root,)
    assert level.height == (PHI-1)/2
    assert _pair(root, h9_field_gram(), root) == 4
    # Compute its action on the x/u plane, rather than just its norm.
    ex=(G(1),*(ZERO for _ in range(9)))
    eu=(*(ZERO for _ in range(9)),G(1))
    reflection=lambda v:tuple(v[i]-_pair(root,h9_field_gram(),v)*root[i]/2 for i in range(10))
    assert reflection(ex) == (PHI,*(ZERO for _ in range(8)),PHI)
    assert reflection(eu) == (-G(1),*(ZERO for _ in range(8)),-PHI)


def test_default_norm_two_enumerator_keeps_its_original_height():
    inverse, schur = _positive_definite_enumeration_data()
    level = _enumerate_height_level(G(1), (), inverse, schur)
    assert len(level.roots) == 2
    assert level.height == PHI
    assert all(_pair(r,h9_field_gram(),r)==2 for r in level.roots)


def test_unsupported_norm_is_rejected_before_search():
    inverse, schur = _positive_definite_enumeration_data()
    with pytest.raises(ValueError,match='unsupported'):
        _enumerate_height_level(G(1), (), inverse, schur, norm=1)


def test_odd_a1_prefix_really_contains_integral_norm_one_reflections():
    path=Path(__file__).resolve().parents[1]/'results/golden_h9_odd_a1_full_norms_h30_20261013.json'
    payload=json.loads(path.read_text())
    assert payload['norms']==[1,2]
    gram=h9_field_gram()
    gram=tuple(tuple(G(1) if i==j==9 else gram[i][j] for j in range(10)) for i in range(10))
    walls=tuple(tuple(G.from_json(v) for v in r) for r in payload['all_walls'])
    positive_norm_one=[r for r in walls if r[0] and _pair(r,gram,r)==1]
    assert len(positive_norm_one)==1
    for root in positive_norm_one:
        assert root[0]==1+PHI
        for i in range(10):
            e=tuple(G(int(i==j)) for j in range(10))
            assert (2*_pair(root,gram,e)).is_integral
