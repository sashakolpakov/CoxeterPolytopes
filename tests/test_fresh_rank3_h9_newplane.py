"""Independent square, array-filter, and cutter-certificate regressions."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import random
import numpy as np
import pytest
from scripts.fresh_rank3_h9_newplane_verify import integral_square_root,integer_sign
from scripts.fresh_h4_blocks_independent_verify import Z,O,P,add,sub,scale,mul,div,sign
from scripts.fresh_rank4_h9_free_search_20261013 import positive
from scripts.fresh_rank3_h9_free_audit import setup,check_cut

def test_integer_square_test_finds_every_small_square():
    for a,b in product(range(-25,26),repeat=2):
        root=integral_square_root(a*a+b*b,2*a*b+b*b)
        assert root is not None
        c,d=root;assert (c*c+d*d,2*c*d+d*d)==(a*a+b*b,2*a*b+b*b)

def test_integer_square_test_rejects_small_nonsquares():
    assert integral_square_root(2,0) is None
    assert integral_square_root(0,1) is None
    assert integral_square_root(-1,0) is None

def test_array_sign_matches_independent_field_arithmetic():
    rng=random.Random(20261013)
    pairs=[(rng.randrange(-100000000,100000001),rng.randrange(-100000000,100000001)) for _ in range(1000)]+[(0,0)]
    a=np.array([v[0] for v in pairs],dtype=np.int64);b=np.array([v[1] for v in pairs],dtype=np.int64)
    assert positive(a,b).tolist()==[sign((F(x),F(y)))>0 for x,y in pairs]

def test_free_slice_array_polynomials_match_direct_field_arithmetic():
    rng=random.Random(20261013)
    for _ in range(300):
        q=rng.choice((1,2));ta,tb,ua,ub,za,zb=[rng.randrange(0,300) for _ in range(6)]
        t,u,z=(F(ta),F(tb)),(F(ua),F(ub)),(F(za),F(zb))
        ea=12*(ta*ta+tb*tb)+18*(ta*ua+tb*ub)+8*(ua*ua+ub*ub)+za*za+zb*zb-q
        eb=12*(2*ta*tb+tb*tb)+18*(ta*ub+tb*ua+tb*ub)+8*(2*ua*ub+ub*ub)+2*za*zb+zb*zb
        la=12*ta+5*ua+3*ub-za+zb;lb=12*tb+3*ua+8*ub+za
        ga=la*la+lb*lb-8*eb;gb=2*la*lb+lb*lb-8*(ea+eb)
        energy=sub(add(add(scale(mul(t,t),12),scale(mul(t,u),18)),add(scale(mul(u,u),8),mul(z,z))),scale(O,q))
        twice=add(add(scale(t,12),mul(add(scale(O,5),scale(P,3)),u)),mul(sub(P,O),z))
        gate=sub(mul(twice,twice),mul(scale(P,8),energy))
        assert energy==(ea,eb) and twice==(la,lb) and gate==(ga,gb)
        assert integer_sign(ea+eb,-eb)==sign((energy[0]+energy[1],-energy[1]))

def test_forged_free_slice_cut_norm_rejected():
    d,gram,einv,inside,escape,restricted,norms=setup();cut=deepcopy(d['cuts'][0]);cut['norm']=3-cut['norm']
    with pytest.raises(AssertionError):check_cut(cut,gram,einv,inside,escape,d['bounds'],restricted,norms)

def test_forged_rank4_projection_identity_rejected():
    from scripts.fresh_rank3_h9_rank4_obstruction import verify,OUTPUT
    import json
    payload=json.loads(OUTPUT.read_text())
    payload['identities'][0]['nonnegative_coefficients'][0]['a']='0'
    with pytest.raises(AssertionError):verify(payload)
