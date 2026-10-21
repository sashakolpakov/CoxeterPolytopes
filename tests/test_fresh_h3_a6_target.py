"""Targeted H3+A6 arithmetic and certificate regressions."""
from copy import deepcopy
from fractions import Fraction as F
import random
import pytest
from scripts.fresh_h3_a6_target_verify import (
    O,P,add,sub,scale,mul,sign,pair,setup,obstruction_payload,read,
)

def test_neighbor_norm_polynomial_matches_independent_host_gram():
    rng=random.Random(20261013)
    for name in ('phi','2phi'):
        _,alpha,B,basis=setup(name);e,h,k,u=basis[:4]
        for _ in range(60):
            x,s,t,z=[(F(rng.randrange(-20,21)),F(rng.randrange(-20,21))) for _ in range(4)]
            r=tuple(add(add(mul(x,e[j]),mul(s,h[j])),add(mul(t,k[j]),mul(z,u[j]))) for j in range(10))
            ss,st,tt,zz=mul(s,s),mul(s,t),mul(t,t),mul(z,z)
            ea=6*(ss[0]+ss[1])+4*(st[0]+2*st[1])+2*(tt[0]+tt[1])+6*zz[0]
            eb=6*(ss[0]+2*ss[1])+4*(2*st[0]+3*st[1])+2*(tt[0]+2*tt[1])+6*zz[1]
            assert pair(r,B,r)==sub((ea,eb),mul(alpha,mul(x,x)))

def test_closed_slice_and_both_unbounded_projection_certificates():
    result=obstruction_payload()
    assert result['status']=='EXACT_MOD3_AND_UNBOUNDED_NEIGHBOR_SLICE_OBSTRUCTIONS'
    assert all(not row['mod3_root_solutions'] for row in result['hosts'])

@pytest.mark.parametrize('name', ['phi','2phi'])
def test_forged_targeted_root_norm_is_rejected(monkeypatch,name):
    import scripts.fresh_h3_a6_target_verify as audit
    path='results/fresh_h3_a6_target_'+name+'_neighbor_20261013.json'
    payload=deepcopy(read(path));payload['cuts'][0]['norm']=3
    original=audit.read
    monkeypatch.setattr(audit,'read',lambda p:payload if p==path else original(p))
    with pytest.raises(AssertionError):audit.verify_host(name)
