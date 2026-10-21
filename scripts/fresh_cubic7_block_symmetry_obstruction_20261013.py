#!/usr/bin/env python3
"""Fraction-only certificate: compatible cuts orthogonal to one swap fail."""
import argparse
import copy
import json
from pathlib import Path
import fresh_cubic7_dihedral_independent_audit_20261013 as f

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/fresh_cubic7_dihedral_h9_2theta_h30_20261013.json'
OUTPUT=ROOT/'results/fresh_cubic7_block_symmetry_obstruction_20261013.json'

def projection(v):
    out=list(v);out[7]=out[8]=f.scale(f.add(v[7],v[8]),f.F(1,2));return tuple(out)

def certificate():
    d=json.loads(SOURCE.read_text());g=tuple(f.vecdecode(r) for r in d['gram']);walls=tuple(f.vecdecode(r) for r in d['all_walls'])
    v=f.vecdecode(d['geometry']['hyperideal_witness']['vector']);p=projection(v)
    b=f.add(f.scale(f.THETA,2),f.mul(f.THETA,f.THETA));lam=f.sub(f.scale(f.THETA,2),f.TWO)
    coeffs=(f.inverse(b),lam,f.divide(lam,f.scale(f.D,2)));ids=(13,4,5)
    rhs=tuple(f.sum_field(f.mul(c,walls[i][j]) for c,i in zip(coeffs,ids)) for j in range(10))
    assert p==rhs and all(f.sign(c)>0 for c in coeffs)
    assert f.sign(f.pair(v,g,v))>0 and all(f.sign(f.pair(v,g,w))<=0 for w in walls)
    assert f.pair(p,g,walls[7])==f.ZERO
    factor=f.divide(f.pair(v,g,walls[7]),f.pair(walls[7],g,walls[7]))
    assert tuple(f.sub(x,y) for x,y in zip(v,p))==tuple(f.mul(factor,x) for x in walls[7])
    for i in ids:assert projection(walls[i])==walls[i]
    return {'status':'ALL_REAL_COMPATIBLE_CUTS_ORTHOGONAL_TO_WALL7_FAIL',
            'source':str(SOURCE.relative_to(ROOT)),'protected_escape':f.vecencode(v),
            'escape_norm':f.encoded(f.pair(v,g,v)),'orthogonal_projection':f.vecencode(p),
            'support_wall_indices_zero_based':list(ids),'positive_coefficients':f.vecencode(coeffs),
            'required_new_cut_property':'Every wall cutting this escape must have unequal coordinates7and8, hence nonzero pairing with old wall7.',
            'scope':'No full-dimensional finite-wall completion using only old-compatible added normals with coordinates7=8 can have finite volume. Other added normals are unrestricted by this argument.'}

def verify(d):
    assert d==certificate()

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=certificate()
    if a.verify:verify(json.loads(OUTPUT.read_text()))
    else:OUTPUT.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    forged=copy.deepcopy(d);forged['positive_coefficients'][0]=f.encoded(f.ONE)
    try:verify(forged)
    except AssertionError:pass
    else:raise AssertionError('forged identity accepted')
    print(json.dumps({'status':'PASS','scope':d['scope'],'forged_identity_rejected':True}))
