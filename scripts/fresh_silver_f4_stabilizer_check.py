#!/usr/bin/env python3
"""Expose the additional eta mirror in the old norm-two stabilizer proof.

Standalone Fraction arithmetic from the independent audit; no project imports.
This checks one explicit involution, not the full group's reflectivity.
"""
import argparse
import json
from pathlib import Path
from scripts.fresh_silver_f4_independent_verify import (
    ROOT,Z,O,P,UNIT,ETA,read,parse,enc,add,neg,sub,scale,mul,div,sign,integral,
    dot,mv,tr,mm,pair,matrix,eye,inverse,oriented_key,
)

def encoded_matrix(a):return [[enc(x) for x in row] for row in a]
def reflection(r,gram):
    cov=mv(gram,r);q=dot(r,cov)
    return tuple(tuple(sub(O if i==j else Z,div(scale(mul(r[i],cov[j]),2),q))
        for j in range(len(r))) for i in range(len(r)))

def verify():
    ids=read('results/fresh_silver_f4_identifications_20261013.json')
    T=matrix(ids['old_carrier_to_F4_squared_scalar1']);B=matrix(ids['target_gram']);oldB=matrix(ids['old_carrier_gram'])
    original=read('data/candidates/silver_k3_norm2_infinite_stabilizer.json')
    old=matrix(read('data/candidates/silver_k3_norm2_vinberg_extension_6_plus_3sqrt2.json')['ordered_prefix_roots'])
    data=read('results/fresh_silver_f4_h9_unit_h30_20261013.json');walls=matrix(data['all_walls'])
    inside=tuple(map(parse,data['strict_interior']));escape=tuple(map(parse,data['geometry']['hyperideal_witness']['vector']))
    delta=tuple(sub(a,b) for a,b in zip(old[4],old[18]));assert pair(delta,oldB,delta)==scale(ETA,2)
    H=reflection(delta,oldB);transported=mm(T,mm(H,inverse(T)))
    r=tuple(neg(div(x,P)) for x in mv(T,delta));q=pair(r,B,r)
    assert q==ETA and all(map(integral,r))
    coefficients=tuple(div(scale(x,2),q) for x in mv(B,r));assert all(map(integral,coefficients))
    assert transported==reflection(r,B)
    assert all(integral(x) for row in transported for x in row)
    assert mm(tr(transported),mm(B,transported))==B and mm(transported,transported)==eye(10)
    fixed=mv(T,tuple(map(parse,original['wall_4_18_swap']['fixed_witness'])))
    assert mv(transported,fixed)==fixed and pair(r,B,fixed)==Z and sign(pair(fixed,B,fixed))<0
    matches=[i for i,w in enumerate(walls) if oriented_key(w)==oriented_key(r)]
    assert matches==[46] and walls[46]==r
    assert sign(pair(inside,B,r))<0 and sign(pair(escape,B,r))<0
    reflected_inside=mv(transported,inside)
    assert pair(reflected_inside,B,r)==neg(pair(inside,B,r)) and sign(pair(reflected_inside,B,r))>0
    assert all(sign(pair(r,B,w))<=0 for i,w in enumerate(walls) if i!=46)
    expected=((3,2),(6,4),(9,6),(8,6),(4,3),(2,2),(3,3),(4,2),(2,1),(1,1))
    assert r==expected
    return {'status':'OLD_NORM2_INVOLUTION_IS_EXISTING_FULL_ALPHABET_REFLECTION',
        'old_stabilizer_source':'data/candidates/silver_k3_norm2_infinite_stabilizer.json',
        'new_chamber_source':'results/fresh_silver_f4_h9_unit_h30_20261013.json',
        'root':[enc(x) for x in r],'norm':enc(q),
        'reflection_coefficients':[enc(x) for x in coefficients],
        'height':enc(div(mul(UNIT,mul(r[0],r[0])),q)),
        'wall_index':46,'interior_pairing':enc(pair(inside,B,r)),
        'escape_pairing':enc(pair(escape,B,r)),
        'old_fixed_witness_on_mirror':True,'reflection_matrix':encoded_matrix(transported),
        'scope':'H exchanges sides of this full-alphabet mirror, already present as wall46. This invalidates transfer of that generator as a chamber stabilizer. No conclusion is drawn about whether the old infinite product preserves some full chamber, or about full-host reflectivity.'}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path)
    a=p.parse_args();d=verify();print(json.dumps({k:v for k,v in d.items() if k not in ('reflection_matrix','reflection_coefficients')},sort_keys=True))
    if a.output:a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
