#!/usr/bin/env python3
"""Stdlib Fraction full census replay, with independent unit balancing."""
import argparse
from functools import lru_cache
import json
import time

import fresh_h9_glued_host_20261013_independent as a
from fresh_h9_glued_host_20261013_independent import *

ORIGINAL_ALPHABET = a.alphabet
UNITS = (sub(T,O),sub(sub(mul(sq(T),T),scale(T,3)),O),
         sub(add(neg(mul(sq(T),T)),scale(T,3)),O))
assert all(integral(u) and integral(inv(u)) and mul(u,inv(u))==O for u in UNITS)


@lru_cache(None)
def balanced(caps):
    if any(sign(c,e)<0 for e,c in enumerate(caps)):
        return ()
    if any(sign(c,e)==0 for e,c in enumerate(caps)):
        return (Z,)
    # Distinct from the producer's least-squares/nearest-lattice heuristic:
    # greedily try each of six unit steps. This only proposes a bijection.
    transformed=caps;unit=O
    steps=tuple(UNITS)+tuple(inv(u)for u in UNITS)
    factors=[tuple(real(sq(u),e)for e in range(4))for u in steps]
    numerical=[real(c,e)for e,c in enumerate(caps)]
    for _ in range(80):
        old=sum(numerical)
        scores=[sum(v/f for v,f in zip(numerical,fs))for fs in factors]
        j=min(range(len(scores)),key=lambda j:scores[j])
        if scores[j]>=old*(1-1e-12):break
        unit=mul(unit,steps[j]);numerical=[v/f for v,f in zip(numerical,factors[j])]
    if sign(unit)<0:unit=neg(unit)
    assert integral(unit) and integral(inv(unit))
    transformed=tuple(divide(c,sq(unit))for c in caps)
    result=tuple(mul(unit,y)for y in ORIGINAL_ALPHABET(transformed))
    assert all(integral(x)and sign(x)>=0 and all(sign(sub(c,sq(x)),e)>=0
               for e,c in enumerate(caps))for x in result)
    return result


def key(r):
    first=next(x for x in r if x!=Z)
    return sign(first),tuple(divide(x,first)for x in r)


def run(path):
    start=time.monotonic();data=json.loads(path.read_text());host=Host(2)
    # Calibration against the older unbalanced independent enumeration.
    for q in POOL:
        assert set(balanced((q,)*4))==set(ORIGINAL_ALPHABET((q,)*4))
    a.alphabet=balanced
    controls=host.controls(data['height_bound'])
    saved={(dec(r['height']),dec(r['norm']),dec(r['time'])):r for r in data['levels']}
    assert len(controls)==len(saved)==data['control_count'] and set(controls)==set(saved)
    base=json.loads((ROOT/data['source']).read_text())
    old=tuple(vec(r)for r in base['all_walls']);oldkeys={key(r)for r in old}
    previous=json.loads((ROOT/data['previous_source']).read_text())
    pwalls=tuple(vec(r)for r in previous['all_walls']);pkeys={key(r)for r in pwalls}
    target=vec(previous['geometry']['vector'])
    dictionary={};count=0;mixed_count=0
    for h,q,x in controls:
        roots,imprimitive,blocks,scalars=host.candidates(x,q)
        row=saved[h,q,x]
        assert set(roots)=={vec(r)for r in row['candidates']}
        assert len(roots)==len(row['candidates'])==row.get('primitive_candidates',0)
        assert imprimitive==row.get('imprimitive',0) and blocks==row['block_vectors']
        count+=len(roots)
        mixed_count+=sum(any(v!=Z for v in r[1:5])and any(v!=Z for v in r[5:9])for r in roots)
        for r in roots:
            k=key(r)
            if k in oldkeys or not host.compatible(r,old):continue
            dictionary[k]=(r,h,q)
    records={key(vec(r['vector'])):r for r in data['compatible_dictionary']}
    assert records.keys()==dictionary.keys()
    for k,(r,h,q)in dictionary.items():
        row=records[k]
        assert (vec(row['vector']),dec(row['height']),dec(row['norm']))==(r,h,q)
        assert row['both_blocks_active']==(any(v!=Z for v in r[1:5])and any(v!=Z for v in r[5:9]))
        assert row['already_in_old22']==(k in pkeys)
        assert row['compatible_with_old22']==(k in pkeys or host.compatible(r,pwalls))
        assert host.pair(r,target)==dec(row['old22_escape_pairing'])
    chosen=[];rejected=[];walls=list(old)
    for i,row in enumerate(data['compatible_dictionary']):
        r=vec(row['vector'])
        if host.compatible(r,walls):chosen.append(i);walls.append(r)
        else:rejected.append(i)
    assert chosen==data['selected_indices'] and rejected==data['rejected_indices']
    assert walls==list(map(vec,data['all_walls']))
    geometry=host.check_geometry(data)
    facepath=ROOT/'results/fresh_h9_quartic_mixed_20261013_h5.json'
    face=json.loads(facepath.read_text())
    expected_face={(*r[:5],r[9])for r in walls if all(x==Z for x in r[5:9])}
    assert expected_face=={vec(r)for r in face['all_walls']}
    face_geometry=Host(1).check_geometry(face)
    return {'status':'PASS_COMPLETE_INDEPENDENT_MIXED_QUARTIC_CENSUS',
            'source':str(path.relative_to(ROOT)), 'controls':len(controls),
            'primitive_candidates':count,'both_active_candidates':mixed_count,
            'old14_compatible_dictionary':len(dictionary),
            'height_ordered_selected':len(chosen),'rejected':len(rejected),
            'h9_geometry':geometry,'h5_geometry':face_geometry,
            'alphabet_calibrations':len(POOL),
            'arithmetic':'Standalone Fraction field. Independent greedy integral-unit balancing, trace-ellipsoid loops, and full Cartesian F4 covector enumeration. No producer imports.',
            'scope':'Full bounded candidate and compatible dictionary replay, not a full-host nonreflectivity theorem.',
            'seconds':round(time.monotonic()-start,3)}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=ROOT/'results/fresh_h9_quartic_mixed_20261013_h1000.json')
    p.add_argument('--output',type=Path,default=ROOT/'results/fresh_h9_quartic_mixed_20261013_independent.json');p.add_argument('--verify',action='store_true');args=p.parse_args()
    result=run(args.source.resolve())
    if args.verify:
        old=json.loads(args.output.read_text());assert {k:v for k,v in old.items()if k!='seconds'}=={k:v for k,v in result.items()if k!='seconds'}
    else:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result)
