#!/usr/bin/env python3
"""Exact bounded search in the new two-fundamental-weight binary slice."""
import argparse
import json
import time
from scripts.fresh_multibasepoint_h9_escape_plane import (
    ROOT,G,PHI,O,Z,read,matrix,mv,dot,enc,local_obstructions,golden_square_root,coxeter_pair_value,
)
from scripts.fresh_rank3_h9_search import alphabet
from scripts.fresh_h4_blocks_h9_search_20261013 import inverse
from scripts.fresh_rank3_h9_cones import compatible,escape


def run(limit):
    start=time.monotonic();source='results/fresh_rank3_h9_cones_20261013.json';base=read(source);case=base['cases'][0]
    gram=matrix(base['gram']);walls=matrix(case['all_walls']);normals=tuple(mv(gram,r) for r in walls)
    norms=tuple(dot(r,n) for r,n in zip(walls,normals));alpha=2*PHI
    e8=tuple(tuple(gram[i][j] for j in range(1,9)) for i in range(1,9));ei=inverse(e8)
    direction=tuple(2*ei[i][4]+2*PHI*ei[i][7] for i in range(8))+(1+PHI,)
    z=(Z,*direction);beta=dot(z,mv(gram,z));assert beta==G(82,107)
    assert direction[-1]*(2-PHI)==1 # The spatial direction is primitive and the slice is saturated.
    v=tuple(map(G.from_json,case['hyperideal_witness']['vector']));inside=tuple(map(G.from_json,case['strict_interior']))
    swap=matrix(base['standard_block_swap_involution']);assert set(mv(swap,r) for r in walls)==set(walls)
    escape_upper=dot(z,mv(gram,v))/(alpha*v[0])
    lower=max((dot(z,n)/(alpha*r[0]) for r,n in zip(walls,normals) if r[0]),default=G(0))
    assert lower<escape_upper and lower>0
    cuts=[];records=[]
    for q in (O,G(2)):
        obstructions=local_obstructions(beta,alpha,q);assert not obstructions
        ts=alphabet(limit,q/beta)
        stats={'norm':q.to_json(),'alphabet_size_including_zero':len(ts),'parity_admissible':0,'integral_square_tests':0,'integral_solutions':0,'compatible_escape_cuts':0,'coxeter_cuts':0,'compatible_swap_orbits':0}
        for t in ts:
            if not t:continue
            a,b=int(t.a),int(t.b)
            if a%2 or b%2!=(1 if q==1 else 0):continue
            stats['parity_admissible']+=1
            square=(beta*t*t-q)/alpha
            if not square.is_integral or square.sign()<=0:continue
            stats['integral_square_tests']+=1;x=golden_square_root(square)
            if x is None or not x.is_integral:continue
            if x.sign()<0:x=-x
            stats['integral_solutions']+=1
            if x<lower*t or not x<escape_upper*t:continue
            r=(x,*(t*y for y in direction));products=tuple(dot(r,n) for n in normals)
            assert all(p.sign()<=0 for p in products) and dot(r,mv(gram,v)).sign()>0
            stats['compatible_escape_cuts']+=1
            if not all(coxeter_pair_value(4*p*p/(q*oldq)) for p,oldq in zip(products,norms)):continue
            stats['coxeter_cuts']+=1
            assert dot(r,mv(gram,r))==q and all(c.is_integral for c in r)
            assert all((2*c/q).is_integral for c in mv(gram,r)) and dot(r,mv(gram,inside)).sign()<0
            mate=mv(swap,r);orbit=tuple(dict.fromkeys((r,mate)));swap_ok=all(compatible(r,s,gram) for i,r in enumerate(orbit) for s in orbit[:i])
            if swap_ok:stats['compatible_swap_orbits']+=1
            cuts.append({'root':enc(r),'norm':q.to_json(),'x':x.to_json(),'t':t.to_json(),
                'escape_pairing':dot(r,mv(gram,v)).to_json(),'wall_pairings':enc(products),
                'height':(alpha*x*x/q).to_json(),'swap_root':enc(mate),'swap_pairing':dot(r,mv(gram,mate)).to_json(),'compatible_swap_orbit':swap_ok})
        records.append(stats);print(json.dumps(stats),flush=True)
    cuts.sort(key=lambda c:G.from_json(c['height']))
    cases=[]
    for i,cut in enumerate(cuts):
        r=tuple(map(G.from_json,cut['root']));mate=tuple(map(G.from_json,cut['swap_root']))
        additions=tuple(dict.fromkeys((r,mate))) if cut['compatible_swap_orbit'] else (r,)
        new=tuple((*walls,*additions));witness=escape(new,gram,20261013+i)
        cases.append({'cut_index':i,'wall_count':len(new),'added_walls':[enc(r) for r in additions],
            'hyperideal_witness':witness,'status':'EXACT_SPACELIKE_ESCAPE' if witness else 'VOLUME_UNDECIDED'})
        print(json.dumps({'cut':i,'walls':len(new),'escape':witness['norm'] if witness else None}),flush=True)
    out={'source':source,'source_case':0,'alpha':alpha.to_json(),'gram':base['gram'],
        'direction':enc(direction),'direction_norm':beta.to_json(),'fundamental_weight_indices':[4,7],
        'positive_t_bound':limit,'completed_parameter_box':True,'equation':'(82+107phi)*t^2-2phi*x^2=q; q in {1,2}',
        'compatibility_lower_ratio':lower.to_json(),'escape_upper_ratio':escape_upper.to_json(),
        'norm_searches':records,'cuts':cuts,'cases':cases,'seconds':round(time.monotonic()-start,3),
        'scope':'Complete stated bounded binary slice; outside the previously obstructed ternary slice. No full-host exhaustion or reflectivity conclusion.'}
    (ROOT/'results/fresh_rank3_h9_newplane_20261013.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'cuts':len(cuts),'seconds':out['seconds']}),flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--limit',type=int,default=1000000);a=p.parse_args();run(a.limit)
