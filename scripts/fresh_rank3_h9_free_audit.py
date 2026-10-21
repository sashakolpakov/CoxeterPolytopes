#!/usr/bin/env python3
"""Independent root and cross-angle replay for the rank-four cutter box."""
import argparse
import json
from pathlib import Path
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,add,neg,sub,scale,mul,div,sign,integral,total,dot,mv,pair,matrix,inverse,
)
from scripts.fresh_multibasepoint_h9_verify import ROOT,read
from scripts.fresh_rank3_h9_newplane_verify import integer_sign

def ip(x):
    assert integral(x)
    return int(x[0]),int(x[1])
def plus(x,y):return x[0]+y[0],x[1]+y[1]
def times(x,y):return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def kscale(x,k):return x[0]*k,x[1]*k
def isum(values):
    a=b=0
    for x,y in values:a+=x;b+=y
    return a,b

def check_cut(c,gram,einv,inside,escape,bounds,restricted,oldnorms):
    q=c['norm'];assert q in (1,2)
    r=tuple(map(parse,c['root']));x,t,u,z=map(parse,c['parameters'])
    assert all(map(integral,r)) and all(map(integral,(x,t,u,z)))
    assert r==(x,*(add(mul(t,einv[k][4]),mul(u,einv[k][7])) for k in range(8)),z)
    assert sign(x)>0 and sign(t)>0 and sign(u)>0 and sign(z)>=0
    for parameter,name in ((t,'t'),(u,'u'),(z,'scalar')):assert sign(sub(scale(O,bounds[name]),parameter))>=0
    xi,ti,ui,zi=map(ip,(x,t,u,z))
    energy=isum((kscale(times(ti,ti),12),kscale(times(ti,ui),18),kscale(times(ui,ui),8),times(zi,zi)))
    assert plus(energy,kscale(times((0,2),times(xi,xi)),-1))==(q,0)
    assert pair(r,gram,r)==scale(O,q)
    assert all(integral(div(scale(cov,2),scale(O,q))) for cov in mv(gram,r))
    assert sign(pair(r,gram,inside))<0
    assert pair(r,gram,escape)==parse(c['escape_pairing']) and sign(pair(r,gram,escape))>0
    assert div(mul(scale(P,2),mul(x,x)),scale(O,q))==parse(c['basepoint_height'])
    parameters=(xi,ti,ui,zi)
    allowed=((0,0),(1,0),(2,0),(3,0),(1,1),(2,1))
    for columns,norm in zip(restricted,oldnorms,strict=True):
        product=isum(times(a,b) for a,b in zip(parameters,columns,strict=True))
        assert integer_sign(*product)<=0
        sq4=kscale(times(product,product),4);den=q*norm
        if integer_sign(sq4[0]-4*den,sq4[1])<0:
            assert sq4 in tuple(kscale(v,den) for v in allowed)
    return r

def setup():
    data=read('results/fresh_rank4_h9_free_search_20261013.json');base=read(data['source'])
    assert data['source_case']==0 and data['fundamental_weight_indices']==[4,7]
    gram=matrix(base['gram']);case=base['cases'][0];walls=matrix(case['all_walls'])
    inside=tuple(map(parse,case['strict_interior']));escape=tuple(map(parse,case['hyperideal_witness']['vector']))
    einv=inverse(tuple(tuple(gram[i][j] for j in range(1,9)) for i in range(1,9)))
    assert (einv[4][4],einv[4][7],einv[7][7])==(scale(O,12),scale(O,9),scale(O,8))
    oldnorms=tuple(pair(r,gram,r) for r in walls);assert all(q in (O,scale(O,2)) for q in oldnorms)
    restricted=tuple(tuple(map(ip,(mul(scale(P,-2),r[0]),r[5],r[8],r[9]))) for r in walls)
    return data,gram,einv,inside,escape,restricted,tuple(int(q[0]) for q in oldnorms)

def verify():
    d,gram,einv,inside,escape,restricted,oldnorms=setup()
    roots=[check_cut(c,gram,einv,inside,escape,d['bounds'],restricted,oldnorms) for c in d['cuts']]
    assert len(set(roots))==len(roots)==335
    counts={q:sum(c['norm']==q for c in d['cuts']) for q in (1,2)}
    assert counts=={1:86,2:249}
    assert counts=={r['norm']:r['compatible_roots'] for r in d['censuses']}
    return {'status':'INDEPENDENT_335_RANK4_CUTTERS_PASS','cuts':335,'norm1':86,'norm2':249,
        'base_walls':len(restricted),'exact_cross_pair_gates':335*len(restricted),'source_bounds':d['bounds'],
        'scope':'All accepted roots, original-cone compatibility, integrality, heights, escape cuts and Coxeter cross-pairs independently replayed. Mutual compatibility among the335 cuts and full parameter-box completeness are separate questions.'}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);a=p.parse_args();result=verify();print(json.dumps(result,sort_keys=True))
    if a.output:a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
