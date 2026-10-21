#!/usr/bin/env python3
"""Independent saved-certificate counting and radical-sign review.

No producer imports, geometric projection reconstruction, orbit enumeration,
or optimizer. Only Python Fraction arithmetic and literal certificate data.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_dominant_cap_compatibility_family_20261017'
TAG=BASE+'_review'
FIRST='fresh_dominant_cap_compatibility_20261017'
WITNESS='fresh_sqrt5_scaled_core_cap_edge_repair_20261017'

def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def mul(x,y):return x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def scale(x,c):return x[0]*c,x[1]*c
def signq(x):
    a,b=x
    if not b:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    if (a>0)==(b>0):return 1 if a>0 else -1
    d=a*a-5*b*b
    return ((d>0)-(d<0))*(1 if a>0 else -1)

TAU=(F(7,2),F(1,2))
ZERO=(F(),F())
ONE=(F(1),F())
def parse(x):return tuple(map(F,x))
def signrad(a,b):
    sa,sb=signq(a),signq(b)
    if not sb:return sa
    if not sa:return sb
    if sa==sb:return sa
    return sa*signq(add(mul(a,a),neg(mul(mul(b,b),TAU))))

def outer_sign(a,b,h):
    """Sign of a+b sqrt(TAU)-sqrt(h), h>0 in the physical embedding."""
    assert signq(h)>0
    s=signrad(a,b)
    if s<=0:return -1
    constant=add(add(mul(a,a),mul(mul(b,b),TAU)),neg(h))
    return signrad(constant,scale(mul(a,b),2))

def read(tag):return json.loads((ROOT/'results'/(tag+'.json')).read_text())
def key(r):return (r['face'],tuple(r['omitted']),r['vertex'],r['provenance'])
def literal_key(s):
    return tuple(sorted(tuple((i,tuple(e))for i,e in q)for q in s['contacts']))
def encoded(x):return json.dumps(x,sort_keys=True,separators=(',',':'))

def run():
    d=read(BASE);first=read(FIRST);witness=read(WITNESS)
    alpha={i:F(a)for i,a in witness['face_multipliers']};mu={i:F(a)for i,a in witness['vertex_multipliers']}
    for v,a,b in first['fixed_witness_repair']['vertex_updates']:
        assert mu[v]==F(a);mu[v]=F(b)
    expected={};by_vertex=Counter()
    for face in d['selected_faces']:
        assert len(face['corners'])==10 and len({tuple(r['omitted'])for r in face['corners']})==10
        assert alpha[face['face']]>0
        for r in face['corners']:
            k=key(r);assert k not in expected;expected[k]=r;by_vertex[r['vertex']]+=1
    assert len(expected)==570 and len(d['selected_faces'])==57
    assert len({f['face']for f in d['selected_faces']})==57
    origins={};counts=Counter();signs={};column_seen=set();new_mu=dict(mu);forests=[]
    for fr in d['forest_records']:
        vi=fr['vertex'];stars={s['index']:s for s in fr['literal_stars']};groups={g['index']:g for g in fr['corner_orbits']}
        assert set(stars)==set(range(len(stars)))and set(groups)==set(range(len(groups)))
        assert len({literal_key(s)for s in stars.values()})==len(stars)
        owner={}
        for g in groups.values():
            ids=g['literal_stars'];assert len(ids)==len(set(ids))
            for i in ids:assert i not in owner;owner[i]=g['index']
            coeff=Counter(r['face']for r in g['origins']);assert len(coeff)==1
            assert sorted(coeff.items())==[tuple(x)for x in g['positive_face_coefficients']]
            for r in g['origins']:
                k=key(r);assert k in expected and k not in origins and r==expected[k]
                assert r['vertex']==vi and r['literal_orientations']==len(ids)
                assert r['active_embeddings']>=len(ids)
                origins[k]=g['index']
        assert set(owner)==set(stars)and sum(len(g['origins'])for g in groups.values())==fr['corner_occurrences']==by_vertex[vi]
        for s in stars.values():
            assert key(s['origin'])in {key(r)for r in groups[owner[s['index']]]['origins']}
            assert len(s['occupied_pairs'])==len({tuple(x)for x in s['occupied_pairs']})==4
            assert all(len(x)==2 and 0<=x[0]<x[1]<13 for x in s['occupied_pairs'])
            mp=s['active_embedding'];assert len(mp)==len({a for a,b in mp})==len({b for a,b in mp})==5
            assert all(signq(parse(h))>0 for h in s['squared_timelike_heights'])
        cert={}
        for r in fr['incompatibility_certificates']:
            a,b=r['stars'];assert a<b and (a,b)not in cert
            sa,sb=stars[a],stars[b];assert literal_key(sa)!=literal_key(sb)
            reason=r['reason']
            if reason=='distinct_full_stars_share_actual_quad_port_pair':
                shared=sorted(set(map(tuple,sa['occupied_pairs']))&set(map(tuple,sb['occupied_pairs'])))
                assert shared and shared==list(map(tuple,r['shared_pairs']))
            else:
                u,v=r['walls'];assert u in(0,1)and v in(0,1)and sa['contacts'][u]!=sb['contacts'][v]
                token=encoded([r['projection'],r['height_product']])
                if token not in signs:
                    x,y=map(parse,r['projection']);h=parse(r['height_product'])
                    signs[token]=(outer_sign(x,y,h),outer_sign(add(x,ONE),y,h))
                sg,sp=signs[token]
                assert sg==r['sign_product']and sp==r['sign_product_plus_one']
                if reason=='positive_actual_wall_product':assert sg>0
                else:assert reason=='actual_wall_product_in_open_minus1_zero'and sg<0<sp
            cert[a,b]=r;counts[reason]+=1
        def rejected(a,b):return tuple(sorted((a,b)))in cert
        for g in groups.values():
            all_bad=all(rejected(a,b)for a,b in combinations(g['literal_stars'],2))
            assert all_bad==g['capacity_one_proved']
        nodes={i for i,g in groups.items()if g['capacity_one_proved']}
        edges={tuple(e)for e in fr['group_incompatibility_edges']}
        assert all(a<b and a in nodes and b in nodes for a,b in edges)
        for a,b in combinations(sorted(nodes),2):
            bad=all(rejected(i,j)for i,j in product(groups[a]['literal_stars'],groups[b]['literal_stars']))
            assert bad==((a,b)in edges)
        local_best=mu.get(vi,F());cc=[]
        for ci in fr['columns']:
            col=d['columns'][ci];assert col['index']==ci and ci not in column_seen;column_seen.add(ci)
            assert col['vertex']==vi and col['vertex_cost']==1
            clique=col['orbit_clique'];assert len(clique)==len(set(clique))and set(clique)<=nodes
            assert all(tuple(sorted(e))in edges for e in combinations(clique,2))
            assert all(any(tuple(sorted((x,y)))not in edges for y in clique)for x in nodes-set(clique))
            coeff=Counter()
            for gi in clique:coeff.update(dict(groups[gi]['positive_face_coefficients']))
            assert sorted(coeff.items())==[tuple(x)for x in col['positive_face_coefficients']]
            load=sum((alpha[i]*m for i,m in coeff.items()),F())
            assert load==F(col['load'])and F(col['baseline_vertex_mass'])==mu.get(vi,F())
            assert col['violated']==(load>mu.get(vi,F()))
            local_best=max(local_best,load);cc.append(ci)
        new_mu[vi]=local_best
        assert len(cert)==fr['tested_literal_pairs']-fr['unknown_literal_pairs']
        counts.update(literal_orientations=len(stars),corner_orbits=len(groups),capacity_one_orbits=len(nodes),
                      group_edges=len(edges),columns=len(cc),corner_occurrences=fr['corner_occurrences'],forests=1)
        forests.append(dict(vertex=vi,orbits=len(groups),literal_stars=len(stars),verified_pair_certificates=len(cert),columns=cc))
    assert set(origins)==set(expected)and column_seen==set(range(len(d['columns'])))
    for k,v in counts.items():assert d['counts'][k]==v,(k,v,d['counts'][k])
    updates=[[v,str(mu.get(v,F())),str(m)]for v,m in sorted(new_mu.items())if m>mu.get(v,F())]
    assert updates==d['fixed_witness_repair']['vertex_updates']
    normalizer=sum(new_mu.values());objective=F(witness['objective_numerator']);ratio=objective/normalizer
    repair=d['fixed_witness_repair']
    assert normalizer==F(repair['normalizer'])and ratio==F(repair['bound'])and repair['bound_below12']==(ratio<12)
    paths=[Path(__file__),ROOT/'scripts'/(BASE+'.py'),ROOT/'results'/(BASE+'.json'),ROOT/'results'/(BASE+'.md'),
           ROOT/'results'/(FIRST+'.json'),ROOT/'results'/(WITNESS+'.json')]
    return dict(status='PASS_INDEPENDENT_CORNER_CLIQUE_COUNT_AND_RADICAL_SIGN_REVIEW',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        counts=dict(counts),independent_radical_expressions=len(signs),forests=forests,
        radical_sign_checks=[dict(expression=json.loads(k),sign_product=v[0],sign_product_plus_one=v[1])for k,v in sorted(signs.items())],
        whole_face_corner_occurrences=len(expected),vertex_updates=updates,normalizer=str(normalizer),bound=str(ratio),
        proof_review='Complete whole-component embeddings give full finite-forest orientation orbits. Identical literal stars are one local occurrence; separate actual prism vertices give global corner multiplicities. Universal pair exclusions and internal capacity one justify every clique charge. Future-side radical signs and identical-normal exceptions are correct.',
        scope='Independent literal certificate consistency and Fraction-only re-evaluation of82 distinct radical expressions used by every saved geometric exclusion. No producer arithmetic imports. Does not independently reconstruct6010 covectors/projections, enumerate raw embeddings, or prove completeness of the underlying raw23 geometric catalogue. Those steps received a source/proof review only. No optimizer or new dimension claim.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(run(),sort_keys=True));dest=ROOT/'results'/(TAG+'.json')
    if a.verify:assert out==json.loads(dest.read_text())
    else:dest.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(out['status'],out['counts'],'radical expressions',out['independent_radical_expressions'],flush=True)
