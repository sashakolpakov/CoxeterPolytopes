#!/usr/bin/env python3
"""Independent endpoint signs, exhaustive covers, and universal escape audit.

Reverses the order of quadratic extensions used by the endpoint constructor:
the endpoint radical is the inner extension and each cap radical the outer.
The exact degree graph is a separately replayed compatibility supergraph.
No producer/search import, floating-point arithmetic, or Boolean solver.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path

import fresh_h9_split_leaf_3040_independent_20261013 as X

ROOT=Path(__file__).resolve().parents[1]
NORMALS=ROOT/'data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json'
ENDPOINTS=ROOT/'data/enumerations/fresh_h9_split_leaf_edge_endpoints_20261013.json'
GRAPH=ROOT/'data/enumerations/fresh_h9_split_leaf_3034_degree_graph_20261013.json'
COVER=ROOT/'results/fresh_h9_split_leaf_endpoint_cover_scout_20261013.json'
OUT=ROOT/'results/fresh_h9_split_leaf_endpoint_independent_20261013.json'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):return json.loads(path.read_text())
def enc(x):return {'a':str(x[0]),'b':str(x[1])}


def outer_sign(field,a,b,d):
    """Sign of a+b*sqrt(d), including dependent and zero radicals."""
    sd=X.sign(d);assert sd>=0
    sa=field.sign(a)
    if not sd:return sa
    sb=field.sign(b)
    if not sa:return sb
    if not sb or sa==sb:return sa
    return sa*field.sign(field.sub(field.mul(a,a),field.mul((d,X.Z),field.mul(b,b))))


def run():
    paths=(NORMALS,ENDPOINTS,GRAPH,COVER)
    hashes={str(p.relative_to(ROOT)):sha(p)for p in paths}
    normal,endpoint,graph,cover=map(load,paths)
    assert endpoint['source_sha256']==graph['source_sha256']==cover['normal_sha256']==sha(NORMALS)
    assert cover['endpoint_sha256']==sha(ENDPOINTS)and cover['graph_sha256']==sha(GRAPH)
    g=tuple(tuple(map(X.parse,row))for row in normal['gram'])
    assert len(g)==16 and all(len(row)==16 for row in g)
    for i in range(16):
        for j in range(16):
            expected=(F(2),F(0))if i==j else X.Z
            if i//8!=j//8:
                same=(i%8)//4==(j%8)//4
                expected=(F(-(3 if same else 8)-(22 if i%8==j%8 else 0),11),F(-(12 if same else 10),11))
            assert g[i][j]==expected
    center=(X.O,)*16;gc=X.mv(g,center)
    assert X.dot(center,gc)==(F(-64),F(-128))
    assert all(X.sign(x)<0 for x in gc)
    ed=X.parse(endpoint['endpoint_radicand']);assert X.sign(ed)>0
    ef=X.Quadratic(ed)
    caps=[]
    for i,raw in enumerate(normal['caps']):
        assert raw['index']==i
        d=X.parse(raw['radicand']);assert X.sign(d)>=0
        cv=tuple(tuple(map(X.parse,v))for v in raw['covector'])
        assert len(cv)==16 and all(len(v)==2 for v in cv)
        f=X.Quadratic(d)
        assert f.sign(f.total(cv))<0
        caps.append((d,cv))
    assert len(caps)==3040
    clauses=[];endpoint_summaries=[];total_signs=0
    for row in endpoint['rows']:
        point=tuple(tuple(map(X.parse,v))for v in row['endpoint_coefficients'])
        assert len(point)==16 and all(len(v)==2 for v in point)
        cov=tuple(ef.dot(point,tuple((a,X.Z)for a in gr))for gr in g)
        assert ef.sign(ef.dot(point,cov))==0
        assert ef.sign(ef.total(cov))<0
        side=row['seed_side'];assert side in(0,1)
        signs=[ef.sign(x)for x in cov]
        assert signs==[0 if i//8==side else-1 for i in range(16)]
        assert signs==row['seed_pairing_signs']
        p=tuple(x[0]for x in point);r=tuple(x[1]for x in point)
        gp,gr=X.mv(g,p),X.mv(g,r)
        pp,rr,pr=X.dot(p,gp),X.dot(r,gr),X.dot(p,gr)
        assert X.sign(pp)<0 and X.sign(rr)>0 and pr==X.Z
        assert X.add(pp,X.mul(ed,rr))==X.Z
        assert all(gp[i]==gr[i]==X.Z for i in range(16)if i//8==side)
        assert X.dot(r,gc)==X.Z and X.sign(X.dot(p,gc))<0
        capsigns=[]
        for d,cv in caps:
            # Cap radical is OUTER here; the constructor uses it INNER.
            a=ef.dot(point,tuple((v[0],X.Z)for v in cv))
            b=ef.dot(point,tuple((v[1],X.Z)for v in cv))
            capsigns.append(outer_sign(ef,a,b,d))
        assert capsigns==row['cap_pairing_signs']
        nonnegative=[i for i,s in enumerate(capsigns)if s>=0]
        zeros=[i for i,s in enumerate(capsigns)if s==0]
        assert nonnegative==row['nonnegative_caps'] and zeros==row['zero_caps']
        assert len(nonnegative)==14 and len(zeros)==2
        clauses.append(nonnegative);total_signs+=len(capsigns)
        endpoint_summaries.append({'side':side,'radical_sign':row['radical_sign'],
                                   'nonnegative_caps':nonnegative,'zero_caps':zeros})
    assert len(clauses)==4 and total_signs==12160
    adj=[int(x,16)for x in graph['adjacency_hex']]
    assert len(adj)==3040
    choices=set()
    for raw in itertools.product(*clauses):
        choice=tuple(sorted(set(raw)))
        if all(adj[i]>>j&1 for i,j in itertools.combinations(choice,2)):
            choices.add(choice)
    assert choices=={(i,j)for i in range(6)for j in range(6,12)}
    branches=[];union=0
    for i,j in sorted(choices):
        allowed=(adj[i]|1<<i)&(adj[j]|1<<j)
        assert allowed.bit_count()==28
        ids=[k for k in range(3040)if allowed>>k&1]
        branches.append({'minimal_cover':[i,j],'all_allowed_entries':ids})
        union|=allowed
    assert branches==cover['branches']
    ids=[i for i in range(3040)if union>>i&1]
    assert len(ids)==332 and ids==cover['all_possible_entries']
    w=[8,8,8,1,-3,-3,-3,-3,-3,-3,-3,5,5,5,5,5]
    assert cover['universal_witness']=={'integer_coefficients':w,'denominator':1}
    w=tuple((F(x),F(0))for x in w);gw=X.mv(g,w)
    square=X.dot(w,gw);time=X.dot(center,gw)
    assert square==(F(7242,11),F(-2800,11)) and X.sign(square)>0
    assert time==(F(-116),F(-232)) and X.sign(time)<0
    assert all(X.sign(x)<0 for x in gw)
    for i in ids:
        d,cv=caps[i];field=X.Quadratic(d)
        value=field.dot(tuple((x,X.Z)for x in w),cv)
        assert field.sign(value)<0,(i,'universal witness fails')
    # Negative controls: equality under a dependent radical must be zero;
    # treating a zero endpoint pairing as a positive-only clause loses a cap.
    sf=X.Quadratic((F(2),F(0)))
    assert outer_sign(sf,(X.Z,X.O),(X.neg(X.O),X.Z),(F(2),F(0)))==0
    assert outer_sign(sf,(X.neg(X.O),X.Z),(X.O,X.O),X.Z)<0
    assert all(set(r['zero_caps'])<=set(r['nonnegative_caps'])for r in endpoint_summaries)
    assert hashes=={str(p.relative_to(ROOT)):sha(p)for p in paths}
    return {'status':'PASS_INDEPENDENT_EXACT_FIXED_BANK_EXCLUSION_REVIEW',
            'source_hashes':hashes,'endpoint_signs_checked_in_reversed_radical_order':total_signs,
            'endpoint_summaries':endpoint_summaries,'cartesian_choices':14**4,
            'compatible_covers':36,'entries_per_cover':28,'union_entries':332,
            'witness_integer_coefficients':[int(x[0])for x in w],
            'witness_square':enc(square),'witness_center_pairing':enc(time),
            'strict_seed_pairings':16,'strict_cap_pairings':332,
            'scope':'The separately certified exact degree compatibility supergraph, four exact necessary endpoint clauses, exhaustive cover enumeration, and strict universal spacelike witness exclude finite-volume Coxeter completions retaining the sixteen seed facets and drawing additional facet mirrors only from this frozen bank of 3034 distinct mirrors. This audit does not exclude new mirrors, removal of seed facets, or arbitrary H9 constructions.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    result=run()
    if args.verify:assert result==load(OUT)
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items()if k not in('source_hashes','endpoint_summaries')})
