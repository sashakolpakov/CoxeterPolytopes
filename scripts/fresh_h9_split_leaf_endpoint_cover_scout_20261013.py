#!/usr/bin/env python3
"""Solver-free endpoint-cover enumeration and an exact universal cone escape.

The endpoint lists and safe compatibility graph have separate exact proofs.
This replay imports neither producer and checks its own geometry from normals.
"""
import argparse
import hashlib
import itertools
import json
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
GRAPH=ROOT/'data/enumerations/fresh_h9_split_leaf_3034_degree_graph_20261013.json'
NORMALS=ROOT/'data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json'
ENDPOINTS=ROOT/'data/enumerations/fresh_h9_split_leaf_edge_endpoints_20261013.json'
OUT=ROOT/'results/fresh_h9_split_leaf_endpoint_cover_scout_20261013.json'
Z=(Q(0),Q(0));O=(Q(1),Q(0))


def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def sub(x,y):return add(x,neg(y))
def mul(x,y):return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def scale(x,c):return x[0]*c,x[1]*c
def sign(x):
    a,b=x;r=2*a+b
    if not b:return (r>0)-(r<0)
    if not r:return (b>0)-(b<0)
    if (r>0)==(b>0):return 1 if r>0 else -1
    return ((r>0)-(r<0))*((r*r>5*b*b)-(r*r<5*b*b))
def parse(x):return Q(x['a']),Q(x['b'])
def enc(x):return {'a':str(x[0]),'b':str(x[1])}
def dot(x,y):
    value=Z
    for a,b in zip(x,y,strict=True):value=add(value,mul(a,b))
    return value
def mv(g,v):return tuple(dot(row,v)for row in g)
def qsign(a,b,d):
    assert sign(d)>=0
    sa,sb=sign(a),sign(b)
    if not sign(d)or not sb:return sa
    if not sa or sa==sb:return sb
    delta=sign(sub(mul(a,a),mul(d,mul(b,b))))
    return 0 if not delta else sa if delta>0 else sb
def bits(m):
    while m:
        bit=m&-m;yield bit.bit_length()-1;m-=bit
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def enumerate_covers(clauses,adj):
    compatible=set()
    for choice in itertools.product(*clauses):
        unique=tuple(sorted(set(choice)))
        if all(adj[i]>>j&1 for i,j in itertools.combinations(unique,2)):
            compatible.add(unique)
    minimal=[]
    for choice in sorted(compatible,key=lambda c:(len(c),c)):
        if not any(set(c)<=set(choice)for c in minimal):minimal.append(choice)
    branches=[];union=0
    for cover in minimal:
        allowed=(1<<len(adj))-1
        for i in cover:allowed&=adj[i]|1<<i
        assert all(allowed>>i&1 for i in cover)
        union|=allowed
        branches.append({'minimal_cover':list(cover),'all_allowed_entries':list(bits(allowed))})
    return sorted(compatible),branches,union


def check_witness(normal,record,allowed):
    gram=tuple(tuple(map(parse,row))for row in normal['gram'])
    # Derive the seed Gram directly rather than trusting archived covectors.
    for i in range(16):
        for j in range(16):
            expected=(Q(2),Q(0))if i==j else Z
            if i//8!=j//8:
                same=(i%8)//4==(j%8)//4
                expected=(Q(-(3 if same else 8)-(22 if i%8==j%8 else 0),11),Q(-(12 if same else 10),11))
            assert gram[i][j]==expected
    a,b=(Q(1,2),Q(1)),(Q(-5,22),Q(1,11))
    assert sign(a)>0 and sign(b)<0 and sign(add(O,scale(a,2)))>0 and sign(add(O,scale(b,2)))>0
    # On symmetric/antisymmetric U,V combinations, the nonzero eigenvalues
    # are -8a,-8b,4+8a,4+8b and six copies of 4. Six directions are radical.
    signature=[9,1,6]
    raw=record['integer_coefficients'];den=record['denominator']
    assert len(raw)==16 and all(type(x)is int for x in raw)and type(den)is int and den>0
    w=tuple((Q(x,den),Q(0))for x in raw);gw=mv(gram,w)
    norm=dot(w,gw);center_pair=dot((O,)*16,gw);center_norm=dot((O,)*16,mv(gram,(O,)*16))
    assert sign(norm)>0 and sign(center_pair)<0 and sign(center_norm)<0
    assert all(sign(x)<0 for x in gw)
    assert norm==(Q(7242,11),Q(-2800,11))
    epsilon=Q(1)
    while True:
        strict_w=tuple(add(x,scale(O,epsilon))for x in w)
        strict_gw=mv(gram,strict_w);strict_norm=dot(strict_w,strict_gw)
        strict_center=dot((O,)*16,strict_gw)
        if sign(strict_norm)>0 and sign(strict_center)<0:break
        epsilon/=2
    assert all(sign(x)<0 for x in strict_gw)
    pairings=[];strict=True
    for index in allowed:
        row=normal['caps'][index];d=parse(row['radicand'])
        assert row['index']==index and sign(d)>=0
        aa,bb=tuple(zip(*(tuple(map(parse,v))for v in row['normal'])))
        ga,gb=mv(gram,aa),mv(gram,bb)
        for actual,saved in zip(zip(ga,gb),row['covector'],strict=True):
            assert qsign(sub(actual[0],parse(saved[0])),sub(actual[1],parse(saved[1])),d)==0
            assert qsign(actual[0],actual[1],d)<=0
        qa=add(dot(aa,ga),mul(d,dot(bb,gb)));qb=scale(dot(aa,gb),2)
        assert qsign(sub(qa,scale(O,2)),qb,d)==0
        a,b=dot(aa,gw),dot(bb,gw);sgn=qsign(a,b,d)
        assert sgn<=0
        strict&=sgn<0
        ca,cb=dot(aa,mv(gram,(O,)*16)),dot(bb,mv(gram,(O,)*16))
        assert qsign(ca,cb,d)<0
        strict_a,strict_b=add(a,scale(ca,epsilon)),add(b,scale(cb,epsilon))
        assert qsign(strict_a,strict_b,d)<0
        pairings.append({'entry':index,'constant':enc(a),'radical_coefficient':enc(b),'sign':sgn,
                         'strict_perturbation_pairing_sign':-1})
    assert qsign(O,neg(O),O)==0 and qsign(neg(O),scale(O,100),Z)<0
    # A strict witness already proves an open ideal-boundary escape patch;
    # if some cap pairing vanished, a small center perturbation would suffice.
    return {'seed_gram_signature':signature,'norm':enc(norm),'center_pairing':enc(center_pair),'center_norm':enc(center_norm),
            'seed_pairings':list(map(enc,gw)),'allowed_cap_pairings':pairings,
            'all_allowed_cap_pairings_strict':strict,'strict_center_perturbation':str(epsilon),
            'strict_coefficients':[str(x[0])for x in strict_w],
            'strict_norm':enc(strict_norm),'strict_center_pairing':enc(strict_center)}


def run(endpoint_source,verify=False):
    graph=json.loads(GRAPH.read_text());normal=json.loads(NORMALS.read_text())
    assert graph['source_sha256']==sha(NORMALS)
    adj=[int(x,16)for x in graph['adjacency_hex']];assert len(adj)==len(normal['caps'])==3040
    assert all(not(m>>i&1)for i,m in enumerate(adj))
    if verify:
        saved=json.loads(OUT.read_text());record=saved['universal_witness']
        assert saved['normal_sha256']==sha(NORMALS)and saved['graph_sha256']==sha(GRAPH)
    else:
        record={'integer_coefficients':[8,8,8,1,-3,-3,-3,-3,-3,-3,-3,5,5,5,5,5],'denominator':1}
    endpoint=json.loads(endpoint_source.read_text())
    assert endpoint['status']=='PASS_EXACT_NECESSARY_IDEAL_EDGE_ENDPOINT_CLAUSES'
    assert endpoint['source_sha256']==sha(NORMALS)
    clauses=[row['nonnegative_caps']for row in endpoint['rows']]
    assert len(clauses)==4 and all(len(c)==14 and len(set(c))==14 for c in clauses)
    compatible,branches,union=enumerate_covers(clauses,adj)
    assert compatible==[(i,j)for i in range(6)for j in range(6,12)]
    assert len(branches)==36 and all(len(row['all_allowed_entries'])==28 for row in branches)
    allowed=list(bits(union));assert len(allowed)==332
    exact=check_witness(normal,record,allowed)
    out={'status':'PASS_EXACT_SOLVER_FREE_FIXED_BANK_FINITE_VOLUME_EXCLUSION',
         'normal_source':str(NORMALS.relative_to(ROOT)),'normal_sha256':sha(NORMALS),
         'graph_source':str(GRAPH.relative_to(ROOT)),'graph_sha256':sha(GRAPH),
         'endpoint_source':str(endpoint_source.relative_to(ROOT)),'endpoint_sha256':sha(endpoint_source),
         'endpoint_clauses':clauses,
         'entry_count':3040,'exact_distinct_mirrors':3034,'endpoint_clause_lengths':list(map(len,clauses)),
         'cartesian_choices':14**4,'compatible_distinct_choices':len(compatible),
         'minimal_cover_count':len(branches),'branches':branches,'all_possible_entries':allowed,
         'discovery':'Simplified by the root agent from saved escape 120; all following checks use the small integer vector directly.',
         'universal_witness':record,'exact_witness_replay':exact,
         'scope':'The separately certified four endpoint clauses and exact compatibility supergraph leave exactly 36 minimal compatible covers. Every cap in any possible selection lies in the archived 332-entry union. One exact spacelike seed-interior witness has nonpositive pairings with every one of those 332 reconstructed normals; its explicit center perturbation is strict. Thus every Coxeter selection from this frozen 3034-distinct-mirror bank has infinite volume. The endpoint-list necessity/completeness and compatibility-supergraph correctness are separate exact premises whose artifacts are hashed here. No LP, numerical cap normals, SAT, or interval sign decision is used in this branch/witness replay. No claim covers new cap mirrors outside the bank or removal of seed facets.'}
    if verify:assert out==saved
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print({'status':out['status'],'minimal_covers':len(branches),'union_entries':len(allowed),
           'strict_cap_pairings':exact['all_allowed_cap_pairings_strict']},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--endpoint-source',type=Path,default=ENDPOINTS);p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.endpoint_source,a.verify)
