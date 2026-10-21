#!/usr/bin/env python3
"""Exact field obstruction to finite completion of a retained deficient cusp.

Uses only Python integers/Fractions. Coordinates in the geometric statement
use diag(Q,1); the inherited cusp file instead records the coefficient of
sqrt(11) in its last coordinate.
"""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
CUSP='results/fresh_h22_lower_new_cusp_20261013.json'
PROPOSAL='data/enumerations/fresh_h22_incidence_field_obstruction_20261013_proposal.json'
CORE='results/fresh_h22_incidence_core_20261013.json'
BINARY='results/fresh_h22_incidence_binary_20261013.json'
OUT='results/fresh_h22_incidence_field_obstruction_20261013.json'
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def mv(M,v):return [dot(r,v) for r in M]
def rank(A):
    A=[list(map(F,r)) for r in A];n=0
    for j in range(len(A[0])):
        i=next((i for i in range(n,len(A)) if A[i][j]),None)
        if i is None:continue
        A[i],A[n]=A[n],A[i];p=A[n][j];A[n]=[x/p for x in A[n]]
        for i in range(n+1,len(A)):
            if A[i][j]:c=A[i][j];A[i]=[x-c*y for x,y in zip(A[i],A[n])]
        n+=1
        if n==len(A):break
    return n
def run():
    raw=read(RAW)['coxeter_input'];Q=[[-int(x) for x in r] for r in raw['basis_gram']]
    W=[[int(x) for x in r] for r in raw['wall_coordinates_as_columns_transposed']]
    data=read(CUSP);assert all(sha(p)==v for p,v in data['sources'].items())
    r=data['representatives'][0];assert r['representative_type']==4
    delta=list(map(F,r['null_ray']));a=delta[:22];assert delta[-1]==2
    ac=mv(Q,a);assert dot(a,ac)==-44
    prop=read(PROPOSAL);assert list(map(F,prop['delta_old']))==a
    eps=F(prop['epsilon']);assert eps>0
    ids=prop['support'];weights=list(map(F,prop['excess_weights']))
    assert len(ids)==len(set(ids))==len(weights)==22 and min(weights)>0
    assert rank([W[i] for i in ids])==22
    rhs=[eps*sum(w[j] for w in W)+sum(t*W[i][j] for i,t in zip(ids,weights)) for j in range(22)]
    assert rhs==a
    # Restore the actual retained cap c4 and its parallel old wall.
    lower=read('results/fresh_h22_new_cusp_completion_20261013_lower.json')
    c=[F(x,lower['common_denominator']) for x in lower['new_cap_old_numerators'][4]]+[2*F(lower['shifts'][4])]
    w=W[26]+[0];assert [x+y for x,y in zip(c,w)]==delta
    Q11=[row+[0] for row in Q]+[[0]*22+[11]]
    assert dot(c,mv(Q11,c))==dot(w,mv(Q11,w))==2
    assert dot(c,mv(Q11,w))==-2
    incident=[W[i]+[0] for i in r['incident_old_wall_indices']]+[c]
    assert len(incident)==18 and rank(incident)==18
    gram=[[dot(u,mv(Q11,v)) for v in incident] for u in incident]
    assert rank(gram)==17
    tangent=r['positive_lineality_basis'][0]
    assert dot(tangent,mv(Q11,tangent))>0
    assert dot(tangent,mv(Q11,delta))==0
    assert all(dot(tangent,mv(Q11,v))==0 for v in incident)
    new=read(CORE)['old_projections']+read(BINARY)['accepted_old_projections']
    vals=[dot(v,ac) for v in new];assert max(vals)==-25
    # Every new physical normal has last coordinate +/-2, hence its cusp
    # pairing is at most -25+4sqrt(11)<0.
    assert 25**2>16*11
    return {'status':'PASS_EXACT_RETAINED_CUSP_FIELD_OBSTRUCTION',
        'sources':{p:sha(p) for p in (RAW,CUSP,PROPOSAL,CORE,BINARY)},
        'physical_null_ray':{'old_part':list(map(str,a)),'last_coordinate':'2 sqrt(11)'},
        'strict_positive_root_identity':{'uniform_weight':str(eps),'excess_indices':ids,'excess_weights':list(map(str,weights)),'spanning_rank':22},
        'retained_parallel_pair':['old26','lower_highest_root_type4'],
        'unchanged_cusp':{'incident_facets':18,'normal_rank':18,'quotient_rank':17,'horosphere_dimension':21,'positive_recession_vector':tangent,'positive_recession_square':str(dot(tangent,mv(Q11,tangent))),'new564_pair_upper':'-25+4 sqrt(11)<0'},
        'field_theorem':'Let F be any real subfield with sqrt(11) not in F. No nonzero normal with a projective representative in F^23, acute to all210 old normals, can be orthogonal to delta=(a,2sqrt(11)). Orthogonality separates over F into r_last=0 and (a,r_old)=0. The strictly positive old-root identity and rank22 force r_old=0.',
        'finite_extension_consequence':'Retaining old26 and cap c4 forces every compatible new normal to pair nonpositively with delta. For each allowed projective field the pairing is strictly negative. Any finite collection of such walls leaves the deficient cusp link unchanged, so it cannot turn the retained chamber into a finite-volume polytope.',
        'scope':'The allowed field may vary between added walls, provided none contains sqrt(11), in the fixed physical coordinates diag(Q,1). This rules out rational-direction completion and many other restricted-field completions while retaining c4 and old210. It does not exclude normals over fields containing sqrt(11), replacements of retained walls, other cusps or other H22 constructions.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();out=run()
    if args.verify:assert out==read(OUT)
    else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
