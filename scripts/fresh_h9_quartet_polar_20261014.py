#!/usr/bin/env python3
"""Exact exterior-ray census and polar-cap test for the nonarithmetic quartet."""
import argparse,json,sys
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import fresh_h9_h8_disjoint_walk_20261014_exact as A
import fresh_h9_h8_disjoint_walk_20261014_next_frontier as D
from coxeter_polytopes.verify import chamber as C
TAG='fresh_h9_quartet_polar_20261014'
BASE='results/fresh_h9_h8_disjoint_walk_20261014_exact.json'

# Adapt the existing field-generic elimination engine locally, without editing
# its source or the frozen quartic implementation. Only embedding0 is used.
old_sign=A.L.sign
def sign(self,embedding=0):
    assert embedding==0
    return old_sign(self)
A.L.sign=sign
C.Golden=A.L;C.ZERO=A.Z

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def relation(g,q,r):
    assert g.sign()<=0
    if not g:return 'right',2,None
    c=g*g/(q*r)
    if(c-1).sign()>=0:return ('parallel'if c==1 else'dotted'),None,c
    allowed={3:A.L(F(1,4)),4:A.L(F(1,2)),5:A.L(A.A.Q(F(3,8),F(1,8))),
             6:A.L(F(3,4)),10:A.L(A.A.Q(F(5,8),F(1,8)))}
    for m,v in allowed.items():
        if c==v:return 'finite',m,c
    return 'forbidden',None,c

def run():
    d=json.loads((ROOT/BASE).read_text())
    for p,h in d['sources'].items():assert sha(p)==h,p
    A.RAD=A.A.Q(*d['field']['new_positive_radical_square'])
    assert A.RAD==A.A.Q(6,F(16,5))
    radical_norm=A.RAD.a*A.RAD.a-5*A.RAD.b*A.RAD.b
    assert radical_norm==F(-76,5) and A.RAD.sign()>0
    assert A.A.Q(A.RAD.a,-A.RAD.b).sign()<0
    Q,N,H,E=D.dm(d['form']),D.dm(d['normals']),D.dv(d['center']),D.dv(d['strict_escape']['point'])
    norms=[A.pair(n,n,Q)for n in N];functionals=[A.mv(Q,n)for n in N]
    census=C.diagnose_chamber_cone(N,Q,H,9)
    names=('finite_vertices','cusps','hyperideal_rays','opposite_sheet_rays')
    counts={k:len(getattr(census,k))for k in names};print('complete cone',counts,flush=True)
    assert set(census.finite_vertices)=={tuple(v['incident'])for v in d['complete_finite_vertex_census']['vertices']}
    records=[];stats=Counter()
    for index,active in enumerate(census.hyperideal_rays):
        support=next(s for s in combinations(active,9)if C._full_rank_kernel([functionals[i]for i in s])is not None)
        r=C._full_rank_kernel([functionals[i]for i in support])
        vals=[A.pair(r,n,Q)for n in N]
        if any(z.sign()>0 for z in vals):r=[-z for z in r];vals=[-z for z in vals]
        assert all(z.sign()<=0 for z in vals)
        assert tuple(i for i,z in enumerate(vals)if not z)==active
        q=A.pair(r,r,Q);assert q.sign()>0
        pairs=[relation(v,q,norms[i])for i,v in enumerate(vals)]
        bad=[i for i,x in enumerate(pairs)if x[0]=='forbidden']
        center=A.pair(r,H,Q);escape=A.pair(r,E,Q)
        compatible=not bad
        stats['polar_candidates']+=1
        stats['compatible'if compatible else'forbidden_pair']+=1
        if compatible:stats['compatible_cuts_escape'if escape.sign()>0 else'compatible_misses_escape']+=1
        records.append(dict(index=index,active=list(active),independent_support=list(support),normal=r,
            norm=q,center_pairing=center,escape_pairing=escape,pairings=vals,
            relations=[dict(kind=t,label=m,invariant=c)for t,m,c in pairs],
            forbidden_indices=bad,coxeter_compatible=compatible,retains_center=center.sign()<0))
    paths=[BASE,f'scripts/{TAG}.py','scripts/fresh_h9_h8_disjoint_walk_20261014_exact.py',
           'scripts/fresh_h9_h8_disjoint_walk_20261014_next_frontier.py',
           'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',
           'src/coxeter_polytopes/verify/chamber.py']
    return A.enc(dict(status='EXACT_QUARTET_EXTERIOR_POLAR_CENSUS',sources={p:sha(p)for p in paths},
        radical=A.RAD,radical_norm=str(radical_norm),field_signature=[2,1],counts=counts,polar_counts=dict(stats),
        complete_census={k:[list(x)for x in getattr(census,k)]for k in names},polar_candidates=records,
        finite_angle_alphabet=[2,3,4,5,6,10],
        completeness='Every exterior extreme ray of the literal19-wall cone is tested as its polar normal. Exact incremental double description is reused through an isolated quartic scalar adapter; finite vertices match the independently verified53-vertex census.',
        angle_scope='All finite Coxeter labels: this quartic field has two real and two complex places, with unique quadratic subfieldQ(sqrt5). A real cyclotomic cosine-square in it is rational or inQ(sqrt5), giving labels2,3,4,5,6,10.',
        scope='Polar hyperplanes of existing exterior vertices only. A compatible polar candidate still needs a common interior, pair compatibility with other new walls and global geometry. Nonpolar caps, changed retained walls and arbitrary completions remain open.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],d['polar_counts'])
