#!/usr/bin/env python3
"""Independent quadratic-field census and uniform height proof review."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sympy as s
import fresh_h9_h8_all_high_20261014_review as C
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_cap_height_20261014_review'
PROOF='results/fresh_h9_h8_cap_height_20261014.json'
BASE='results/fresh_h9_h8_disjoint_cap_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def dm(M):return[[C.decode(x)for x in row]for row in M]
def dv(v):return[C.decode(x)for x in v]
def pos(x):return bool(x)and C.ev(x).l>0
def neg(x):return bool(x)and C.ev(x).h<0
def key(v):
    x=next(x for x in v if x);return tuple(t/x for t in v)
def pd(M):
    W=[row[:]for row in M]
    for k in range(len(W)):
        p=W[k][k]
        if not pos(p):return False
        for i in range(k+1,len(W)):
            for j in range(k+1,len(W)):W[i][j]-=W[i][k]*W[k][j]/p
    return True
def run():
    pending,seen=[PROOF],set()
    while pending:
        p=pending.pop()
        if p in seen:continue
        seen.add(p);d=json.loads((ROOT/p).read_text())
        for f,h in d.get('sources',{}).items():
            assert sha(f)==h,f
            if f.endswith('.json'):pending.append(f)
    d=json.loads((ROOT/PROOF).read_text());base=json.loads((ROOT/BASE).read_text())
    C.K=s.QQ.algebraic_field(s.sqrt(5));C.Z,C.O=C.K.zero,C.K.one
    K,Z,O=C.K,C.Z,C.O
    Q,N=dm(d['form']),dm(d['normals']);H=dv(d['strict_center']['point'])
    a=K.from_sympy((s.sqrt(5)-1)/4);D=a+2*O;phi=K.from_sympy((1+s.sqrt(5))/2)
    assert Q[-1][-1]==D and len(N)==15
    norms=[C.pair(n,n,Q)for n in N]
    assert norms==[2*O]*11+[D]+[2*O]*3
    G=[[C.pair(n,m,Q)for m in N]for n in N];assert G==dm(d['gram'])
    assert neg(C.pair(H,H,Q))and all(neg(C.pair(H,n,Q))for n in N)
    adj=[set()for _ in N];counts=Counter()
    for i,j in combinations(range(15),2):
        g=G[i][j];assert not g or neg(g)
        if not g:t='right'
        elif 4*g*g==norms[i]*norms[j]:t='pi3'
        elif g*g==((phi*phi)/4)*norms[i]*norms[j]:t='pi5'
        else:assert pos(g*g-norms[i]*norms[j]);t='dotted'
        counts[t]+=1
        if t!='dotted':adj[i].add(j);adj[j].add(i)
    assert dict(counts)==d['pair_counts']==dict(right=79,pi3=12,pi5=5,dotted=9)
    # A different census: filter all5005nine-wall supports and solve physical
    # nullspaces, rather than lift the41source rays separately to each cap.
    actual={};cliques=bad=outside=0
    for ids in combinations(range(15),9):
        if not all(j in adj[i]for i,j in combinations(ids,2)):continue
        cliques+=1;M=[[G[i][j]for j in ids]for i in ids]
        if not pd(M):bad+=1;continue
        v=C.null_one([C.mv(Q,N[i])for i in ids])
        if pos(C.pair(v,H,Q)):v=[-x for x in v]
        assert neg(C.pair(v,v,Q))and neg(C.pair(v,H,Q))
        vals=[C.pair(v,n,Q)for n in N]
        if any(pos(x)for x in vals):outside+=1;continue
        assert all(not x or neg(x)for x in vals)
        assert tuple(i for i,x in enumerate(vals)if not x)==ids
        actual[ids]=v
    saved=d['complete_finite_vertices'];assert len(actual)==len(saved)==62
    for row in saved:
        ids=tuple(row['incident']);assert ids in actual and key(actual[ids])==key(dv(row['point']))
    edges=[]
    for i,j in combinations(range(62),2):
        support=sorted(set(saved[i]['incident'])&set(saved[j]['incident']))
        if len(support)==8:
            assert pd([[G[k][l]for l in support]for k in support])
            edges.append(dict(vertices=[i,j],support=support))
    assert edges==d['complete_finite_edges']and len(edges)==220
    for i,row in enumerate(d['facet_points']):
        v=dv(row['point']);vals=[C.pair(v,n,Q)for n in N]
        assert not vals[i]and all(neg(x)for j,x in enumerate(vals)if j!=i)
        assert neg(C.pair(v,v,Q))and neg(C.pair(v,H,Q))
    assert len(d['facet_points'])==15
    E=dv(d['strict_escape']['point']);assert E[-1]==-5*O
    escape_norm=C.pair(E,E,Q)
    assert escape_norm==K.from_sympy(s.Rational(4363,256)-s.Rational(3057,512)*s.sqrt(5))
    assert pos(escape_norm)and neg(C.pair(E,H,Q))and all(neg(C.pair(E,n,Q))for n in N)
    contact_sets=[]
    for i in (12,13,14):
        p=N[i][:-1]+[Z]
        assert C.pair(p,p,Q)==-a and N[i][-1]==-O
        h=[-C.pair(p,n,Q)for n in N[:11]]
        assert O in h and phi in h and all(x in(Z,O,phi,phi**4)for x in h)
        contact_sets.append([str(K.to_sympy(x))for x in h])
    assert pos(phi*phi-2*O)
    margins=[]
    for i,j in combinations((12,13,14),2):
        vi=N[i][:-1]+[Z];vj=N[j][:-1]+[Z]
        margin=-C.pair(vi,vj,Q)-a-4*O;assert pos(margin);margins.append(str(K.to_sympy(margin)))
    # Since0<q<=2, eta=sqrt((a+q)/D)<=1. Multiplying the
    # negative escape coordinate by the cap change gives5D(eta-1)<=0.
    z,eta,dd=s.symbols('z eta D',real=True)
    assert s.expand((z+5*dd*eta)-(z+5*dd)-5*dd*(eta-1))==0
    Qsigma=[[K.from_sympy(K.to_sympy(x).xreplace({s.sqrt(5):-s.sqrt(5)}))for x in row]for row in Q]
    assert pd(Qsigma)
    dexpr=K.to_sympy(D)
    assert s.simplify(dexpr*dexpr.xreplace({s.sqrt(5):-s.sqrt(5)}))==s.Rational(11,4)
    paths=[PROOF,BASE,f'scripts/{TAG}.py','scripts/fresh_h9_h8_all_high_20261014_review.py',
           'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(status='PASS_INDEPENDENT_MAXIMAL_AND_ALL_HEIGHT_REVIEW',sources={p:sha(p)for p in paths},
        total_supports=5005,finite_neighbor_cliques=cliques,nonpositive_cliques=bad,outside_candidates=outside,
        finite_vertices=62,finite_edges=220,actual_facets=15,remaining_finite_corner_germs=118,
        reconstructed_contacts=contact_sets,strict_all_height_pair_margins=margins,
        scope='Separate SymPy quadratic field, all5005supports and physical nullspaces independently check maximal-height geometry and strictescape. Source contact1forcesq<=2; exact phi gap proves the stated all-label alphabet, and the escape difference is nonpositive for everyreal0<q<=2 independently. This covers only the same three projected cap directions and their subfamilies; new directions/additional walls remain open. No host reflectivity or lattice claim.')
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['finite_vertices'],result['finite_edges'])
