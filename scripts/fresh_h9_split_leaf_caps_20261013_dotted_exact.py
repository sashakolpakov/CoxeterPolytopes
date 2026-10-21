#!/usr/bin/env python3
"""Exact 36 dotted caps and constructive completion of the combined228 bank.

Discovery uses linear programming only to propose rational escapes. Verification
uses exact arithmetic; a final failed-bank claim requires coverage of every
maximal compatible clique, not the absence of a numerical survivor.
"""
import argparse,hashlib,json
from fractions import Fraction as Q
from itertools import combinations,product
from pathlib import Path
from scripts.fresh_h9_split_leaf_caps_20261013_exact import maximal_cliques
from scripts.fresh_h9_quartic_completion_20261013_arithmetic import K,Quadratic,dot,matvec
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data/enumerations/fresh_h9_split_leaf_caps_20261013_exact.json'
OUT=ROOT/'data/enumerations/fresh_h9_split_leaf_caps_20261013_dotted_exact.json'
P=(Q(0),Q(1));L=Quadratic(K,P);Z=L.zero;O=L.one

def gparse(a):return tuple(Q(t) for t in a)
def total(xs):
    r=Z
    for x in xs:r=L.add(r,x)
    return r
def times(x,c):return L.mul(x,L.rational(c))
def numerical(x):
    phi=(1+5**.5)/2
    return float(x[0][0])+float(x[0][1])*phi+(float(x[1][0])+float(x[1][1])*phi)*phi**.5

def setup():
    source=json.loads(BASE.read_text());G=tuple(tuple(L.lift(gparse(x)) for x in r) for r in source['gram'])
    a=L.lift(gparse(source['a']));b=L.lift(gparse(source['b']));S=(O,)*16
    assert a==L.lift((Q(1,2),Q(1))) and b==L.lift((Q(-5,22),Q(1,11)))
    eps=(1,)*4+(-1,)*4
    for i in range(16):
        for j in range(16):
            expected=L.rational(2) if i==j else (Z if i//8==j//8 else
                L.neg(L.add(L.add(a,times(b,eps[i%8]*eps[j%8])),L.rational(2 if i%8==j%8 else 0))))
            assert G[i][j]==expected
    rows=[tuple(L.lift(gparse(x)) for x in row['normal']) for row in source['caps']]
    d=(P,P) # phi + phi*sqrt(phi)
    assert L.sign(L.sub(d,L.rational(2)))>0
    patterns=[]
    for ii,jj in product(combinations(range(4),2),combinations(range(4,8),2)):
        pos=set(ii+jj);eta=tuple(1 if i in pos else -1 for i in range(8))
        c=L.div(L.add(d,O),times(a,16));w=times(L.sub(d,O),Q(1,8))
        # With this sign, x_i=1 when eta_i=+1 and x_i=d otherwise.
        r=tuple(L.add(c,times(w,e)) for e in eta)+tuple(L.sub(c,times(w,e)) for e in eta)
        gr=matvec(L,G,r)
        expected=tuple(L.neg(O if e==1 else d) for e in eta)+tuple(L.neg(d if e==1 else O) for e in eta)
        assert gr==expected and dot(L,r,gr)==L.rational(2)
        rows.append(r);patterns.append(sorted(pos))
    cov=[matvec(L,G,r) for r in rows]
    assert all(dot(L,r,gr)==L.rational(2) for r,gr in zip(rows,cov))
    for row,gr in zip(source['caps'],cov):
        assert gr==tuple(L.neg(L.lift(gparse(v))) for v in row['x']+row['y'])
    assert all(L.sign(dot(L,S,gr))<0 for gr in cov)
    adj=[0]*len(rows);edges=[]
    for i,j in combinations(range(len(rows)),2):
        val=dot(L,rows[i],cov[j]);sg=L.sign(val)
        # Any element u+v sqrt(phi) with v!=0 has a nonreal conjugate,
        # so cannot be a finite Coxeter value, all of whose conjugates are real.
        finite=val[1]==K.zero and val[0] in (K.zero,K.rational(-1),K.neg(P))
        okay=finite or L.sign(L.add(val,L.rational(2)))<=0
        if okay:
            assert sg<=0;adj[i]|=1<<j;adj[j]|=1<<i;edges.append([i,j,L.encode(val)])
    cliques=maximal_cliques(adj)
    masks=[sum(1<<i for i in c) for c in cliques]
    metadata={'field':'Q(phi,sqrt(phi)), phi^2=phi+1, both chosen positive',
              'source':str(BASE.relative_to(ROOT)),'dotted_magnitude':L.encode(d),
              'new_cap_patterns':patterns,'new_cap_normals':[[L.encode(x) for x in r] for r in rows[192:]],
              'cap_count':len(rows),'compatibility_edges':len(edges),
              'edge_digest':hashlib.sha256(json.dumps(edges,separators=(',',':')).encode()).hexdigest(),
              'maximal_clique_count':len(cliques),'maximum_caps':max(map(len,cliques)),
              'selected_fourteen_caps':list(next(c for c in cliques if len(c)==14)),
              'clique_digest':hashlib.sha256(json.dumps(cliques,separators=(',',':')).encode()).hexdigest()}
    assert len(rows)==228 and len(edges)==8460 and len(cliques)==156990 and max(map(len,cliques))==14
    return G,rows,cov,S,cliques,masks,metadata

def exact_witness(G,rows,S,raw):
    w=tuple(L.rational(Q(x)) for x in raw);gw=matvec(L,G,w)
    assert L.sign(dot(L,w,gw))>0 and L.sign(dot(L,S,gw))<0
    assert all(L.sign(x)<0 for x in gw)
    signs=[L.sign(dot(L,r,gw)) for r in rows]
    assert 0 not in signs
    cuts=sum(1<<i for i,x in enumerate(signs) if x>0)
    return cuts,{'coefficients':[str(Q(x)) for x in raw],
                 'square':L.encode(dot(L,w,gw)),
                 'center_pairing':L.encode(dot(L,S,gw)),
                 'cut_caps':[i for i,x in enumerate(signs) if x>0]}

def propose(G,rows,S,clique,index):
    import numpy as np
    from scipy.optimize import linprog
    gn=np.array([[numerical(x) for x in r] for r in G]);rs=np.array([[numerical(x) for x in r] for r in rows])
    eig,Qm=np.linalg.eigh(gn);inds=list(np.where(eig< -1e-7)[0])+list(np.where(eig>1e-7)[0]);N=Qm[:,inds]*np.sqrt(abs(eig[inds]))
    if sum(N[:,0])<0:N[:,0]*=-1
    wall=np.vstack([N,rs[list(clique)]@N]);rng=np.random.default_rng(20261013+index)
    candidates=list(np.eye(9))+list(-np.eye(9))+list(rng.normal(size=(30,9)))
    for direction in candidates:
        ans=linprog(-direction,A_ub=wall[:,1:],b_ub=wall[:,0],bounds=[(None,None)]*9,method='highs')
        assert ans.success
        if ans.x@ans.x<=1.0001:continue
        c=np.linalg.pinv(N.T)@np.r_[1,ans.x]
        for exp in [8,10,12]:
            raw=[Q(str(float(x))).limit_denominator(10**exp)+Q(1,10**(exp//2)) for x in c]
            try:cuts,w=exact_witness(G,rows,S,raw)
            except AssertionError:continue
            if any((cuts>>i)&1 for i in clique):continue
            return cuts,w
    raise RuntimeError('No exact escape found; candidate requires full geometry investigation')

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    G,rows,cov,S,cliques,masks,meta=setup();print({k:meta[k] for k in ['cap_count','compatibility_edges','maximal_clique_count','maximum_caps']},flush=True)
    if args.verify:
        saved=json.loads(OUT.read_text());assert saved['metadata']==meta;remaining=masks[:]
        assert saved['status']=='PASS_EXACT_COMBINED228_BANK' and saved['remaining_maximal_cliques']==0
        for old in saved['witnesses']:
            cut,checked=exact_witness(G,rows,S,old['coefficients']);assert checked==old
            remaining=[m for m in remaining if m&cut]
        assert not remaining
        bad=json.loads(json.dumps(saved['witnesses'][0]));bad['coefficients']=['0']*16
        try:exact_witness(G,rows,S,bad['coefficients'])
        except AssertionError:pass
        else:raise AssertionError('Forged escape accepted')
        print({'status':'PASS_EXACT_COMBINED228_BANK','witnesses':len(saved['witnesses']),'forged_escape_rejected':True});return
    remaining=list(zip(cliques,masks));witnesses=[]
    for iteration in range(1000):
        if not remaining:break
        clique,_=max(remaining,key=lambda x:len(x[0]));cut,w=propose(G,rows,S,clique,iteration);witnesses.append(w)
        remaining=[(c,m) for c,m in remaining if m&cut]
        out={'metadata':meta,'witnesses':witnesses,'remaining_maximal_cliques':len(remaining),
             'status':'PASS_EXACT_COMBINED228_BANK' if not remaining else 'PARTIAL',
             'scope':'Every compatible subset of this fixed228-cap bank has infinite volume only when remaining_maximal_cliques is zero. Other caps are not excluded.'}
        OUT.write_text(json.dumps(out,indent=2)+'\n')
        print({'witnesses':len(witnesses),'remaining':len(remaining)},flush=True)
    assert not remaining,'Search budget reached without a full-bank conclusion'
if __name__=='__main__':main()
