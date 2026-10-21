#!/usr/bin/env python3
"""Exact 560-wall family obtained from binary PG(2,4) angle profiles."""
import argparse, json
from collections import Counter
from fractions import Fraction as F
import numpy as np
from fresh_h22_incidence_core_20261013 import core,full_bank,read,sha,ROOT,classify

OUT='results/fresh_h22_incidence_binary_20261013.json'
SYMS='data/enumerations/fresh_h20_global_interior_sections_20261013.json'
def build():
    cert,Q,W,V=core(); Q=np.array(Q,dtype=np.int64); W=np.array(W,dtype=np.int64)
    P,L=cert['incidence_parts']; v,vp=np.array(V,dtype=np.int64)
    A=-(W[P]@Q@W[L].T)
    equations=[sum(int(x)<<i for i,x in enumerate((A[:,j]+A[:,0])%2)) for j in range(1,21)]
    piv={}
    for row in equations:
        while row:
            b=row.bit_length()-1
            if b in piv: row^=piv[b]
            else: piv[b]=row; break
    basis=[]
    for i in range(21):
        if i in piv: continue
        x=1<<i
        for p,row in sorted(piv.items()):
            if (row&x).bit_count()%2: x^=1<<p
        assert all((row&x).bit_count()%2==0 for row in equations)
        basis.append(x)
    assert len(basis)==12
    WQ=W@Q; wn=np.einsum('ij,ij->i',W,WQ); rows=[]; masks=[]; candidates=0
    for b in range(1<<len(basis)):
        mask=0
        for i in range(len(basis)):
            if b>>i&1: mask^=basis[i]
        if mask.bit_count()!=12: continue
        u=np.array([mask>>i&1 for i in range(21)],dtype=np.int64); c=u@A
        if np.any(c%2) or max(c)>4: continue
        candidates+=1; nx=4*vp-u@W[P]
        assert np.all(nx%2==0); x=nx//2
        assert int(x@Q@x)==-2 and np.array_equal(x@Q@np.array([v,vp]).T,[-4,-4])
        pp=x@WQ.T
        try:
            for p,n in zip(pp,wn): classify(int(p),2,int(n))
        except AssertionError: continue
        rows.append(x); masks.append(mask)
    C=np.array(rows,dtype=np.int64); assert C.shape==(280,22)
    assert 22**2*max(int(abs(a).max()) for a in (Q,W,C))**3<2**63
    cc=C@Q@C.T; assert np.all(np.diag(cc)==-2)
    upper=(cc+4)[np.triu_indices(280,1)]
    assert set(upper.tolist())=={-2,-1,0}
    counts=Counter()
    for p in upper: counts[classify(int(p),2,2)]+=2
    for p in cc.ravel(): counts[classify(int(p)-4,2,2)]+=1
    assert sum(counts.values())==560*559//2
    h=np.array(cert['center']); assert np.all(C@Q@h==-8)
    oldcounts=Counter()
    for row in C@WQ.T:
        for p,n in zip(row,wn): oldcounts[classify(int(p),2,int(n))]+=2
    # All 560 walls are right to two same-sign incidence caps and dotted to
    # the two opposite-sign incidence caps.
    assert np.all(C@Q@np.array(V).T==-4)
    meta=read(SYMS)['ambient_symmetry_generators']; lookup={tuple(r):i for i,r in enumerate(C)}
    perms=[]
    for m in meta:
        g=np.array(m['matrix'],dtype=np.int64)
        assert np.array_equal(g.T@Q@g,Q) and np.array_equal(g@h,h)
        assert np.array_equal(W@g.T,W[m['wall_permutation']])
        perm=[lookup[tuple(r)] for r in C@g.T]; assert len(set(perm))==280; perms.append(perm)
    seen={0}; todo=[0]
    while todo:
        i=todo.pop()
        for p in perms:
            j=p[i]
            if j not in seen: seen.add(j); todo.append(j)
    assert len(seen)==280
    # These are exactly the symmetry generators in the inherited complete
    # banks. Their signed versions and lower repair placements are invariant.
    bank=full_bank(Q.tolist(),[C[0].tolist()],new_core_facets=774)
    oldescape=read('results/fresh_h22_new_cusp_completion_20261013_lower_escape.json')
    X=np.array(oldescape['numerator'],dtype=np.int64); values=C@(Q@X[:22])+2*abs(int(X[22]))
    ids=np.flatnonzero(values>0); assert len(ids)==8
    return {'status':'PASS_EXACT_560_INCIDENCE_PROFILE_WALLS',
        'sources':{p:sha(p) for p in (SYMS,'results/fresh_h22_incidence_core_20261013_full.json')},
        'binary_kernel_dimension':12,'profile_candidates':candidates,'accepted_old_projections':C.tolist(),'binary_masks':masks,
        'recipe':'u has twelve 1s; each line sum is even and at most4; x=(4 vprime-sum(u_i P_i))/2. Keep the all210 Coxeter-compatible projections. Add (x,+/-2).',
        'projection_norm':-2,'new_wall_norm':2,'strict_center_pair':-8,
        'new_new_pair_counts':dict(counts),'new_old210_pair_counts':dict(oldcounts),
        'new_core4_pair_counts':{'right':1120,'dotted':1120},
        'symmetry_orbit':{'size':280,'generator_permutations':perms,'representative_index':0},
        'retained_bank_representative_gate':bank,
        'prior_escape_cuts':[{'index':int(i),'pairing':str(F(int(values[i]),oldescape['denominator']))} for i in ids],
        'facets':5398110,
        'scope':'A coherent Coxeter chamber retaining the preceding full bank and four incidence caps, adding560 actual walls. Representative gates extend by the verified transitive orbit and inherited bank invariance. Previous strict escape is cut. Global volume requires a separate certificate.'},Q,W,C
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out,*_=build()
    if a.verify: assert out==read(OUT)
    else: (ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['facets'])
