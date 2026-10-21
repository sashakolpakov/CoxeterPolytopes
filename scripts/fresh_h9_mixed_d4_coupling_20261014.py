#!/usr/bin/env python3
"""Coupled D4 lift of a mixed compact H5 prism, with exact polar closure.

All predicates are exact in L=Q(sqrt5,sqrt((16+4sqrt5)/11)).  Cone clipping
enumerates every projective extreme ray, rather than only finite vertices.
No finite-volume or new-dimension claim follows from a local face lift.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from fresh_h9_compactification_20261013 import (
    L,K,R,PHI,Z,O,sm,dot,mv,pair,eye,inv,null1,inertia,geometry,sqrat)

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_mixed_d4_coupling_20261014'

def enc(x):
    if isinstance(x,(L,K)):return x.enc()
    if isinstance(x,(list,tuple)):return [enc(y) for y in x]
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    return x

def rank(A):
    if not A:return 0
    M=[r[:] for r in A];k=0
    for j in range(len(M[0])):
        p=next((i for i in range(k,len(M)) if M[i][j]),None)
        if p is None:continue
        M[k],M[p]=M[p],M[k];a=M[k][j];M[k]=[x/a for x in M[k]]
        for i in range(k+1,len(M)):
            if M[i][j]:
                a=M[i][j];M[i]=[x-a*y for x,y in zip(M[i],M[k])]
        k+=1
        if k==len(M):break
    return k

def normalize(v):
    a=next(x for x in v if x)
    if a.sign()<0:a=-a
    return tuple(x/a for x in v)

def cone(Q,N):
    n=len(Q);rows=[mv(Q,w) for w in N];base=[]
    for i,row in enumerate(rows):
        if rank([rows[j] for j in base]+[row])>len(base):base.append(i)
        if len(base)==n:break
    assert len(base)==n
    Ai=inv([rows[i] for i in base])
    rays={normalize([-Ai[j][i] for j in range(n)]) for i in range(n)}
    used=base[:];steps=[];cached={}
    for k in [i for i in range(len(N)) if i not in base]:
        active={v:frozenset(i for i in used if not dot(rows[i],v)) for v in rays}
        values={v:dot(rows[k],v) for v in rays}
        pos=[v for v in rays if values[v].sign()>0]
        neg=[v for v in rays if values[v].sign()<0]
        keep={v for v in rays if values[v].sign()<=0};made=0
        for a in pos:
            for b in neg:
                common=tuple(sorted(active[a]&active[b]))
                if common not in cached:cached[common]=rank([rows[i] for i in common])
                if cached[common]!=n-2:continue
                v=normalize([values[a]*y-values[b]*x for x,y in zip(a,b)])
                assert not dot(rows[k],v)
                assert all(dot(rows[i],v).sign()<=0 for i in used)
                keep.add(v);made+=1
        steps.append(dict(wall=k,before=len(rays),positive=len(pos),negative=len(neg),
                          adjacent_crossings=made,after=len(keep)))
        rays=keep;used.append(k)
    result=[]
    for v in sorted(rays,key=lambda x:str(enc(x))):
        q=pair(v,v,Q);active=[i for i,row in enumerate(rows) if not dot(row,v)]
        assert rank([rows[i] for i in active])==n-1
        assert all(dot(row,v).sign()<=0 for row in rows)
        result.append(dict(vector=v,norm=q,kind={-1:'finite',0:'ideal',1:'spacelike'}[q.sign()],incident=active))
    dual=[sm(v[j] for v in N) for j in range(n)]
    assert pair(dual,dual,Q).sign()<0
    assert all(pair(row['vector'],dual,Q).sign()<0 for row in result)
    return dict(counts=dict(Counter(row['kind'] for row in result)),rays=result,
                initial_independent_walls=base,clipping_steps=steps,
                positive_dual_sum=dual,positive_dual_sum_norm=pair(dual,dual,Q))

COS={2:Z,3:L(F(1,4)),4:L(F(1,2)),5:(L(PHI)+1)/4,6:L(F(3,4)),10:(L(PHI)+2)/4}
def relation(g,q,r):
    if not g:return dict(type='right',label=2)
    if g.sign()>0:return dict(type='obtuse',product=g)
    c=g*g/(q*r)
    if (c-1).sign()>=0:return dict(type='parallel' if c==O else 'dotted',squared_weight=c)
    for m,w in COS.items():
        if c==w:return dict(type='finite',label=m,squared_weight=c)
    return dict(type='forbidden_finite',squared_weight=c)

def polar_bank(Q,N,h,census):
    result=[];norms=[pair(n,n,Q) for n in N]
    for ri,row in enumerate(census['rays']):
        if row['kind']!='spacelike':continue
        v=row['vector'];q=row['norm']
        rels=[relation(pair(v,n,Q),q,s) for n,s in zip(N,norms)]
        bad=[i for i,r in enumerate(rels) if r['type'] in ('obtuse','forbidden_finite')]
        result.append(dict(ray=ri,normal=v,norm=q,relations=rels,bad_old_walls=bad,
                           center_pairing=pair(v,h,Q),passes=not bad))
    return result

def escape(Q,N,h,census):
    for row in census['rays']:
        if row['kind']!='spacelike':continue
        for k in range(1,10):
            v=[x+y/F(10**k) for x,y in zip(row['vector'],h)]
            if pair(v,v,Q).sign()<=0 or pair(v,h,Q).sign()>=0:continue
            v=[L(F(round(x.approx()*10**9),10**9)) for x in v]
            vals=[pair(v,n,Q) for n in N]
            if pair(v,v,Q).sign()>0 and pair(v,h,Q).sign()<0 and all(x.sign()<0 for x in vals):
                return dict(vector=v,norm=pair(v,v,Q),center_pairing=pair(v,h,Q),all_pairings=vals)
    raise AssertionError('No strict rational escape certified')

def build():
    G,Q5,_,W,V,inds,h5,cert,_,_,_,_=geometry()
    D=eye(4)
    for i in range(4):D[i][i]=L(2)
    for j in (1,2,3):D[0][j]=D[j][0]=-O
    Di=inv(D);omega=[Di[i][1] for i in range(4)]
    assert pair(omega,omega,D)==O
    Q=[row+[Z]*4 for row in Q5]+[[Z]*6+row for row in D]
    N=[w+[Z]*4 for w in W]
    N[6]=[x/2 for x in W[6]]+[-x for x in omega]
    N += [[Z]*6+v for v in eye(4)]
    assert inertia(Q)==[9,1] and all(pair(n,n,Q)==L(2) for n in N)
    h=h5+[-sm(row)/100 for row in Di]
    assert pair(h,h,Q).sign()<0 and all(pair(h,n,Q).sign()<0 for n in N)
    source=cone(Q5,W);assert source['counts']=={'finite':10}
    seed=cone(Q,N);print('seed',seed['counts'],flush=True)
    # Independent support enumeration for the small seed, using the same field.
    enumerated=set();dual=seed['positive_dual_sum']
    for ids in combinations(range(11),9):
        v=null1([mv(Q,N[i]) for i in ids])
        if v is None:continue
        if pair(v,dual,Q).sign()>0:v=[-x for x in v]
        if all(pair(v,n,Q).sign()<=0 for n in N):enumerated.add(normalize(v))
    assert enumerated=={tuple(row['vector']) for row in seed['rays']}
    initial=N[:];census=seed;rounds=[];complete_polar_closure=False
    for step in range(4):
        bank=polar_bank(Q,N,h,census);passed=[i for i,r in enumerate(bank) if r['passes']]
        mutual=[]
        for i,j in combinations(passed,2):
            a,b=bank[i],bank[j]
            mutual.append(dict(candidates=[i,j],**relation(pair(a['normal'],b['normal'],Q),a['norm'],b['norm'])))
        assert all(r['type'] in ('right','finite','parallel','dotted') for r in mutual)
        rounds.append(dict(wall_count=len(N),complete_ray_census=census,polar_candidates=bank,
                           passing_polar_indices=passed,mutual_relations=mutual))
        print('round',step,'walls',len(N),'rays',census['counts'],'passing',passed,flush=True)
        if not passed:complete_polar_closure=True;break
        assert all(bank[i]['center_pairing'].sign()<0 for i in passed)
        N += [list(bank[i]['normal']) for i in passed]
        census=cone(Q,N)
    facets=[];pairs=[];kinds=Counter()
    for i,n in enumerate(N):
        q=pair(n,n,Q);p=[a-pair(h,n,Q)*b/q for a,b in zip(h,n)]
        assert pair(p,p,Q).sign()<0
        assert all(pair(p,m,Q).sign()==(0 if i==j else -1) for j,m in enumerate(N))
        facets.append(p)
        for j in range(i):
            rel=relation(pair(n,N[j],Q),q,pair(N[j],N[j],Q))
            assert rel['type'] in ('right','finite','parallel','dotted')
            pairs.append(dict(pair=[j,i],**rel));kinds[rel['type']+str(rel.get('label',''))]+=1
    normR=R.a*R.a-5*R.b*R.b
    assert normR==F(16,11) and not sqrat(normR) and not sqrat(normR/5)
    end_product=pair(initial[5],initial[6],Q);cycle=end_product*end_product
    assert cycle.b and inertia(G)==[4,1]
    files=[f'scripts/{TAG}.py','scripts/fresh_h9_compactification_20261013.py',
           'scripts/fresh_h9_makarov_face_extension_20261013.py']
    return enc(dict(status='PASS_EXACT_MIXED_COMPACT_H5_COUPLED_D4_AND_POLAR_DIAGNOSTIC',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
        field=dict(base='Q(sqrt5)',extension='positive sqrt(R)',R=R,norm_R=str(normR),nonnormal_quartic=True),
        source_form=Q5,source_normals=W,source_complete_ray_census=source,
        support_gram=D,support_weight=omega,form=Q,signature=[9,1],
        seed_normals=initial,center=h,center_norm=pair(h,h,Q),seed_complete_ray_census=seed,
        seed_support_enumeration_count=len(enumerated),polar_rounds=rounds,
        exhausted_iterated_polar_closure=complete_polar_closure,polar_round_limit=4,
        final_normals=N,final_complete_ray_census=census,final_actual_facet_points=facets,
        final_pairs=pairs,final_pair_counts=dict(kinds),strict_escape=escape(Q,N,h,census),
        nonquasiarithmetic_obstruction=dict(end_cartan_two_cycle=cycle,
            fixed_principal_gram=G,fixed_principal_signature=[4,1],
            embedding='fix sqrt5 and negate sqrt(R)',
            scope='Any finite-volume completion retaining these seed facets fails the conjugate semidefiniteness condition. This does not establish finite volume.'),
        scope='Exact coupled D4 lift of a known mixed compact H5 source and the recorded complete polar rounds. Complete cone rays and a strict rational escape prove the resulting H9 chamber has infinite volume. Exhaustion of extreme-ray polar candidates does not exclude arbitrary new walls.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=build();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['final_normals']),d['final_complete_ray_census']['counts'],
          'polar_closure_exhausted',d['exhausted_iterated_polar_closure'],flush=True)
