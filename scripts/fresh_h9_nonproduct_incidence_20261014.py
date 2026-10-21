#!/usr/bin/env python3
"""Bounded 16-facet incidence search with overlapping Lanner4 carriers.

Saved models are only incidence proposals. Exact source labels and
spherical predicates drive proposal generation; metric realization is
tested separately. No previous frozen output is modified.
"""
from pathlib import Path
from itertools import combinations, permutations
from collections import Counter
from fractions import Fraction as F
import argparse, hashlib, json
import z3
from fresh_h9_incidence_repair_20261014 import source, family, finite_labels
from fresh_h9_incidence_repair_20261014_verify import spherical

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fresh_h9_nonproduct_incidence_20261014.json'

def start(pattern):
    row,G,ps=source(); H,fixed,_=family(row,G,ps[pattern]);L=finite_labels(H)
    for i,j in fixed:L[i][j]=L[j][i]=2
    for r in L:r.append(2)
    L.append([2]*15+[1])
    return L

def cycles(L,require_new):
    out={}
    for B in combinations(range(16),4):
        if (15 in B)!=require_new:continue
        for tail in permutations(B[1:]):
            C=(B[0],)+tail
            if C[1]>C[-1]:continue
            edges={tuple(sorted((C[k],C[(k+1)%4])))for k in range(4)}
            if any(L[i][j]!=2 for i,j in combinations(B,2)if(i,j)not in edges):continue
            for high in edges:
                target={p:4 if p==high else 3 for p in edges}
                if any(L[i][j]not in(2,m)for(i,j),m in target.items()):continue
                additions=tuple(sorted((i,j,m)for(i,j),m in target.items()if L[i][j]==2))
                if additions:out[(B,tuple(sorted(target.items())))]=additions
    return [(B,dict(es),adds)for(B,es),adds in out.items()]

def recipes(pattern,limit):
    L=start(pattern);a=cycles(L,True);b=cycles(L,False);out={};examined=0
    a.sort(key=lambda r:(len(r[2]),r[2],r[0]));b.sort(key=lambda r:(len(r[2]),r[2],r[0]))
    for B,E,A in a:
        for C,F,Q in b:
            examined+=1
            if len(set(B)&set(C))not in(1,2):continue
            if any(p in E and E[p]!=m for p,m in F.items()):continue
            add=tuple(sorted(set(A)|set(Q)))
            # Explicit changes in existing old finite labels are not used in
            # this bounded family; all contact changes occur at current rights.
            if len(add)>6:continue
            T=[r[:]for r in L]
            for i,j,m in add:T[i][j]=T[j][i]=m
            if any(not spherical(T,D)for X in(B,C)for D in combinations(X,3)):continue
            # Reject an unrepairable affine triangle, since all three entries
            # are fixed finite contacts and no dotted choice can break it.
            if any(all(T[i][j]==3 for i,j in combinations(D,2))for D in combinations(range(16),3)):continue
            out.setdefault(add,[list(B),list(C)])
            if len(out)>=limit:return L,list(out.items()),examined
    ordered=sorted(out.items(),key=lambda z:(len(z[0]),z[0]))
    return L,ordered[:limit],examined

def levels(L):
    ans={0:[()]}
    for k in range(1,11):
        ans[k]=[]
        for B in ans[k-1]:
            for j in range(B[-1]+1 if B else 0,16):
                C=B+(j,)
                if spherical(L,C):ans[k].append(C)
    return ans

def missing(V,n=16):
    faces={k:set(B for W in V for B in combinations(W,k))for k in range(1,11)}
    out=[]
    for k in range(2,11):
        for B in combinations(range(n),k):
            if B not in faces[k]and all(C in faces[k-1]for C in combinations(B,k-1)):out.append(B)
    return out

def product_obstruction(M):
    # A compact Coxeter induced subdiagram with precisely these disjoint
    # minimal nonfaces is a forbidden times_0 product by Alexandrov Thm A.
    for R in combinations(M,4):
        if not any(len(A)==2 for A in R)or not any(len(A)>=3 for A in R):continue
        if any(set(A)&set(B)for A,B in combinations(R,2)):continue
        U=set().union(*map(set,R))
        if all(A in R or not set(A)<=U for A in M):return [list(A)for A in R]
    return None

def search(L,timeout=5000,rounds=12,pin_vertex=False):
    S=levels(L);rights=[p for p in combinations(range(16),2)if L[p[0]][p[1]]==2]
    pinned=next((B for B in S[9]if 15 in B),None)if pin_vertex else None
    pinned_rights=set(combinations(pinned,2))if pinned is not None else set()
    rights=[p for p in rights if p not in pinned_rights]
    X={p:z3.Bool('d_%d_%d'%p)for p in rights};solver=z3.Tactic('sat').solver();solver.set(timeout=timeout)
    NX={p:z3.Not(x)for p,x in X.items()}
    alive={}
    def face(B):
        if B not in alive:alive[B]=z3.And([NX[p]for p in combinations(B,2)if p in X])
        return alive[B]
    print('  spherical',len(S[9]),'PD9',len(S[8]),'PD8','variables',len(X),'pin',pinned,flush=True)
    contain={B:[]for k in range(1,9)for B in S[k]}
    for V in S[9]:
        for k in range(1,9):
            for B in combinations(V,k):contain[B].append(V)
    for B in S[10]:solver.add(z3.Not(face(B)))
    pending=[]
    for B,Vs in contain.items():
        if len(B)==8:solver.add(z3.Implies(face(B),z3.PbEq([(face(V),1)for V in Vs],2)if Vs else False))
        elif len(B)==1:solver.add(z3.Implies(face(B),z3.Or([face(V)for V in Vs])))
        else:pending.append((B,Vs))
    history=[];models=[]
    for step in range(rounds):
        status=solver.check();hist={'status':str(status)};history.append(hist)
        if status!=z3.sat:break
        model=solver.model();D={p for p,x in X.items()if z3.is_true(model.eval(x,model_completion=True))}
        V=[W for W in S[9]if not any(p in D for p in combinations(W,2))]
        bad=[(B,Ws)for B,Ws in pending if not any(p in D for p in combinations(B,2))and not any(not any(p in D for p in combinations(W,2))for W in Ws)]
        if bad:
            for B,Ws in bad:solver.add(z3.Implies(face(B),z3.Or([face(W)for W in Ws])))
            hist['unextendible_faces']=len(bad);continue
        M=missing(V);ob=product_obstruction(M)
        hist.update(vertices=len(V),missing_counts=dict(Counter(map(len,M))),product_obstruction=ob)
        if not ob and any(len(B)in(4,5)for B in M):
            models.append({'dotted_pairs':sorted(D),'vertices':V,'minimal_nonfaces':M});break
        solver.add(z3.Or([z3.Not(x)if p in D else x for p,x in X.items()]))
    return {'spherical_counts':{str(k):len(v)for k,v in S.items()},'pinned_vertex':pinned,'history':history,'models':models}

def main():
    p=argparse.ArgumentParser();p.add_argument('--patterns',type=int,default=2);p.add_argument('--per-pattern',type=int,default=8);p.add_argument('--timeout-ms',type=int,default=5000);p.add_argument('--pin-vertex',action='store_true');a=p.parse_args();rows=[]
    for pattern in range(a.patterns):
        L,R,total=recipes(pattern,a.per_pattern);print('recipes',pattern,total,flush=True)
        for ri,(adds,carriers)in enumerate(R):
            H=[r[:]for r in L]
            for i,j,m in adds:H[i][j]=H[j][i]=m
            r={'pattern':pattern,'recipe':ri,'recipe_pairs_examined':total,'additions':adds,'carriers':carriers,'labels':H,**search(H,a.timeout_ms,pin_vertex=a.pin_vertex)};rows.append(r)
            print('done',pattern,ri,'history',r['history'],'models',len(r['models']),flush=True)
            paths=[Path(__file__).resolve(),ROOT/'scripts/fresh_h9_incidence_repair_20261014.py',ROOT/'scripts/fresh_h9_incidence_repair_20261014_verify.py',ROOT/'results/fresh_h9_incidence_surgery_20261013.json',ROOT/'results/fresh_h9_incidence_surgery_20261013_dotted.json']
            OUT.write_text(json.dumps({'status':'BOUNDED_OVERLAPPING_LANNER4_INCIDENCE_SEARCH','settings':vars(a),'sources':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in paths},'cases':rows,'scope':'Sixteen facets, two prescribed overlapping Lanner4 cycles and additional right/dotted incidences. Satisfying models are abstract proposals only, not metric or compactness certificates.'},indent=2)+'\n')
if __name__=='__main__':main()
