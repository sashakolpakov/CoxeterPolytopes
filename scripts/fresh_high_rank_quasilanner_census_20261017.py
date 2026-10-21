#!/usr/bin/env python3
"""Exact finite/affine graph census of minimal indefinite diagrams of orders6..10."""
import argparse, hashlib, json
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import sympy as s
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C
import fresh_nikulin_two_budget_frontier_20261015_review as A
import fresh_sqrt5_rank14_core_spectator_domain_20261017 as I

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_high_rank_quasilanner_census_20261017'

def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def restrict(g, ids):
    ids=set(ids)
    return {i:{j:m for j,m in g[i].items() if j in ids} for i in sorted(ids)}
def edges(g): return [[i,j,m] for i,row in g.items() for j,m in sorted(row.items()) if i<j]
def gram(g):
    ids=sorted(g)
    return s.Matrix([[2 if i==j else -({2:0,3:1,4:s.sqrt(2)}[g[i].get(j,2)]) for j in ids] for i in ids])

def run():
    finite={C.code(g):name for name,g in C.models() if len(g)<=9}
    affine={}; affine_models=[]
    for r in range(2,9):
        names=['A','B','C']+(['D'] if r>=4 else [])
        names += ['E'+str(r)] if r in (6,7,8) else []
        names += ['F4'] if r==4 else []
        for name in names:
            m=A.model(name,r);G=m['G'];g={i:{} for i in range(r+1)}
            for i,j in combinations(g,2):
                if not G[i,j]: continue
                label=3 if G[i,j]==-1 else 4
                assert G[i,j] in (-1,-s.sqrt(2))
                g[i][j]=g[j][i]=label
            full=name if name.startswith(('E','F')) else name+str(r)
            if sum(map(len,g.values()))==2*len(g)-2:
                key=C.code(g)
                if key in affine: assert affine[key]['rank']==r
                else: affine[key]=dict(name='affine_'+full,rank=r)
            if r>=4: affine_models.append((full,g))

    def component(g):
        n=len(g);e=sum(map(len,g.values()))//2
        if e==n-1:
            key=C.code(g)
            if key in finite: return dict(kind='finite',name=finite[key],rank=n)
            if key in affine: return dict(kind='affine',**affine[key])
        if e==n and n>=3 and all(len(row)==2 and set(row.values())=={3} for row in g.values()):
            return dict(kind='affine',name='affine_A'+str(n-1),rank=n-1)
        return None

    def deletion(g,v):
        parts=C.pieces(restrict(g,set(g)-{v}))
        out=[component(p) for p in parts]
        if None in out: return None
        return sorted(out,key=lambda z:(z['kind'],z['name']))

    diagrams=[];graphs=[];colors=[];buckets=defaultdict(list)
    provenance=[];counts=Counter();censuses=[];trace=hashlib.sha256()
    for name,base in affine_models:
        n=len(base)+1;v=n-1;local=Counter()
        for contacts in product((2,3,4),repeat=n-1):
            if all(x==2 for x in contacts): continue
            counts['attachment_candidates']+=1;local['candidates']+=1
            g={i:dict(row) for i,row in base.items()};g[v]={}
            for i,m in enumerate(contacts):
                if m!=2:g[i][v]=g[v][i]=m
            classified=[];failure=None
            for removed in range(n):
                row=deletion(g,removed)
                if row is None: failure=removed;break
                classified.append(row)
            trace.update((name+':'+','.join(map(str,contacts))+':'+str(failure)+'\n').encode())
            if failure is not None: continue
            counts['retained_attachments']+=1;local['retained']+=1
            cg,key=I.refinement(g);index=None;mapping=None
            for ci in buckets[key]:
                mp=I.isomorphism(g,graphs[ci],cg,colors[ci])
                if mp is not None:index=ci;mapping=mp;break
            if index is None:
                index=len(diagrams);mapping={i:i for i in g}
                graphs.append(g);colors.append(cg);buckets[key].append(index)
                ideal=[i for i,row in enumerate(classified) if any(p['kind']=='affine' for p in row)]
                assert ideal and all(len(classified[i])==1 for i in ideal)
                G=gram(g);det=s.simplify(G.det());assert det<0,(n,edges(g),det)
                # Proper principal matrices are PSD by exact component types.
                # Interlacing plus negative determinant gives precisely one negative eigenvalue.
                thresholds=[]
                if len(ideal)==1:
                    j=ideal[0]
                    for b in range(n):
                        if b==j:continue
                        ids=[i for i in range(n) if i!=b]
                        inv=G.extract(ids,ids).inv();q=s.simplify(1/inv[ids.index(j),ids.index(j)])
                        assert q>0 and q.is_Rational
                        thresholds.append([b,str(q)])
                diagrams.append(dict(index=index,order=n,edges=edges(g),determinant=str(det),
                    deletions=classified,ideal_deletions=ideal,finite_deletions=[i for i in range(n) if i not in ideal],
                    sole_ideal_diagonal_decrease_thresholds=thresholds))
            provenance.append(dict(affine_source=name,contacts=list(contacts),diagram=index,
                isomorphism=[[i,mapping[i]] for i in range(n)]))
        censuses.append(dict(affine_source=name,order=n,**dict(local)))
        print('AFFINE_ATTACHMENT_CENSUS',name,n,dict(local),'unique',len(diagrams),flush=True)
    by_order=Counter(d['order'] for d in diagrams)
    # Independent recorded comparison to Chein's printed tables; these counts
    # are not used to prune or accept any graph.
    assert dict(sorted(by_order.items()))=={6:12,7:3,8:4,9:4,10:3}
    counts.update(unique_diagrams=len(diagrams),affine_source_models=len(affine_models))
    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'scripts/{A.TAG}.py',f'scripts/{I.TAG}.py']
    return dict(status='PASS_EXACT_HIGH_RANK_QUASILANNER_ATTACHMENT_CENSUS',sources={p:sha(p) for p in paths},
        counts=dict(counts),orders=[[n,c] for n,c in sorted(by_order.items())],affine_censuses=censuses,
        diagrams=diagrams,attachment_provenances=provenance,complete_rejection_trace_sha256=trace.hexdigest(),
        primary_reference=dict(title='M. Chein (1969), Recherche des graphes des matrices de Coxeter hyperboliques',
            url='https://www.numdam.org/item/M2AN_1969__3_3_3_0.pdf',table_pages='printed15–16'),
        completeness='For orders6..10 every minimal indefinite diagram is noncompact, contains an affine connected subdiagram, and minimality makes it an affine diagram plus one adjacent node. Every proper connected5-node graph is finite/affine, excluding labels5 or larger; parallel and dotted bonds are likewise excluded. Enumerate every nonzero attachment over labels2,3,4 to every affine source of the required order, then classify all node deletions exactly by finite/affine components. Finite/affine classification and the order10 ceiling are stated published inputs.',
        scope='Exact higher-order diagram census and vertex-deletion/diagonal-threshold data. It does not identify a general induced simplex face with a Coxeter simplex, settle external white attachments for arbitrary black subsets, or bound ambient Coxeter dimension.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],d['orders'],flush=True)
