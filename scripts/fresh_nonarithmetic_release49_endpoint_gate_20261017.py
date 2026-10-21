#!/usr/bin/env python3
"""All-label mandatory endpoint gates for the release49 hexagonal patch.

The parent supplies exact residual coordinates. All finite rank<=5 diagrams
are tested combinatorially, then each real norm branch has a strict obtuse
retained-wall certificate. No numerical branch acceptance is used.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

import fresh_nonarithmetic_replace47_close12_boundary_20261017 as E

ROOT,P,U,I,D=E.ROOT,E.P,E.U,E.I,E.D
TAG='fresh_nonarithmetic_release49_endpoint_gate_20261017'
CONTROL=ROOT/'results'/f'{E.TAG}.json'
PARENT=ROOT/'results'/f'{E.X.TAG}.json'
FRAMES=[[60,58,7],[58,7,33]]
def decode(v):return [I(Q(a),Q(b))for a,b in v]
def encmat(m):return [P.enc(v)for v in m]
def overlap(a,b):return a.l<=b.h and b.l<=a.h
def white(a,b):return U.rnd((2*a[0]*b[0]+a[0]*b[1]+a[1]*b[0]+2*a[1]*b[1])/3)
def combine(c,vs):return [U.idot(c,[v[j]for v in vs])for j in range(4)]
def profile(i):
    if i==0:return [I(),I()]
    a=[I(),I()];a[0 if i<=3 else 1]=-D.theta(3+(i-1)%3);return a
def profile_labels(i):
    if i==0:return [2,2]
    a=[2,2];a[0 if i<=3 else 1]=3+(i-1)%3;return a


def canonical(g):
    n=len(g)
    return min(tuple(g[p[i]][p[j]]for i in range(n)for j in range(i+1,n))for p in permutations(range(n)))


def finite_catalogue():
    out={}
    def add(name,n,edges):
        g=[[2]*n for _ in range(n)]
        for i in range(n):g[i][i]=1
        for a,b,m in edges:g[a][b]=g[b][a]=m
        out[n,canonical(g)]=name
    for n in range(1,6):add(f'A{n}',n,[(i,i+1,3)for i in range(n-1)])
    for n in range(2,6):add(f'B{n}',n,[(i,i+1,4 if i==0 else 3)for i in range(n-1)])
    for n in (4,5):add(f'D{n}',n,[(0,1,3),(0,2,3),(0,3,3)]+([(3,4,3)]if n==5 else []))
    add('F4',4,[(0,1,3),(1,2,4),(2,3,3)])
    add('I2(5)',2,[(0,1,5)])
    add('H3',3,[(0,1,5),(1,2,3)])
    add('H4',4,[(0,1,5),(1,2,3),(2,3,3)])
    return out


CAT=finite_catalogue()
def finite_components(g):
    left=set(range(len(g)));answer=[]
    while left:
        comp=set();todo=[min(left)]
        while todo:
            a=todo.pop()
            if a in comp:continue
            comp.add(a);todo.extend(b for b in left if g[a][b]>2 and b not in comp)
        left-=comp;ids=sorted(comp)
        key=(len(ids),canonical([[g[i][j]for j in ids]for i in ids]))
        if key not in CAT:return None
        answer.append(CAT[key])
    return sorted(answer)


class Space:
    def __init__(self,control):
        self.G=[decode(r)for r in control['residual_gram']]
        z=E.X.geometry();self.nodes={};self.pair_cache={}
        for row in control['projections']:
            i=row['wall'];cp=[E.X.literal(z,i,j)for j in (5,6)]
            self.nodes[i]=dict(v=decode(row['coordinates']),q=I(*map(Q,row['square'])),cp=cp)
        # Exact contacts are inherited from the frozen native-gate packages.
        for i,pi in [(7,4),(33,0),(58,0),(60,6),(8,0),(30,0),(44,0),(52,0)]:
            cp=profile(pi);assert all(overlap(a,b)for a,b in zip(cp,self.nodes[i]['cp']))
            self.nodes[i]['cp']=cp;self.nodes[i]['q']=U.rnd(2-white(cp,cp))
        self.literal={frozenset((60,58)):I(),frozenset((58,7)):I(-1),frozenset((7,33)):I(-1)}
    def pair(self,a,b):return U.idot(a,U.imv(self.G,b))
    def rp(self,a,b):
        if a==b:return self.nodes[a]['q']
        key=frozenset((a,b))
        if key in self.pair_cache:return self.pair_cache[key]
        direct=self.pair(self.nodes[a]['v'],self.nodes[b]['v'])
        if key in self.literal:
            value=U.rnd(self.literal[key]-white(self.nodes[a]['cp'],self.nodes[b]['cp']))
            assert overlap(value,direct);direct=value
        self.pair_cache[key]=direct;return direct
    def recipe(self,ids,pi,ms):
        cp=profile(pi);square=U.rnd(2-white(cp,cp))
        gram=[[self.rp(a,b)for b in ids]for a in ids];inv,pivots=P.C53.inverse(gram)
        contacts=[U.rnd(-D.theta(m)-white(cp,self.nodes[i]['cp']))for i,m in zip(ids,ms)]
        coeff=U.imv(inv,contacts);part=combine(coeff,[self.nodes[i]['v']for i in ids])
        residual=U.rnd(square-U.idot(contacts,coeff))
        seeds=[]
        for seed in (8,30,44,52):
            if seed in ids:continue
            ss=[self.rp(seed,i)for i in ids];cc=U.imv(inv,ss)
            ks=U.rnd(self.nodes[seed]['q']-U.idot(ss,cc))
            if ks.l>0:seeds.append((ks,seed,cc))
        assert seeds,'frame has one positive orthogonal dimension'
        ks,seed,kc=min(seeds,key=lambda t:float((t[0].h-t[0].l)/t[0].l))
        kernel=[U.rnd(a-b)for a,b in zip(self.nodes[seed]['v'],combine(kc,[self.nodes[i]['v']for i in ids]))]
        proof=dict(frame=ids,profile=pi,labels=list(ms),gram=encmat(gram),inverse_pivots=pivots,
            contacts=P.enc(contacts),particular=P.enc(part),kernel_seed=seed,kernel=P.enc(kernel),
            kernel_square=ks.enc(),residual_square=residual.enc())
        if residual.h<0:return proof,[]
        assert residual.l>0,('unresolved tangent',proof)
        scale=U.sqrt(residual/ks);proof['positive_scale']=scale.enc();branches=[]
        for sign in (-1,1):
            root=[U.rnd(a+sign*scale*b)for a,b in zip(part,kernel)]
            assert overlap(self.pair(root,root),square)
            # Factored pair products prevent cancellation in root coordinates.
            products={}
            for target in sorted(self.nodes):
                if target==49:continue
                if target in ids:val=-D.theta(ms[ids.index(target)])
                else:
                    tt=[self.rp(target,i)for i in ids]
                    kk=ks if target==seed else U.rnd(self.rp(seed,target)-U.idot(kc,tt))
                    val=U.rnd(U.idot(coeff,tt)+sign*scale*kk+white(cp,self.nodes[target]['cp']))
                products[target]=val
            bad=[(i,v)for i,v in products.items()if v.l>0]
            assert bad,('acute_or_unresolved_branch',ids,pi,ms,sign)
            target,val=max(bad,key=lambda item:item[1].l)
            branches.append(dict(sign=sign,coordinates=P.enc(root),obtuse_retained_wall=target,
                strict_positive_product=val.enc(),positive_retained_products=[[i,v.enc()]for i,v in bad]))
        return proof,branches


def corner_graph(ids,pi,ms,pair_indices):
    a,b=pair_indices;anchors=[ids[a],ids[b]]
    old_profiles={7:4,33:0,58:0,60:6}
    g=[[2]*5 for _ in range(5)]
    for i in range(5):g[i][i]=1
    g[0][1]=g[1][0]=3
    for j,node in enumerate(anchors):
        for i,m in enumerate(profile_labels(old_profiles[node])):g[i][2+j]=g[2+j][i]=m
    for i,m in enumerate(profile_labels(pi)):g[i][4]=g[4][i]=m
    old=2 if set(anchors)=={60,58}else 3
    g[2][3]=g[3][2]=old
    for j,idx in enumerate((a,b)):g[2+j][4]=g[4][2+j]=ms[idx]
    return g


def run():
    control=json.loads(CONTROL.read_text())['after64'];sp=Space(control)
    records=[];tot=Counter()
    for parity,ids in enumerate(FRAMES):
        trials=[];recipes=[];counts=Counter()
        for pi in range(7):
            for ms in product(range(2,6),repeat=3):
                types=[finite_components(corner_graph(ids,pi,ms,ab))for ab in ((0,1),(1,2))]
                if any(t is None for t in types):
                    counts['nonfinite_contact_recipes']+=1
                    trials.append(dict(profile=pi,labels=list(ms),finite_corner_types=types,status='nonfinite_corner'))
                    continue
                proof,branches=sp.recipe(ids,pi,ms)
                counts['finite_corner_recipes']+=1;counts['obtuse_real_branches']+=len(branches)
                if not branches:counts['negative_norm_residual']+=1
                trials.append(dict(profile=pi,labels=list(ms),finite_corner_types=types,status='obtuse_branches'if branches else'negative_norm_residual',recipe=len(recipes)))
                recipes.append(dict(**proof,branches=branches))
        assert sum(counts[k]for k in ('nonfinite_contact_recipes','finite_corner_recipes'))==448
        assert counts['finite_corner_recipes']==[5,13][parity]
        assert counts['obtuse_real_branches']==[10,26][parity]
        assert counts['negative_norm_residual']==0
        records.append(dict(parity=parity,frame=ids,required_finite_endpoint_pairs=[ids[:2],ids[1:]],
            high_contact_blockers=[dict(anchor=i,old_neighbor=6 if i in (7,60)else 7,label=5 if i==60 else 3)for i in ids],
            counts=dict(counts),finite_contact_trials=trials,exact_recipes=recipes))
        tot.update(counts);print('ENDPOINT_GATE',ids,dict(counts),flush=True)
    paths=[Path(__file__),Path(E.__file__),Path(E.X.__file__),Path(E.X.C.__file__),CONTROL,PARENT]
    return dict(status='PASS_EXACT_ALL_LABEL_RELEASE49_HEXAGON_MANDATORY_ENDPOINT_GATES',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        white_support=E.S,deleted_old_wall=49,retained_old_walls=[i for i in control['labels']if i!=49],
        old_cycle=[7,33,30,59,60,58],counts=dict(tot),records=records,
        finite_connected_catalogue=[dict(rank=n,code=list(code),type=name)for(n,code),name in sorted(CAT.items())],
        completeness='Every actual S6 facet is right to H4 and has one of the seven finite A2 profiles. Each frame anchor has a fixed nonright neighbor in a required finite corner. A new bond m>=6 gives a3x3 principal determinant at most8-2(3+1)=0, including its nonpositive cycle term, so all three labels belong to2..5. Each of448 profiles is tested against the complete finite rank<=5 classification. The3-contact frame plus norm determines at most two real roots, both tested with exact retained-wall obstructions.',
        scope='No acute actual S6 facet can make both indicated old edges finite after deleting49 while retaining the other63 old walls. This excludes both cyclic placements of the specified three-pentagon hexagon patch. No obstruction to a different patch topology, more deletions, or caps missing S6; no new chamber or finite-volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
