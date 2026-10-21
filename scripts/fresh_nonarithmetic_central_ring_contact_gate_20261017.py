#!/usr/bin/env python3
"""Exact finite-contact obstruction for the twenty literal central rings.

All low contacts use the unique four-contact inverse-Gram norm equation.
Every possible high contact is isolated and solved without a label cutoff.
The rings are exact recipe objects; this does not promote their ambient gates.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_replace47_close12_boundary_20261017 as E
import fresh_nonarithmetic_release47_replace12_frontier_20261017 as F

P,U,I,D=E.P,E.U,E.I,E.D
ROOT=E.ROOT
TAG='fresh_nonarithmetic_central_ring_contact_gate_20261017'
CONTROL=ROOT/'results'/f'{E.TAG}.json'
FRONTIER=ROOT/'results'/f'{F.TAG}.json'
RINGS=ROOT/'results/fresh_nonarithmetic_release47_nonzero_a2_bridge_20261017.json'


def decode(a):return [I(Q(lo),Q(hi)) for lo,hi in a]
def encmat(a):return [P.enc(row) for row in a]
def mid(a):return float((a.l+a.h)/2)
def overlap(a,b):return a.l<=b.h and b.l<=a.h
def white(a,b):return U.rnd((2*a[0]*b[0]+a[0]*b[1]+a[1]*b[0]+2*a[1]*b[1])/3)
def combine(c,vec):return [U.idot(c,[v[j] for v in vec]) for j in range(4)]


class Space:
    def __init__(self,control):
        self.G=[decode(r) for r in control['residual_gram']]
        _,_,num,_,_=F.A.setup()
        self.num_basis=np.array([num[i] for i in E.FRAME])
        self.nodes={};self.literal={};self.pair_cache={};self.recipes=[]
        for r in control['projections']:
            i=r['wall']
            if i not in [7,8,30,44,52]:continue
            cp=[I(),I(-1)] if i==7 else [I(),I()]
            vec=decode(r['coordinates'])
            self.nodes[str(i)]=dict(name=str(i),v=vec,cp=cp,q=U.rnd(2-white(cp,cp)))
        for a,b,m in [(7,8,5),(7,30,5),(8,52,2),(30,52,4)]:self.literal[frozenset([str(a),str(b)])]=-D.theta(m)
    def pair(self,a,b):return U.idot(a,U.imv(self.G,b))
    def rp(self,a,b):
        if a['name']==b['name']:return a['q']
        x=self.literal.get(frozenset([a['name'],b['name']]))
        key=frozenset([a['name'],b['name']])
        direct=self.pair_cache.get(key)
        if direct is None:
            direct=self.pair(a['v'],b['v']);self.pair_cache[key]=direct
        if x is None:return direct
        val=U.rnd(x-white(a['cp'],b['cp']))
        assert overlap(val,direct),(a['name'],b['name'],val.enc(),direct.enc())
        return val
    def profile(self,index):
        if index==0:return [I(),I()]
        node=0 if index<=3 else 1;m=3+(index-1)%3
        out=[I(),I()];out[node]=-D.theta(m);return out
    def solutions(self,frame,profile,labels):
        cp=self.profile(profile);q=U.rnd(2-white(cp,cp))
        gram=[[self.rp(a,b) for b in frame] for a in frame]
        inv,pivots=P.C53.inverse(gram)
        contacts=[U.rnd(-D.theta(m)-white(cp,n['cp'])) for m,n in zip(labels,frame)]
        coeff=U.imv(inv,contacts);part=combine(coeff,[a['v'] for a in frame])
        residual=U.rnd(q-U.idot(contacts,coeff))
        kernels=[]
        for key,seed in self.nodes.items():
            if key in {a['name'] for a in frame}:continue
            ss=[self.rp(seed,a) for a in frame]
            c=U.imv(inv,ss);v=[U.rnd(a-b) for a,b in zip(seed['v'],combine(c,[a['v'] for a in frame]))]
            sq=U.rnd(seed['q']-U.idot(ss,c))
            if sq.l>0:kernels.append((sq,key,v))
        assert kernels,'positive one-dimensional kernel required'
        ks,ki,kv=min(kernels,key=lambda t:float((t[0].h-t[0].l)/t[0].l))
        base=dict(frame=[a['name'] for a in frame],profile=profile,labels=list(labels),
            gram=encmat(gram),inverse_pivots=pivots,contacts=P.enc(contacts),
            particular=P.enc(part),kernel_seed=ki,kernel=P.enc(kv),
            kernel_square=ks.enc(),residual_square=residual.enc())
        if residual.h<0:return [],base
        assert residual.l>0,('zero/unresolved norm branch',base)
        scale=U.sqrt(residual/ks);out=[]
        for sign in [-1,1]:
            v=[U.rnd(a+sign*scale*b) for a,b in zip(part,kv)]
            assert overlap(self.pair(v,v),q)
            for a,c in zip(frame,contacts):assert overlap(self.pair(v,a['v']),c)
            out.append((sign,v))
        base['positive_scale']=scale.enc()
        return out,base
    def solution_pair(self,frame,recipe,sign,target):
        inv,_=P.C53.inverse([[self.rp(a,b) for b in frame] for a in frame])
        contacts=decode(recipe['contacts']);coeff=U.imv(inv,contacts)
        seed=self.nodes[recipe['kernel_seed']]
        seed_products=[self.rp(seed,a) for a in frame];seed_coeff=U.imv(inv,seed_products)
        target_products=[self.rp(target,a) for a in frame]
        particular=U.idot(coeff,target_products)
        if target['name']==seed['name']:kernel=I(*map(Q,recipe['kernel_square']))
        elif target['name'] in {a['name'] for a in frame}:kernel=I()
        else:kernel=U.rnd(self.rp(seed,target)-U.idot(seed_coeff,target_products))
        scale=I(*map(Q,recipe['positive_scale']))
        return U.rnd(particular+sign*scale*kernel)
    def add(self,name,frame,profile,labels,target):
        sols,proof=self.solutions(frame,profile,labels);assert len(sols)==2
        errors=[float(np.max(abs(np.array([mid(q) for q in v])@self.num_basis-np.array(target)))) for s,v in sols]
        selected=min(range(2),key=lambda i:errors[i]);assert errors[selected]<1e-5 and errors[1-selected]>1e-3,(name,errors)
        sign,vec=sols[selected];cp=self.profile(profile)
        obj=dict(name=name,v=vec,cp=cp,q=U.rnd(2-white(cp,cp)))
        scalar_products={key:self.solution_pair(frame,proof,sign,target) for key,target in self.nodes.items()}
        for key,value in scalar_products.items():self.pair_cache[frozenset([name,key])]=value
        self.nodes[name]=obj
        for n,m in zip(frame,labels):self.literal[frozenset([name,n['name']])]=-D.theta(m)
        proof.update(name=name,chosen_sign=sign,branch_discovery_errors=errors,coordinates=P.enc(vec),white_contacts=P.enc(cp),
            factored_projected_products={k:v.enc() for k,v in scalar_products.items()})
        self.recipes.append(proof);return obj


def run():
    control=json.loads(CONTROL.read_text())['after64']
    front=json.loads(FRONTIER.read_text());rings=json.loads(RINGS.read_text())
    records=[];total=Counter()
    for index,st in enumerate(rings['stages']):
        if st['best_coherent'] is None:continue
        sp=Space(control)
        group=next(g for g in front['changed_seed_pools'] if front['replacement_trials'][g['replacement']]['label_to52']==st['replacement_label'])
        replacement=front['replacement_trials'][group['replacement']]
        n12=sp.add('12new',[sp.nodes[str(i)] for i in [7,8,52]],0,[2,2,st['replacement_label']],replacement['projection'])
        seed=group['seeds'][st['seed']];patch=group['patches'][st['seed']]
        p=patch['P'][st['P']];q=patch['Q'][st['Q']];r=st['best_coherent']['R']
        assert seed['rotation']==0
        L=sp.add('L',[sp.nodes['7'],n12,sp.nodes['52']],seed['profile'],seed['labels'],seed['projection'])
        V=sp.add('P',[sp.nodes['52'],sp.nodes['30'],L],p['profile'],p['labels'],p['projection'])
        W=sp.add('Q',[sp.nodes['30'],sp.nodes['7'],L],q['profile'],q['labels'],q['projection'])
        R=sp.add('R',[sp.nodes['30'],V,W],r['profile'],r['labels'],r['projection'])
        ring=[L,V,R,W];rawring=[seed,p,r,q]
        gram=[[sp.rp(a,b) for b in ring] for a in ring];inv,pivots=P.C53.inverse(gram)
        low=[]
        for profile in range(7):
            cp=sp.profile(profile);square=U.rnd(2-white(cp,cp));counts=Counter();minimum=None;digest=sha256()
            shifts=[white(cp,n['cp']) for n in ring]
            for labels in product(range(2,6),repeat=4):
                cc=[U.rnd(-D.theta(m)-s) for m,s in zip(labels,shifts)]
                diff=U.rnd(U.idot(cc,U.imv(inv,cc))-square)
                assert diff.l>0 or diff.h<0,('unresolved low contact norm',index,profile,labels,diff.enc())
                sign=1 if diff.l>0 else -1;counts[str(sign)]+=1
                lower=diff.l if sign==1 else -diff.h
                minimum=lower if minimum is None else min(minimum,lower)
                digest.update(json.dumps([labels,diff.enc()],separators=(',',':')).encode())
            assert sum(counts.values())==256
            low.append(dict(profile=profile,contact_quadruples=256,sign_counts=dict(counts),
                minimum_absolute_norm_error=str(minimum),ordered_interval_digest=digest.hexdigest()))
        edge=[p['labels'][2],r['labels'][1],r['labels'][2],q['labels'][2]]
        eligible=[i for i in range(4) if rawring[i]['profile']==0 and edge[i-1]==edge[i]==2]
        assert not any(i in eligible and (i+2)%4 in eligible for i in range(4))
        high=[]
        for j in eligible:
            other=[i for i in range(4) if i!=j]
            for m in [2,3,4,5]:
                labels=[m if i==(j+2)%4 else 2 for i in other]
                sols,recipe=sp.solutions([ring[i] for i in other],0,labels)
                for sign,v in sols:
                    val=sp.solution_pair([ring[i] for i in other],recipe,sign,ring[j])
                    gap=U.rnd(val+D.theta(6))
                    assert gap.l>0,('high contact range unresolved',index,j,m,sign,val.enc())
                    high.append(dict(high_vertex=j,opposite_label=m,branch=sign,
                        recipe=recipe,product=val.enc(),gap_above_minus_sqrt3=gap.enc()))
        total['rings']+=1;total['low_contact_quadruples']+=7*256;total['high_real_branches']+=len(high)
        records.append(dict(parent_stage=index,replacement_label=st['replacement_label'],
            seed=st['seed'],P=st['P'],Q=st['Q'],recipes=sp.recipes,
            residual_ring_gram=encmat(gram),inverse_pivots=pivots,low_contact_exclusions=low,
            ring_edge_labels=edge,eligible_high_vertices=eligible,high_contact_exclusions=high))
        print('EXACT_RING',index,'low',1792,'high',len(high),flush=True)
    assert total==Counter(rings=20,low_contact_quadruples=35840,high_real_branches=64)
    paths=[Path(__file__),Path(E.__file__),Path(E.X.__file__),CONTROL,FRONTIER,RINGS,
           Path(F.__file__),Path(F.A.__file__)]
    return dict(status='PASS_EXACT_ALL_LABEL_CENTRAL_CAP_OBSTRUCTION_FOR20_LITERAL_RINGS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        counts=dict(total),records=records,
        finite_profile_completeness='An actual S6 facet is right to H4 and has one of the seven A2 profiles. A bond m>=6 is isolated in every required finite corner: with any other incident nonright bond its3x3 determinant is at most8-2(3+1)=0. The literal eligible high vertices contain no opposite pair, so at most one central high bond is possible. The other three contacts are two right and one label2..5. All-low assignments are exhausted by the4x4 inverse-Gram norm equation; all high assignments are solved by their other three contacts, without a bound on the high label.',
        scope='Exact obstruction to one actual S6 cap making all four consecutive ring-edge endpoints finite, for the20 specified exact recipe rings. No claim that these numerical parent rings pass every ambient native gate, and no obstruction to more walls, changed rings, caps missing S6 or a different patch topology. No old frozen geometry is modified.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
