#!/usr/bin/env python3
"""Exact mandatory A--R endpoint gate for the new four-cap pentagon disk.

The frozen physical interval engine is shared.  This successor reconstructs
both successive norm frames and covers every symbolic high-label interval.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import isqrt
from pathlib import Path
import fresh_new_cavity_construction_20261017 as N
import fresh_three_cap_pentagon_changed_seed_gate_20261017 as E

S,G,P,U,I,D,ROOT=N.S,N.G,E.P,E.U,E.I,E.D,N.ROOT
TAG='fresh_new_cavity_construction_gate_20261017'
RF=[8,44,63]
TARGETS=[4,30,53,9,10,7,8,44,63]

@lru_cache(None)
def theta(m):return D.theta(m)

def exclusion(value):
    if value.l>0:return dict(kind='strict_obtuse',product=value.enc())
    if value.h<=-2:return None
    for m in range(2,65):
        lower,upper=-theta(m+1),-theta(m)
        if lower.h<value.l and value.h<upper.l:
            return dict(kind='strict_consecutive_angle_gap',product=value.enc(),labels=[m,m+1],
                        lower_threshold=lower.enc(),upper_threshold=upper.enc())
    return None

def sqrt_closed(x):
    assert x.l>=0 and x.h>=x.l
    def lo(a):return F(isqrt(a.numerator*(1<<512)//a.denominator),1<<256)
    return I(lo(x.l),lo(x.h)+F(1,1<<256))

class Engine:
    def __init__(self):
        self.sp=G.Space(json.loads(G.CONTROL.read_text())['after64'])
        self.z=G.E.X.geometry();self.cache={}
        self.gram=E.exact_gram(self.sp,self.z,RF)
        self.inv,self.pivots=P.C53.inverse(self.gram)
        self.seed=30;cc=[self.pair(self.seed,i)for i in RF]
        self.kc=U.imv(self.inv,cc)
        self.kq=U.rnd(self.sp.nodes[self.seed]['q']-U.idot(cc,self.kc))
        assert self.kq.l>0
        self.kprod={i:U.rnd(self.pair(self.seed,i)-U.idot(self.kc,[self.pair(j,i)for j in RF]))for i in TARGETS}
        self.kprod[self.seed]=self.kq
        self.counts=Counter()

    def pair(self,a,b):
        if a==b:return self.sp.nodes[a]['q']
        key=tuple(sorted((a,b)))
        if key not in self.cache:
            v=U.rnd(G.E.X.literal(self.z,a,b)-G.white(self.sp.nodes[a]['cp'],self.sp.nodes[b]['cp']))
            assert G.overlap(v,self.sp.rp(a,b));self.cache[key]=v
        return self.cache[key]

    def evaluate(self,recipe,arlabel,box):
        """A leaf proves all branches impossible, or requests a finer box."""
        for parameter,value in box.items():
            gap=exclusion(-value)
            if gap:return dict(reason='assigned_high_contact_in_strict_angle_gap',parameter=parameter,rejection=gap)
        cp=G.profile(recipe['profile']);rq=U.rnd(2-G.white(cp,cp))
        rc=[-box['R'+str(j)]if j in recipe['high_positions']else-theta(m)for j,m in enumerate(recipe['labels'])]
        rc=[U.rnd(v-G.white(cp,self.sp.nodes[i]['cp']))for v,i in zip(rc,RF)]
        coeff=U.imv(self.inv,rc);remainder=U.rnd(rq-U.idot(rc,coeff))
        if remainder.h<0:return dict(reason='R_no_real_norm_branch',R_remaining_norm=remainder.enc())
        scale=sqrt_closed(I(max(F(0),remainder.l),remainder.h)/self.kq)
        ar=box['AR']if arlabel==6 else theta(arlabel)
        out=[]
        for rs in[-1,1]:
            rp={i:U.rnd(U.idot(coeff,[self.pair(j,i)for j in RF])+rs*scale*self.kprod[i])for i in TARGETS}
            for i,v in zip(RF,rc):rp[i]=v
            rbad=[]
            for i in TARGETS:
                if i in RF:continue
                v=U.rnd(rp[i]+G.white(cp,self.sp.nodes[i]['cp']))
                ex=exclusion(v)
                if ex:rbad.append(dict(wall=i,**ex))
            if rbad:
                out.append(dict(R_sign=rs,reason='R_native_obstruction',rejections=rbad));continue
            gram=[[self.sp.nodes[7]['q'],self.pair(7,8),rp[7]],
                  [self.pair(7,8),self.sp.nodes[8]['q'],rp[8]],
                  [rp[7],rp[8],rq]]
            try:inv,piv=P.C53.inverse(gram)
            except AssertionError:return None
            cov={i:[self.pair(7,i),self.pair(8,i),rp[i]]for i in TARGETS}
            seeds=[]
            for seed in[30,44,53,4,63]:
                kc=U.imv(inv,cov[seed]);kq=U.rnd(self.sp.nodes[seed]['q']-U.idot(cov[seed],kc))
                if kq.l>0:seeds.append((seed,kc,kq))
            if not seeds:return None
            seed,kc,kq=min(seeds,key=lambda r:float((r[2].h-r[2].l)/r[2].l))
            # The positive complementary line proves that this nondegenerate
            # three-frame is Lorentzian.  Only then does a negative remaining
            # norm rule out real roots.
            arem=U.rnd(2-ar*ar*inv[2][2])
            if arem.h<0:
                out.append(dict(R_sign=rs,reason='A_no_real_norm_branch',A_remaining_norm=arem.enc(),
                                A_frame_gram=G.encmat(gram),inverse_pivots=piv,kernel_seed=seed,kernel_square=kq.enc()));continue
            asc=sqrt_closed(I(max(F(0),arem.l),arem.h)/kq);abranches=[]
            for az in[-1,1]:
                bad=[]
                for i in TARGETS:
                    if i in[7,8]:continue
                    kb=kq if i==seed else U.rnd(self.pair(seed,i)-U.idot(kc,cov[i]))
                    value=U.rnd(-ar*U.idot(inv[2],cov[i])+az*asc*kb)
                    ex=exclusion(value)
                    if ex:bad.append(dict(wall=i,**ex))
                if not bad:return None
                abranches.append(dict(sign=az,rejections=bad))
            out.append(dict(R_sign=rs,reason='all_A_branches_excluded',A_frame_gram=G.encmat(gram),
                inverse_pivots=piv,A_remaining_norm=arem.enc(),kernel_seed=seed,kernel_square=kq.enc(),
                branches=abranches))
        return dict(R_remaining_norm=remainder.enc(),R_scale=scale.enc(),branches=out)

    def cover(self,recipe,arlabel,box,depth=0):
        self.counts['tested_boxes']+=1
        if self.counts['tested_boxes']%100==0:print('INTERVAL_BOXES',dict(self.counts),'depth',depth,flush=True)
        leaf=self.evaluate(recipe,arlabel,box)
        if leaf is not None:
            self.counts['certified_boxes']+=1
            self.counts['maximum_depth']=max(self.counts['maximum_depth'],depth)
            return dict(box={k:v.enc()for k,v in box.items()},certificate=leaf)
        assert box and depth<24,('unresolved_box',recipe,arlabel,{k:v.enc()for k,v in box.items()},depth)
        key=max(box,key=lambda k:box[k].h-box[k].l);v=box[key];mid=(v.l+v.h)/2
        children=[]
        for w in[I(v.l,mid),I(mid,v.h)]:
            new=dict(box);new[key]=w;children.append(self.cover(recipe,arlabel,new,depth+1))
        return dict(split_parameter=key,split_value=str(mid),children=children)

def run():
    engine=Engine();dom=S.domain(RF);assert len(dom['low'])==96 and len(dom['tails'])==10
    assert [(p,a,b)for p in range(7)for a in range(2,7)for b in range(2,7)if S.finite_corner(7,8,p,a,b)]==[(0,2,2)]
    assert not set(TARGETS)&set(S.DELETE+[49,62])
    rows=[];high=I(theta(6).l,F(2))
    for at,r in enumerate(dom['low']+dom['tails']):
        if at>=96:print('HIGH_R_FAMILY',at-96,r,flush=True)
        labels=[m for m in range(2,7)if N.finite_at8({'recipe':r},m)is not None]
        cases=[]
        for m in labels:
            box={'R'+str(j):high for j in r['high_positions']}
            if m==6:box['AR']=high
            proof=engine.cover(r,m,box)
            cases.append(dict(AR_label=m,AR_high_symbolic=m==6,finite_corner_type=N.finite_at8({'recipe':r},m),proof=proof))
        rows.append(dict(R_recipe=r,complete_AR_cases=cases))
        if(at+1)%10==0:print('FOUR_CENTER_EXACT',at+1,'of106',dict(engine.counts),flush=True)
    parent=json.loads((ROOT/'results'/f'{S.TAG}.json').read_text())
    paths=[Path(__file__),Path(N.__file__),Path(S.__file__),Path(E.__file__),Path(G.__file__),Path(G.E.X.__file__),G.CONTROL,G.PARENT,S.H.PARENT,S.H.BOUNDARY,ROOT/'results'/f'{S.TAG}.json']
    return dict(status='PASS_EXACT_ALL_LABEL_FOUR_CENTER_CHANGED_PENTAGON_ENDPOINT_GATE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        boundary_cycle=S.CYCLE,deleted_walls=S.DELETE,exact_old_deletion_cost=parent['exact_old_deletion_cost'],
        patch=dict(cap_arcs={'A':[7,8],'R':RF,'B':[63,30],'C':[30,7]},interior_cycle=['A','R','B','C'],
                   central_diagonals=[['A','B'],['R','C']],mandatory_corner=[8,'A','R']),
        R_frame_gram=G.encmat(engine.gram),R_inverse_pivots=engine.pivots,R_kernel_seed=engine.seed,R_kernel_square=engine.kq.enc(),
        counts=dict(engine.counts,R_low_recipes=96,R_high_families=10,AR_cases=sum(len(r['complete_AR_cases'])for r in rows)),records=rows,
        completeness='The old finite H4 component on A2+7+8 forces A right to all S6,7,8. R must give both finite corners S6+8+44+R and S6+44+63+R. The complete rank<=5 finite classification gives96 low R recipes and10 isolated high-label families. The mandatory corner S6+8+A+R gives every possible A-R label. Each symbolic high cosine is relaxed to the closed real interval[sqrt3,2]. Every R norm branch and then every A norm branch is covered by the recorded interval boxes, with strict native obtuse or consecutive-angle-gap witnesses.',
        scope='No compatible A,R pair exists for this specified four-center disk with the retained walls. Therefore neither central diagonal can realize the disk. No exclusion of other four-cap incidence templates, caps with zero boundary arcs, further wall deletion/deformation, or general finite-volume completion. The physical interval engine and prior finite classification are shared; no new chamber or global volume claim.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(run()));dest=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert out==json.loads(dest.read_text())
    else:dest.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
