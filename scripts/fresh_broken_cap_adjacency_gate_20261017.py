#!/usr/bin/env python3
"""All-label A--X gate for the two quadrilateral disks sharing X.

Reuses the frozen nested norm/interval engine with a different literal frame
and a new complete finite-corner domain. No numerical scout is an acceptance
or exclusion oracle.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import fresh_new_cavity_construction_gate_20261017 as Q
import fresh_broken_cap_adjacency_20261017 as T

ROOT,G,P,U,I=Q.ROOT,Q.G,Q.P,Q.U,Q.I
TAG='fresh_broken_cap_adjacency_gate_20261017'
# This module instance is private to this process. Frozen source bytes are
# unchanged. The inherited evaluator uses RF as its first literal frame.
Q.RF=[8,30,63]

class Engine(Q.Engine):
    def __init__(self):
        self.sp=G.Space(json.loads(G.CONTROL.read_text())['after64'])
        self.z=G.E.X.geometry();self.cache={}
        self.gram=Q.E.exact_gram(self.sp,self.z,Q.RF)
        self.inv,self.pivots=P.C53.inverse(self.gram)
        self.seed=44;cc=[self.pair(self.seed,i)for i in Q.RF]
        self.kc=U.imv(self.inv,cc)
        self.kq=U.rnd(self.sp.nodes[self.seed]['q']-U.idot(cc,self.kc));assert self.kq.l>0
        self.kprod={i:U.rnd(self.pair(self.seed,i)-U.idot(self.kc,[self.pair(j,i)for j in Q.RF]))for i in Q.TARGETS}
        self.kprod[self.seed]=self.kq;self.counts=Counter()

def domain():
    rows={}
    for pi,m8,m30,m63,ma in product(range(7),range(2,7),range(2,7),range(2,7),range(2,7)):
        c1=T.corner(pi,m30,m63,0,2);c2=T.corner(pi,m8,ma,0,2)
        if c1 is None or c2 is None:continue
        key=pi,m8,m30,m63
        if key not in rows:rows[key]=dict(profile=pi,labels=[m8,m30,m63],
            high_positions=[j for j,m in enumerate([m8,m30,m63])if m==6],old_corner_type=c1,A_cases=[])
        rows[key]['A_cases'].append(dict(label=ma,corner_type=c2))
    return list(rows.values())

def run():
    engine=Engine();dom=domain();assert len(dom)==148 and sum(len(r['A_cases'])for r in dom)==344
    assert [(p,a,b)for p in range(7)for a in range(2,7)for b in range(2,7)if T.S.finite_corner(7,8,p,a,b)]==[(0,2,2)]
    assert not set(Q.TARGETS)&set(T.S.DELETE+[49,62])
    high=I(Q.theta(6).l,F(2));records=[]
    for index,r in enumerate(dom):
        cases=[]
        for ar in r['A_cases']:
            m=ar['label'];box={'R'+str(j):high for j in r['high_positions']}
            if m==6:box['AR']=high
            assert len(box)<=2
            proof=engine.cover(r,m,box)
            cases.append(dict(AX_label=m,symbolic_high=m==6,finite_corner_type=ar['corner_type'],proof=proof))
        records.append(dict(X_recipe={k:v for k,v in r.items()if k!='A_cases'},A_cases=cases))
        if(index+1)%10==0:print('BROKEN_EXACT',index+1,'of148',dict(engine.counts),flush=True)
    parent=json.loads((ROOT/'results'/f'{T.S.TAG}.json').read_text())
    paths=[Path(__file__),Path(Q.__file__),ROOT/'results'/f'{Q.TAG}.json',Path(T.__file__),Path(T.S.__file__),Path(G.__file__),
           Path(G.E.X.__file__),G.CONTROL,G.PARENT,T.H.PARENT,T.H.BOUNDARY,ROOT/'results'/f'{T.S.TAG}.json']
    return dict(status='PASS_EXACT_ALL_LABEL_SHARED_X_DISCONNECTED_BOUNDARY_GATE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        deleted_walls=T.S.DELETE,boundary_cycle=T.S.CYCLE,exact_old_deletion_cost=parent['exact_old_deletion_cost'],
        mandatory_corners=[[0,1,2,3,5,6,7,8,'A'],[0,1,2,3,5,6,8,'A','X'],[0,1,2,3,5,6,30,63,'X']],
        X_frame=Q.RF,X_frame_gram=G.encmat(engine.gram),X_inverse_pivots=engine.pivots,
        X_kernel_seed=engine.seed,X_kernel_square=engine.kq.enc(),
        counts=dict(engine.counts,X_recipe_shapes=148,AX_cases=344,maximum_high_parameters=2),records=records,
        completeness='Finite H4 on A2+7+8 forces A right to all S6,7,8. X has one of the seven A2 profiles. Complete finite rank-five classification of the two required corners A2+30+63+X and A2+8+A+X gives148 X shapes/344 cases. Symbolic6 means every isolated finite I2 label>=6, never a cutoff. Exact closed high-cosine boxes cover both norm branches for X in8,30,63, then both A branches in7,8,X. Every leaf has a necessary native or prescribed-angle obstruction; unresolved inversions/signs fail rather than reject.',
        topology='X covers oldedge30-63 and also meets old8. The disks[7,8,X,30] and[8,44,63,X] require the three saved corners, so their proposed five-cap union is impossible while the retained oldwalls are fixed.',
        scope='An all-label obstruction for this disconnected-boundary X placement, with all seven A2 profiles. It does not exclude a zero-boundary-arc X dotted to30/63, other incidences, changed oldwalls or arbitrary finite-volume completion. The frozen exact interval engine is reused under an explicitly changed literal frame; no independent arithmetic implementation and no new chamber.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],flush=True)
