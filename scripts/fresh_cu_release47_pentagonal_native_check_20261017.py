#!/usr/bin/env python3
"""Exact native-contact check of the six provisional pentagonal improvements.

The physical interval engine is inherited. Numerical vectors only select
explicit norm branches; they never establish a native contact equality.
The moving C wall is unnecessary for this local K/R obstruction.
"""
import argparse,json,math
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import numpy as np
import fresh_dotted_c10_domain_20261017 as D
import fresh_cu_release47_pentagonal_cap_20261017 as N

ROOT,P,U,I,Q,G=D.ROOT,D.P,D.U,D.I,D.Q,D.G
TAG='fresh_cu_release47_pentagonal_native_check_20261017'
TARGETS=[7,8,10,30,44,47,52,53]

class Gram:
    def __init__(self,d,ns):
        self.d,self.ns=d,dict(ns);self.recipes=[]
        self.cp={i:[D.C.literal(d['z'],i,j)for j in [5,6]]for i in TARGETS}
        self.g={(i,j):d['pair'](i,j)for i in TARGETS for j in TARGETS}
        self.q={i:self.g[i,i]for i in TARGETS}

    def add(self,name,frame,profile,labels,numeric):
        cp=G.profile(profile);q=U.rnd(2-G.white(cp,cp))
        gram=[[self.g[i,j]for j in frame]for i in frame]
        inv,piv=P.C53.inverse(gram)
        rhs=[U.rnd(-Q.theta(m)-G.white(cp,self.cp[i]))for i,m in zip(frame,labels)]
        coeff=U.imv(inv,rhs);rem=U.rnd(q-U.idot(rhs,coeff));assert rem.l>0
        choices=[]
        for seed in TARGETS:
            if seed in frame:continue
            cov=[self.g[seed,i]for i in frame];kc=U.imv(inv,cov)
            kq=U.rnd(self.q[seed]-U.idot(cov,kc))
            if kq.l>0:choices.append((float((kq.h-kq.l)/kq.l),seed,kc,kq))
        assert choices,('no certified positive complement',name)
        _,seed,kc,kq=min(choices);scale=Q.sqrt_closed(rem/kq)
        branches=[]
        for sign in [-1,1]:
            row={}
            for j in self.cp:
                col=[self.g[i,j]for i in frame]
                kb=kq if j==seed else U.rnd(self.g[seed,j]-U.idot(kc,col))
                row[j]=U.rnd(U.idot(coeff,col)+sign*scale*kb)
            for i,v in zip(frame,rhs):row[i]=v
            error=max(abs(float((row[i].l+row[i].h)/2)-float(numeric@N.M@self.ns[i]))/(1+abs(float(numeric@N.M@self.ns[i])))for i in TARGETS)
            branches.append((error,sign,row))
        branches.sort(key=lambda r:r[0]);error,sign,row=branches[0]
        assert error<1e-6 and branches[1][0]>1e-4,('ambiguous exact branch',name,[(a,b)for a,b,_ in branches])
        for j,v in row.items():self.g[name,j]=self.g[j,name]=v
        self.g[name,name]=q;self.q[name]=q;self.cp[name]=cp;self.ns[name]=numeric
        self.recipes.append(dict(name=name,frame=frame,profile=profile,labels=labels,
            gram=G.encmat(gram),inverse_pivots=piv,projected_contacts=P.enc(rhs),
            kernel_seed=seed,kernel_square=kq.enc(),remaining_square=rem.enc(),
            selected_sign=sign,positive_scale=scale.enc(),numerical_branch_error=error))

def run():
    p=ROOT/'results'/f'{N.TAG}.json';doc=json.loads(p.read_text())
    assert len(doc['improvements'])==6 and not doc['unresolved']
    d=D.data();_,_,ns,_,_,_,*_=N.U.setup();records=[]
    for imp in doc['improvements']:
        assert imp['stage']=='R'
        r=doc['R_caps'][imp['index']];e=doc['events'][r['event']]
        l=doc['L_seeds'][e['L']];k=doc['K_seeds'][e['K']];pseed=doc['P_seeds'][e['P']]
        assert pseed['profile']==0 and pseed['labels']==[3,2,2]
        assert np.max(np.abs(np.array(pseed['normal'])-ns[47]))<1e-8
        assert e['family']in[25,27,29] and e['profiles']['68']==0 and e['AK']==2
        assert l['profile']in[1,4] and l['labels']==[3,2,2]
        assert k['profile']==0 and k['labels']==[4,2,4]
        assert r['profile']==0 and r['labels']==[2,2,2]
        gr=Gram(d,ns)
        gr.add('L',[30,52,47],l['profile'],[3,2,2],np.array(l['normal']))
        gr.add('K',[30,'L',47],0,[4,2,4],np.array(k['normal']))
        # All three moving-A families prescribe the same fixed A7/A30/AK
        # contacts. This norm frame defines their common local branch.
        gr.add('A',[7,30,'K'],0,[3,2,2],np.array(e['normals']['68']))
        gr.add('R',[7,'A',47],0,[2,2,2],np.array(r['normal']))
        kr=gr.g['K','R'];assert kr.l>-2 and kr.h<0
        nominal=round(math.pi/math.acos(-float((kr.l+kr.h)/2)/2))
        gaps=[]
        for m in range(max(2,nominal-3),nominal+4):
            lower,upper=-Q.theta(m+1),-Q.theta(m)
            if lower.h<kr.l and kr.h<upper.l:gaps.append((m,lower,upper))
        assert len(gaps)==1,('not a certified angle gap',kr.enc(),nominal)
        m,lower,upper=gaps[0]
        records.append(dict(R_index=r['index'],event=e['index'],family=e['family'],
            u=e['u'],L_profile=l['profile'],numerical_KR_label=r['KR_relation'],
            norm_recipes=gr.recipes,native_KR=kr.enc(),strict_angle_gap=[m,m+1],
            lower_threshold=lower.enc(),upper_threshold=upper.enc(),
            conclusion='This literal branch is not a Coxeter chamber; the unprescribed K/R product lies strictly between consecutive native finite angles.'))
        print('PENTAGON_EXACT_KR',r['index'],float(kr.l),float(kr.h),'gap',m,m+1,flush=True)
    paths=[Path(__file__),Path(N.__file__),p,Path(D.__file__),Path(D.C.__file__),Path(Q.__file__),Path(G.__file__),
           ROOT/'results'/f'{D.TAG}.json',ROOT/'results'/f'{D.C.TAG}.json']
    return dict(status='PASS_EXACT_PENTAGONAL_PROVISIONAL_KR_REJECTIONS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        records=records,exact_improvements=0,
        scope='All six saved provisional nine-end records, each reconstructed by exact fixed L/K/A/R norm recipes. P3 is old47. Numerical data choose explicit signs only. The inherited physical interval engine is reused. This local native obstruction does not certify the remaining scout chambers or exclude other pentagonal geometry, labels or moving parameters. The exact63-wall/199-finite/ten-end infinite-volume control is unchanged.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],len(out['records']),flush=True)
