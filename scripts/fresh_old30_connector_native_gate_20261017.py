#!/usr/bin/env python3
"""Exact all-label Z gate over twenty selected literal prefix branches.

Uses inherited outward rational physical intervals. Numerical vectors select
explicit algebraic branches, not angle equalities or acceptance decisions.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_old30_connector_flip_20261017 as SC
import fresh_new_cavity_construction_gate_20261017 as Q
P,U,I,G,ROOT=Q.P,Q.U,Q.I,Q.G,Q.ROOT
TAG='fresh_old30_connector_native_gate_20261017'
TARGETS=[4,7,8,9,10,11,30,44,53,63]

class Gram:
    def __init__(self,base,ns):
        self.base,self.ns=base,ns;self.cp={i:base.sp.nodes[i]['cp']for i in TARGETS};self.recipes=[]
        # Auxiliary old normals may meet H4. Project BOTH arguments off
        # the entire S6, rather than applying the A2-only shortcut valid
        # when one argument is already rightH4.
        white=[0,1,2,3,5,6]
        wg=[[I(2)if i==j else I()for j in range(6)]for i in range(6)]
        for i,j,m in[(0,1,5),(1,2,3),(2,3,3),(4,5,3)]:wg[i][j]=wg[j][i]=-Q.theta(m)
        for i in range(6):
            for j in range(6):assert G.overlap(wg[i][j],G.E.X.literal(base.z,white[i],white[j]))
        wi,_=P.C53.inverse(wg)
        cov={i:[G.E.X.literal(base.z,i,j)for j in white]for i in TARGETS}
        co={i:U.imv(wi,cov[i])for i in TARGETS}
        self.g={(i,j):U.rnd(G.E.X.literal(base.z,i,j)-U.idot(cov[i],co[j]))for i in TARGETS for j in TARGETS}
        self.q={i:self.g[i,i]for i in TARGETS}
        for i in TARGETS:
            for j in TARGETS:assert G.overlap(self.g[i,j],base.sp.rp(i,j))
    def pair(self,i,j):return self.g[i,j]
    def norm_frame(self,frame,profile,contacts):
        cp=G.profile(profile);q=U.rnd(2-G.white(cp,cp));g=[[self.pair(i,j)for j in frame]for i in frame]
        inv,piv=P.C53.inverse(g);rhs=[U.rnd(v-G.white(cp,self.cp[i]))for v,i in zip(contacts,frame)]
        coeff=U.imv(inv,rhs);rem=U.rnd(q-U.idot(rhs,coeff));choices=[]
        for seed in TARGETS:
            if seed in frame:continue
            sc=[self.pair(seed,i)for i in frame];kc=U.imv(inv,sc);kq=U.rnd(self.q[seed]-U.idot(sc,kc))
            if kq.l>0:choices.append((float((kq.h-kq.l)/kq.l),seed,kc,kq))
        assert choices,('no_positive_complement',frame)
        _,seed,kc,kq=min(choices)
        # Positive complementary line + ambient Lorentz signature proves
        # the nondegenerate three-frame Lorentzian before testing rem.
        info=dict(frame=frame,profile=profile,native_contacts=P.enc(contacts),frame_gram=G.encmat(g),inverse_pivots=piv,kernel_seed=seed,kernel_square=kq.enc(),remaining_square=rem.enc())
        if rem.h<0:return [],info
        sc=Q.sqrt_closed(I(max(F(0),rem.l),rem.h)/kq);out=[]
        for sg in[-1,1]:
            row={}
            for j in self.cp:
                col=[self.pair(i,j)for i in frame];kb=kq if j==seed else U.rnd(self.pair(seed,j)-U.idot(kc,col))
                row[j]=U.rnd(U.idot(coeff,col)+sg*sc*kb)
            for i,v in zip(frame,rhs):row[i]=v
            out.append((sg,row))
        return out,info
    def add(self,name,frame,profile,labels,numeric):
        contacts=[-Q.theta(m)for m in labels];rows,info=self.norm_frame(frame,profile,contacts);assert len(rows)==2
        errors=[]
        for sg,row in rows:
            errs=[abs(float((row[i].l+row[i].h)/2)-float(numeric@SC.M@self.ns[i]))/(1+abs(float(numeric@SC.M@self.ns[i])))for i in TARGETS]
            errors.append((max(errs),sg,row))
        errors.sort(key=lambda x:x[0]);error,sign,row=errors[0]
        assert error<1e-6 and errors[1][0]>1e-4,('ambiguous_branch_selection',name,[(v,s)for v,s,r in errors])
        cp=G.profile(profile);q=U.rnd(2-G.white(cp,cp))
        for j,v in row.items():self.g[name,j]=self.g[j,name]=v
        self.g[name,name]=q;self.q[name]=q;self.cp[name]=cp
        info.update(name=name,selected_sign=sign,selection='The stated sign defines the exact branch; matching the saved numerical products is only a discovery sanity check.',numerical_relative_error=error)
        self.recipes.append(info)
    def reject(self,profile,labels,box):
        contacts=[-box if m==6 else-Q.theta(m)for m in labels];rows,info=self.norm_frame(['X','Y','R'],profile,contacts)
        if not rows:return dict(reason='negative_remaining_square',frame_certificate=info)
        cp=G.profile(profile);branches=[]
        for sg,row in rows:
            bad=[]
            for i in TARGETS:
                v=U.rnd(row[i]+G.white(cp,self.cp[i]));ex=Q.exclusion(v)
                if ex:bad.append(dict(wall=i,**ex))
            if not bad:return None
            branches.append(dict(sign=sg,rejection=bad[0]))
        return dict(reason='both_native_branches_excluded',frame_certificate=info,branches=branches)

def run():
    parent=ROOT/'results/fresh_interior_connector_disk_20261017.json';scout=ROOT/'results'/f'{SC.TAG}.json';doc=json.loads(parent.read_text());probe=json.loads(scout.read_text());base=Q.Engine()
    _,_,ns,_,_,_=SC.H.setup();records=[];counts=Counter();high=I(Q.theta(6).l,F(2))
    assert len(doc['seeds'])==20
    def cover(g,pi,ms,box=None,depth=0):
        counts['tested_boxes']+=1;out=g.reject(pi,ms,box)
        if out is not None:
            counts['certified_leaves']+=1;counts['maximum_depth']=max(counts['maximum_depth'],depth)
            return dict(cosine_box=box.enc()if box is not None else None,certificate=out)
        assert box is not None and depth<24,('unresolved_exact_native_gate',pi,ms,depth,box.enc()if box else None)
        mid=(box.l+box.h)/2
        return dict(split=str(mid),children=[cover(g,pi,ms,I(box.l,mid),depth+1),cover(g,pi,ms,I(mid,box.h),depth+1)])
    profile_counts=Counter()
    for i,rec in enumerate(doc['seeds']):
        a,c,x,y,r=map(np.array,rec['normals']);xp,yp=rec['X_profile'],rec['Y_profile'];profile_counts[xp,yp]+=1
        m8,maxx,mcx=rec['X_labels'];my8,myr,myx=rec['Y_labels'];g=Gram(base,ns)
        g.add('A',[7,8,53],0,[2,2,3],a)
        g.add('C',[7,30,'A'],0,[2,2,rec['AC']],c)
        g.add('R',[8,44,63],0,[2,rec['R44'],2],r)
        g.add('X',[8,'A','C'],xp,[m8,maxx,mcx],x)
        g.add('Y',[8,'R','X'],yp,[my8,myr,myx],y)
        cases=[]
        for pi,mx,my,mr in product(range(7),range(2,7),range(2,7),range(2,7)):
            if SC.P.corner3(xp,yp,pi,myx,mx,my)is None or SC.P.corner3(yp,0,pi,myr,my,mr)is None:continue
            ms=[mx,my,mr];nh=ms.count(6);assert nh<=1
            counts['high_families'if nh else'low_recipes']+=1
            proof=cover(g,pi,ms,high if nh else None)
            cases.append(dict(profile=pi,labels=ms,high_positions=[j for j,m in enumerate(ms)if m==6],proof=proof))
        records.append(dict(seed=i,exact_prefix_recipes=g.recipes,complete_Z_cases=cases))
        print('OLD30_NATIVE_GATE',i+1,dict(counts),flush=True)
    assert counts['low_recipes']==295 and counts['high_families']==2
    assert profile_counts=={(0,0):14,(0,1):2,(0,4):2,(1,0):2}
    paths=[Path(__file__),Path(SC.__file__),Path(SC.P.__file__),Path(Q.__file__),Path(Q.E.__file__),Path(G.__file__),Path(G.E.X.__file__),G.CONTROL,G.PARENT,SC.H.PARENT,SC.H.BOUNDARY,parent,scout]
    return dict(status='PASS_EXACT_ALL_LABEL_Z_NATIVE_GATE_ON_TWENTY_LITERAL_PREFIXES',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},counts=dict(counts),records=records,
      mandatory_corners=['S6+X+Y+Z','S6+Y+R+Z'],targets=TARGETS,predecessor_profile_counts={f'{a},{b}':c for(a,b),c in sorted(profile_counts.items())},
      completeness='For an actual S6 facet Z, finite extension of H4+A2 forces rightH4 and one of the seven A2 profiles. Both mandatory finite corners are classified as finite diagrams on A2 plus three roots. Any bond>=6 in a spherical diagram is an isolated I2 component, so replacing each such label by6 preserves the finite-domain test. The enumerated domain consists of295 low recipes and2 single high families. Each high magnitude is relaxed to the complete closed interval[sqrt3,2]; both norm branches are certified impossible by retained native contacts, or the remaining norm is negative after certifying a positive complementary line. Assigned-label gap search through64 only supplies necessary rejection witnesses and is not a candidate-label cutoff.',
      scope='Exact obstruction to the new Z duty on the twenty specified algebraic prefix branches. This does not certify all native contacts or completeness of the parent numerical prefix search, exclude continuously varied prefixes or other disk flips, or establish an ambient volume result. Physical interval and finite-diagram helpers are inherited; no new chamber is promoted.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],flush=True)
